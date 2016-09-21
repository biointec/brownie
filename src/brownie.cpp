/***************************************************************************
 *   Copyright (C) 2014 - 2016 Jan Fostier (jan.fostier@intec.ugent.be)    *
 *   This file is part of Brownie                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "brownie.h"
#include "settings.h"
#include "kmeroverlap.h"
#include "kmeroverlaptable.h"
#include "kmertable.h"
#include "readcorrection.h"
#include "refcomp.h"

using namespace std;

Brownie::Brownie(int argc, char** args)
{
        settings.parseCommandLineArguments(argc, args, libraries);
        Kmer::setWordSize(settings.getK());
        RKmer::setWordSize(settings.getK() - KMERBYTEREDUCTION * 4);
}

void Brownie::stageOne()
{
        // ============================================================
        // STAGE 1 : PARSE THE READS
        // ============================================================

        cout << "Entering stage 1" << endl;
        cout << "================" << endl;
        if (!stageOneNecessary()) {
                cout << "Files produced by this stage appear to be present, "
                "skipping stage 1..." << endl << endl;
                return;
        }
        KmerTable *readParser = new KmerTable(settings);
        cout << "Generating kmers with k = " << Kmer::getK()
             << " from input files..." << endl;
        Util::startChrono();
        readParser->parseInputFiles(libraries);
        size_t kmerGOne = readParser->getNumKmersCovGTOne();
        size_t allKmers = readParser->getNumKmers() ;
        cout << "Parsed input files (" << Util::stopChronoStr() << ")" << endl;
        cout << "Total number of unique kmers in table: "
             << allKmers << " (" << kmerGOne << " with coverage > 1)" << endl;

#ifdef DEBUG
        readParser->validateStage1();
#endif

        // write kmers file containing all kmers with cov > 1
        cout << "Writing kmer file...";
        cout.flush();
        Util::startChrono();
        //readParser->writeAllKmers(getKmerFilename());
        readParser->writeKmersWithCovGTOne(getKmerFilename());
        cout << "done (" << Util::stopChronoStr() << ")" << endl;

        delete readParser;

        // write metadata for all libraries
        libraries.writeMetadata(settings.getTempDirectory());

        cout << "Stage 1 finished.\n" << endl;
}

void Brownie::stageTwo()
{
        // ============================================================
        // STAGE 2 : KMER OVERLAP TABLE
        // ============================================================

        cout << "Entering stage 2" << endl;
        cout << "================" << endl;

        if (!stageTwoNecessary()) {
                cout << "Files produced by this stage appear to be present, "
                "skipping stage 2..." << endl << endl;
                return;
        }

        // create a kmer table from the reads
        KmerOverlapTable overlapTable(settings);
        Util::startChrono();
        cout << "Building kmer overlap table...";
        overlapTable.loadKmersFromDisc(getKmerFilename());
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << "Number of kmers loaded: " << overlapTable.size() << endl;

        // find the overlap between kmers
        Util::startChrono();
        cout << "Finding overlaps between kmers..." << endl;
        overlapTable.parseInputFiles(libraries);
        cout << "Done building overlap table ("
             << Util::stopChronoStr() << ")" << endl;
        cout << "Overlap table contains " << overlapTable.size()
             << " nodes" << endl;

#ifdef DEBUG
        overlapTable.validateStage2();
#endif

        // extract nodes and arcs from the kmer table
        overlapTable.extractNodes(getNodeFilename(2),
                                  getArcFilename(2),
                                  getMetaDataFilename(2));

        overlapTable.clear();   // clear memory !
        cout << "Stage 2 finished.\n" << endl;
}

void Brownie::stageThree()
{
        // ============================================================
        // STAGE 3 : MULTIPLICITY COUNTING - KMER SPECTRUM FITTING
        // ============================================================

        cout << "Entering stage 3" << endl;
        cout << "================" << endl;

        if (!stageThreeNecessary()) {
                cout << "Files produced by this stage appear to be present, "
                "skipping stage 3..." << endl << endl;
                return;
        }

        // build pre graph and simplify it
        DBGraph graph(settings);

        Util::startChrono();
        cout << "Loading graph... ";
        cout.flush();
        graph.loadFromFile(getNodeFilename(2),
                             getArcFilename(2),
                             getMetaDataFilename(2));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << graph.getGraphStats() << endl;

        Util::startChrono();
        graph.generateKmerSpectrum(settings.getTempDirectory(), libraries);
        cout << "Done generating spectrum (" << Util::stopChronoStr() << ")" << endl;

        cout << "Writing graph..." << endl;
        graph.writeGraphBin(getBinNodeFilename(3),
                            getBinArcFilename(3),
                            getMetaDataFilename(3));

#ifdef DEBUG
        graph.sanityCheck();
#endif
        graph.clear();
        cout << "Stage 3 finished.\n" << endl;
}

void Brownie::stageFour()
{
        // ============================================================
        // STAGE 4 : GRAPH SIMPLIFICATION
        // ============================================================

        cout << "Entering stage 4" << endl;
        cout << "================" << endl;

        if (!stageFourNecessary()) {
                cout << "Files produced by this stage appear to be present, "
                        "skipping stage 4..." << endl << endl;
                return;
        }

        DBGraph graph(settings);
        Util::startChrono();
        cout << "Creating graph... "; cout.flush();
        graph.loadGraphBin(getBinNodeFilename(3),
                           getBinArcFilename(3),
                           getMetaDataFilename(3));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << graph.getGraphStats() << endl;

        graph.loadKmerSpectrumFit(getSpectrumFitFilename());
        double cutoff = graph.getCovCutoff();

#ifdef DEBUG
{
        Util::startChrono();
        cout << "Building kmer - node/position index... "; cout.flush();
        graph.buildKmerNPPTable();      // build kmer-NPP index
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        RefComp refComp("genome.fasta");
        refComp.validateGraph(graph);
        vector<size_t> trueMult;
        refComp.getNodeMultiplicity(graph, trueMult);
        graph.setTrueNodeMultiplicity(trueMult);
}
#endif

        // TIP CLIPPING
        Util::startChrono();
        cout << "Cleaning graph (tips, cov-cutoff = " << cutoff
             << ", lmax = " << libraries.getAvgReadLength() << ")\n";
        while (graph.clipTips(cutoff, libraries.getAvgReadLength())) {
                graph.concatenateNodes();
                cout << "\tGraph contains " << graph.getNumValidNodes() << " nodes" << endl;
        }
        cout << "Done (" << Util::stopChronoStr() << ")\n" << endl;

        // BUBBLE DETECTION
        Util::startChrono();
        cout << "Cleaning graph (bubbles, cov-cutoff = " << cutoff
             << ", lmax = " << libraries.getAvgReadLength() << ", maxvisits = "
             << settings.getBubbleDFSNodeLimit() << ", threads = "
             << settings.getNumThreads() << ")\n";
        while (graph.bubbleDetection(cutoff, libraries.getAvgReadLength())) {
                graph.concatenateNodes();
                cout << "\tGraph contains " << graph.getNumValidNodes() << " nodes" <<  endl;
        }
        cout << "Done (" << Util::stopChronoStr() << ")\n" << endl;

        // FLOW CORRECTION
        Util::startChrono();
        cout << "Cleaning graph (flow correction)\n";
        while (graph.flowCorrection()) {
                graph.concatenateNodes();
                cout << "\tGraph contains " << graph.getNumValidNodes() << " nodes" <<  endl;
        }
        cout << "Done (" << Util::stopChronoStr() << ")\n" << endl;

        //graph.writeCytoscapeGraph(settings.getTempDirectory() + "tip", 24450, 3);

#ifdef DEBUG
        Util::startChrono();
        cout << "Building kmer - node/position index... "; cout.flush();
        graph.buildKmerNPPTable();      // build kmer-NPP index
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        RefComp refComp("genome.fasta");
        refComp.validateGraph(graph);
        vector<size_t> trueMult;
        refComp.getNodeMultiplicity(graph, trueMult);
        graph.setTrueNodeMultiplicity(trueMult);

        for (int i = 1; i < graph.getNumNodes(); i++) {
                SSNode node = graph.getSSNode(i);
                if (!node.isValid())
                        continue;
                if (trueMult[i] == 0)
                        cout << "Node " << i << " with avgKmerCov " << node.getAvgKmerCov() << " is false." << endl;
        }
#endif

        cout << graph.getGraphStats() << endl;
        cout << "Writing graph..." << endl;
        graph.writeGraph(getNodeFilename(4),
                         getArcFilename(4),
                         getMetaDataFilename(4));

#ifdef DEBUG
        graph.sanityCheck();
#endif
        graph.clear();
        cout << "Stage 4 finished.\n" << endl;
}

void Brownie::stageFive()
{
        cout << "Entering stage 5" << endl;
        cout << "================" << endl;

        if (!stageFiveNecessary()) {
                cout << "Files produced by this stage appear to be present, "
                        "skipping stage 5..." << endl << endl;
                return;
        }

        // Build a DBG from stage 4 files on disk
        DBGraph graph(settings);
        Util::startChrono();
        cout << "Creating graph... "; cout.flush();
        graph.loadFromFile(getNodeFilename(4),
                             getArcFilename(4),
                             getMetaDataFilename(4));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << "Graph contains " << graph.getNumNodes() << " nodes and "
             << graph.getNumArcs() << " arcs" << endl;

        cout << "Writing cytoscape graph: " << endl;
        graph.writeCytoscapeGraph(settings.getTempDirectory() + "stage4", 1, 5);

        /*Util::startChrono();
        ReadCorrectionHandler rcHandler(graph, settings);
        rcHandler.doErrorCorrection(libraries);

        cout << "Error correction completed in " << Util::stopChronoStr() << endl;*/
        cout << "Stage 5 finished\n" << endl;
        graph.clear();
}

void Brownie::stageSix()
{
        cout << "Entering stage 6" << endl;
        cout << "================" << endl;

        // Build a DBG from stage 4 files on disk
        DBGraph graph(settings);
        Util::startChrono();
        cout << "Creating graph... "; cout.flush();
        graph.loadFromFile(getNodeFilename(4),
                             getArcFilename(4),
                             getMetaDataFilename(4));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << graph.getGraphStats() << endl;

        vector<NodeChain> trueNodeChain;

#ifdef DEBUG
{
        Util::startChrono();
        cout << "Building kmer - node/position index... "; cout.flush();
        graph.buildKmerNPPTable();      // build kmer-NPP index
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        RefComp refComp("genome.fasta");
        refComp.validateGraph(graph);
        refComp.getTrueNodeChain(graph, trueNodeChain);
        vector<size_t> trueMult;
        refComp.getNodeMultiplicity(graph, trueMult);
        graph.setTrueNodeMultiplicity(trueMult);
}
#endif

        Util::startChrono();

        graph.loadNodeChainContainer(libraries, trueNodeChain);
        graph.writeCytoscapeGraph(settings.getTempDirectory() + "cytRed");

#ifdef DEBUG
        graph.sanityCheck();
#endif
        cout << graph.getGraphStats() << endl;

        graph.writeGraphFasta();

        cout << "Repeat resolution completed in " << Util::stopChronoStr() << endl;
        cout << "Stage 6 finished\n" << endl;
        graph.clear();
}

void Brownie::writeGraphFasta()
{
        DBGraph graph(settings);
        if (settings.getSkipStage4()) {
                graph.loadGraphBin(getBinNodeFilename(3),
                        getBinArcFilename(3),
                        getMetaDataFilename(3));
                graph.writeGraph(getNodeFilename(4),
                                 getArcFilename(4),
                                 getMetaDataFilename(4));
                graph.writeGraphFasta();
        } else if (settings.getSkipStage5()) {
                graph.loadFromFile(getNodeFilename(4),
                        getArcFilename(4),
                        getMetaDataFilename(4));
                graph.writeGraphFasta();
        }
        graph.clear();
}

int main(int argc, char** args)
{
        try {
                Brownie brownie(argc, args);

                cout << "Welcome to Brownie v." << BROWNIE_MAJOR_VERSION << "."
                     << BROWNIE_MINOR_VERSION << "." << BROWNIE_PATCH_LEVEL;
#ifdef DEBUG
                cout << " (debug mode)" << endl;
#else
                cout << " (release mode)" << endl;
#endif
                cout << "Today is " << Util::getDateTime() << endl;

                brownie.stageOne();
                brownie.stageTwo();
                brownie.stageThree();
                brownie.stageFour();
                brownie.stageFive();
                exit(0);
                brownie.stageSix();

                brownie.writeGraphFasta();
        } catch (exception &e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }
        cout << "Exiting... bye!" << endl << endl;
        return EXIT_SUCCESS;
}
