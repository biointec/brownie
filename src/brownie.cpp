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
#include "findGap.h"
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

        cout << "\nEntering stage 1" << endl;
        cout << "================" << endl;

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
        if (Kmer::getK()>31)
                readParser->writeAllKmers(getKmerFilename());
        else
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

        cout << "\nEntering stage 2" << endl;
        cout << "================" << endl;

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

        cout << "\nEntering stage 3" << endl;
        cout << "================" << endl;

        // build pre graph and simplify it
        DBGraph graph(settings);

        Util::startChrono();
        cout << "Loading graph... ";
        cout.flush();
        graph.loadGraph(getNodeFilename(2),
                        getArcFilename(2),
                        getMetaDataFilename(2));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << graph.getGraphStats() << endl;

        Util::startChrono();
        graph.generateKmerSpectrum(settings.getTempDirectory(), libraries);
        cout << "Done generating spectrum (" << Util::stopChronoStr() << ")" << endl;

        cout << "Writing graph..." << endl;
        graph.writeGraph(getNodeFilename(3),
                         getArcFilename(3),
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

        cout << "\nEntering stage 4" << endl;
        cout << "================" << endl;

        DBGraph graph(settings);
        Util::startChrono();
        cout << "Creating graph... "; cout.flush();
        graph.loadGraph(getNodeFilename(3),
                        getArcFilename(3),
                        getMetaDataFilename(3));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << graph.getGraphStats() << endl;
        double cutoff = settings.getCutOffValue();
        graph.loadKmerSpectrumFit(getSpectrumFitFilename());
        if (cutoff == 0)
                cutoff = graph.getCovCutoff();
        FindGap  findGap (libraries, settings, graph);

        #ifdef DEBUG

        Util::startChrono();
        cout << "Building kmer - node/position index... "; cout.flush();
        graph.buildKmerNPPTable();      // build kmer-NPP index
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        RefComp refComp("genome.fasta");
        // refComp.validateGraph(graph);
        vector<size_t> trueMult;
        refComp.getNodeMultiplicity(graph, trueMult);
        graph.setTrueNodeMultiplicity(trueMult);
        #endif
        graph.buildKmerNPPTable();      // build kmer-NPP index
        graph.findbreakpoints("breakpoints.fasta");
        bool change = true;
        cutoff = cutoff*2/3;

        graph.buildKmerNPPTable();      // build kmer-NPP index
        graph.findbreakpoints("breakpoints.fasta");
        while (change){
                change = false;
                // TIP CLIPPING
                Util::startChrono();
                size_t lmax = libraries.getAvgReadLength();
                if (lmax < settings.getK()*2)
                        lmax = settings.getK()*2;
                cout << "Cleaning graph (tips, cov-cutoff = " << cutoff
                << ", lmax = " << lmax << ")\n";
                while (graph.clipTips(cutoff, lmax)) {
                        graph.concatenateNodes();
                        cout << "\tGraph contains " << graph.getNumValidNodes() << " nodes" << endl;
                        change = true;
                }
                cout << "Done (" << Util::stopChronoStr() << ")\n" << endl;
                // BUBBLE DETECTION
                Util::startChrono();
                lmax = libraries.getAvgReadLength();
                if (lmax < settings.getK()*2 )
                        lmax = settings.getK()*2+5;
                bool bubbleChange =true;
                while (bubbleChange){
                        size_t len = settings.getK() +5;
                        bubbleChange = false;
                        while (len <=lmax){
                                cout << "Cleaning graph (bubbles, cov-cutoff = " << cutoff
                                << ", lmax = " << len << ", maxvisits = "
                                << settings.getBubbleDFSNodeLimit() << ", threads = "
                                << settings.getNumThreads() << ")\n";
                                if (graph.bubbleDetection(cutoff, len)){
                                        graph.concatenateNodes();
                                        cout << "\tGraph contains " << graph.getNumValidNodes() << " nodes" <<  endl;
                                        change = true;
                                        bubbleChange = true;
                                }
                                else
                                        break;
                                len = len + settings.getK();
                        }
                }
                cout << "Done (" << Util::stopChronoStr() << ")\n" << endl;
                graph.extractStatistic();
                graph.removeChimericLinksByFlow(cutoff, libraries.getAvgReadLength());
                graph.concatenateNodes();
        }

        graph.buildKmerNPPTable();      // build kmer-NPP index
        graph.findbreakpoints("breakpoints.fasta");
        #ifdef DEBUG
        findGap.extractBreakpointSubgraph("breakpoints.fasta", settings.getTempDirectory()+"Stage4_");
        Util::startChrono();
        cout << "Building kmer - node/position index... "; cout.flush();
        graph.buildKmerNPPTable();      // build kmer-NPP index
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        refComp.getNodeMultiplicity(graph, trueMult);
        graph.setTrueNodeMultiplicity(trueMult);
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

        cout << "\nEntering stage 5" << endl;
        cout << "================" << endl;

        // Build a DBG from stage 4 files on disk
        DBGraph graph(settings);
        Util::startChrono();
        cout << "Creating graph... "; cout.flush();
        graph.loadGraph(getNodeFilename(4),
                        getArcFilename(4),
                        getMetaDataFilename(4));

        cout << "Graph contains " << graph.getNumNodes() << " nodes and "
             << graph.getNumArcs() << " arcs" << endl;
        graph.buildKmerNPPTable();
        ReadCorrectionHandler rcHandler(graph, settings);
        rcHandler.doErrorCorrection(libraries);

        cout << "Error correction completed in " << Util::stopChronoStr() << endl;
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
        graph.loadGraph(getNodeFilename(4),
                        getArcFilename(4),
                        getMetaDataFilename(4));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << graph.getGraphStats() << endl;
        graph.loadKmerSpectrumFit(getSpectrumFitFilename());

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

void Brownie::assembleModule()
{
        if (stageOneNecessary())
                stageOne();
        else
                cout << "Files produced by this stage appear to"
                        " be present, skipping stage 1...\n";
        if (stageTwoNecessary())
                stageTwo();
        else
                cout << "Files produced by this stage appear to"
                        " be present, skipping stage 2...\n";
        if (stageThreeNecessary())
                stageThree();
        else
                cout << "Files produced by this stage appear to"
                        " be present, skipping stage 3...\n";
        if (stageFourNecessary())
                stageFour();
        else
               cout << "Files produced by this stage appear to"
                        " be present, skipping stage 4...\n";
        //if (stageFiveNecessary())
                stageFive();
        //else
                cout << "Files produced by this stage appear to"
                        " be present, skipping stage 5...\n";
        /*if (stageSixNecessary())
                stageSix();
        else
                cout << "Files produced by this stage appear to"
                        " be present, skipping stage 6...\n";*/
}

void Brownie::visualizeModule()
{
        // Build a DBG from stage 4 files on disk
        DBGraph graph(settings);
        Util::startChrono();
        cout << "Creating graph... "; cout.flush();
        graph.loadGraphBin(getBinNodeFilename(3),
                           getBinArcFilename(3),
                           getMetaDataFilename(3));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << graph.getGraphStats() << endl;

        graph.writeCytoscapeGraph(settings.getTempDirectory() + "graph", 1, 10);

        graph.clear();
}

void Brownie::compareModule()
{
        int stage = settings.getRunSpecificStage();

        // load the DBG
        DBGraph graph(settings);
        Util::startChrono();
        cout << "Loading stage " << stage << " graph... "; cout.flush();
        graph.loadGraph(getNodeFilename(stage),
                        getArcFilename(stage),
                        getMetaDataFilename(stage));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << graph.getGraphStats() << "\n\n";

        // load the reference filename
        Util::startChrono();
        cout << "Loading reference file " << settings.getReferenceFilename() << "..."; cout.flush();
        RefComp refComp(settings.getReferenceFilename());
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << "Reference file contains " << refComp.getNumContigs()
             << " contigs with a total size of " << refComp.getSize() << " bp.\n";

        Util::startChrono();
        cout << "Building kmer - node/position index... "; cout.flush();
        graph.buildKmerNPPTable();      // build kmer-NPP index
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        refComp.validateGraph(graph);
        vector<size_t> trueMult;
        refComp.getNodeMultiplicity(graph, trueMult);

        cout << "Compare module finished.\n" << endl;
}

void Brownie::run()
{
        switch (settings.getCommand()) {
                case Command::assemble:
                        assembleModule();
                        break;
                case Command::compare:
                        compareModule();
                        break;
                case Command::visualize:
                        visualizeModule();
                        break;
                case Command::none:
                        cerr << "brownie: no command specified\n";
                        cerr << "Try 'brownie --help' for more information" << endl;
                        exit(EXIT_FAILURE);
                        break;
        }
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
                brownie.run();
        } catch (exception &e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }

        cout << "Exiting... bye!" << endl << endl;
        return EXIT_SUCCESS;
}
