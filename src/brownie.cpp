/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
 *   Copyright (C) 2014, 2015 Mahdi Heydari (mahdi.heydari@intec.ugent.be) *
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
#include "graph.h"
#include "kmertable.h"

#include <string>

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
        cout << "Parsed input files (" << Util::stopChronoStr() << ")" << endl;

        cout << "Total number of unique kmers in table: "
        << readParser->getNumKmers() << " ("
        << readParser->getNumKmersCovGTOne() << " with coverage > 1)" << endl;

        #ifdef DEBUG
        readParser->validateStage1();
        #endif

        // write kmers file containing all kmers with cov > 1
        cout << "Writing kmer file...";
        cout.flush();
        Util::startChrono();
        readParser->writeKmersWithCovGTOne(getKmerFilename());
        cout << "done (" << Util::stopChronoStr() << ")" << endl;

        delete readParser;
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
        // STAGE 3 : MULTIPLICITY COUNTING
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
        cout << "Creating graph... ";
        cout.flush();
        graph.createFromFile(getNodeFilename(2),
                             getArcFilename(2),
                             getMetaDataFilename(2));
        cout << "done (" << graph.getNumNodes() << " nodes, "
        << graph.getNumArcs() << " arcs)" << endl;

        Util::startChrono();
        graph.countNodeandArcFrequency(libraries);
        cout << "Done counting multiplicity (" << Util::stopChronoStr() << ")" << endl;

        cout << "Extracting graph..." << endl;
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

        cout << "Entering stage 4" << endl;
        cout << "================" << endl;

       /* if (!stageFourNecessary()) {
                cout << "Files produced by this stage appear to be present, "
                "skipping stage 4..." << endl << endl;
                return;
        }*/
        DBGraph testgraph(settings);
        testgraph.createFromFile(getNodeFilename(3),
                                 getArcFilename(3),
                                 getMetaDataFilename(3));
        cout<<"initial kmerCoverage : "<<testgraph.estimatedKmerCoverage<<endl;
        string command="rm cov/*.pdf && rm cov/*.dat";
        system(command.c_str());
        testgraph.clipTips(0);
        testgraph.mergeSingleNodes();
        testgraph.filterCoverage(0);
        testgraph.mergeSingleNodes();
        testgraph.extractStatistic(1);

        DBGraph graph(settings);
        graph.estimatedArcCoverageMean= testgraph.estimatedArcCoverageMean;
        graph.estimatedArcCoverageSTD=sqrt(graph.estimatedArcCoverageMean);
        graph.estimatedKmerCoverage=testgraph.estimatedKmerCoverage;
        graph.estimatedMKmerCoverageSTD=sqrt(graph.estimatedKmerCoverage);
        cout<<"estimated Arc Coverage Mean: "<<graph.estimatedArcCoverageMean<<endl;
        cout<<"estimated Arc Coverage STD: "<<graph.estimatedArcCoverageSTD<<endl;
        cout<<"estimated Kmer Coverage Mean: "<<graph.estimatedKmerCoverage<<endl;
        cout<<"estimated Kmer Coverage STD: "<<graph.estimatedMKmerCoverageSTD<<endl;
        testgraph.clear();

        Util::startChrono();
        cout << "Creating graph... ";
        cout.flush();
        graph.createFromFile(getNodeFilename(3),
                             getArcFilename(3),
                             getMetaDataFilename(3));
        string nodeFileName=getNodeFilename(4);
        string arcFileName=getArcFilename(4);
        string metaDataFileName=getMetaDataFilename(4);

        cout << "done (" << graph.getNumNodes() << " nodes, "
        << graph.getNumArcs() << " arcs)" << endl;
        cout << "Created graph in "
        << Util::stopChrono() << "s." << endl;
        Util::startChrono();
        #ifdef DEBUG
        graph.compareToSolution();
        #endif
        //variables
        bool simplified = true;
        size_t bigestN50=0;
        int round=1;
        size_t minN50=100;
        int maxSizeNodeForDel=100;


        while (simplified) {
                //*******************************************************
                bool tips=graph.clipTips(round);
                if (tips) {
                        graph.mergeSingleNodes();
                        if(! graph.continueEdit(bigestN50, nodeFileName,arcFileName,metaDataFileName))
                                break;
                        #ifdef DEBUG
                        graph.compareToSolution();
                        #endif
                }
                //*******************************************************
                bool coverage = graph.filterCoverage(round);

                if (coverage) {
                        graph.mergeSingleNodes();
                        if (!graph.continueEdit(bigestN50, nodeFileName,arcFileName,metaDataFileName))
                                break;
                        #ifdef DEBUG
                        graph.compareToSolution();
                        #endif
                }
                //*******************************************************
                bool chimeric = graph.filterChimeric( round);

                if (chimeric) {
                        graph.mergeSingleNodes();
                        if (!graph.continueEdit(bigestN50, nodeFileName,arcFileName,metaDataFileName))
                                break;
                        #ifdef DEBUG
                        graph.compareToSolution();
                        #endif

                }
                //*******************************************************
                bool bubble=false;
                bubble= graph.bubbleDetection(round);

                if (bubble) {
                        graph.mergeSingleNodes();
                        if (!graph.continueEdit(bigestN50, nodeFileName,arcFileName,metaDataFileName))
                                break;
                        #ifdef DEBUG
                        graph.compareToSolution();
                        #endif

                }
                //*******************************************************

                graph.extractStatistic(1);
                #ifdef DEBUG
                cout<<"estimated Arc Coverage Mean: "<<graph.estimatedArcCoverageMean<<endl;
                cout<<"estimated Arc Coverage STD: "<<graph.estimatedArcCoverageSTD<<endl;
                cout<<"estimated Kmer Coverage Mean: "<<graph.estimatedKmerCoverage<<endl;
                cout<<"estimated Kmer Coverage STD: "<<graph.estimatedMKmerCoverageSTD<<endl;
                #endif
                //*******************************************************
                bool link=false;
                if (bigestN50>minN50) {
                        link =graph.removeIncorrectLink();
                        #ifdef DEBUG
                        graph.compareToSolution();
                        #endif
                }

                if (link) {
                        graph.mergeSingleNodes();
                        if (!graph.continueEdit(bigestN50, nodeFileName,arcFileName,metaDataFileName))
                                break;
                        #ifdef DEBUG
                        graph.compareToSolution();
                        #endif

                }
                //*******************************************************
                bool deleted=false;
                if (bigestN50>minN50) {
                        deleted=graph.deleteUnreliableNodes(maxSizeNodeForDel, round);
                        #ifdef DEBUG
                        graph.compareToSolution();
                        #endif
                }

                if (deleted) {
                        graph.mergeSingleNodes();
                        if (!graph.continueEdit(bigestN50, nodeFileName,arcFileName,metaDataFileName))
                                break;
                        #ifdef DEBUG
                        graph.compareToSolution();
                        #endif
                }
                //*******************************************************
                simplified = tips  || chimeric  || coverage || link|| deleted || bubble;
                round++;
        }

        cout << " Ghraph correction completed in "
        << Util::stopChrono() << "s." << endl;
        Util::startChrono();
        #ifdef DEBUG
        graph.sanityCheck();
        command="pdftk cov/*.pdf cat output allpdfFiles.pdf";
        system(command.c_str());
        #endif
        graph.clear();
        cout << "Stage 4 finished.\n" << endl;
}

void Brownie::stageFive()
{
        cout << "Entering stage 5" << endl;
        cout << "================" << endl;
        DBGraph graph(settings);
        Util::startChrono();
        cout << "Creating graph... ";
        cout.flush();
        graph.createFromFile(getNodeFilename(4),
                             getArcFilename(4),
                             getMetaDataFilename(4));
        cout << "done (" << graph.getNumNodes() << " nodes, "
        << graph.getNumArcs() << " arcs)" << endl;
        cout << "Created graph in "
        << Util::stopChrono() << "s." << endl;
        #ifdef DEBUG
        graph.compareToSolution();
        #endif
        //graph.writeCytoscapeGraph(0);
        Util::startChrono();
        graph.errorCorrection(libraries);
        cout << "Error correction completed in "
        << Util::stopChrono() << "s." << endl;
        cout << "Stage 5 finished.\n" << endl;
        graph.clear();
}

int main(int argc, char** args)
{
        try {
                Brownie brownie(argc, args);

                cout << "Welcome to Brownie v." << BROWNIE_MAJOR_VERSION << "."
                << BROWNIE_MINOR_VERSION << "." << BROWNIE_PATCH_LEVEL << endl;
                cout << "Today is " << Util::getTime() << endl;
                brownie.stageOne();
                brownie.stageTwo();
                brownie.stageThree();
                brownie.stageFour();
                brownie.stageFive();
        } catch (exception &e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }

        cout << "Exiting... bye!" << endl << endl;
        return EXIT_SUCCESS;
}

