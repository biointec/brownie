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
#include "correction.h"
#include "kmertable.h"
#include "alignment.h"
#include <string>

using namespace std;

Brownie::Brownie(int argc, char** args)
{
        settings.parseCommandLineArguments(argc, args, libraries);

        Kmer::setWordSize(settings.getK());
        RKmer::setWordSize(settings.getK() - KMERBYTEREDUCTION * 4);
}

void Brownie::printInFile()
{
        string konsole = settings.getTempDirectory()+ "/konsole.txt";
        freopen(konsole.c_str(), "a", stdout);
        cout<<"-------------------------------------------------------------------------------"<<endl;
        time_t t;
        time(&t);
        cout << ctime(&t) << endl; // e.g., Fri May 02 17:57:14 2003
        cout<<"-------------------------------------------------------------------------------"<<endl;
        cout << "Welcome to Brownie\n" << endl;
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
        size_t kmerGOne= readParser->getNumKmersCovGTOne();
        size_t allKmers=readParser->getNumKmers() ;
        Util::startChrono();
        readParser->parseInputFiles(libraries);
        cout << "Parsed input files (" << Util::stopChronoStr() << ")" << endl;
        cout << "Total number of unique kmers in table: "
        << allKmers << " ("
        << kmerGOne << " with coverage > 1)" << endl;

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
        graph.writeGraphBin(getBinNodeFilename(3),
                            getBinArcFilename(3),
                            getMetaDataFilename(3));

#ifdef DEBUG
        graph.sanityCheck();
#endif
        graph.clear();
        cout << "Stage 3 finished.\n" << endl;
        if (settings.getSkipStage4()) {
                writeGraphExplicit(3);
        }
}

void Brownie::stageFour()
{
        // ============================================================
        // STAGE 4 : GRAPH SIMPLIFICATION
        // ============================================================

        cout << "Entering stage 4" << endl;
        cout << "================" << endl;
//#ifndef DEBUG
        if (!stageFourNecessary()) {
                cout << "Files produced by this stage appear to be present, "
                        "skipping stage 4..." << endl << endl;
                return;
        }
//#endif
        Util::startChrono();
        cout << "Creating graph... ";
        cout.flush();
        double readLength=libraries.getReadLength();
        if (readLength<=settings.getK() ||readLength>500)
                readLength=250;
        settings.setReadLenght(readLength);
        DBGraph testgraph(settings);
        testgraph.loadGraphBin(getBinNodeFilename(3),
                               getBinArcFilename(3),
                               getMetaDataFilename(3));
        cout << "done (" << testgraph.getNumNodes() << " nodes, "
             << testgraph.getNumArcs() << " arcs), ("
             << Util::stopChronoStr() << ")" << endl;
        cout << "Initial kmerCoverage : "<<testgraph.estimatedKmerCoverage<<endl;
        cout<<endl<<"Maximum node size to delete is::     "<<testgraph.maxNodeSizeToDel<<endl;
        string command = "mkdir " + settings.getTempDirectory() + "cov";
        system(command.c_str());
        command="rm "+settings.getTempDirectory() + "cov/* && rm "+settings.getTempDirectory() + "cov/*.dat";
        cout<<command<<endl;
        system(command.c_str());
        testgraph.clipTips(0);
        testgraph.mergeSingleNodes(true);
        testgraph.filterCoverage(testgraph.cutOffvalue);
        testgraph.mergeSingleNodes(true);
        testgraph.extractStatistic(0);
        DBGraph graph(settings);
        graph.loadGraphBin(getBinNodeFilename(3),
                               getBinArcFilename(3),
                               getMetaDataFilename(3));
        graph.estimatedKmerCoverage=testgraph.estimatedKmerCoverage;
        graph.estimatedMKmerCoverageSTD=testgraph.estimatedMKmerCoverageSTD;
        testgraph.clear();
        Util::startChrono();
        cout << "Creating graph... ";
        cout.flush();
        cout << "done (" << graph.getNumNodes() << " nodes, "
        << graph.getNumArcs() << " arcs)" << endl;
        cout << "Created graph in "
        << Util::stopChrono() << "s." << endl;
        Util::startChrono();
        #ifdef DEBUG
        cout<<"estimated Kmer Coverage Mean: "<<graph.estimatedKmerCoverage<<endl;
        cout<<"estimated Kmer Coverage STD: "<<graph.estimatedMKmerCoverageSTD<<endl;
        graph.compareToSolution(getTrueMultFilename(3), true);
        #endif
        //variables
        int round=1;
        graph.updateGraphSize();
        double coverageCutOff=1;
        do{
                graph.filterCoverage(coverageCutOff);
                bool simplified = true;
                while (simplified ) {// &&
                        //
                        //*******************************************************
                        graph.updateCutOffValue(round);
                        bool tips=graph.clipTips(round);
                        if(tips) {
                                graph.mergeSingleNodes(true);
                                #ifdef DEBUG
                                graph.compareToSolution(getTrueMultFilename(3), false);
                                #endif
                        }
                        //*******************************************************
                        bool bubble=false;
                        size_t depth=settings.getK()+1;
                        #ifdef DEBUG
                        graph.updateCutOffValue(round);
                        graph.compareToSolution(getTrueMultFilename(3), false);
                        #endif
                        bubble= graph.bubbleDetection(depth);
                        graph.mergeSingleNodes(true);
                        bool continuEdit=false;
                        size_t maxDepth=(round)*100>1000?1000:(round)*100;
                        while(depth<maxDepth){
                                depth=depth+150;
                                cout<<"bubble depth:"<<depth <<endl;
                                continuEdit= graph.bubbleDetection(depth);
                                if (continuEdit)
                                        graph.mergeSingleNodes(true);

                                bubble=false?continuEdit:bubble;
                        }
                        #ifdef DEBUG
                        graph.compareToSolution(getTrueMultFilename(3),false);
                        graph.updateCutOffValue(round);
                        #endif
                        graph.extractStatistic(round);
                        bool deleted=graph.deleteUnreliableNodes( round);
                        continuEdit=deleted;
                        while(continuEdit){
                                continuEdit=graph.deleteUnreliableNodes( round);
                                graph.mergeSingleNodes(false);
                                graph.extractStatistic(round);
                        }
                        while(graph.mergeSingleNodes(true));
                        #ifdef DEBUG
                        graph.compareToSolution(getTrueMultFilename(3),false);
                        graph.updateCutOffValue(round);
                        cout<<"estimated Kmer Coverage Mean: "<<graph.estimatedKmerCoverage<<endl;
                        cout<<"estimated Kmer Coverage STD: "<<graph.estimatedMKmerCoverageSTD<<endl;
                        #endif
                        simplified = tips    || deleted || bubble;
                        graph.updateGraphSize();
                        round++;
                }
                coverageCutOff++;
                if (coverageCutOff>graph.certainVlueCov)
                        break;
                graph.updateGraphSize();
                cout<<"coverage cut off ::"<<coverageCutOff<<endl;
                cout<<"certainVlueCov   ::"<<graph.certainVlueCov<<endl;;
        }
        while(graph.sizeOfGraph>settings.getGenomeSize());
        #ifdef DEBUG
        graph.compareToSolution(getTrueMultFilename(3), false);
        #endif
        cout<<"graph size:"<<graph.sizeOfGraph<<endl;
        graph.writeGraph(getNodeFilename(4),getArcFilename(4),getMetaDataFilename(4));
        cout << " Ghraph correction completed in "
        << Util::stopChrono() << "s." << endl;
        Util::startChrono();

#ifdef DEBUG
        graph.writeCytoscapeGraph(0);
        graph.sanityCheck();
        command="pdftk "+settings.getTempDirectory()+ "cov/*.pdf cat output allpdfFiles.pdf";
        system(command.c_str());
#endif

        graph.clear();
        cout << "Stage 4 finished.\n" << endl;
}

void Brownie::stageFive()
{
        if (!stageFiveNecessary()) {
                         cout << "Files produced by this stage appear to be present, "
                         "skipping stage 5..." << endl << endl;
                         return;
        }
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
        graph.compareToSolution(getTrueMultFilename(4),true);
        graph.updateGraphSize();
        graph.writeCytoscapeGraph(0);
        #endif
        Util::startChrono();
        ReadCorrection rc(graph, settings);
        rc.errorCorrection(libraries);
        cout << "Error correction completed in "
        << Util::stopChrono() << "s." << endl;
        cout << "Stage 5 finished.\n" << endl;
        graph.clear();
}

void Brownie::writeGraphExplicit(int stage)
{
        // build pre-graph and simplify it
        DBGraph graph(settings);
        graph.createFromFile(getNodeFilename(stage),
                             getArcFilename(stage),
                             getMetaDataFilename(stage));
        graph.writeGraphExplicit();
        graph.clear();
}


int main(int argc, char** args)
{
        try {
                Brownie brownie(argc, args);
#ifndef DEBUG
                cout<<"running in Release mode"<<endl;
               // brownie.printInFile();
#endif
#ifdef DEBUG
                cout<<"In DEBUG mode"<<endl;
#endif

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

