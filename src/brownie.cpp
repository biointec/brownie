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
void Brownie::printInFile() {
        string konsole=settings.getTempDirectory()+ "/konsole.txt";
        freopen(konsole.c_str(),"a",stdout);
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
        readParser->writeAllKmers(getKmerFilename());
        //readParser->writeKmersWithCovGTOne(getKmerFilename());
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

        graph.writeGraphExplicit(3);

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

       /*if (!stageFourNecessary()) {
                         cout << "Files produced by this stage appear to be present, "
                         "skipping stage 4..." << endl << endl;
                         return;
        }*/
        DBGraph testgraph(settings);
        testgraph.createFromFile(getNodeFilename(3),
                                 getArcFilename(3),
                                 getMetaDataFilename(3));
        cout<<"initial kmerCoverage : "<<testgraph.estimatedKmerCoverage<<endl;
        string command="mkdir "+settings.getTempDirectory()+"cov";
        system(command.c_str());
        command="rm "+settings.getTempDirectory() + "cov/* && rm "+settings.getTempDirectory() + "cov/*.dat";
        cout<<command<<endl;
        system(command.c_str());
        testgraph.clipTips(0);
        testgraph.mergeSingleNodes(true);
        testgraph.filterCoverage(0);
        testgraph.mergeSingleNodes(true);
        testgraph.extractStatistic(0);
        DBGraph graph(settings);
        graph.estimatedKmerCoverage=testgraph.estimatedKmerCoverage;
        graph.estimatedMKmerCoverageSTD=testgraph.estimatedMKmerCoverageSTD;
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
        //#ifdef DEBUG
        graph.compareToSolution();
        //#endif
        //variables


        bool simplified = true;
        size_t bigestN50=0;
        int round=1;
        size_t minN50=100;

        while (simplified ) {
                //*******************************************************
                graph.updateCutOffValue(round);
                bool tips=graph.clipTips(round);
                if (tips) {
                        graph.mergeSingleNodes(false);
                        graph.compareToSolution();
                }
                if (graph.sizeOfGraph<settings.getGenomeSize()){
                         while(graph.mergeSingleNodes(true));
                         break;
                }
                //*******************************************************
                graph.updateCutOffValue(round);
                bool bubble=false;
                bubble= graph.bubbleDetection(round);
                if (bubble) {
                        graph.mergeSingleNodes(false);
                        graph.compareToSolution();

                }

                graph.extractStatistic(round);
                bool deleted=graph.deleteUnreliableNodes( round);
                while(graph.mergeSingleNodes(true));
                graph.compareToSolution();
                graph.updateCutOffValue(round);
                cout<<"estimated Kmer Coverage Mean: "<<graph.estimatedKmerCoverage<<endl;
                cout<<"estimated Kmer Coverage STD: "<<graph.estimatedMKmerCoverageSTD<<endl;
                simplified = tips    || deleted || bubble; //link |||| coverage
                round++;

        }
        graph.writeCytoscapeGraph(0);
        graph.clipTips(0);
        graph.deleteUnreliableNodes(0);
        graph.mergeSingleNodes(true);
        graph.bubbleDetection(0);
        graph.mergeSingleNodes(true);
        graph.writeGraph( nodeFileName,arcFileName,metaDataFileName);
        cout << " Ghraph correction completed in "
        << Util::stopChrono() << "s." << endl;
        Util::startChrono();
        #ifdef DEBUG
        graph.sanityCheck();
        command="pdftk "+settings.getTempDirectory()+ "cov/*.pdf cat output allpdfFiles.pdf";
        system(command.c_str());
        #endif

        graph.writeGraphExplicit(4);

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
        graph.compareToSolution();
        graph.writeCytoscapeGraph(0);
        #endif

        cout<<"N50 size is: " <<graph.getN50()<<endl;
        Util::startChrono();
        ReadCorrection rc(graph, settings);
        rc.errorCorrection(libraries);
        cout << "Error correction completed in "
        << Util::stopChrono() << "s." << endl;
        cout << "Stage 5 finished.\n" << endl;
        graph.clear();
}

int main(int argc, char** args)
{
        try {
                Brownie brownie(argc, args);
                bool debug=false;
                #ifdef DEBUG
                debug=true;
                #endif
                //if(!debug)
                     //  brownie.printInFile();
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

