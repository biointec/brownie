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
#include "kmeroverlaptable.h"
#include "kmertable.h"
#include "readcorrection.h"

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
        cout << "-------------------------------------------------------------------------------" << endl;
        time_t t;
        time(&t);
        cout << ctime(&t) << endl; // e.g., Fri May 02 17:57:14 2003
        cout << "-------------------------------------------------------------------------------" << endl;
        cout << "Welcome to Brownie\n" << endl;
}

void Brownie::parameterEstimationInStage4(DBGraph &graph){
        cout << endl << " ================ Parameter Estimation ===============" << endl;
        double  estimatedKmerCoverage = 0,
                estimatedMKmerCoverageSTD = 0,
                cutOffvalue = 0,
                readLength = 0;
        readLength = libraries.getAvgReadLength();
        if (readLength <= settings.getK() || readLength > 500)
                readLength = 150;
        cout << "Loading test graph for initial parameter estimation" << endl;
        DBGraph testgraph(settings);
        testgraph.loadGraphBin(getBinNodeFilename(3),
                               getBinArcFilename(3),
                               getMetaDataFilename(3));
        testgraph.readLength = readLength;
        #ifdef DEBUG
        testgraph.compareToSolution(getTrueMultFilename(3), true);
        #endif
        testgraph.clipTips(0);
        testgraph.mergeSingleNodes();
        testgraph.filterCoverage(testgraph.cutOffvalue);
        testgraph.mergeSingleNodes();
        testgraph.extractStatistic(0);
        cout << "Estimated Kmer coverage mean: " << testgraph.estimatedKmerCoverage << endl;
        cout << "Estimated Kmer coverage std:  " << testgraph.estimatedMKmerCoverageSTD << endl;

        estimatedKmerCoverage = testgraph.estimatedKmerCoverage;
        estimatedMKmerCoverageSTD = testgraph.estimatedMKmerCoverageSTD;
        double estimatedErroneousKmerCoverage = 1+ estimatedKmerCoverage/100;
        double e = 2.718281;
        double c = estimatedErroneousKmerCoverage/estimatedKmerCoverage;

        if (settings.getCutOffValue() > 0)
                cutOffvalue = settings.getCutOffValue();
        else
                cutOffvalue = (estimatedErroneousKmerCoverage-estimatedKmerCoverage)* (log(e)/log(c));
        testgraph.updateGraphSize();
        testgraph.clear();
           //initialize values for graph parameter based on test graph.
        graph.estimatedKmerCoverage = estimatedKmerCoverage;
        graph.estimatedMKmerCoverageSTD = estimatedMKmerCoverageSTD;
        graph.cutOffvalue = cutOffvalue;
        graph.readLength = readLength;
        graph.maxNodeSizeToDel = readLength*4;
        graph.redLineValueCov = cutOffvalue;
        graph.certainVlueCov = cutOffvalue*.3;
        graph.safeValueCov = cutOffvalue*.7;

        cout << "Estimated Kmer coverage:        " << graph.estimatedKmerCoverage << endl;
        cout << "Estimated MKmer coverage std:   " << graph.estimatedMKmerCoverageSTD << endl;
        cout << "Cut-off value:                  " << graph.cutOffvalue << endl;
        cout << "Read length:                    " << graph.readLength << endl;
        cout << "Max node size to delete:        " << graph.maxNodeSizeToDel << endl;
        cout << "Red line value cov:             " << graph.redLineValueCov << endl;
        cout << "Certain value cov:              " << graph.certainVlueCov << endl;
        cout << "Safe value cov:                 " << graph.safeValueCov << endl;
        cout << "End of parameter estimation ... " << endl;
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
        parameterEstimationInStage4( graph );
        Util::startChrono();
        cout << "Creating graph... ";
        graph.loadGraphBin(getBinNodeFilename(3),
                           getBinArcFilename(3),
                           getMetaDataFilename(3));

        cout.flush();
        cout << "done (" << graph.getNumNodes() << " nodes, "
             << graph.getNumArcs() << " arcs)" << endl;
        cout << "Created graph in "
             << Util::stopChrono() << "s." << endl;

#ifdef DEBUG
        graph.compareToSolution(getTrueMultFilename(3), true);
#endif
        Util::startChrono();
        graph.graphPurification(getTrueMultFilename(3), libraries);
#ifdef DEBUG
        graph.compareToSolution(getTrueMultFilename(3), false);
#endif
        cout << "Graph size: " << graph.sizeOfGraph << " bp" << endl;
        graph.writeGraph(getNodeFilename(4),getArcFilename(4),getMetaDataFilename(4));
        cout<<"N50 is: "<<graph.n50<<endl;
        cout << "Graph correction completed in "
             << Util::stopChrono() << "s." << endl;
        Util::startChrono();

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
        graph.createFromFile(getNodeFilename(4),
                             getArcFilename(4),
                             getMetaDataFilename(4));
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        cout << "Graph contains " << graph.getNumNodes() << " nodes and "
             << graph.getNumArcs() << " arcs" << endl;

#ifdef DEBUG
        graph.compareToSolution(getTrueMultFilename(4),true);
        graph.updateGraphSize();

#endif

        Util::startChrono();
        ReadCorrectionHandler rcHandler(graph, settings);
        rcHandler.doErrorCorrection(libraries);

        cout << "Error correction completed in " << Util::stopChronoStr() << endl;
        cout << "Stage 5 finished\n" << endl;
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
                graph.createFromFile(getNodeFilename(4),
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
                //brownie.printInFile();
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
                brownie.writeGraphFasta();
        } catch (exception &e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }
        cout << "Exiting... bye!" << endl << endl;
        return EXIT_SUCCESS;
}

