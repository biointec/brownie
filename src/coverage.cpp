/***************************************************************************
 *   Copyright (C) 2015 - 2016 Jan Fostier (jan.fostier@intec.ugent.be)    *
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

#include "graph.h"
#include "settings.h"
#include "library.h"
#include <cmath>

using namespace std;

// ============================================================================
// PRIVATE COVERAGE.CPP (STAGE 3)
// ============================================================================

void DBGraph::parseReads(size_t thisThread, const vector<string>& readBuffer,
                         KmerCountTable& table)
{
        for (size_t i = 0; i < readBuffer.size(); i++) {
                const string& read = readBuffer[i];

                KmerIt it(read);
                if (!it.isValid())
                        continue;

                // process the first kmer
                NodeID nodeID = table.incrKmerCount(it.getKmer()); // atomic
                if (nodeID != 0) {
                        SSNode node = getSSNode(nodeID);
                        node.incReadStartCov(); // atomic
                        node.incKmerCov();      // atomic
                }

                // process the remainder of the kmers
                NodeID prevID = it.hasRightOverlap() ? nodeID : 0;
                for (it++; it.isValid(); it++ ) {
                        nodeID = table.incrKmerCount(it.getKmer());    // atomic
                        if (nodeID == 0) {
                                prevID = 0;
                                continue;
                        }

                        // we've found the kmer, increase the node coverage
                        SSNode node = getSSNode(nodeID);
                        node.incKmerCov();      // atomic

                        // if the previous node was valid and different, increase the arc coverage
                        if ((prevID != 0) && (prevID != nodeID)) {
                                getSSNode(prevID).getRightArc(nodeID)->incReadCov();
                                getSSNode(nodeID).getLeftArc(prevID)->incReadCov();
                        }

                        prevID = it.hasRightOverlap() ? nodeID : 0;
                }
        }
}

void DBGraph::kmerCountWorkerThread(size_t thisThread, LibraryContainer* inputs,
                                    KmerCountTable* table)
{
        // local storage of reads
        vector<string> myReadBuf;

        size_t blockID, recordOffset;
        while (inputs->getReadChunk(myReadBuf, blockID, recordOffset))
                parseReads(thisThread, myReadBuf, *table);
}

double DBGraph::getInitialKmerCovEstimate(double errLambda, double p) const
{
        // sanity checks
        assert(errLambda > 0.0);
        assert(p > 0.0);
        assert(p < 1.0);

        // Given a Poisson distribution for the error model, find a cutoff
        // value for the coverage for which the probability of observing
        // a coverage is less than p under this error model
        double cutoff = ceil(errLambda);
        for ( ; cutoff < 10.0 * errLambda; cutoff++)
                if (Util::poissonPDF((unsigned int)cutoff, errLambda) < p)
                        break;

        size_t totCoverage = 0, totSize = 0;
        for ( NodeID id = 1; id <= numNodes; id++ ) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                if (node.getAvgKmerCov() < cutoff)
                        continue;
                totCoverage += node.getKmerCov();
                totSize += node.getMarginalLength();
        }

        if (totSize > 0)
                return (double)totCoverage / (double)totSize;
        else
                return 2.0 * errLambda;      // pathological case
}

// ============================================================================
// PUBLIC COVERAGE.CPP (STAGE 3)
// ============================================================================

void DBGraph::generateKmerSpectrum(const string& tempdir, LibraryContainer &inputs)
{
        KmerCountTable kmerCountTable;
        const unsigned int& numThreads = settings.getNumThreads();

        cout << "Building kmer-count table... "; cout.flush();
        buildKmerCountTable(kmerCountTable);
        cout << "done (size: " << kmerCountTable.size() << " kmers)" << endl;

        cout << "Number of threads: " << numThreads << endl;
        cout << "Generating k-mer spectrum: " << endl;

        inputs.startIOThreads(settings.getThreadWorkSize(),
                              settings.getThreadWorkSize() * settings.getNumThreads());

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&DBGraph::kmerCountWorkerThread,
                                          this, i, &inputs, &kmerCountTable);

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
        inputs.joinIOThreads();

        // build a kmer-spectrum from the table
        kmerCountTable.getKmerSpectrum(kmerSpectrum);
        kmerCountTable.clear();

        double avgKmerCov = getInitialKmerCovEstimate(2.0, 0.01);
        cout << "Initial coverage estimate: " << avgKmerCov << endl;
        kmerSpectrum.fitKmerSpectrum(avgKmerCov);

        kmerSpectrum.writeSpectrum(tempdir + "spectrum.txt");
        kmerSpectrum.writeGNUPlotFile(tempdir + "spectrum.gnu");
        kmerSpectrum.writeSpectrumFit(tempdir + "spectrum.fit");

        cout << kmerSpectrum << endl;
}

int DBGraph::getExpMult(double obsKmerCov) const
{
        if (obsKmerCov < getCovCutoff())
                return 0;
        int expMult = round(obsKmerCov / getAvgKmerCov());
        return (expMult == 0) ? 1 : expMult;
}

double DBGraph::getObsProb(unsigned int obsCov, unsigned int mult) const
{
        return kmerSpectrum.evalSpec(obsCov, mult);
}

double DBGraph::getObsProbLog(double obsCov, int ML, unsigned int mult) const
{
        return kmerSpectrum.evalSpecLog(obsCov, mult);
}

void DBGraph::writeNodeFile(const std::string& filename) const
{
        ofstream ofs(filename.c_str());
        for ( NodeID id = 1; id <= numNodes; id++ ) {
                        SSNode node = getSSNode(id);
                        if (!node.isValid())
                                continue;
                        ofs << id << "\t" << node.getMarginalLength() << "\t" << node.getAvgKmerCov() << "\t" << node.getReadStartCov() << endl;
        }
        ofs.close();
}
