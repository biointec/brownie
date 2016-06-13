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

using namespace std;

// ============================================================================
// PRIVATE COVERAGE.CPP (STAGE 3)
// ============================================================================

void DBGraph::parseReads(size_t thisThread,
                         vector<string>& readBuffer)
{
        for (size_t i = 0; i < readBuffer.size(); i++) {
                const string& read = readBuffer[i];

                KmerIt it(read);
                if (!it.isValid())
                        continue;

                // increase the read start coverage (only for the first valid kmer)
                NodePosPair result = findNPP(it.getKmer());
                if (result.getNodeID() != 0) {
                        SSNode node = getSSNode(result.getNodeID());
                        node.setReadStartCov(node.getReadStartCov()+1);
                }

                NodeID prevID = 0;
                for (KmerIt it(read); it.isValid(); it++ ) {
                        Kmer kmer = it.getKmer();
                        NodePosPair result = findNPP(kmer);
                        if (!result.isValid()) {
                                prevID = 0;
                                continue;
                        }

                        // we've found the kmer, increase the node coverage
                        NodeID thisID = result.getNodeID();
                        SSNode node = getSSNode(thisID);
                        node.incKmerCov();

                        // if the previous node was valid and different, increase the arc coverage
                        if ((prevID != 0) && (prevID != thisID)) {
                                getSSNode(prevID).getRightArc(thisID)->incReadCov();
                                getSSNode(thisID).getLeftArc(prevID)->incReadCov();
                        }

                        if (it.hasRightOverlap())
                                prevID = thisID;
                        else
                                prevID = 0;
                }
        }
}

void DBGraph::workerThread(size_t thisThread, LibraryContainer* inputs)
{
        // local storage of reads
        vector<string> myReadBuf;

        size_t blockID, recordOffset;
        while (inputs->getReadChunk(myReadBuf, blockID, recordOffset))
                parseReads(thisThread, myReadBuf);
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

void DBGraph::countNodeandArcFrequency(LibraryContainer &inputs)
{
        const unsigned int& numThreads = settings.getNumThreads();

        cout << "Building kmer-node table... "; cout.flush();
        buildKmerNPPTable();
        cout << "done" << endl;

        cout << "Number of threads: " << numThreads << endl;
        cout << "Counting node and arc frequency: " << endl;

        inputs.startIOThreads(settings.getThreadWorkSize(),
                              settings.getThreadWorkSize() * settings.getNumThreads());

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&DBGraph::workerThread, this, i, &inputs);

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        inputs.joinIOThreads();

        destroyKmerNPPTable();
}

void DBGraph::fitKmerSpectrum(const string& tempdir)
{
        double avgKmerCov = getInitialKmerCovEstimate(2.0, 0.01);
        cout << "Initial coverage estimate: " << avgKmerCov << endl;

        size_t numComponents = 3;
        double xMax = numComponents * avgKmerCov;

        map<unsigned int, double> data;
        for ( NodeID id = 1; id <= numNodes; id++ ) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                if (node.getAvgKmerCov() > xMax)
                        continue;

                data[round(node.getAvgKmerCov())] += node.getMarginalLength();
        }

        // also add the number of unique k-mers to the model
        double temp;
        ifstream ifs((tempdir + "numkmers.stage1").c_str());
        ifs >> temp;
        ifs.close();
        data[1] = temp;

        vector<double> mu(numComponents), var(numComponents), MC(numComponents);

        for (size_t i = 0; i < numComponents; i++) {
                mu[i] = i * avgKmerCov + 2.0;
                var[i] = 1.01 * mu[i];
                MC[i] = 1.0;
        }

        Util::binomialMixtureEM(data, mu, var, MC, 50);

        double yMax = 1.2 * MC[1] * Util::negbinomialPDF(mu[1], mu[1], var[1]);

        for (int j = 0; j < numComponents; j++) {
                cout << "Average of component " << j << ": " << mu[j] << endl;
                cout << "Mixing coefficient of component " << j << ": " << MC[j] << endl;
                cout << "Variance of component " << j << ": " << var[j] << endl;
        }

        // estimate the genome size
        size_t genomeSize = 0;
        for (int j = 1; j <= numComponents; j++)
                genomeSize += j * MC[j];
        cout << "Estimated genome size: " << genomeSize << endl;

        ofstream ofs((tempdir + "spectrum.txt").c_str());
        for (const auto& it : data) {
                double fit = MC[0] * Util::geometricPDF(it.first, mu[0]);
                for (size_t i = 1; i < numComponents; i++)
                        fit += MC[i] * Util::negbinomialPDF(it.first, mu[i], var[i]);
                ofs << it.first << "\t" << it.second << "\t" << fit << endl;
        }
        ofs.close();

        ofs.open((tempdir + "spectrum.gnu").c_str());
        ofs << "set output \"spectrum.ps\"" << endl;
        ofs << "set terminal postscript landscape" << endl;
        ofs << "set xrange [0:" << (size_t)xMax << "]" << endl;
        ofs << "set yrange [0:" << (size_t)yMax << "]" << endl;
        //ofs << "set style line 1 lt 2 lc rgb \"black\" lw 3" << endl;
        //ofs << "set style line 2 lt 2 lc rgb \"red\" lw 3" << endl;
        ofs << "plot \"spectrum.txt\" using 1:2 title \'k-mer spectrum\' with lines, ";
        ofs << "\"spectrum.txt\" using 1:3 title \'fitted model\' with lines lt 3 lw 3 lc 3" << endl;
        ofs.close();

        // figure out how many of the nodes are correctly flagged as erroneous
       /* double fp = 0, tp = 0, tn = 0, fn = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                bool actual = (trueMult[id] >= 1);
                bool predicted = (node.getExpMult(avgKmerCov) >= 1);

                //cout << actual << " " << predicted << endl;
                if (actual && predicted)
                        tp++;
                if (!actual && predicted)
                        fp++;
                if (actual && !predicted)
                        fn++;
                if (!actual && !predicted)
                        tn++;
        }

        double prec = tp / (tp + fp);
        double rec = tp / (tp + fn);
        double F1 = 2 * prec * rec / (prec + rec);

        cout << "Precision: " << prec << endl;
        cout << "Recall: " << rec << endl;
        cout << "F1 score: " << F1 << endl;*/
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
