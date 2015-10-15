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

#include "graph.h"
#include "kmernode.h"
#include "readfile/fastafile.h"
#include "settings.h"
#include "mappedread.h"
#include <cmath>
#include <fstream>

using namespace std;

// ============================================================================
// READ MAPPING PRIVATE
// ============================================================================

void DBGraph::detectPairedReadDirection(ReadLibrary& input,
                                        const KmerNodeTable& table)
{
        TString read1, read2;

        const double threshold = 0.05;
        const size_t minOcc = 100;

        ifstream ifs(input.getPairedFilename().c_str(),
                     std::ios::in | std::ios::binary);

        vector<size_t> dir(4, 0);

        while (ifs.good()) {
                read1.read(ifs);
                read2.read(ifs);

                if (read1.getLength() == 0)
                        continue;
                if (read2.getLength() == 0)
                        continue;

                MappedRead mRead1(read1, table, *this);
                MappedRead mRead2(read2, table, *this);

                // if one the reads can not be mapped, get out
                if ((!mRead1.isMapped()) || (!mRead2.isMapped()))
                        continue;

                // if one of the reads maps to multiple nodes, get out
                if ((!mRead1.mapsToSingle()) || (!mRead2.mapsToSingle()))
                        continue;

                // if the reads to not map to the same double stranded node, get out
                if (mRead1.getSingleNodeID() != abs(mRead2.getSingleNodeID()))
                        continue;

                bool sameDir = mRead1.getSingleNodeID() == mRead2.getSingleNodeID();
                if (sameDir) {
                        if (mRead1.getInsertSize(mRead2) > 0)
                                dir[0]++; // R-R case
                        else
                                dir[1]++; // L-L case
                } else {
                        read2.reverseComplement();
                        MappedRead mRead2a(read2, table, *this);

                        if (mRead1.getFirstNodeReadStartPos() <
                            mRead2a.getFirstNodeReadStartPos())
                                dir[2]++; // R-L case
                        else
                                dir[3]++; // L-R case
                }

                // find the maximum and the second best in dir
                size_t max = 0, smax = 0, maxi = 0;
                for (size_t i = 0; i < 4; i++) {
                        if (dir[i] > max) {
                                smax = max;
                                max = dir[i];
                                maxi = i;
                        } else if (dir[i] > smax) {
                                smax = dir[i];
                        }
                }

                // check for exit criterion
                if ((max >= minOcc) && (((double)smax / max) <= threshold)) {
                        switch (maxi) {
                                case 0:
                                        input.setReadDirType(RR);
                                        break;
                                case 1:
                                        input.setReadDirType(LL);
                                        break;
                                case 2:
                                        input.setReadDirType(RL);
                                        break;
                                case 3:
                                        input.setReadDirType(LR);
                                        break;
                                default:
                                        break;
                        }
                        break;
                }
        }

        ifs.close();
}

void DBGraph::mapIsolatedReadsToGraph(ReadLibrary &input,
                                      const KmerNodeTable& table)
{
        size_t numReads = input.getNumIsolatedReads();
        input.setNumMappedIsolated(0);

        ofstream ofs(input.getMappedIsolatedFilename().c_str(),
                     std::ios::out | std::ios::binary);

        if (numReads == 0) {
                ofs.close();
                return;
        }

        ifstream ifs(input.getIsolatedFilename().c_str(),
                     std::ios::in | std::ios::binary);

        TString read;
        size_t readCount = 0, notMapped = 0, mapped = 0;

        while (ifs.good()) {
                read.read(ifs);
                if (read.getLength() == 0)
                        continue;

                if (readCount++ % OUTPUT_FREQUENCY == 0) {
                        cout << "\tProcessing read " << readCount
                             << "/" << numReads << "\r";
                        cout.flush();
                }

                MappedRead mRead(read, table, *this);

                if (!mRead.isMapped()) {
                        notMapped++;
                        continue;
                }

                mRead.write(ofs);
                mapped++;
        }

        cout << "\tProcessing read " << numReads  << "/" << numReads << endl;
        cout << "\tNumber of reads that can be mapped to graph: "
             << mapped << " (" << 100.0*mapped/(double)numReads << "%)" << endl;
        cout << "\tNumber of pairs that can not be mapped to graph: "
             << notMapped << " (" << 100.0*notMapped/(double)numReads << "%)" << endl;

        input.setNumIsolatedReads(mapped);
        ofs.close();
}

void DBGraph::mapPairedReadsToGraph(ReadLibrary &input,
                                    const KmerNodeTable& table)
{
        size_t numPairs = input.getNumPairedReads() / 2;
        input.setNumMappedPaired(0);

        if (numPairs == 0)
                return;

        ifstream ifs(input.getPairedFilename().c_str(),
                     std::ios::in | std::ios::binary);
        ofstream ofsp(input.getMappedPairedFilename().c_str(),
                      std::ios::out | std::ios::binary);
        ofstream ofsi(input.getMappedIsolatedFilename().c_str(),
                      std::ios::out | std::ios::binary | std::ios::app);

        TString read1, read2;
        size_t pairCount = 0, mapToSame = 0, notMapped = 0,
               mappedIsolated = 0, mappedPaired = 0;


        map<size_t, size_t> distribution;
        Observations ilEst;

        while (ifs.good()) {
                read1.read(ifs);
                read2.read(ifs);

                if (read1.getLength() == 0)
                        continue;
                if (read2.getLength() == 0)
                        continue;

                pairCount++;

                if (pairCount % OUTPUT_FREQUENCY == 0) {
                        cout << "\tProcessing pair " << pairCount
                             << "/" << numPairs << "\r";
                        cout.flush();
                }

                if (input.needToRCRead1())
                        read1.reverseComplement();
                if (input.needToRCRead2())
                        read2.reverseComplement();

                MappedRead mRead1(read1, table, *this);
                MappedRead mRead2(read2, table, *this);

                // if both reads can not be mapped, get out
                if ((!mRead1.isMapped()) && (!mRead2.isMapped())) {
                        notMapped++;
                        continue;
                }

                if (!mRead1.isMapped()) {
                        mRead2.write(ofsi);
                        mappedIsolated++;
                        continue;
                }

                if (!mRead2.isMapped()) {
                        mRead1.write(ofsi);
                        mappedIsolated++;
                        continue;
                }

                // if the reads map to connected nodes, get out
                /*if (mRead1.mapsToConnected(mRead2)) {
                        mapToConnected++;
                        continue;
                }*/

               /* if (mRead1.mapsToSingle() && mRead2.mapsToSingle())
                        if (mRead1.getSingleNodeID() == 50 && mRead2.getSingleNodeID() == 46){
                                if (mRead1.getFirstNodeReadStartPos() < getSSNode(mRead1.getSingleNodeID()).getMarginalLength() - 3800 + 3.0 * 300 +  62 )
                                        continue;
                                double il = getSSNode(50).getMarginalLength() + 62 + mRead2.getFirstNodeReadStartPos() - mRead1.getFirstNodeReadStartPos() + mRead2.getReadLength();
                                const_cast<Input&>(input).updateInsertSize(il);
                                mapToSame++;

                            //    cout << read1 << endl;
                            //    cout << read2 << endl;
                                cout << il << endl;
                            //    exit(0);
                                continue;
                        }*/

             /*   size_t k = Kmer::getWordSize();
                if (mRead1.mapsToSameNode(mRead2) && getSSNode(mRead1.getSingleNodeID()).getMarginalLength() > 1500) {
                        if (mRead1.getFirstNodeReadStartPos() + k < 1000)
                                for (size_t i = 0; i < 10; i++) {
                                        size_t pos = mRead2.getFirstNodeReadStartPos() + i;
                                        if (pos >= 1000)
                                                if (distribution.find(pos) == distribution.end())
                                                        distribution[pos] = 1;
                                                else
                                                        distribution[pos]++;
                                }
                }*/

                // if both kmers map to the same node
                if (mRead1.mapsToSameNode(mRead2)) {
                        size_t il = mRead1.getInsertSize(mRead2);

                        // make an unbiased estimation of the insert length
                        if (ilEst.getNumObs() > 50) {
                                SSNode node = getSSNode(mRead1.getSingleNodeID());
                                if (mRead1.getFirstNodeReadStartPos() + ilEst.getAverage() + 5.0 * ilEst.getStdev() > node.getMarginalLength())
                                        continue;
                        }

                        ilEst.addObsOnline(il);
                        mapToSame++;
                        if (distribution.find(il) == distribution.end())
                                distribution[il] = 1;
                        else
                                distribution[il]++;

                        continue;
                }

                mRead1.write(ofsp);
                mRead2.write(ofsp);

                read1.write(ofsp);
                read2.write(ofsp);

                mappedPaired++;
        }

        cout << "\tProcessing pair " << numPairs  << "/" << numPairs << endl;

        cout << "Estimator average : " << ilEst.getAverage() << endl;
        cout << "Estimator stddev : " << ilEst.getStdev() << endl;

        size_t sum = 0;
        for (map<size_t, size_t>::iterator it = distribution.begin(); it != distribution.end(); it++) {
                sum += it->first * it->second;
        }

        double Q1, Q2;
        size_t currSum = 0;

        // Estimate the left sigma point
        double target = 0.158655254;
        typedef map<size_t, size_t>::iterator DistIt;
        for (DistIt it = distribution.begin(); it != distribution.end(); it++) {
                currSum += it->first * it->second;
                if (currSum > target*sum) {
                        if (it == distribution.begin()) {
                                Q1 = it->first;
                                break;
                        }

                        double y2 = currSum;
                        double y1 = currSum - it->first * it->second;

                        double p = (target*sum - y2) / (y1 - y2);

                        double x2 = it->first;
                        double x1 = (--it)->first;

                        Q1 = p*x1 + (1.0-p)*x2;
                        break;
                }
        }

        // Estimate the right sigma point
        target = 1.0-0.158655254;
        currSum = 0;
        for (DistIt it = distribution.begin(); it != distribution.end(); it++) {
                currSum += it->first * it->second;
                if (currSum > target*sum) {
                        if (it == distribution.begin()) {
                                Q2 = it->first;
                                break;
                        }

                        double y2 = currSum;
                        double y1 = currSum - it->first * it->second;

                        double p = (target*sum - y2) / (y1 - y2);

                        double x2 = it->first;
                        double x1 = (--it)->first;

                        Q2 = p*x1 + (1.0-p)*x2;
                        break;
                }
        }

        cout << "CENTRAL AVERAGE: " << 0.5*(Q2+Q1) << endl;
        cout << "CENTRAL SIGMA : " << 0.5*(Q2-Q1) << endl;

        const_cast<ReadLibrary&>(input).setInsertSize(ilEst.getAverage());
        const_cast<ReadLibrary&>(input).setStddev(ilEst.getStdev());
        const_cast<ReadLibrary&>(input).setCentralInsertSize(0.5*(Q2+Q1));
        const_cast<ReadLibrary&>(input).setCentralStddev(0.5*(Q2-Q1));

        /*ofstream ofs("cliff.txt", std::ios::out);
        for (map<size_t, size_t>::iterator it = distribution.begin(); it != distribution.end(); it++)
                ofs << it->first << "\t" << it->second << endl;
        ofs.close();*/

        cout.precision(4);

        cout << "\tNumber of pairs that map to the same node: " << mapToSame
             << " (" << 100.0*mapToSame/(double)numPairs << "%) -- est. insert size: "
             << input.getFragmentSize() << " +/- " << input.getFragmentStd() << endl;
        cout << "\tNumber of pairs that map to different nodes: " << mappedPaired
             << " (" << 100.0*mappedPaired/(double)numPairs << "%)" << endl;
        cout << "\tNumber of pairs that can be partially mapped to graph: "
             << mappedIsolated << " (" << 100.0*mappedIsolated/(double)numPairs << "%)" << endl;
        cout << "\tNumber of pairs that can not be mapped to graph: "
             << notMapped << " (" << 100.0*notMapped/(double)numPairs << "%)" << endl;

        input.setNumMappedPaired(mappedPaired);
        input.setNumMappedIsolated(input.getNumMappedIsolated() + mappedIsolated);

        ifs.close();
        ofsi.close();
        ofsp.close();
}

// ============================================================================
// READ MAPPING PUBLIC
// ============================================================================

void DBGraph::mapReads() {

      /*  table = new KmerNodeTable(settings, numNodes);
        table->populateTable(nodes);

        // first of all, detect the read direction type
        cout << "Trying to determine the read direction for "
                "paired-end read file: " << endl;
        for (size_t i = 0; i < settings.getInputCount(); i++) {
                Input &input = const_cast<Input&>(settings.getInput(i));
                if (!input.isPaired())
                        continue;

                cout << "\t" << input.getFilename() << "... "; cout.flush();
                detectPairedReadDirection(input, *table);
                cout << input.getReadDirType() << endl;
        }

        // first of all, map reads to the graph
        Util::startChrono();
        vector<vector<MappedRead> > mappedRead(settings.getInputCount());
        for (size_t i = 0; i < settings.getInputCount(); i++) {
                Input &input = const_cast<Input&>(settings.getInput(i));

                cout << "Mapping reads from file " << i+1 << "/"
                     << settings.getInputCount() << ": "
                     << input.getFilename()
                     << " (type: " << input.getFileType()
                     << ", reads: " << input.getReadType() << ")" << endl;

                mapIsolatedReadsToGraph(input, *table);
                mapPairedReadsToGraph(input, *table);

                input.writeMetadata();
        }

        cout << "Time for mapping reads: " << Util::stopChrono() << endl;*/
}
