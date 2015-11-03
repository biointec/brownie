/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
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

#include <thread>
#include "library.h"
#include "readcorrection.h"

using namespace std;

void ReadCorrectionJan::workerThread(size_t myID, LibraryContainer& libraries)
{
        // local storage of reads
        vector<ReadRecord> myReadBuf;

        while (true) {
                size_t blockID, recordID;
                bool result = libraries.getRecordChunk(myReadBuf, blockID, recordID);

                if (result)
                        libraries.commitRecordChunk(myReadBuf, blockID, recordID);

                if (!result)
                        break;
        }
}

void ReadCorrectionJan::doErrorCorrection(LibraryContainer& libraries)
{
        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        libraries.startIOThreads(settings.getThreadWorkSize(),
                                 settings.getThreadWorkSize() * settings.getNumThreads(),
                                 true);

        // read all input data
        //libraries.threadReads();

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&ReadCorrectionJan::workerThread,
                                          this, i, ref(libraries));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        libraries.joinIOThreads();
}
