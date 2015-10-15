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

#include "library.h"

#include "readfile/fastafile.h"
#include "readfile/fastqfile.h"
#include "readfile/rawfile.h"
#include "readfile/sequencefile.h"
#include "readfile/samfile.h"

#include <iostream>
#include <algorithm>
#include <functional>

using namespace std;

// ============================================================================
// FILETYPE ENUM
// ============================================================================

std::ostream &operator<<(std::ostream &out, const FileType &fileType) {
        switch (fileType) {
                case FASTA :
                        out << "fastA";
                        break;
                case FASTA_GZ :
                        out << "fastA.gz";
                        break;
                case FASTQ :
                        out << "fastQ";
                        break;
                case FASTQ_GZ :
                        out << "fastQ.gz";
                        break;
                case SAM:
                        out << "sam";
                        break;
                case SAM_GZ:
                        out << "sam.gz";
                        break;
                case RAW:
                        out << "raw";
                        break;
                case RAW_GZ:
                        out << "raw.gz";
                        break;
                case UNKNOWN_FT:
                        out << "unknown";
                        break;
        }

        return out;
}

// ============================================================================
// READTYPE ENUM
// ============================================================================

std::ostream &operator<<(std::ostream &out, const ReadType &readType) {
        switch (readType) {
                case SHORTREAD :
                        out << "short";
                        break;
                case SHORTPAIRED :
                        out << "short paired";
                        break;
                case LONGREAD :
                        out << "long";
                        break;
                case LONGPAIRED :
                        out << "long paired";
                        break;
                case REFERENCEREAD :
                        out << "reference";
                        break;
        }

        return out;
}

// ============================================================================
// READ DIRECTION TYPE ENUM
// ============================================================================

std::ostream &operator<<(std::ostream &out, const ReadDirType &ReadDirType) {
        switch (ReadDirType) {
                case UNKNOWN :
                        out << "unknown";
                        break;
                case (RR) :
                        out << "[---> ... --->]";
                        break;
                case (RL) :
                        out << "[---> ... <---]";
                        break;
                case (LR) :
                        out << "[<--- ... --->]";
                        break;
                case (LL) :
                        out << "[<--- ... <---]";
                        break;
        }

        return out;
}

// ============================================================================
// READLIBRARY CLASS
// ============================================================================

ReadLibrary::ReadLibrary(const std::string& inputFilename_,
                         const std::string& outputFilename_) :
        inputFilename(inputFilename_), outputFilename(outputFilename_),
        fileType(UNKNOWN_FT)
{
        // try to figure out the file format based on the extension
        string extension;

        if (inputFilename.length() >= 4)
                extension = inputFilename.substr(inputFilename.length() - 4);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == ".SAM") {
                fileType = SAM;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 4);
        } else if (extension == ".RAW") {
                fileType = RAW;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 4);
        }

        if (inputFilename.length() >= 6)
                extension = inputFilename.substr(inputFilename.length() - 6);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == ".FASTA") {
                fileType = FASTA;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 6);
        } else if (extension == ".FASTQ") {
                fileType = FASTQ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 6);
        }

        if (inputFilename.length() >= 7)
                extension = inputFilename.substr(inputFilename.length() - 7);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == ".RAW.GZ") {
                fileType = RAW_GZ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 7);
        } else if (extension == ".SAM.GZ") {
                fileType = SAM_GZ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 7);
        }

        if (inputFilename.length() >= 9)
                extension = inputFilename.substr(inputFilename.length() - 9);
        transform(extension.begin(), extension.end(), extension.begin(), ::toupper);
        if (extension == ".FASTA.GZ") {
                fileType = FASTA_GZ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 9);
        } else if (extension == ".FASTQ.GZ") {
                fileType = FASTQ_GZ;
                baseFilename = inputFilename.substr(0, inputFilename.length() - 9);
        }

        if (fileType == UNKNOWN_FT) {
                cerr << "Brownie: don't know how to open file: '" << inputFilename << "'\n";
                cerr << "Expected one of the following extensions: .fastq, .fasta, .sam, .raw (or .gz variants thereof)\n";
                exit(EXIT_FAILURE);
        }

        if (!Util::fileExists(inputFilename)) {
                cerr << "Brownie: cannot open input read file: '" << inputFilename << "'\n";
                cerr << "Exiting... " << endl;
                exit(EXIT_FAILURE);
        }

        // set the outputFilename only if not specified by the user
        if (outputFilename.empty()) {
                ostringstream oss;
                oss << baseFilename << "." << fileType;
                outputFilename = oss.str();
        }
}

ReadFile* ReadLibrary::allocateReadFile() const
{
        switch (getFileType()) {
                case FASTA :
                        return new FastAFile ( false );
                case FASTA_GZ :
                        return new FastAFile ( true );
                case FASTQ :
                        return new FastQFile ( false );
                case FASTQ_GZ :
                        return new FastQFile ( true );
                case SAM :
                        return new SamFile ( false );
                case SAM_GZ :
                        return new SamFile ( true );
                case RAW :
                        return new RawFile ( false );
                case RAW_GZ :
                        return new RawFile ( true );
                case UNKNOWN_FT:
                        assert(false);
                        return NULL;
        }

        assert ( false );
        return NULL;
}

// ============================================================================
// LIBRARY CONTAINER CLASS
// ============================================================================

bool LibraryContainer::getReadChunk(vector<string>& localBuffer,
                                    size_t targetNumKmers)
{
        // clear the buffer
        localBuffer.clear();

        // wait until reads become available
        std::unique_lock<std::mutex> lock(inputReadMutex);
        readBufFull.wait(lock, [this]{return ((!actReadBuf->empty()) ||
                                              (!readThreadWorking));});

        size_t currBlockSize = 0;
        while (currBlockSize < targetNumKmers && (!actReadBuf->empty())) {
                size_t numKmers = actReadBuf->back().size() + 1 - Kmer::getK();
                localBuffer.push_back(actReadBuf->back());
                actReadBuf->pop_back();
                currBlockSize += numKmers;
        }

        // if you took the last reads, notify the input thread
        if (actReadBuf->empty())
                readBufEmpty.notify_one();

        lock.unlock();

        return !localBuffer.empty();
}

void LibraryContainer::countReadFrequency( ReadLibrary& input)
{
        const size_t targetBufferSize = 100000; /*settings.getNumThreads() *
                                        settings.getThreadWorkSize();*/

        // read counters
        size_t totNumReads = 0, numTooShort = 0, totReadLength = 0;

        // open the read file
        ReadFile *readFile = input.allocateReadFile();
        readFile->open(input.getFilename());

        // aux variables
        string read, description;

        while (true) {
                // fill up the idle read buffer (only this thread has access)
                idlReadBuf->clear();
                size_t thisBufferSize = 0;
                while (thisBufferSize < targetBufferSize) {
                        readFile->getNextRead(read, description);
                        if (!readFile->good() && read.empty())
                                break;

                        totNumReads++;
                        totReadLength += read.size();

                        // if the read is short than k, no need to process it
                        if (read.size() < Kmer::getK()) {
                                numTooShort++;
                                continue;
                        }

                        idlReadBuf->push_back(read);
                        thisBufferSize += read.size() + 1 - Kmer::getK();
                }

                cout << "Number of reads processed: " << totNumReads << "\r";

                // wait until active buffer is empty
                std::unique_lock<std::mutex> lock(inputReadMutex);
                readBufEmpty.wait(lock, [this]{return actReadBuf->empty();});

                // swap active buffer and idle buffer
                std::swap(actReadBuf, idlReadBuf);

                // notify workers that more work is available
                readBufFull.notify_all();
                lock.unlock();

                // file has completely been read
                if (!readFile->good())
                        break;
        }

        cout << "Number of reads processed: " << totNumReads << endl;
        cout << "Average read length: " << totReadLength/totNumReads << endl;
	input.setReadLength(totReadLength/totNumReads);
	input.setNumOfReads(totNumReads);

        readFile->close();

        // free temporary memory
        delete readFile;

        // issue a warning about skipped reads
        if (numTooShort > 0) {
                double perc = 100.0 * (double)numTooShort/(double)totNumReads;
                streamsize prev = cout.precision (2);

                cout << "\n\tWARNING: Skipped " << numTooShort << " reads ("
                     << perc << "% of total) because they were shorter than "
                     << "the kmer size" << endl;

                cout.precision(prev);
        }
}

void LibraryContainer::initiateReadThreading()
{
        readThreadWorking = true;
        actReadBuf = new vector<string>;
        idlReadBuf = new vector<string>;
}

void LibraryContainer::threadReads()
{
        // read all input data
        for (size_t i = 0; i < getSize(); i++) {
                ReadLibrary &input = getInput(i);

                cout << "Processing file " << i+1 << "/" << getSize() << ": "
                     << input.getFilename() << ", type: " << input.getFileType()
                     << endl;

                countReadFrequency(input);

        }

        // wait until active buffer is empty
        unique_lock<std::mutex> lock(inputReadMutex);
        readThreadWorking = false;
        readBufFull.notify_all();
        lock.unlock();
}

void LibraryContainer::finalizeReadThreading()
{
        delete idlReadBuf; idlReadBuf = NULL;
        delete actReadBuf; actReadBuf = NULL;
}
