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

#ifndef READLIBRARY_H
#define READLIBRARY_H

#include "global.h"
#include "tkmer.h"
#include "gaussval.h"

#include <string>
#include <sstream>
#include <fstream>

#include <mutex>
#include <condition_variable>
#include <atomic>

// ============================================================================
// FILETYPE ENUM
// ============================================================================

typedef enum { FASTQ, FASTA, FASTQ_GZ, FASTA_GZ,
               SAM, SAM_GZ, RAW, RAW_GZ, UNKNOWN_FT} FileType;

std::ostream &operator<<(std::ostream &out, const FileType &fileType);

// ============================================================================
// READTYPE ENUM
// ============================================================================

typedef enum { SHORTREAD, SHORTPAIRED, LONGREAD,
               LONGPAIRED, REFERENCEREAD } ReadType;

std::ostream &operator<<(std::ostream &out, const ReadType &readType);

// ============================================================================
// READ DIRECTION ENUM
// ============================================================================

typedef enum { UNKNOWN = 0, LL = 1, LR = 2, RR = 3, RL = 4 } ReadDirType;

std::ostream &operator<<(std::ostream &out, const ReadDirType &readType);

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class ReadFile;

// ============================================================================
// READ LIBRARY CLASS
// ============================================================================

class ReadLibrary
{
private:
        FileType fileType;              // file type (FASTQ, FASTA, etc.)

        // TODO below: not set not initialized yet
        ReadType readType;

        std::string inputFilename;      // name of the input file
        std::string outputFilename;     // name of the output file
        std::string baseFilename;       // base filename

        std::string filenamePairedReads; // for paired reads in separate files
        std::string directory;

        // ==== begin of metadata ====
        size_t readLength;              // size of the reads
        size_t numIsolatedReads;        // number of reads that are isolated
        size_t numPairedReads;          // number of reads that are paired

        ReadDirType readDirType;        // direction of the reads

        size_t numMappedIsolated;
        size_t numMappedPaired;

        double fragLen;                 // fragment length
        double fragStd;                 // fragment standard deviation

        double centralInsertSize;
        double centralStddev;
        double genomeSize;
        size_t numOfReads;
        double coverage;                // coverage of this library
                                        // == number of times a nucleotide is contained in a read
        double NOM(double lA, double y) const {
                double Sa = fragLen - 2.0 * readLength + Kmer::getK() - 0.5;
                double Ea = lA + fragLen - readLength - 0.5;
                double oosts = 0.70710678118654752440 / fragStd; // 1/sqrt(2)/sigma

                double argS = oosts * (y-Sa);
                double argE = oosts * (y-Ea);
                return 0.5 * fragStd * fragStd * getFragStartCov() * (erf(argS) - erf(argE));
        }

        double DEN(double lA, double x) const {
                double Sa = fragLen - 2.0 * readLength + Kmer::getK() - 0.5;
                double Ea = lA + fragLen - readLength - 0.5;
                double oosts = 0.70710678118654752440 / fragStd;   // 1/sqrt(2)/sigma
                double stops = 0.7978845608028653558 * fragStd;    // sqrt (2/Pi)*sigma

                double Sax = Sa - x;
                double Eax = Ea - x;
                double Saxs = Sax * oosts;
                double Eaxs = Eax * oosts;

                return 0.5 * getFragStartCov() * ( erf(Saxs)*Sax - erf(Eaxs)*Eax + stops*(exp(-Saxs*Saxs)-exp(-Eaxs*Eaxs)) );
        }

        double DENZ(double lA, double x) const {
                double Sa = fragLen - 2.0 * readLength + Kmer::getK() - 0.5;
                double Ea = lA + fragLen - readLength - 0.5;

                if (x > Sa) {
                        if (x > Ea)
                                return getFragStartCov() * (Ea - Sa);
                        else
                                return getFragStartCov() * (x - Sa);
                }
                else
                        return 0;
        }

        // ==== end of metadata ====

		/**
		 * Parses the Basename from the filename
		 */
		void parseBasename() {
			size_t location = inputFilename.find_last_of('/');
			if(location != inputFilename.npos) {
				// was found so get substring from there!
				baseFilename = inputFilename.substr(location + 1);
			} else {
				//not found, relative filename, just copy to basename
				baseFilename = inputFilename;
			}
		}

public:
        /**
         * Default constructor
         * @param inputFilename Input filename
         * @param outputFilename Output filename
         */
        ReadLibrary(const std::string& inputFilename,
                    const std::string& outputFilename);

        /**
         * Default constructor
         * @param fileType The file type (fastA, fastQ, etc)
         * @param readType The read type (short, long, reference, etc)
         * @param filename The filename
         * @param directory Directory to store temporary output files
         */
        ReadLibrary(FileType fileType_, ReadType readType_, std::string filename_,
              std::string directory_) : fileType(fileType_),
              readType(readType_), inputFilename(filename_), filenamePairedReads(""),
              directory(directory_), readLength(0),
              numIsolatedReads(0), numPairedReads(0), readDirType(UNKNOWN),
              numMappedIsolated(0), numMappedPaired(0), fragLen(0.0),
              fragStd(0.0) { parseBasename(); }

        /**
         * Default constructor (paired files)
         * @param fileType The file type (fastA, fastQ, etc)
         * @param readType The read type (short, long, reference, etc)
         * @param filename The filename
         * @param pairedFilename The filename of the paired reads
         * @param directory Directory to store temporary output files
         */
        ReadLibrary(FileType fileType_, ReadType readType_, std::string filename_,
              std::string filenamePairedReads_, std::string directory_) :
              fileType(fileType_), readType(readType_), inputFilename(filename_),
              filenamePairedReads(filenamePairedReads_), directory(directory_),
              readLength(0), numIsolatedReads(0),
              numPairedReads(0), readDirType(UNKNOWN), numMappedIsolated(0),
              numMappedPaired(0), fragLen(0.0), fragStd(0.0) { parseBasename(); }

        /**
         * Allocate and return the correct read file for a certain input
         * @return An allocated readfile
         */
        ReadFile* allocateReadFile() const;

        //coment by mahdi
        size_t getGenomeSize() const {
                return genomeSize;
        }

        void setGenomeSize(size_t target) {
                genomeSize = target;
        }

        size_t getNumOfReads() const {
                return numOfReads;
        }

        void setNumOfReads(size_t target) {
                numOfReads = target;
        }
        std::string getOutputFileName(){
	  return this->outputFilename;
	}





        /**
         * Get the file type
         * @return The file type
	 *
	 *
         */

        FileType getFileType() const {
                return fileType;
        }

        /**
         * Get the read type
         * @return The read type
         */
        ReadType getReadType() const {
                return readType;
        }

        /**
         * Get the read direction type
         * @return The read direction type
         */
        ReadDirType getReadDirType() const {
                return readDirType;
        }

        /**
         * Get the read length
         * @return The read length
         */
        size_t getReadLength() const {
                return readLength;
        }

        /**
         * Check whether read 1 needs to be RC'ed
         * @return true of false
         */
        bool needToRCRead1() const {
                return ((readDirType == LR) || (readDirType == LL));
        }

        /**
         * Check whether read 2 needs to be RC'ed
         * @return true of false
         */
        bool needToRCRead2() const {
                return ((readDirType == RL) || (readDirType == LL));
        }

        /**
         * Set the read direction type
         * @param Read direction type
         */
        void setReadDirType(ReadDirType readDirType_) const {
                const_cast<ReadDirType&>(readDirType) = readDirType_;
        }

        /**
         * Get the Basename
         * @return The basename
         */
        std::string getBasename() const {
                return baseFilename;
        }

        /**
         * Get the filename
         * @return The filename
         */
        std::string getFilename() const {
                return inputFilename;
        }

        /**
         * Get the filename of the paired reads
         * @return The filename
         */
        std::string getFilenamePairedReads() const {
                return filenamePairedReads;
        }

        /**
         * Check whether we deal with paired-end reads
         * @return True or false
         */
        bool isPaired() const {
                return (readType == SHORTPAIRED || readType == LONGPAIRED);
        }

        /**
         * Get the isolated paired filename
         * @return The isolated filename
         */
        std::string getIsolatedFilename() const {
                std::ostringstream oss;
                oss << directory << "/" << baseFilename << "_i.stage1";
                return oss.str();
        }

        /**
         * Get the internal paired filename
         * @return The paired filename
         */
        std::string getPairedFilename() const {
                std::ostringstream oss;
                oss << directory << "/" << baseFilename << "_p.stage1";
                return oss.str();
        }

        /**
         * Get the metadata filename
         * @return The metadata filename
         */
        std::string getMetadataFilename() const {
                std::ostringstream oss;
                oss << directory << "/" << baseFilename << "_md.txt";
                return oss.str();
        }

        /**
         * Get the mapped isolated read filename
         * @return The mapped paired readfilename
         */
        std::string getMappedIsolatedFilename() const {
                std::ostringstream oss;
                oss << directory << "/" << baseFilename << "_mi.stage4";
                return oss.str();
        }

        /**
         * Get the mapped paired read filename
         * @return The mapped paired readfilename
         */
        std::string getMappedPairedFilename() const {
                std::ostringstream oss;
                oss << directory << "/" << baseFilename << "_mp.stage4";
                return oss.str();
        }

        /**
         * Set the insert size
         * @param target Target insert size
         */
        void setInsertSize(double target) {
                fragLen = target;
        }

        /**
         * Get the insert size
         * @return The insert size
         */
        double getFragmentSize() const {
                return fragLen;
        }

        /**
         * Get the marginal insert size
         * @return The marginal insert size
         */
        double getMarginalInsertSize() const {
                return fragLen - Kmer::getK() + 1;
        }

        /**
         * Get the distance from which the estimates are unbiased
         * @return The unbiased start position
         */
        double getUnbiasedStartPos() const {
                return getMarginalInsertSize() + 3.0 * getFragmentStd();
        }

        /**
         * Set the standard deviation
         * @param target Target standard deviation
         */
        void setStddev(double target) {
                fragStd = target;
        }

        /**
         * Get the standard variation of the insert size
         * @return The standard variation
         */
        double getFragmentStd() const {
                return fragStd;
        }

        /**
         * Set the insert size
         * @param target Target insert size
         */
        void setCentralInsertSize(double target) {
                centralInsertSize = target;
        }

        /**
         * Get the insert size
         * @return The insert size
         */
        double getCentralInsertSize() const {
                return centralInsertSize;
        }

        /**
         * Set the standard deviation
         * @param target Target standard deviation
         */
        void setCentralStddev(double target) {
                centralStddev = target;
        }

        /**
         * Get the standard variation of the insert size
         * @return The standard variation
         */
        double getCentralStdDev() const {
                return centralStddev;
        }

        /**
         * Read the metadata from file
         */
        void writeMetadata() const {
                std::ofstream ofs(getMetadataFilename().c_str());
                ofs << readLength << std::endl
                    << numIsolatedReads << "\t" << numPairedReads << std::endl
                    << numMappedIsolated << "\t" << numMappedPaired << std::endl
                    << coverage << std::endl
                    << (int)readDirType << std::endl
                    << fragLen << "\t" << fragStd << std::endl
                    << centralInsertSize << "\t" << centralStddev << std::endl;
                ofs.close();
        }

        /**
         * Read the metadata from file
         */
        void readMetadata() {
                std::ifstream ifs(getMetadataFilename().c_str());
                if (!ifs)
                        return;
                int readDirTypeInt;
                ifs >> readLength
                    >> numIsolatedReads >> numPairedReads
                    >> numMappedIsolated >> numMappedPaired
                    >> coverage
                    >> readDirTypeInt
                    >> fragLen >> fragStd
                    >> centralInsertSize >> centralStddev;
                ifs.close();
                readDirType = (ReadDirType)readDirTypeInt;
        }

        /**
         * Set the read length
         * @param target The read length
         */
        void setReadLength(size_t target) {
                readLength = target;
        }

        /**
         * Set the number of isolated reads
         * @param target The number of isolated reads
         */
        void setNumIsolatedReads(size_t target) {
                numIsolatedReads = target;
        }

        /**
         * Set the number of paired-end reads
         * @param target The number of paired-end reads
         */
        void setNumPairedReads(size_t target) {
                numPairedReads = target;
        }

        /**
         * Get the number of isolated reads
         * @return The number of isolated reads
         */
        size_t getNumIsolatedReads() const {
                return numIsolatedReads;
        }

        /**
         * Get the number of paired-end reads
         * @return The number of paired-end reads
         */
        size_t getNumPairedReads() const {
                return numPairedReads;
        }

        /**
         * Get the total number of reads
         * @return The total number of reads
         */
        size_t getNumReads() const {
                return numPairedReads + numIsolatedReads;
        }

        /**
         * Set the number of mapped isolated reads
         * @param target The number of mapped isolated reads
         */
        void setNumMappedIsolated(size_t target) {
                numMappedIsolated = target;
        }

        /**
         * Set the number of mapped paired-end reads
         * @param target The number of mapped paired-end reads
         */
        void setNumMappedPaired(size_t target) {
                numMappedPaired = target;
        }

        /**
         * Get the number of mapped isolated reads
         * @return The number of mapped isolated reads
         */
        size_t getNumMappedIsolated() const {
                return numMappedIsolated;
        }

        /**
         * Get the number of mapped paired-end reads
         * @return The number of mapped paired-end reads
         */
        size_t getNumMappedPaired() const {
                return numMappedPaired;
        }

        /**
         * Get the closest k-mer for which paired-end reads will contain info
         * @param MNL Marginal node length
         * @return The closest k-mer position
         */
        double getPERKmerBegin(double MNL) const {
                return std::max<double>(fragLen - 2.0 * readLength +
                        Kmer::getK() - 0.5, MNL - 0.5);
        }

        /**
         * Get the furthest k-mer for which paired-end reads will contain info
         * @param MNL Marginal node length
         * @return The furthest k-mer position
         */
        double getPERKmerEnd(double MNL) const {
                return MNL + fragLen - Kmer::getK() - 0.5 ;
        }

        /**
         * Get the closest k-mer for which paired-end reads will contain info within 3 sigma
         * @param MNL Marginal node length
         * @return The closest k-mer position
         */
        GaussVal getPERKmerBeginGV(double MNL) const {
                return GaussVal(getPERKmerBegin(MNL), fragStd);
        }

        /**
         * Get the furthest k-mer for which paired-end reads will contain info within 3 sigma
         * @param MNL Marginal node length
         * @return The furthest k-mer position
         */
        GaussVal getPERKmerEndGV(double MNL) const {
                return GaussVal(getPERKmerEnd(MNL), fragStd);
        }

        /**
         * Get the closest k-mer for which paired-end reads will contain unbiased info
         * @param MNL Marginal node length
         * @return The closest k-mer position
         */
        double getPERKmerUnbiasedBegin(double MNL) const {
                // Note the important difference with getPERKmerUnbiasedEnd:
                // we do not use getPERKmerBeginGV because the we want to omit
                // the use of the max() function.  The unbiased region
                // can extent within the "source" node !!
                return fragLen - 2.0 * readLength + Kmer::getK() - 0.5 + 1.5 * getFragmentStd();
        }

        /**
         * Get the furthest k-mer for which paired-end reads will contain unbiased info
         * @param MNL Marginal node length
         * @return The furthest k-mer position
         */
        double getPERKmerUnbiasedEnd(double MNL) const {
                return getPERKmerEndGV(MNL).getAverage() -1.5 * getFragmentStd();
        }

        /**
         * Get the smallest read start position where paired-end reads will land
         * @param MNL Marginal node length
         * @return The closest read start position
         */
        double getPERReadStartBegin(double MNL) const {
                return std::max<double>(fragLen - 2.0 * readLength +
                        Kmer::getK() - 0.5, MNL - 0.5);
        }

        /**
         * Get the furthest read start position where paired-end reads will land
         * @param MNL Marginal node length
         * @return The furthest read start position
         */
        double getPERReadStartEnd(double MNL) const {
                return MNL + fragLen - readLength - 0.5;
        }

        /**
         * Get the smallest read start position where paired-end reads will land
         * @param MNL Marginal node length
         * @return The closest read start position
         */
        double getReadKmerBegin(double MNL) const {
                return MNL - 0.5;
        }

        /**
         * Get the furthest read start position where paired-end reads will land
         * @param MNL Marginal node length
         * @return The furthest read start position
         */
        double getReadKmerEnd(double MNL) const {
                return MNL + readLength - Kmer::getK() - 0.5;
        }

        /**
         * Get the smallest read start position where paired-end reads will land
         * @param MNL Marginal node length
         * @return The closest read start position
         */
        GaussVal getReadKmerBeginGV(double MNL) const {
                return GaussVal(getReadKmerBegin(MNL), 0.0);
        }

        /**
         * Get the furthest read start position where paired-end reads will land
         * @param MNL Marginal node length
         * @return The furthest read start position
         */
        GaussVal getReadKmerEndGV(double MNL) const {
                return GaussVal(getReadKmerEnd(MNL), 0.0);
        }

        /**
         * Set the coverage of this library (= number of times a nucleotide is covered)
         * @return The coverage of this library
         */
        void setCoverage(double targetCoverage) {
                coverage = targetCoverage;
        }

        /**
         * Set the coverage of this library (= number of times a nucleotide is covered)
         * @return The coverage of this library
         */
        void setReadStartCoverage(double targetCoverage) {
                coverage = targetCoverage * (getReadLength() - Kmer::getK() + 1.0);
        }

        /**
         * Get the coverage of this library (= number of times a nucleotide is covered)
         * @return The coverage of this library
         */
        double getCoverage() const {
                return coverage;
        }

        /**
         * Get the read start coverage of this library
         * @return The read start coverage of this library
         */
        double getReadStartCov() const {
                return getCoverage() / (getReadLength() - Kmer::getK() + 1.0);
        }

        /**
         * Get the fragment start coverage of this library
         * @return The fragment start coverage of this library
         */
        double getFragStartCov() const {
                return 0.5 * getReadStartCov();
        }

        double getNumExpectedReads(double lA, double PB, double lB) const {
                double numReads = getReadStartCov() * (readLength - Kmer::getK() - (PB - lA));
                if (numReads < 0)
                        numReads = 0;
                return numReads;
        }

        double getNumExpectedPER(double lA, double PB, double lB) const {
                if (fragStd < 0.001)
                        return DENZ(lA, PB+lB-0.5)-DENZ(lA, PB-readLength+Kmer::getK()-0.5);
                else
                        return DEN(lA, PB+lB-0.5)-DEN(lA, PB-readLength+Kmer::getK()-0.5);
        }

        double getExpOffset(double lA, double PB, double lB) const {
                if (fragStd < 0.001)
                        return 0;

                double intBE = PB + lB - 0.5;
                double intBS = PB-readLength+Kmer::getK()-0.5;
                double A = NOM(lA, intBE);
                double B = NOM(lA, intBS);
                double C = DEN(lA, intBE);
                double D = DEN(lA, intBS);

                return (A - B) / (C - D);
        }
};

// ============================================================================
// LIBRARYCONTAINER CLASS
// ============================================================================

class LibraryContainer
{
private:
        std::vector<ReadLibrary> container;     // all library files

        bool readThreadWorking;                 // read thread still working
        std::mutex inputReadMutex;              // read buffer mutex
        std::condition_variable readBufFull;    // read buffer full condition
        std::condition_variable readBufEmpty;   // read buffer empty condition

        std::vector<std::string> *actReadBuf;   // active read buffer
        std::vector<std::string> *idlReadBuf;   // idle read buffer

        void countReadFrequency( ReadLibrary& input);

public:
        /**
         * Add a read library to the container
         * @arg library Library to add
         */
        void insert(const ReadLibrary& library) {
                container.push_back(library);
        }

        /**
         * Get the number of inputs
         * @return The number of inputs
         */
        size_t getSize() const {
                return container.size();
        }

        /**
         * Get a specified input
         * @param index Input identifier
         * @return A reference to the input
         */
        ReadLibrary& getInput(size_t index) {
                assert(index < container.size());
                return container[index];
        }

        /**
         * Get a specified input
         * @param index Input identifier
         * @return A reference to the input
         */
        const ReadLibrary& getInput(size_t index) const {
                assert(index < container.size());
                return container[index];
        }

        /**
         * Write the metadata for all input files
         */
        void writeMetadata() const {
                for (size_t i = 0; i < getSize(); i++) {
                        const ReadLibrary &input = getInput(i);
                        input.writeMetadata();
                }
        }

        /**
         * Load the metadata for all input files
         */
        void readMetadata() {
                for (size_t i = 0; i < getSize(); i++) {
                        ReadLibrary &input = getInput(i);
                        input.readMetadata();
                }
        }

        /**
         * Get next read chunk from the input
         * @param buffer Buffer in which to store the reads
         * @param targetNumKmers Target number of targetKmers
         * @return False if no more reads are available, true otherwise
         */
        bool getReadChunk(std::vector<std::string>& localBuffer,
                          size_t targetNumKmers);

        /**
         * Initialize read threading
         */
        void initiateReadThreading();

        /**
         * Get all the reads in chunks
         */
        void threadReads();

        /**
         * Finalize read threading
         */
        void finalizeReadThreading();
        double getReadLength(){
                size_t readLengthAvg=0;
                size_t totalNumOfReads=0;
                for (auto it : container){
                        readLengthAvg=((readLengthAvg *totalNumOfReads)+it.getReadLength()*it.getNumReads())/(totalNumOfReads+it.getNumReads());
                        totalNumOfReads=totalNumOfReads+it.getNumReads();
                }
                return readLengthAvg;

        }
};

#endif
