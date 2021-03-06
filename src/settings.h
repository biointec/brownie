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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class LibraryContainer;

// ============================================================================
// SETTINGS CLASS
// ============================================================================

class Settings
{
private:
        /**
         * Print Brownie program info
         */
        void printProgramInfo() const;

        /**
         * Print Brownie usage instructions
         */
        void printUsage() const;

        unsigned int kmerSize;          // user specified kmer size
        size_t numThreads;              // number of threads
        bool doubleStranded;            // double stranded sequences
        std::string pathtotemp;         // directory specified by user

        int essaMEMSparsenessFactor;    // sparseness factor for essaMEM
        int bubbleDFSNodeLimit;         // maximum number of visited nodes during bubble detection
        int readCorrDFSNodeLimit;       // maximal number of visited nodes during read mapping
        double covCutoff;               // coverage cutoff value to separate true and false nodes based on their node-kmer-coverage
        bool skipStage4;                // true if stage 4 should be skipped
        bool skipStage5;                // true if stage 5 should be skipped

public:
        /**
         * Default constructor
         */
        Settings();

        /**
         * Parse the command line arguments
         * @param argc Command line argument count
         * @param args Command line arguments
         * @param libCont Library container (output)
         */
        void parseCommandLineArguments(int argc, char **args,
                                       LibraryContainer& libCont);

        /**
         * Get the user-specified kmer hash length
         * @return The hash length
         */
        unsigned int getK() const {
                return kmerSize;
        }

        /**
         * Get the number of threads
         * @return The number of threads
         */
        size_t getNumThreads() const {
                return numThreads;
        }

        /**
         * Check of if the reads are double stranded
         * @return true of false
         */
        bool isDoubleStranded() const {
                return doubleStranded;
        }

        /**
         * Get the temporary working directory
         * @return The temporary working directory
         */
        std::string getTempDirectory() const {
                return pathtotemp;
        }

        /**
         * Prepend the temporary working directory to a file
         * @return The expanded directory
         */
        std::string addTempDirectory(const std::string filename) const {
                return pathtotemp + filename;
        }

        /**
         * Get the IO block size in number of kmers
         * @return The IO block size in number of kmers
         */
        size_t getThreadWorkSize() const {
                return 100000;
        }

        /**
         * Get the essaMEM sparseness factor
         * @return The essaMEM sparseness factor
         */
        int getEssaMEMSparsenessFactor() const {
                return essaMEMSparsenessFactor;
        }

        /**
         * Get the maximum number of nodes visited during a DFS during bubble detection
         * @return The maximum number of nodes visited during bubble detection
         */
        int getBubbleDFSNodeLimit() const {
                return bubbleDFSNodeLimit;
        }

        /**
         * Get the maximum number of nodes visited during a DFS during read correction
         * @return The maximum number of nodes visited during read correction
         */
        int getReadCorrDFSNodeLimit() const {
                return readCorrDFSNodeLimit;
        }

        /**
         * True if stage 4 should be skipped
         * @return True if stage 4 should be skipped
         */
        bool getSkipStage4() const {
                return skipStage4;
        }

        /**
         * True if stage 5 should be skipped
         * @return True if stage 5 should be skipped
         */
        bool getSkipStage5() const {
                return skipStage5;
        }

        /**
         * Get the coverage cutoff value for a node
         * @return The coverage cutoff value for a node
         */
        double getCutOffValue(){
                return covCutoff;
        }
};

#endif
