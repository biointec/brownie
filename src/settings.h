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
        size_t genomeSize;              // genome size
        bool doubleStranded;            // double stranded sequences
        std::string pathtotemp;         // directory specified by user

        int ESSASparsenessFactor;       // sparseness factor for essamem
        int bubbleDFSNodeLimit;         // maximal number of visited nodes for bubble detection
        int readCorrDFSNodeLimit;       // maximal search depth for stage 5
        int essa_factor;             // sparseness factor for essamem
        int max_visits;                 // maximal number of visited nodes for bubble detection
        int max_depth;                  // maximal search depth for stage 5
        double cutoff;                   //cutoff value to separate true and false nodes based on their node-kmer-coverage
        bool skip_stage_4;
        bool skip_stage_5;

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

         void parseCommandLineArgumentsMain(int argc, char** args,LibraryContainer& libCont);
         /**
          * @param libCont Library container (output)
          */
         void checkInputArguments(LibraryContainer& libCont);

        /**
         * Parse the command line arguments
         * @param argc Command line argument count
         * @param args Command line arguments
         * @param libCont Library container (output)
         */
        void parseCommandLineArgumentsEC(int argc, char **args,
                                       LibraryContainer& libCont, size_t startArg);



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

        size_t getGenomeSize() const{
                return this->genomeSize;
        }

        void setGenomeSize(){
                this->genomeSize=genomeSize;
        }

        int getESSASparsenessFactor() const {
                return ESSASparsenessFactor;
        }

        int getBubbleDFSNodeLimit() const {
                return bubbleDFSNodeLimit;
        }

        int getReadCorrDFSNodeLimit() const {
                return readCorrDFSNodeLimit;
        }
        double getCutOffValue(){
                return cutoff;
        }
        bool getSkipStage4() const {
                return skip_stage_4;
        }

        bool getSkipStage5() const {
                return skip_stage_5;
        }
};

#endif
