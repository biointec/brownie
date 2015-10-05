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

#ifndef BROWNIE_H
#define BROWNIE_H

#include "settings.h"
#include "library.h"
#include <vector>

class Brownie {

private:
        Settings settings;              // settings object
        LibraryContainer libraries;     // read libraries

public:
        /**
         * Constructor
         * @param argc Argument count
         * @param args Argument strings
         */
        Brownie(int argc, char ** args);

        /**
         * Execute stage one
         */
        void stageOne();

        /**
         * Execute stage two
         */
        void stageTwo();

        /**
         * Execute stage three
         */
        void stageThree();

        /**
         * Execute stage four
         */
        void stageFour();

        /**
         * Execute stage four
         */
        void stageFive();

        /**
         * Get the node filename
         * @return String containing the node filename
         */
        std::string getNodeFilename(int filestage) const {
                char stageStr[4];
                sprintf(stageStr, "%d", filestage);
                return settings.addTempDirectory("nodes.stage") + stageStr;
        }

        /**
         * Get the arc filename
         * @return String containing the arc filename
         */
        std::string getArcFilename(int filestage) const {
                char stageStr[4];
                sprintf(stageStr, "%d", filestage);
                return settings.addTempDirectory("arcs.stage") + stageStr;
        }

        /**
         * Get the stage 2 metadata filename
         * @return String containing the metadata filename
         */
        std::string getMetaDataFilename(int filestage) const {
                char stageStr[4];
                sprintf(stageStr, "%d", filestage);
                return settings.addTempDirectory("metadata.stage") + stageStr;
        }

        /**
         * Get the kmer filename
         * @return The kmer filename
         */
        std::string getKmerFilename() const {
                return settings.addTempDirectory("kmers.stage1");
        }

        /**
         * Check if it is necessary to perform stage one
         * @return True of false
         */
        bool stageOneNecessary() const {
                /*for (size_t i = 0; i < container.size(); i++) {
                        const ReadLibrary &input = container[i];
                        if (!ReadLibrary::fileExists(input.getMetadataFilename()))
                                return true;
                        if (!ReadLibrary::fileExists(input.getIsolatedFilename()))
                                return true;
                        if (!ReadLibrary::fileExists(input.getPairedFilename()))
                                return true;
                }*/

                return !ReadLibrary::fileExists(getKmerFilename());
        }

        /**
         * Check if it is necessary to perform stage two
         * @return True of false
         */
        bool stageTwoNecessary() const {
                if (!ReadLibrary::fileExists(getNodeFilename(2)))
                        return true;
                if (!ReadLibrary::fileExists(getArcFilename(2)))
                        return true;
                return !ReadLibrary::fileExists(getMetaDataFilename(2));
        }

        /**
         * Check if it is necessary to perform stage three
         * @return True or false
         */
        bool stageThreeNecessary() const {
                if (!ReadLibrary::fileExists(getNodeFilename(3)))
                        return true;
                if (!ReadLibrary::fileExists(getArcFilename(3)))
                        return true;
                return !ReadLibrary::fileExists(getMetaDataFilename(3));
        }

        /**
         * Check if it is necessary to perform stage four
         * @return True or false
         */
        bool stageFourNecessary() const {
                if (settings.getSkipStage4())
                        return false;
                if (!ReadLibrary::fileExists(getNodeFilename(4)))
                        return true;
                if (!ReadLibrary::fileExists(getArcFilename(4)))
                        return true;
                return !ReadLibrary::fileExists(getMetaDataFilename(4));
        }

        /**
         * Check if it is necessary to perform stage five
         * @return True of false
         */
        bool stageFiveNecessary() const {
                if (settings.getSkipStage5())
                        return false;
                return true;
        }
        void printInFile();
        void writeGraphExplicit(int stage);
};

#endif
