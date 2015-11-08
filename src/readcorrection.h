/***************************************************************************
 *   Copyright (C) 2015 Jan Fostier (jan.fostier@intec.ugent.be)           *
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

#ifndef READCORRECTION_H
#define READCORRECTION_H

#include "settings.h"
#include "graph.h"

class ReadCorrectionJan
{
private:
        const DBGraph &dbg;
        const Settings &settings;

        /**
         * Entry routine for worker thread
         * @param myID Unique threadID
         * @param libaries Library container with libraries to be corrected
         */
        void workerThread(size_t myID, LibraryContainer& libraries);

public:
        /**
         * Default constructor
         */
        ReadCorrectionJan(const DBGraph& g, const Settings& s) :
                dbg(g), settings(s) {}

        /**
         * Perform error correction in the libaries
         * @param libraries Library container with libraries to be corrected
         */
        void doErrorCorrection(LibraryContainer &libraries);
};

#endif
