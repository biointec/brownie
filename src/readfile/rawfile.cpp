/***************************************************************************
 *   Copyright (C) 2010 Jan Fostier (jan.fostier@intec.ugent.be)           *
 *   Original Velvet code by Daniel Zerbino (zerbino@ebi.ac.uk)            *
 *                                                                         *
 *   This file is part of Velvet 2.0                                       *
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

#include "rawfile.h"

using namespace std;

// ============================================================================
// FASTA FILE
// ============================================================================

bool RawFile::getNextRead(std::string &read, std::string &description)
{
        // empty output strings
        read.clear();
        description.clear();

        // read until something non-empty is encountered
        while (rfHandler->good() && read.empty()) {
                const char * buffer = rfHandler->getLine();
                read.append(buffer);
                if (!read.empty() && read[read.size() - 1] == '\n')
                        read.erase(read.size() - 1);
        }

        return !read.empty();
}

void RawFile::writeRead(const std::string &read, const std::string &description)
{
        // in the raw format, we simply write the reads, and nothing else
        rfHandler->writeLine(read);
}
