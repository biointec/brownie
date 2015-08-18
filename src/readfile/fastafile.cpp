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

#include "fastafile.h"
#include <cstdlib>

using namespace std;

// ============================================================================
// FASTA FILE
// ============================================================================

bool FastAFile::getNextRead(string &read, string &description)
{
        // empty output strings
        read.clear();
        description.clear();

        // read lines until '>' encountered (skip ; lines)
        string dummy;
        while (rfHandler->good() && rfHandler->peekCharacter() != '>')
                rfHandler->getLine();

        // read description lines
        if (rfHandler->good()) {
                const char *result = rfHandler->getLine();
                description.append(result + 1);
        }
        if (!description.empty() && description[description.size() - 1] == '\n')
                description.erase(description.size() - 1);

        // read the actual read
        while (rfHandler->good() && rfHandler->peekCharacter() != '>') {
                const char *result = rfHandler->getLine();
                read.append(result);
                if (!read.empty() && read[read.size() - 1] == '\n')
                        read.erase(read.size() - 1);
        }

        return !read.empty();
}

void FastAFile::writeRead(const string& read, const string& description)
{
        rfHandler->writeChar('>');
        rfHandler->writeLine(description);
        rfHandler->writeLine(read);
}
