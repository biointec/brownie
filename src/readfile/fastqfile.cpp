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

#include "fastqfile.h"
#include <cstdlib>

using namespace std;

// ============================================================================
// FASTQ FILE
// ============================================================================

bool FastQFile::getNextRead(string &read, string &description)
{
        // empty output strings
        read.clear();
        description.clear();

        // read sequence identifier
        char c = rfHandler->getCharacter();

        // end of file might be reached
        if (!good())
                return false;
        if (c != '@')
                throw ios::failure("File doesn't appear to be in fastQ format");
        description = rfHandler->getLine();
        if (!description.empty() && description[description.size() - 1] == '\n')
                description.erase(description.size() - 1);

        // read the actual read
        read = rfHandler->getLine();
        if (!read.empty() && read[read.size() - 1] == '\n')
                read.erase(read.size() - 1);

        // read the + line
        rfHandler->getLine();
        // read the quality scores
        rfHandler->getLine();

        return !read.empty();
}

void FastQFile::writeRead(const string& read, const string& description)
{
        cerr << "Writing of FastQ files is yet implemented" << endl;
        exit(EXIT_FAILURE);
}
