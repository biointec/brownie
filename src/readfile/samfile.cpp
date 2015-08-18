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

#include "samfile.h"

using namespace std;

// ============================================================================
// SAM FILE
// ============================================================================

bool SamFile::getNextRead(string &read, string &description)
{
        // empty output strings
        read.clear();
        description.clear();

        // read past header files
        while (rfHandler->good() && rfHandler->peekCharacter() == '@')
                rfHandler->getLine();

        // read sequence identifier
        const char *result = rfHandler->getLine();
        istringstream oss(result);

        string qname, flag, rname, position, dummy;
        oss >> description >> flag >> rname >> position >> dummy >> dummy >>
                        dummy >> dummy >> dummy >> read;

        if (read.empty())       // end of file might be reached
                return false;

        return true;
}

void SamFile::writeRead(const std::string& read, const std::string& description)
{
        // not yet implemented
}

