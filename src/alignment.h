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

#include "string.h"
#include "global.h"
#include "tstring.h"
#include <iostream>
using namespace std;



class NW_Alignment {
private:

    double matchScore;
    double mismatchPenalty;
    double gapPenalty;
    double maxGap=3;
    double max(double x, double y);
    double max(double x, double y, double z);
    void    traceback(string& s1,string& s2,char **traceback );
    void init();
public:
    NW_Alignment();
    double alignment(string &s1, string &s2);
    double enhancedAlignment(string &s1, string &s2);
    double get_similarity_per(string s1,string s2);
    double get_similarity_perEnhanced(string s1,string s2);
    int findDirectSim(string const &a, string const &b);
};

