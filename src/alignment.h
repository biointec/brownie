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

class alignment
{
public:
    void init();
    void cleanup();

    int nw(string  seq_1,std::string   seq_2, string&  seq_1_al, string&  seq_2_al, bool prm);
    void  dpm_init( int ** F, char ** traceback, int L1, int L2, int d );
    int nw_align( int ** F, char ** traceback, string  seq_1, string  seq_2, string&  seq_1_al, string&  seq_2_al,int  d );
    int  max( int f1, int f2, int f3, char * ptr );
    void  print_matrix( int ** F, string seq_1, string seq_2 );
    void  print_traceback( char ** traceback, string seq_1, string seq_2 );
    void  print_al( string& seq_1_al, string& seq_2_al );

};


class NW_Alignment {
private:

    double matchScore;
    double mismatchPenalty;
    double gapPenalty;
    double max(double x, double y);
    double max(double x, double y, double z);
    void    traceback(string& s1,string& s2,char **traceback );
    void init();
public:
    NW_Alignment();
    double alignment(string &s1, string &s2);
    double get_similarity_per(string s1,string s2);
};

