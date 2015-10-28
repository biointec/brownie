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

#include "alignment.h"
#include <string>
#include <stdlib.h>
using namespace std;

void alignment::init()
{
    // Called before each testfunction is executed
}

void alignment::cleanup()
{
    // Called after every testfunction
}

int alignment::nw(
        string       seq_1,          /*  Needleman-Wunsch   */
        string       seq_2,          /*  algorithm for      */
        string&      seq_1_al,       /*  global alignment   */
        string&      seq_2_al,       /*  of nt sequence.    */
        bool         prm
      )
{
        int  d = 2 ;                 /* gap penalty */

        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

        // Dynamic programming matrix
        int ** F = new int * [ L2+1 ];
        for( int i = 0; i <= L2; i++ )  F[ i ] = new int [ L1 ];

        // Traceback matrix
        char ** traceback = new char * [ L2+1 ];
        for( int i = 0; i <= L2; i++ )  traceback[ i ] = new char [ L1 ];

        // Initialize traceback and F matrix (fill in first row and column)
        dpm_init( F, traceback, L1, L2, d );

        // Create alignment
        nw_align( F, traceback, seq_1, seq_2, seq_1_al, seq_2_al, d );

        #if DEBUG
            int  L_al = seq_1_al.length();
            cout << "Length after alignment: " << L_al << endl;
        #endif

        if( prm )
        {
                cout << "\nDynamic programming matrix: " << "\n\n";
                print_matrix( F, seq_1, seq_2 );

                cout << "\nTraceback matrix: " << "\n\n";
                print_traceback( traceback, seq_1, seq_2 );

                cout << endl;
        }

        //for( int i = 0; i <= L2; i++ )  delete F[ i ];
        //delete[] F;
        //for( int i = 0; i <= L2; i++ )  delete traceback[ i ];
        //delete[] traceback;

        return  0 ;
}
void  alignment::dpm_init( int ** F, char ** traceback, int L1, int L2, int d )
{
        F[ 0 ][ 0 ] =  0 ;
        traceback[ 0 ][ 0 ] = 'n' ;

        int i=0, j=0;

        for( j = 1; j <= L1; j++ )
        {
                F[ 0 ][ j ] =  -j * d ;
                traceback[ 0 ][ j ] =  '-' ;
        }
        for( i = 1; i <= L2; i++ )
        {
                F[ i ][ 0 ] =  -i * d ;
                traceback[ i ][ 0 ] =  '|' ;
        }
}

int  alignment::max( int f1, int f2, int f3, char * ptr )
{
        int  max = 0 ;

        if( f1 >= f2 && f1 >= f3 )
        {
                max = f1 ;
                *ptr = '|' ;
        }
        else if( f2 > f3 )
        {
                max = f2 ;
                *ptr = '\\' ;
        }
        else
        {
                max = f3 ;
                *ptr = '-' ;
        }

        return  max ;
}


int alignment::nw_align(                  // Needleman-Wunsch algorithm
              int **     F,
              char **    traceback,
              string     seq_1,
              string     seq_2,
              string&    seq_1_al,
              string&    seq_2_al,
              int        d         // Gap penalty
            )
{
        int        k = 0, x = 0, y = 0;
        int        fU, fD, fL ;
        char       ptr, nuc ;
        int        i = 0, j = 0;

        const int  a =  2;   // Match
        const int  b = -1;   // Mismatch

        const int  s[ 4 ][ 4 ] = { { a, b, b, b },    /* substitution matrix */
                                   { b, a, b, b },
                                   { b, b, a, b },
                                   { b, b, b, a } } ;

        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

        for( i = 1; i <= L2; i++ )
        {
                for( j = 1; j <= L1; j++ )
                {
                        nuc = seq_1[ j-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  x = 0 ;  break ;
                                case 'C':  x = 1 ;  break ;
                                case 'G':  x = 2 ;  break ;
                                case 'T':  x = 3 ;
                        }

                        nuc = seq_2[ i-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  y = 0 ;  break ;
                                case 'C':  y = 1 ;  break ;
                                case 'G':  y = 2 ;  break ;
                                case 'T':  y = 3 ;
                        }

                        fU = F[ i-1 ][ j ] - d ;
                        fD = F[ i-1 ][ j-1 ] + s[ x ][ y ] ;
                        fL = F[ i ][ j-1 ] - d ;

                        F[ i ][ j ] = max( fU, fD, fL, &ptr ) ;

                        traceback[ i ][ j ] =  ptr ;
                }
        }
        i-- ; j-- ;

        while( i > 0 || j > 0 )
        {
                switch( traceback[ i ][ j ] )
                {
                        case '|' :      seq_1_al += '-' ;
                                        seq_2_al += seq_2[ i-1 ] ;
                                        i-- ;
                                        break ;

                        case '\\':      seq_1_al += seq_1[ j-1 ] ;
                                        seq_2_al += seq_2[ i-1 ] ;
                                        i-- ;  j-- ;
                                        break ;

                        case '-' :      seq_1_al += seq_1[ j-1 ] ;
                                        seq_2_al += '-' ;
                                        j-- ;
                }
                k++ ;
        }

        reverse( seq_1_al.begin(), seq_1_al.end() );
        reverse( seq_2_al.begin(), seq_2_al.end() );

        return  0 ;
}

void  alignment::print_matrix( int ** F, string seq_1, string seq_2 )
{
        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << "   ";
        }
        cout << "\n  ";

        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout.width( 3 );
                        cout << F[ i ][ j ] << " ";
                }
                cout << endl;
        }
}


void  alignment::print_traceback( char ** traceback, string seq_1, string seq_2 )
{
        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

        cout << "    ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << " ";
        }
        cout << "\n  ";

        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout << traceback[ i ][ j ] << " ";
                }
                cout << endl;
        }
}


void  alignment::print_al( string& seq_1_al, string& seq_2_al )
{
        cout << seq_1_al << endl;
        cout << seq_2_al << endl;
}
NW_Alignment::NW_Alignment() {
    init();
}
void NW_Alignment::init()
{
    matchScore=1;
    mismatchPenalty=1;
    gapPenalty=5;
}

double NW_Alignment::max(double x, double y)
{
    return x > y ? x : y;
}

double NW_Alignment::max(double x, double y, double z)
{
    return x > y ? max(x, z) : max(y, z);
}
double NW_Alignment::get_similarity_per(string s1,string s2){
  double fullScore= alignment(s1,s1);
  double simScore=alignment(s1,s2);
  return (simScore/fullScore)*100;
}
//size of two string should be same to avoid doing alignment in the first try
int NW_Alignment::findDirectSim(string const &a, string const &b) {
        int d = 0;
        if (a.length()!=b.length())
                return 0;
        if (a==b)
                return a.length();
        for (unsigned int i = 0; i < a.length() && i < b.length(); i++) {
                if (a[i] == b[i] || a[i] == 'N' || b[i] == 'N') {
                        d++;
                }
        }
        return d;

}
double NW_Alignment::get_similarity_perEnhanced(string s1,string s2){
        if (s1==s2)
                return 100;
        if (abs(s1.length()-s2.length())>maxGap )
                return 0;
        double fullScore=100;
        double simScore=0;
        fullScore= enhancedAlignment(s1,s1);
        simScore=enhancedAlignment(s1,s2);
        if (simScore<0)
                return 0;

        return (simScore/fullScore)*100;
}
void NW_Alignment::traceback(string& s1,string& s2, char **traceback ){
        string news1="";
        string news2="";
        int i=s1.length();
        int j=s2.length();
        while( i > 0 || j > 0 )
        {
                if (i>0 && j>0 && traceback[i][j]=='\\') {
                        news1 += s1[ i-1 ] ;
                        news2 += s2[ j-1 ] ;
                        i-- ;
                        j-- ;

                } else {
                        if (i>0 && traceback[i][j]=='|') {
                                news2 += '-' ;
                                news1 += s1[ i-1 ] ;
                                i-- ;

                        } else {
                                if (j>0 ) {
                                        news2 += s2[ j-1 ] ;
                                        news1 += '-' ;
                                        j-- ;

                                } else {
                                        if (i>0) {
                                                news2 += '-' ;
                                                news1 += s1[ i-1 ] ;
                                                i-- ;
                                        }
                                }
                        }
                }

        }

        reverse( news1.begin(), news1.end() );
        reverse( news2.begin(), news2.end() );
        s1=news1;
        s2=news2;

}
double NW_Alignment::alignment(string &s1, string &s2)
{
        int n = s1.length() + 1, m = s2.length() + 1, i, j;

        //int s[n][m];

        //char traceback[n][m];
        char **tracebackArr=new char*[n];
        for(int i = 0; i < n; i++)
        {
                tracebackArr[i] = new char[m];
        }

        int **s = new int*[n];
        for(int i = 0; i < n; i++)
        {
                s[i] = new int[m];
        }
        for (int i=0; i<n; i++) {
                s[i][0]=i*-1;
        }
        for (int j = 0; j < m; j++)
        {
                s[0][ j] = j*-1;
        }
        for (int i = 1; i <= n-1; i++)
        {
                for (int j = 1; j <= m-1; j++) // for (int j = jIndexMin; j <= i+(d); j++)//
                {

                        int scroeDiag = 0;

                        if(s1[i-1]==s2[j-1]|| s1[i-1]=='N'||s2[i-1]=='N')
                                //if (s1.substr(i - 1, 1) == s2.substr(j - 1, 1))
                                scroeDiag = s[i - 1][ j - 1] + matchScore;   //match
                                else
                                        scroeDiag = s[i - 1][ j - 1]  -mismatchPenalty; //substitution
                                        int scroeLeft = s[i][ j - 1] - gapPenalty; //insert
                                        int scroeUp = s[i - 1][ j] - gapPenalty;  //delete
                                        int maxScore =  max(scroeDiag,scroeLeft,scroeUp);//  Math.Max(Math.Max(scroeDiag, scroeLeft), scroeUp);
                                s[i][ j] = maxScore;
                        if (scroeDiag==maxScore) {
                                tracebackArr[i][j]='\\';
                        } else {
                                if (maxScore==scroeLeft) {
                                        tracebackArr[i][j]='-';
                                } else {
                                        if(maxScore==scroeUp) {
                                                tracebackArr[i][j]='|';
                                        }
                                }
                        }

                }
        }
        traceback(s1, s2,tracebackArr);
        int result=s[n-1][m-1];
        for(int i = 0; i < n; i++)
        {
                delete[] s[i];
        }
        delete[] s;

        for(int i = 0; i < n; i++)
        {
                delete[] tracebackArr[i];
        }
        delete[] tracebackArr;
        //return 0;
        return result;

}
double NW_Alignment::enhancedAlignment(string &s1, string &s2){
        if (s1==s2)
                return (s1.length()*matchScore);

        int n = s1.length() + 1, m = s2.length() + 1, i, j;
        if (abs(m-n)>maxGap || m<maxGap || n<maxGap)
            return 0;
        char **tracebackArr=new char*[n];
        for(int i = 0; i < n; i++)
        {
                tracebackArr[i] = new char[m];
        }

        int **s = new int*[n];
        for(int i = 0; i < n; i++)
        {
                s[i] = new int[m];
        }
        int p=0;
        for (int i=0; i<n; i++) {
                s[i][p]=i*-1;
                if (i>maxGap)
                p++;
        }
        p=0;
        for (int j = 0; j < m; j++)
        {
                s[p][ j] = j*-1;
                if (j>maxGap)
                p++;
        }

        for (int i = 1; i <= n-1; i++)
        {
                int jIndexMin=i-maxGap >0?i-maxGap:1;
                int jIndexMax=i+maxGap<m?i+maxGap:m-1;
            for (int j = jIndexMin; j <= jIndexMax; j++)//for (int j = 1; j <= m-1; j++) // for (int j = jIndexMin; j <= i+(d); j++)//
                {
                        int scroeDiag = 0;
                        if(s1[i-1]==s2[j-1]|| s1[i-1]=='N'||s2[i-1]=='N')
                                scroeDiag = s[i - 1][ j - 1] + matchScore;   //match
                        else
                                scroeDiag = s[i - 1][ j - 1]  -mismatchPenalty; //substitution
                        int scroeLeft = s[i][ j - 1] - gapPenalty; //insert
                        int scroeUp = s[i - 1][ j] - gapPenalty;  //delete
                        int maxScore =  max(scroeDiag,scroeLeft,scroeUp);//  Math.Max(Math.Max(scroeDiag, scroeLeft), scroeUp);
                        s[i][ j] = maxScore;
                        if (scroeDiag==maxScore) {
                                tracebackArr[i][j]='\\';
                        } else {
                                if (maxScore==scroeLeft) {
                                        tracebackArr[i][j]='-';
                                } else {
                                        if(maxScore==scroeUp) {
                                                tracebackArr[i][j]='|';
                                        }
                                }
                        }

                }
        }
        traceback(s1, s2,tracebackArr);
        int result=s[n-1][m-1];
        for(int i = 0; i < n; i++)
        {
                delete[] s[i];
        }
        delete[] s;
        for(int i = 0; i < n; i++)
        {
                delete[] tracebackArr[i];
        }
        delete[] tracebackArr;
        return result;
}
