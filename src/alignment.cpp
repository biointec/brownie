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
                for (int j = 1; j <= m-1; j++)
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






















