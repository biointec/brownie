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

NW_Alignment::NW_Alignment():maxSize(100),maxGap(3) {
        matchScore=1;
        mismatchPenalty=1;
        gapPenalty=5;
        alocateMemory();
}

NW_Alignment::~NW_Alignment()
{

        deAlocateMemory();
}

void NW_Alignment::deAlocateMemory()
{
        for(int i = 0; i < maxSize; i++) {
                delete[] s[i];
        }
        delete[] s;
        for(int i = 0; i < maxSize; i++) {
                delete[] tracebackArr[i];
        }
        delete[] tracebackArr;}
void NW_Alignment::alocateMemory()
{
        tracebackArr=new char*[maxSize];
        for(int i = 0; i < maxSize; i++)
        {
                tracebackArr[i] = new char[maxSize];
        }
        s = new int*[maxSize];
        for(int i = 0; i < maxSize; i++)
        {
                s[i] = new int[maxSize];
        }
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
        double fullScore= enhancedAlignment(s1,s1);
        cout<<fullScore<<endl;
        double simScore=enhancedAlignment(s1,s2);
        if (s1==s2)
                return 100;
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


double NW_Alignment::enhancedAlignment(string &s1, string &s2){
       if (s1.length()>=maxSize || s2.length()>=maxSize){
                deAlocateMemory();
                maxSize=max(s1.length(),s2.length())+2;
                alocateMemory();
        }
        if (s1==s2)
                return s1.length()*matchScore;

        int n = s1.length() + 1, m = s2.length() + 1;
        int d=3;
        if (abs(m-n)>d || m<d || n<d)
                return 0;

        int p=0;
        for (int i=0; i<n; i++) {
                s[i][p]=i*-1;
                if (i>d)
                        p++;
        }
        p=0;
        for (int j = 0; j < m; j++)
        {
                s[p][ j] = j*-1;
                if (j>d)
                        p++;
        }
        for (int i = 1; i <= n-1; i++)
        {
                int jIndexMin=i-d >0?i-d:1;
                int jIndexMax=i+d<m?i+d:m-1;
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
                                if(maxScore==scroeUp) {
                                        tracebackArr[i][j]='|';
                                }
                                else {
                                        if (maxScore==scroeLeft) {
                                                tracebackArr[i][j]='-';
                                        }
                                }
                        }
                }
        }
        traceback(s1, s2,tracebackArr);
        int result=s[n-1][m-1];

        return result;
}




int NW_Alignment::findQualityDistance(string a, string b,string q) {
        int i=0;
        int j=0;
        int d=0;
        while(i<a.length()&&j<b.length()) {
                if (a[i]!=b[j] ) {
                        d=d+int(q[i]);
                }
                i++;
                j++;
        }
        return d;

}


/*
 * compare three string lines base per base if any char is different in any of these
 * lines put star in the mapped line
 *
 */
string NW_Alignment::getaligned(string a,string b,string c) {
        int i=0;
        string result="";
        while (i<a.length()) {
                if ((a[i]==b[i])&& (a[i]==c[i])) {
                        result=result+'*';
                }
                else {
                        result=result+a[i];
                }
                i++;
        }
        return result;
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


// ============================================================================
// ALIGNMENT CLASS
// ============================================================================

int AlignmentJan::align(const string& s1, const string& s2)
{
        // reallocate memory if necessary
        int thisMaxDim = max(s1.length(), s2.length());
        if (thisMaxDim > maxDim) {
                maxDim = thisMaxDim;
                delete [] M;
                M = new int[(2*maxIndel+1) * (maxDim+1)];
        }

        // initialize the borders of the matrix
        for (int i = 0; i <= maxIndel; i++)
                (*this)(i, 0) = i * gap;

        for (int i = 0; i <= maxIndel; i++)
                (*this)(0, i) = i * gap;

        // initialize the rest of the bulk of the matrix
        for (int i = 1; i <= (int)s1.length(); i++) {
                for (int j = max(1, i - maxIndel); j <= min((int)s2.length(), i + maxIndel); j++) {

                        bool hit = (s1[i-1] == s2[j-1]);
                        if ((s1[i-1] == 'N') || (s2[j-1] == 'N'))
                                hit = true;

                        int thisMatch = (*this)(i-1, j-1) + ((hit) ? match : mismatch);
                        int thisDel = (j < i + maxIndel) ? (*this)(i-1, j) + gap : thisMatch-1;
                        int thisIns = (j > i - maxIndel) ? (*this)(i, j-1) + gap : thisMatch-1;

                        int score = max(thisMatch, max(thisDel, thisIns));

                        (*this)(i, j) = score;
                }
        }

        return (*this)(s1.length(), s2.length());
}

AlignmentJan::AlignmentJan(int maxDim_, int maxIndel_, int match_,
                           int mismatch_, int gap_) : maxDim(maxDim_),
                           maxIndel(maxIndel_), match(match_),
                           mismatch(mismatch_), gap(gap_)
{
        M = new int[(2*maxIndel+1) * (maxDim+1)];
}

void AlignmentJan::printMatrix() const
{
        for (int l = 0; l < 2*maxIndel+1; l++) {
                for (int k = 0; k < maxDim+1; k++)
                        cout << M[k * (2 * maxIndel + 1) + l] << "\t";
                cout << endl;
        }
}

void AlignmentJan::printAlignment(const string& s1, const string& s2) const
{
        string al1, al2;

        int i = s1.size();
        int j = s2.size();
        while (i > 0 || j > 0) {
                bool hit = (s1[i-1] == s2[j-1]);
                if ((s1[i-1] == 'N') || (s2[j-1] == 'N'))
                        hit = true;

                if ((i > 0) && (j > 0) && ((*this)(i, j) == ((*this)(i-1, j-1) + ((hit) ? match : mismatch)))) {
                        al1.push_back(s1[i-1]);
                        al2.push_back(s2[j-1]);
                        i--;
                        j--;
                } else if ((j < i + maxIndel) && ((*this)(i, j) == (*this)(i-1, j) + gap)) {
                        al1.push_back(s1[i-1]);
                        al2.push_back('-');
                        i--;
                } else {
                        al1.push_back('-');
                        al2.push_back(s2[j-1]);
                        j--;
                }
        }

        reverse(al1.begin(), al1.end());
        reverse(al2.begin(), al2.end());

        cout << al1 << endl;
        for (int i = 0; i < max((int)al1.size(), (int)al2.size()); i++)
                if (al1[i] == al2[i] || al1[i] == 'N' || al2[i] == 'N')
                        cout << "|";
                else
                        cout << "*";
        cout << "\n" << al2 << endl;
}
