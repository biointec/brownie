/***************************************************************************
 *   Copyright (C) 2014 - 2016 Jan Fostier (jan.fostier@intec.ugent.be)    *
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
bool AlignmentJan::detectLocalHammering(const string& s1, const string& s2)const{

       if (s1==s2)
               return false;
       string al1, al2;
       int i = s1.size();
       int j = s2.size();
       size_t numOfErrors=0;
       vector<int> errorPosList;
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

        //find the location of erros in the read and save it in a list
        for (int i = 0; i < max((int)al1.size(), (int)al2.size()); i++){
                 if (al1[i] != al2[i] && al1[i] != 'N' && al2[i] != 'N'){
                         numOfErrors++;
                         errorPosList.push_back(i);
                }
        }

        int currentBinStart = 0;
        //if the error rate in the read is less than 4% its normal
        if ( numOfErrors < 4*s1.length()/100 )
                return false;

        double patternLessNess = 1;
        double patternLessNessCutOff = 1000000000;
        //calculate patternLessNess value which shows the degree in which errors are spread with the same distance from each other in the read.
        //This value is 1 if the distance between errors is exactly the same. The value increases if the errors are locally located in some part of the reads.

        for ( size_t i = 1 ;i <errorPosList.size(); i++ ){
                //cout<<"double(1/(errorPosList[i]-errorPosList[i-1])) " <<(1/double(errorPosList[i]-errorPosList[i-1])) <<endl ;
                //cout<<"double(max((int)al1.size(), (int)al2.size())/(double)numOfErrors)" <<double(max((int)al1.size(), (int)al2.size())/(double)numOfErrors)<<endl;
                patternLessNess = patternLessNess* double(1/double(errorPosList[i]-errorPosList[i-1])) * double(max((int)al1.size(), (int)al2.size())/(double)numOfErrors);
        }
        int binSize = 20;
        int NumbOfbin = max((int)al1.size(), (int)al2.size()) / binSize ;
        float expectedErrorInBin = (float)numOfErrors/ (float) NumbOfbin;
        int numOfErrorIncurrBin = 0;
        i = 0 ;
        while ( currentBinStart + binSize < max((int)al1.size(), (int)al2.size())){
                numOfErrorIncurrBin = 0;
                for (size_t i =0; i <errorPosList.size();i++){
                        if ( errorPosList[i]> currentBinStart && errorPosList[i] < currentBinStart + binSize ){
                                numOfErrorIncurrBin++;
                                if ( numOfErrorIncurrBin  >expectedErrorInBin *7 && patternLessNess > patternLessNessCutOff ){
                                        cout<<"From "<<numOfErrors<< " erros In read " <<numOfErrorIncurrBin << " of Errors are in bin starts from " << currentBinStart
                                        <<" to "<< currentBinStart + binSize<<endl;
                                        printAlignment(s1,s2);
                                        cout << "patternLessNess: " << patternLessNess <<endl;
                                        return true;
                                }
                        }
                        if ( errorPosList[i] > currentBinStart + binSize)
                                break;
                }
                //to have the overlap bins we go 1/2 of bin size further in each step
                currentBinStart = currentBinStart + binSize /2 ;
        }
        return false;
}
