
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <cstring>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <omp.h>
#include <assert.h>
#include <iomanip>      // std::setprecision
#include <algorithm>
using namespace std;


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
    double get_similarity_per(string s1,string s2);
    void printMatrix(char **s, int m, int n);
    double enhancedAlignment(string &s1, string &s2);
};

NW_Alignment::NW_Alignment() {
    init();
}
void NW_Alignment::init()
{
    matchScore=1;
    mismatchPenalty=1;
    gapPenalty=20;
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

void NW_Alignment::printMatrix(char **s, int m, int n){
    for (int i=0;i<m;i++ ){
	for (int j=0;j<n;j++)
	    cout<< s[i][j]<<"      ";
	cout<<"\n";
    }
}

double NW_Alignment::enhancedAlignment(string &s1, string &s2){
	if (s1==s2)
	    return s1.length()*matchScore;
        int n = s1.length() + 1, m = s2.length() + 1, i, j;
	int d=3;
	if (abs(m-n)>d || m<d || n<d)
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

string erroneousReadsFileName;
string perfectReadsFileName;
string correctedReadsFileName;
string notCorrectedReadFileName;



int findDifference(string a, string b) {
    NW_Alignment al;
    al.enhancedAlignment(a,b);
    int i=0;
    int j=0;
    int d=0;
    while(i<a.length()&&j<b.length()) {
        if (a[i]!=b[j]) {
            d++;
        }
        i++;
        j++;
    }
    return d;

}
int findQualityDistance(string a, string b,string q) {
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
bool isEqual(string a, string b){
    int i=0;
    int j=0;
    int d=0;
    bool equal=true;
    while(i<a.length()&&j<b.length()) {
        if (a[i]!=b[j]&& a[i]!='N') {
	    equal=false;
        }
       i++;
       j++;
    }
    return equal;
}
string getaligned(string a,string b,string c) {
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


void validateCorrectionResult() {

    NW_Alignment a;
    ifstream perfectReadsStream;
    perfectReadsStream.open(perfectReadsFileName.c_str());
    ifstream correctedReadsStream;
    correctedReadsStream.open(correctedReadsFileName.c_str());
    ifstream  erroneousReadsFileStream;
    erroneousReadsFileStream.open(erroneousReadsFileName.c_str());

    ofstream notCorrectedReadF;
    notCorrectedReadF.open(notCorrectedReadFileName.c_str(), ios::out);

    ofstream testStream;
    string WorseName=notCorrectedReadFileName+"Worse.fastq";
    testStream.open(WorseName.c_str(), ios::out);

    
    ofstream notCorrectedFastq;
    string fastqName=notCorrectedReadFileName+".fastq";
    notCorrectedFastq.open(fastqName.c_str(), ios::out);

    double fullyRecoveredReads=0;
    double allErrorNum=0;
    double sumOfQualityDistances=0;
    double existErrorNum=0;
    long allReads=0;
    string erroneousRead, erroneousReadDesc, qualityProfile,erroneousReadName ;
    string perfectRead, perfectDesc;
    string correctedRead, correctedDesc;
    double truePositive=0;
    double trueNegative=0;
    double falsePositive=0;
    double falseNegative=0;

    int truePositiveFullRec=0;
    int trueNegativeFullRec=0;
    int falsePositiveFullRec=0;
    int falseNegativeFullRec=0;
    double maximumNumOfError=0;
    double maximumQualityDistance=0;
    double avgOfErrorsInRead=0;
    double sumOfAlignmentDist=0;
    double numOfAllChanges=0;
    int numOfReads=0;
    double numberOfAllErrors=0;
    int i=0;
    while (     getline(erroneousReadsFileStream, erroneousReadName)
                &&getline(perfectReadsStream, perfectDesc)&&  getline(correctedReadsStream, correctedDesc) ) {
        getline(erroneousReadsFileStream,erroneousRead);
        getline(perfectReadsStream,perfectRead);
        getline(correctedReadsStream,correctedRead);

        getline(erroneousReadsFileStream, erroneousReadDesc);
        getline(erroneousReadsFileStream, qualityProfile);

        getline(correctedReadsStream, correctedDesc);
        getline(correctedReadsStream, correctedDesc);
        numOfReads=numOfReads+1;
	if ( (int)numOfReads % 10000 == 0 ) {
                        cout <<" read Number:: " <<numOfReads <<"\r";
                        cout.flush();
        }	
        int allErrorINRead=findDifference(perfectRead,erroneousRead);;
        int exsistErrorInRead=0;
        int qualityDistance=0;
        int numOfChangesInRead=0;
        double initialQualityDistance=0;
        
	double numberOfErrorInRead=0;
        if(allErrorINRead>maximumNumOfError) {
            maximumNumOfError=allErrorINRead;
        }
        initialQualityDistance=findQualityDistance(perfectRead, erroneousRead,qualityProfile);
        if (initialQualityDistance>maximumQualityDistance) {
            maximumQualityDistance=initialQualityDistance;
        }
        if (!isEqual( perfectRead,correctedRead)) {
            double maxSim =a.enhancedAlignment(perfectRead,perfectRead);
            sumOfAlignmentDist +=maxSim - a.enhancedAlignment(perfectRead,correctedRead);
            string copyOfcorrectedRead=correctedRead;
	    string copyOfperfecRead=perfectRead;

            exsistErrorInRead=findDifference(perfectRead,correctedRead);
            a.enhancedAlignment(copyOfcorrectedRead,erroneousRead);
            //modifyQualityBasedOnAlignment(qualityProfile,erroneousRead);
            qualityDistance= findQualityDistance(erroneousRead,copyOfcorrectedRead,qualityProfile);
            numOfChangesInRead=findDifference(erroneousRead,copyOfcorrectedRead);
             
	    a.enhancedAlignment(erroneousRead,copyOfperfecRead);
	    numberOfErrorInRead=findDifference(erroneousRead,copyOfperfecRead);
	    numberOfAllErrors+=numberOfErrorInRead;
	    
            if (allErrorINRead==0) {
                falsePositiveFullRec++;
            } else {
                falseNegativeFullRec++;
            }

            notCorrectedReadF<<erroneousReadName <<"	|"<< "read number:	"<<i <<"	num of all Errors:	"<<allErrorINRead<<"	number of remaining Errors:	"<< exsistErrorInRead<<endl;
            string alignedError="";
            alignedError=getaligned(erroneousRead,perfectRead,correctedRead);
	    
	   

            notCorrectedReadF<<"E:"<<erroneousRead<<endl;
            notCorrectedReadF<<"P:"<<perfectRead<<endl;
            notCorrectedReadF<<"C:"<<correctedRead<<endl;
            notCorrectedReadF<<"A:"<<alignedError<<endl;
            if ( exsistErrorInRead>allErrorINRead) {
                testStream<<erroneousReadName <<"	|"<< "read number:	"<<i <<"	num of all Errors:	"<<allErrorINRead<<"	number of remaining Errors:	"<< exsistErrorInRead<<endl;
                testStream<<"E:"<<erroneousRead<<endl;
                testStream<<"P:"<<perfectRead<<endl;
                testStream<<"C:"<<correctedRead<<endl;
                testStream<<"A:"<<alignedError<<endl;
            }
            notCorrectedReadF<<"A:"<<qualityProfile<<endl;

            notCorrectedFastq<<erroneousReadName<<endl;
            notCorrectedFastq<<erroneousRead<<endl;
            notCorrectedFastq<<'+'+erroneousReadName<<endl;
            notCorrectedFastq<<qualityProfile<<endl;

        } else {
            fullyRecoveredReads++;
            if (allErrorINRead>0) {
                truePositiveFullRec++;
            } else {
                trueNegativeFullRec++;
            }
        }
        char correctChar;
        char modifiedChar;
        char erroneousChar;
        for (int i=0; i<correctedRead.length(); i++) {
            modifiedChar=correctedRead[i];
            erroneousChar=erroneousRead[i];
            correctChar=perfectRead[i];
            if (erroneousChar!=correctChar) {
                if (modifiedChar==correctChar) {
                    truePositive++;
                }
                else {
                    falseNegative++;
                }

            } else {
                if(modifiedChar==correctChar) {
                    trueNegative++;
                }
                else {
                    falsePositive++;
                }

            }
        }
        numOfAllChanges=numOfAllChanges+numOfChangesInRead;
        sumOfQualityDistances=sumOfQualityDistances+qualityDistance;
        existErrorNum=existErrorNum+exsistErrorInRead;
        allErrorNum=allErrorNum+allErrorINRead;
        allReads++;
	i++;
    }
    double gain=(double)(truePositive-falsePositive)/(double)(truePositive+falseNegative);
    double gainFR=(double)(truePositiveFullRec-falsePositiveFullRec)/(double)(truePositiveFullRec+falseNegativeFullRec);
    erroneousReadsFileStream.close();
    perfectReadsStream.close();
    correctedReadsStream.close();
    double corrected=allErrorNum- existErrorNum;
    double correctionAVG=((double) truePositive/(double)(truePositive+falseNegative))*100; 
    std::cout << std::fixed;
    //std::cout << std::setprecision(2) << f << '\n'

    cout <<endl<< "<<<Report for reads>>>" << endl;
    cout << "----------------------------------------------------\n" << endl;
    cout<<"Number of reads for comparison: "<<(int)numOfReads<<endl;
    cout<<"Maximum number Of Errors in one read is: "<< (int)maximumNumOfError <<", and maximum quality distances is: "<<(int)maximumQualityDistance <<endl;

    
    cout <<endl<< "<<<Alignment based report>>>" << endl;
    cout << "----------------------------------------------------" << endl;
    cout<<"Sum of  errors in the Data Set: "<<(int)numberOfAllErrors<<", (" <<(double)numberOfAllErrors/(double)numOfReads<<" per read)"<<endl;
    cout<<"Sum of alignment distances between corrected and perfect read is: "<<(int)sumOfAlignmentDist<<", (" <<(double)sumOfAlignmentDist/(double)numOfReads<<" per read)"<<endl;
    
    cout <<endl<< "<<<The evaluation report based on base pairs>>>"<<endl;
    cout << "----------------------------------------------------" << endl;
    cout<<"Among "<<(int)(truePositive+falseNegative)<<" of Errors "<<(int)truePositive<<" number of them are corrected which is: ("<<std::setprecision(2)<<correctionAVG<<"%)"<<endl ;
    cout<<"	TP:"<<(int)truePositive<<"	TN:"<<(int)trueNegative<<"	FP:"<<(int)falsePositive <<"	FN:"<<(int)falseNegative<<endl;
    cout <<"    The Gain value percentage is: ("<<std::setprecision(2)<<gain*100<<"%)" <<endl;
    
    cout<<endl<<"<<<The evaluation report based full read recovery>>>"<<endl;
    cout << "----------------------------------------------------" << endl;
    cout<<" Among "<<(int)allReads<<" of reads "<<(int)fullyRecoveredReads<<", number of them are fully recovered which is: ("<<((double)fullyRecoveredReads/(double)allReads)*100<<")%"<<endl ;
    cout<<"	TP: "<<(int)truePositiveFullRec<<"	TN: "<<(int)trueNegativeFullRec<<"	FP: "<<(int)falsePositiveFullRec<<"	FN: "<<(int)falseNegativeFullRec<<endl;
    cout <<"    The Gain value percentage for full recovery of reads is: ("<<gainFR*100<<"%)" <<endl;

    cout<<endl<< "<<<Quality based reports>>>"<<endl;
    cout << "----------------------------------------------------" << endl;
    cout<<"The sum of quality distances between original reads and corrected reads is: "<<(int)sumOfQualityDistances<<" ,("<<double( sumOfQualityDistances)/(double)(numOfReads)<<"   per read) ,and the number of changes is:   "<<numOfAllChanges<<", (" <<(double)numOfAllChanges/(double) numOfReads <<" per read)" <<endl;
    
}
int main(int argc, char ** args)
{
    cout << "Welcome to comparison software\n" << endl;
    cout << "----------------------------------------------------\n" << endl;
    if (argc < 4) {
        cout << "ERROR: Wrong number of arguments!\n";
        cout << "Usage: ./comparison <erroneousReads.fq> <perfectReads.fa> <correctedReads.fq> (opt: notCorrectedReadFileName)" << endl << flush;
        cout << "Order of reads in all input files should be the same.\nPlease look at the readme file for more information." << endl << flush;
        return 0;
    }

    erroneousReadsFileName=args[1];
    perfectReadsFileName=args[2];
    correctedReadsFileName=args[3];
    if (argc==5) {
        notCorrectedReadFileName=args[4];
    }
    validateCorrectionResult();
    cout << "\n----------------------------------------------------\n" << endl;
    cout << "Exiting... bye!" << endl << endl;
    return 0;
}
