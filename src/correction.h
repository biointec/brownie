#ifndef CORRECTION_H
#define CORRECTION_H

#include <string>
#include "graph.h"

using namespace std;

struct readStructStr
{
        int intID;
        string strID;
        string corrctReadContent;
        string erroneousReadContent;
        string qualityProfile ;
        string orientation;

};

class ReadCorrection 
{
private:
        DBGraph &dbg;
        ReadLibrary *library;
        Settings &settings;
        int kmerSize;
        ifstream readsFile;
        ofstream outFastq;
        vector<string> references;
        vector<sparseSA*> saVec;
public:
        enum readCorrectionStatus {

                kmerfound,
                fullHealing,
                parHealing,
                kmerNotfound,
                graphIsMissing,
                anotherKmer

        };
        ReadCorrection(DBGraph &g, Settings &s);
        void CorrectErrorsInLibrary(ReadLibrary *input);
        void readInputReads(vector<readStructStr> &reads, double &numOfReads);
        void correctRead(readStructStr &readInfo, int &numOfSupportedReads);
        void errorCorrection(LibraryContainer &libraries);
        bool cheakForAnswer(Kmer kmer, int startOfRead, string & correctRead,
                            string & erroneousRead, string & guessedRead,
                            string &qualityProfile,
                            readCorrectionStatus &status);
        bool recKmerCorrection(string &kmerStr, const string & qualityProfile,
                               int kmerStart, int round);
        int lowQualityPos(string quality, int startOfRead, string kmer,
                          int round);
        string applyINDchanges(string reference, string read);
        bool checkForIndels(string ref, string query, const int maxError, string& qualityProfile, string& newRead);
        bool recursiveCompare(SSNode leftNode,string nodeContent,int  startOfNode,int startOfRead, string & correctRead, string & erroneousRead,string & guessedRead , string &qualityProfile,readCorrectionStatus &status);
        bool findBestMatch(vector<string>& results, string erroneousRead, string qualityProfile,bool rightDir , string& bestrightMatch, SSNode &leftNode,  int readLength);
        int findDifference(string guessedRead, string originalRead, string &qualityProfile, int startOfRead);
        int findDifference(string a, string b);
        void getAllRightSolutions(SSNode rootNode, string readPart,string qualityProfile ,unsigned int depth, std::vector<std::string> & results);
        void getAllLeftSolutions(SSNode rootNode,string readPart,string qualityProfile , unsigned int depth, std::vector<std::string> & results);
};

#endif
