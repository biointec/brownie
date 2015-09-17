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
        enum readCorrectionStatus {

                kmerfound,
                fullHealing,
                parHealing,
                kmerNotfound,
                graphIsMissing,
                anotherKmer

        };
        DBGraph &dbg;
        ReadLibrary *library;
        Settings &settings;
        int kmerSize;
        vector<string> references;
        vector<sparseSA*> saVec;
        
        ifstream readsFile;
        ofstream outFastq;
        
        void readInputReads(vector<readStructStr> &reads, double &numOfReads);
        void correctReads(vector<readStructStr> &reads,
                          int &numOfSupportedReads);
        void writeOutputReads(vector<readStructStr> const &reads);
        void printProgress(clock_t const &begin, double const &numOfReads,
                           double const &numOfAllReads,
                           double const &numOfSupportedReads);
        void correctRead(readStructStr &readInfo, int &numOfSupportedReads);
        bool checkForAnswer(Kmer kmer, int startOfRead, string & correctRead,
                            string & erroneousRead, string & guessedRead,
                            string &qualityProfile,
                            readCorrectionStatus &status);
        bool recKmerCorrection(string &kmerStr, const string & qualityProfile,
                               int kmerStart, int round);
        int lowQualityPos(string quality, int startOfRead, string kmer,
                          int round);
        string applyINDchanges(string reference, string read);
        bool checkForIndels(string ref, string query, int const maxError,
                            string& qualityProfile, string& newRead);
        bool recursiveCompare(SSNode leftNode,string nodeContent,
                              int startOfNode,int startOfRead,
                              string &correctRead, string &erroneousRead,
                              string &guessedRead, string &qualityProfile,
                              readCorrectionStatus &status);
        bool findBestMatch(vector<string>& results, string erroneousRead,
                           string qualityProfile, bool rightDir,
                           string &bestrightMatch, SSNode &leftNode,
                           int readLength);
        int findDifference(string guessedRead, string originalRead,
                           string &qualityProfile, int startOfRead);
        int findDifference(string a, string b);
        void getAllRightSolutions(SSNode rootNode, string readPart,
                                  string qualityProfile, unsigned int depth,
                                  std::vector<std::string> &results);
        void getAllLeftSolutions(SSNode rootNode, string readPart,
                                 string qualityProfile, unsigned int depth,
                                 std::vector<std::string> & results);
public:
        ReadCorrection(DBGraph &g, Settings &s);
        void errorCorrection(LibraryContainer &libraries);
        void CorrectErrorsInLibrary(ReadLibrary *input);
};

#endif
