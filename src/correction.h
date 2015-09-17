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
        void printProgress(clock_t const &begin, double numOfReads,
                           double numOfAllReads,
                           double numOfSupportedReads);
        void correctRead(readStructStr &readInfo, int &numOfSupportedReads);
        bool checkForAnswer(Kmer const &kmer, int startOfRead,
                            string &correctRead, string &erroneousRead,
                            string &guessedRead, string &qualityProfile,
                            readCorrectionStatus &status);
        bool recKmerCorrection(string &kmerStr, string const &qualityProfile,
                               int kmerStart, int round);
        int lowQualityPos(string quality, int startOfRead,
                          string const &kmer, int round);
        string applyINDchanges(string const &reference, string const &read);
        bool checkForIndels(string const &ref, string query,
                            int maxError, string &qualityProfile,
                            string &newRead);
        bool recursiveCompare(SSNode const &leftNode, string const &nodeContent,
                              int startOfNode, int startOfRead,
                              string &correctRead, string &erroneousRead,
                              string &guessedRead, string &qualityProfile,
                              readCorrectionStatus &status);
        bool findBestMatch(vector<string> &results, string erroneousRead,
                           string qualityProfile, bool rightDir,
                           string &bestrightMatch, int readLength);
        int findDifference(string const &guessedRead, string const &originalRead,
                           string const &qualityProfile, int startOfRead);
        int findDifference(string const &a, string const &b);
        void getAllRightSolutions(SSNode const &rootNode, string const &readPart,
                                  string const &qualityProfile,
                                  unsigned int depth,
                                  std::vector<std::string> &results);
        void getAllLeftSolutions(SSNode const &rootNode, string const &readPart,
                                 string const &qualityProfile,
                                 unsigned int depth,
                                 std::vector<std::string> & results);
        void getAllSolutions(SSNode const &rootNode, string const &readPart,
                             string const &qualityProfile, unsigned int depth,
                             std::vector<std::string> &results, bool forward);
public:
        ReadCorrection(DBGraph &g, Settings &s);
        void errorCorrection(LibraryContainer &libraries);
        void CorrectErrorsInLibrary(ReadLibrary *input);
};

#endif
