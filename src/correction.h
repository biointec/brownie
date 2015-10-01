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
        int numOfReads;
        int numOfAllReads;
        int numOfSupportedReads;
        
        ifstream readsFile;
        ofstream outFastq;
        
        void readInputReads(vector<readStructStr> &reads);
        void correctReads(vector<readStructStr> &reads);
        void writeOutputReads(vector<readStructStr> const &reads);
        void printProgress(clock_t const &begin);
        void correctRead(readStructStr &readInfo, int &supportedReads);
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
                            int maxError, string const &qualityProfile,
                            string &newRead);
        bool recursiveCompare(SSNode const &leftNode, string const &nodeContent,
                              int startOfNode, int startOfRead,
                              string &correctRead, string &erroneousRead,
                              string &guessedRead, string const &qualityProfile,
                              readCorrectionStatus &status);
        bool findBestMatch(vector<string> const &results, string &erroneousRead,
                           string &qualityProfile, bool rightDir,
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
        bool correctionByKmer(readCorrectionStatus &status, string &erroneousRead, string &guessedRead, string &correctRead, string &qualityProfile);
        bool correctionByMEM(readCorrectionStatus &status, string &erroneousRead, string &guessedRead, string &correctRead, string &qualityProfile);
        bool correctionByMEM(vector<match_t> &matches, string &reference, readCorrectionStatus &status, string &erroneousRead, string &guessedRead, string &correctRead, string &qualityProfile);
public:
        ReadCorrection(DBGraph &g, Settings &s);
        void errorCorrection(LibraryContainer &libraries);
        void CorrectErrorsInLibrary(ReadLibrary *input);
};

#endif
