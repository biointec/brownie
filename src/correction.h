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
        string originalContent;
        string qProfile ;
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
        bool correctRead(readStructStr &readInfo);
        bool checkForAnswer(Kmer const &kmer, int startOfRead,
                            string const &original, string &guess,
                            string const &qProfile,
                            readCorrectionStatus &status);
        bool recKmerCorrection(string &kmerStr, string const &qProfile,
                               int kmerStart, int round);
        int lowQualityPos(string quality, int startOfRead,
                          string const &kmer, int round);
        string applyINDchanges(string const &reference, string const &read);
        bool checkForIndels(string const &ref, string query,
                            int maxError, string const &qProfile,
                            string &newRead);
        bool recursiveCompare(SSNode const &leftNode, string const &nodeContent,
                              int startOfNode, int startOfRead, string const &original,
                              string &guess, string const &qProfile,
                              readCorrectionStatus &status);
        bool expand(SSNode const &node, string const &original, string const &qProfile, string &guess, pair<int, int> bounds, bool forward, readCorrectionStatus &status);
        bool findBestMatch(vector<string> const &results, string &original,
                           bool rightDir, string &bestMatch, int readLength);
        
        int findDifference(string const &guess, string const &original,
                           string const &qProfile, int start);
        int findDifference(string const &a, string const &b);
        vector<string> getAllSolutions(SSNode const &rootNode, string const &readPart,
                             string const &qProfile, bool forward);
        
        bool correctionByKmer(readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
        bool findKmer(readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
        bool findSimilarKmer(readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
        
        bool correctionByMEM(readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
        bool correctionByMEM(vector<match_t> &matches, string &reference, readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
public:
        ReadCorrection(DBGraph &g, Settings &s);
        void errorCorrection(LibraryContainer &libraries);
        void CorrectErrorsInLibrary(ReadLibrary *input);
};

#endif
