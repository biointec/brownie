#ifndef CORRECTION_H
#define CORRECTION_H

#include <string>
#include "graph.h"
#include "dijkstra.h"
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

        const float maxErrorRate=.3;
        const size_t avgQualityError=40;
        DBGraph &dbg;
        ReadLibrary *library;
        Settings &settings;
        size_t kmerSize;

        int numOfReads;
        int numOfAllReads;
        int numOfSupportedReads;
        double minSimPer;
        double maxTimePerRead;
        vector<string> references;
        vector<sparseSA*> saVec;
        NW_Alignment Nw;

        ifstream readsFile;
        ofstream outFastq;

        Dijkstra dijk;

        void readInputReads(vector<readStructStr> &reads);
        void correctReads(vector<readStructStr> &reads);
        void writeOutputReads(vector<readStructStr> const &reads);
        void printProgress(clock_t const &begin);
        bool correctRead(readStructStr &readInfo);
        bool  firstHitErrorCorrection(readStructStr &readInfo);
        bool  makeBridgeErrorCorrection(readStructStr &readInfo);
        bool checkForAnswer(Kmer const &kmer, int startOfRead,
                            string const &original, string &guess,
                            string const &qProfile,
                            readCorrectionStatus &status);
        bool recKmerCorrection(string &kmerStr, string const &qProfile,
                               int kmerStart, int round);
        int lowQualityPos(string quality, int startOfRead,
                          string const &kmer, int round);
        string applyINDchanges(string const &reference, string const &read);
        bool recursiveCompare(SSNode &leftNode, string const &nodeContent,
                              int startOfNode, int startOfRead, string const &original,
                              string &guess, string const &qProfile,
                              readCorrectionStatus &status);
        bool expand(SSNode &node, string const &original, string const &qProfile, string &guess, pair<int, int> bounds, bool forward, readCorrectionStatus &status);
        bool findBestMatch(vector<string> const &results, string &original,
                           bool rightDir, string &bestMatch);
        int findDifference(string const &a, string const &b);
        vector<string> getAllSolutions(SSNode  &rootNode, string const &readPart, bool forward);
        int findRecSolutionsForward(vector<string> &results, SSNode const rootNode,
                                      string const &readPart,
                                      string currentPath,clock_t& start);
        int findRecSolutionsBackward(vector<string> &results, SSNode const rootNode,
                                      string const &readPart,
                                      string currentPath,clock_t& start);
        bool findRecSolutionsRec(vector<string> &results, SSNode const rootNode,
                                      string const &readPart,
                                      string currentPath,bool forward);
        bool correctionByKmer(readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
        bool findKmer(readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
        bool findSimilarKmer(readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
        bool findSimilarKmer2(readStructStr &readInfo,readCorrectionStatus &status);
        void updateExtremeValues (readStructStr &readInfo,NodePosPair& result,int& curReadLeftExtreme,
                                          int &curReadRightExtreme , int& startOfRead, int &nodeLeftExtreme ,int &nodeRightExtreme  );
        void checkReadSize(string &guess,readStructStr &readInfo);
        bool correctionByMEM(readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
        bool correctionByMEM(vector<match_t> &matches, string const &reference, readCorrectionStatus &status, string const &original, string &guess, string const &qProfile);
        int shiftSize(match_t const &m, string const &reference);
        bool seedIsContainedInKmer(match_t const &m, string const &reference, int shift);
        bool correctRead(readStructStr &readInfo,readCorrectionStatus &status);
        bool newCorrectReadByKmer(readStructStr &readInfo,readCorrectionStatus &status, int  startOfRead , Kmer& kmer, NodePosPair &result);
        bool findMiddlePart2( NodePosPair& result,int& startOfRead,int& curReadLeftExtreme, int &curReadRightExtreme , int &nodeLeftExtreme,
                                      int &nodeRightExtreme, readStructStr &readInfo,string& guess );
        void findMiddlePart(SSNode &currNode, SSNode &prevNode,SSNode firstNode , NodePosPair& result,
                                    int& startOfRead,int& curReadLeftExtreme, int &curReadRightExtreme , int &prevReadRightExtreme
                                    ,int &maxRightInRead,int &minLeftInRead,readStructStr& readInfo, string& guess );
        bool fillGap(readStructStr &readInfo,NodePosPair& result,int& curReadLeftExtreme,
                             int &curReadRightExtreme , int& startOfRead, int &prevReadRightExtreme,
                             SSNode &prevNode, SSNode &currNode ,SSNode firstNode ,string &guess, int maxRightInRead, int minLeftInRead,
                     int &nodeLeftExtreme,int &nodeRightExtreme);
        void findBridge(vector<string> &results  ,SSNode startNode,SSNode &endNode,string& readPart,string currentPath);
        string checkForIndels(string const &ref, string query);
public:
        ReadCorrection(DBGraph &g, Settings &s);
        void errorCorrection(LibraryContainer &libraries);
        void CorrectErrorsInLibrary(ReadLibrary *input);
};

#endif
