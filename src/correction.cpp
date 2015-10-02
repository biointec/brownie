#include "alignment.h"
#include "graph.h"
#include "kmernode.h"
#include "settings.h"
#include <list>
#include <queue>
#include <string.h>
#include <stack>
#include "library.h"

#include "correction.h"
#include "essamem.h"

/**
 * ctor
 */
ReadCorrection::ReadCorrection(DBGraph &g, Settings &s) :
        dbg(g), library(NULL), settings(s), kmerSize(s.getK()),
        numOfSupportedReads(0), numOfAllReads(0), numOfReads(0)
{
        cout << endl << "welcome to error correction part" << endl;
        //*******************************************************
        cout << "populateTable" << endl;
        dbg.populateTable();
        //*******************************************************
        cout << "making reference for essaMEM" << endl;
        references = makeRefForessaMEM(dbg);
        cout << "using essaMEM" << endl;
        string meta = "noise";
        for (unsigned int i = 0; i < references.size(); i++) {
                sparseSA* sa1 = init_essaMEM(references[i], meta);
                saVec.push_back(sa1);
        }
}

/**
 * read the reads from the inputfile
 */
void ReadCorrection::readInputReads(vector<readStructStr> &reads) {
        int batchSize = OUTPUT_FREQUENCY;
        for (int i = 0; i < batchSize; ++i) {
                readStructStr readInfo;
                if (getline(readsFile, readInfo.strID)
                                && getline(readsFile, readInfo.erroneousReadContent)) {
                        if (FASTQ == library->getFileType()) {
                                getline(readsFile,readInfo.orientation );
                                getline(readsFile,readInfo.qualityProfile );
                        } else {
                                readInfo.orientation='+';
                                int qvLength = readInfo.qualityProfile.length();
                                for (unsigned int i = 0; i < qvLength; i++) {
                                        readInfo.qualityProfile[i] = '#';
                                }
                        }
                        numOfReads++;
                        readInfo.strID.erase(0, 1);
                        readInfo.strID = "@" + readInfo.strID;
                        readInfo.intID = numOfReads;
                        reads.push_back(readInfo);
                } else {
                        readsFile.close();
                }
        }
}

/**
 * Correct a batch of reads
 */
void ReadCorrection::correctReads(vector<readStructStr> &reads) {
        int supportedReads = 0;
#pragma omp parallel for reduction(+:supportedReads)
        for (int i = 0; i < reads.size(); ++i) {
                correctRead(reads[i], supportedReads);
        }
        numOfSupportedReads += supportedReads;
}

/**
 * write the reads to the outputfile
 */
void ReadCorrection::writeOutputReads(vector<readStructStr> const &reads) {
        for (int i = 0; i < reads.size(); ++i) {
                readStructStr readInfo = reads[i];
                //write ID
                outFastq << readInfo.strID << endl;
                //write correction
                outFastq << readInfo.corrctReadContent << endl;
                //write orientation
                outFastq << readInfo.orientation << endl;
                //write qualityProfile
                outFastq << readInfo.qualityProfile << endl;
        }
}

/**
 * Try to find a kmer from the read in the graph and correct the read
 * @return true if correction is found
 */
bool ReadCorrection::findKmer(readCorrectionStatus &status, string const &erroneousRead,
                string &guessedRead, string &qualityProfile) {
        bool found = false;
        TString read = erroneousRead;
        int startOfRead = 0;
        for (TStringIt it = read.begin(); it != read.end(); it++ ) {
                Kmer kmer = *it;
                if (!checkForAnswer(kmer, startOfRead, erroneousRead, guessedRead, qualityProfile, status)) {
                        startOfRead++;
                        if (status == anotherKmer
                                        || status == kmerNotfound) {
                                continue;
                        }
                } else {
                        found = true;
                        break;
                }
        }
        return found;
}

/**
 * Try to find a kmer from the read that is similar to a kmer in the graph and correct the read
 * @return true if correction is found
 */
bool ReadCorrection::findSimilarKmer(readCorrectionStatus &status,
                string const &erroneousRead, string &guessedRead, string &qualityProfile) {
        bool found = false;
        int readLength = erroneousRead.length();
        for (int kmerStart = 0; kmerStart + kmerSize < readLength
                        && kmerSize <= readLength; kmerStart += 5) {
                string tempstr = erroneousRead.substr(kmerStart, kmerSize);
                if (!dbg.kmerExistsInGraph(tempstr)) {
                        if(recKmerCorrection(tempstr, qualityProfile, kmerStart, 1)) {
                                Kmer kmer = tempstr;
                                if (checkForAnswer(kmer, kmerStart, erroneousRead, guessedRead, qualityProfile, status)) {
                                        found = true;
                                        break;
                                }
                        }
                }
        }
        return found;
}

/**
 * Try to correct read using Kmers
 * @return true if correction is found
 */
bool ReadCorrection::correctionByKmer(readCorrectionStatus &status, string const &erroneousRead, string &guessedRead, string &correctRead, string &qualityProfile) {
        bool found = false;
        if (erroneousRead.length() >= kmerSize) {
                found = findKmer(status, erroneousRead, guessedRead, qualityProfile);
        }
        if (status == kmerNotfound) {
                found = findSimilarKmer(status, erroneousRead, guessedRead, qualityProfile);
        }
        return found;
}

/**
 * Try to correct read using MEMs
 * @return true if correction is found
 */
bool ReadCorrection::correctionByMEM(vector<match_t> &matches, string &reference, readCorrectionStatus &status, string &erroneousRead, string &guessedRead, string &correctRead, string &qualityProfile) {
        if (matches.size() > 0) {
                match_t bestMatch;
                string bestRefstr = "";
                int bestShiftLeft = 0;
                unsigned int cutLength = 0;
                for (unsigned int i = 0; i < matches.size(); i++) {
                        match_t m = matches[i];
                        unsigned int shiftLeft = kmerSize - m.len < m.query ? kmerSize - m.len : m.query;
                        int j = m.ref;
                        unsigned int k = 0;
                        int lastValue = 0;
                        if (m.len == kmerSize) {
                                continue; //we already in previous step checked this kmer
                        }
                        while(j > 0 && reference[j-1] != '>' && k < shiftLeft) {
                                lastValue = k;
                                k++;
                                j--;
                        }
                        shiftLeft = k;
                        lastValue = 0;
                        j = m.ref-shiftLeft;
                        k = 0;
                        while(reference[j] != '<'
                                        && k <= kmerSize
                                        && (reference[j] == 'A'
                                                || reference[j] == 'T'
                                                || reference[j] == 'C'
                                                || reference[j] == 'G'
                                                || reference[j] == 'a'
                                                || reference[j] == 't'
                                                || reference[j] == 'c'
                                                || reference[j] == 'g')) {
                                lastValue = k;
                                j++;
                                k++;
                        }
                        cutLength = lastValue;
                        string refstr = reference.substr(m.ref - shiftLeft, cutLength);
                        string foundKmerStr = erroneousRead.substr(m.query - shiftLeft,cutLength);
                        float d=findDifference(refstr, foundKmerStr);
                        if (refstr.length()/(d+1) > 6 &&  cutLength==kmerSize && foundKmerStr.length()==kmerSize && refstr.length()==kmerSize) {
                                bestMatch=m;
                                bestRefstr=refstr;
                                bestShiftLeft=shiftLeft;
                                Kmer kmer = bestRefstr;
                                int startOfRead = bestMatch.query - bestShiftLeft;
                                if (checkForAnswer(kmer, startOfRead, erroneousRead, guessedRead, qualityProfile, status)) {
                                        return true;
                                }
                        }
                }
        }
        return false;
}

/**
 * for every reference sequence, find MEMs and process them
 * @return true if correcction is found
 */
bool ReadCorrection::correctionByMEM(readCorrectionStatus &status, string &erroneousRead, string &guessedRead, string &correctRead, string &qualityProfile) {
        for(unsigned int r = 0; r < saVec.size(); ++r) {
                int min_len = kmerSize / 3;
                bool print = 0; //not sure what it prints if set to 1
                long memCounter = 0; //this does not really work
                string query = erroneousRead;
                vector<match_t> matches; //will contain the matches
                saVec[r]->MEM(query, matches, min_len, print, memCounter, true, 1);
                sort(matches.begin(), matches.end(), sortMemResultbySize);
                if(correctionByMEM(matches, references[r], status, erroneousRead, guessedRead, correctRead, qualityProfile)) {
                        return true;
                }
        }
        return false;
}

/**
 * correct a single read
 */
void ReadCorrection::correctRead(readStructStr &readInfo, int &supportedReads) {
        unsigned int readLength=readInfo.erroneousReadContent.length();
        if (readLength < 50) {
                //TODO GILES > should short reads be skipped?
                int stop = 0;
                stop++;
        }
        readCorrectionStatus status = kmerNotfound;
        string erroneousRead =  readInfo.erroneousReadContent;
        string guessedRead = erroneousRead;
        string correctRead = readInfo.corrctReadContent;
        string qualityProfile = readInfo.qualityProfile;
        //find in regular way
        bool found = correctionByKmer(status, erroneousRead, guessedRead, correctRead, qualityProfile);
        //find with essaMEM
        if (status == kmerNotfound) {
                found = correctionByMEM(status, erroneousRead, guessedRead, correctRead, qualityProfile);
        }
        //all attempt for error correction finished
        //now write the guessed Read into the file
        if (found) {
                supportedReads += 1;
        }
        if (guessedRead.length() > readLength) {
                guessedRead = guessedRead.substr(0, readLength);
        }
        if (guessedRead.length() < readLength) {
                cout << readInfo.intID << endl;
                guessedRead = readInfo.erroneousReadContent;
        }
        if (guessedRead.length() != readInfo.qualityProfile.length()) {
                if (guessedRead.length() > readInfo.qualityProfile.length()) {
                        guessedRead = guessedRead.substr(0, readInfo.qualityProfile.length());
                }
                if (guessedRead.length() < readInfo.qualityProfile.length()) {
                        readInfo.qualityProfile = readInfo.qualityProfile.substr(0, guessedRead.length());
                }
        }
        //save the correction so we can write it to file later
        readInfo.corrctReadContent = guessedRead;
}

/**
 * estimate the remaining time
 */
void ReadCorrection::printProgress(clock_t const &begin) {
        clock_t endt=clock();
        double passedTime=double(endt-begin)/(double) (CLOCKS_PER_SEC*60);
        double progress = ((double)numOfReads/(double)numOfAllReads)*100;
        double remainingTime=((100-progress)*passedTime)/progress;
        cout << "Processing read number " << numOfReads << "/" << numOfAllReads
                << " which  means: " << progress
                << "%  progress. Approximate remaining time is "
                << remainingTime
                << " minutes. Error correction is less than "
                << (numOfSupportedReads / numOfReads) * 100 << "%"
                << "\r";
        cout.flush();
}

/**
 * process one read file for error correction
 */
void ReadCorrection::CorrectErrorsInLibrary(ReadLibrary *input) {
        library = input;
        numOfAllReads = library->getNumOfReads();
        //opening file for input and output
        readsFile.open(library->getFilename().c_str(), ios::in);
        outFastq.open(library->getOutputFileName(), ios::out);
        clock_t begin=clock();
        while (readsFile.is_open()) {
                vector<readStructStr> reads;
                readInputReads(reads);
                correctReads(reads);
                writeOutputReads(reads);
                printProgress(begin);
        }
        outFastq.close();
        cout << endl << "Number of reads which are found in graph: "
                << numOfSupportedReads << " out of " << numOfReads << " which is "
                << (100 * numOfSupportedReads / numOfReads) << "%" << endl;
}

/**
 * process all read files for error correction
 */
void ReadCorrection::errorCorrection(LibraryContainer &libraries) {
        for (size_t i = 0; i <libraries.getSize(); i++) {
                ReadLibrary &input = libraries.getInput(i);

                cout << "Processing file " << i+1 << "/" << libraries.getSize()
                        << ": " << input.getFilename() << ", type: "
                        << input.getFileType() << endl;
                CorrectErrorsInLibrary(&input);
        }
}

/**
 * 
 * @return 
 */
bool ReadCorrection::checkForAnswer(Kmer const &kmer, int startOfRead,
                string const &erroneousRead,
                string &guessedRead, string const &qualityProfile,
                readCorrectionStatus &status) {
        NodePosPair result = dbg.getNodePosPair(kmer);
        if (!result.isValid()) {
                return false;
        }
        if (!dbg.kmerExistsInGraph(kmer)) {
                return false;
        }
        status = kmerfound;
        int startOfNode = result.getPosition();
        SSNode leftNode = dbg.getSSNode(result.getNodeID());
        string nodeContent = leftNode.getSequence();
        return recursiveCompare(leftNode, nodeContent, startOfNode, startOfRead,
                        erroneousRead, guessedRead, qualityProfile, status);
}

/**
 * 
 * @return 
 */
bool ReadCorrection::recKmerCorrection(string &kmerStr,
                string const &qualityProfile,
                int kmerStart, int round) {
        if(round > 5) {
                return false;
        }
        string oriKmer = kmerStr;
        int pos = lowQualityPos(qualityProfile, kmerStart, oriKmer, round);
        char currentBase = kmerStr[pos - kmerStart];
        char bases[4] = {'A', 'C', 'G', 'T'};
        for(char const &x : bases) {
                if (currentBase != x) {
                        kmerStr[pos - kmerStart] = x;
                        if (dbg.kmerExistsInGraph(kmerStr)) {
                                return true;
                        }
                }
        }
        for(char const &x : bases) {
                if (currentBase != x) {
                        kmerStr[pos - kmerStart] = x;
                        if (recKmerCorrection(kmerStr, qualityProfile,
                                                kmerStart, round + 1)) {
                                return true;
                        }
                }
        }
        //no corrections found
        return false;
}

/**
 * 
 * @return 
 */
int ReadCorrection::lowQualityPos(string quality, int startOfRead,
                string const &kmer, int round) {
        int lowIndex = startOfRead;
        int preLow = 0;
        int preIndex = -1;
        for(int j = 0; j < round; ++j) {
                char lowValue = '~';
                int k = 0;
                for(unsigned int i = startOfRead;
                                i < (startOfRead + kmer.length()); ++i) {
                        if (lowValue > quality[i]
                                        && quality[i] >= preLow
                                        && preIndex != i ) {
                                lowIndex = i;
                                lowValue = quality[i];
                                if(kmer[k] == 'N') {
                                        quality[i] = '!';
                                }
                        }
                        k++;
                }
                preLow = lowValue;
                preIndex = lowIndex;
        }
        return lowIndex;
}

/**
 * 
 * @return 
 */
string ReadCorrection::applyINDchanges(string const &reference,
                string const &read) {
        string newRead;
        for (unsigned int i = 0; i < reference.length() && i < read.length();
                        ++i) {
                char cInRef = reference[i];
                char cInRead = read[i];
                if (cInRead != '-' && cInRef != '-') {
                        newRead = newRead + cInRead;
                }
                else {
                        if(cInRef == '-') {
                                //insertion occured delete char from read
                                //do nothing
                        } else {
                                if(cInRead == '-') {
                                        //deletion occured insert char to read
                                        newRead = newRead + cInRef;
                                }
                        }
                }
        }
        return newRead;
}

/**
 * 
 * @return 
 */
bool ReadCorrection::checkForIndels(string const &ref, string query,
                int const maxError, string const &qualityProfile,
                string &newRead) {
        NW_Alignment Nw;
        string refCopy1 = ref;

        double fullsim = Nw.alignment(refCopy1, refCopy1);
        double sim = Nw.alignment(refCopy1, query);
        if ((sim / fullsim) > .8 ) {
                newRead = applyINDchanges(refCopy1, query);
                double dif = findDifference(ref, newRead, qualityProfile, 0);
                if (dif < maxError) {
                        return true;
                }
        }
        return false;
}

bool ReadCorrection::expand(SSNode const &leftNode, string const &erroneousRead, string const &qualityProfile, string &guessedRead, pair<int, int> bounds, bool forward, readCorrectionStatus &status) {
        bool found = false;
        int readLength = erroneousRead.length() < qualityProfile.length() ? erroneousRead.length() : qualityProfile.length();
        if (leftNode.getNumRightArcs() > 0) {
                string remainingInRead = erroneousRead.substr(bounds.first, bounds.second);
                string remainingInQuality = qualityProfile.substr(bounds.first, bounds.second);
                vector<string> rightResults = getAllSolutions(leftNode, remainingInRead, remainingInQuality, forward);
                if (rightResults.size() > 0) {
                        string bestMatch;
                        if (findBestMatch(rightResults, remainingInRead, true, bestMatch, readLength)) {
                                guessedRead.replace(bounds.first, bestMatch.length(), bestMatch);
                                found = true;
                        }
                }
                else {
                        status = graphIsMissing;
                }
        } else {
                status = graphIsMissing;
        }
        return found;
}

/**
 * 
 * @return 
 */
bool ReadCorrection::recursiveCompare(SSNode const &leftNode,
                string const &nodeContent, 
                int startOfNode, int startOfRead,
                string const &erroneousRead,
                string &guessedRead,
                string const &qualityProfile,
                readCorrectionStatus &status) {

        string initialGuessedRead = guessedRead;
        int readLength = erroneousRead.length() < qualityProfile.length() ? erroneousRead.length() : qualityProfile.length();
        int nodeLength = nodeContent.length();

        bool expandFromRight = nodeLength - startOfNode >= readLength - startOfRead ? false : true;
        bool expandFromLeft = startOfNode - startOfRead >= 0 ? false : true;
        bool rightFound = false;
        bool leftFound = false;
        bool middleFound = false;

        int readRightExtreme = nodeLength - startOfNode >= readLength - startOfRead ? readLength : nodeLength - startOfNode + startOfRead;
        int readLeftExtreme = startOfNode - startOfRead >= 0 ? 0 : startOfRead - startOfNode;
        string commonInRead = erroneousRead.substr(readLeftExtreme, readRightExtreme - readLeftExtreme);

        int nodeRightExtreme = nodeLength - startOfNode >= readLength - startOfRead ? readLength - startOfRead + startOfNode : nodeLength;
        int nodeLeftExtreme = startOfNode - startOfRead >= 0 ? startOfNode - startOfRead : 0;
        string commonInNode = nodeContent.substr(nodeLeftExtreme, nodeRightExtreme - nodeLeftExtreme );

        int errorCommonThreshold = sqrt(commonInRead.length()) * 33;
        int maxError = dbg.getReadLength() * 33 * .2;

        if (expandFromLeft || expandFromRight) {
                if (findDifference(commonInNode, erroneousRead, qualityProfile, readLeftExtreme) < errorCommonThreshold) {
                        //TODO GILES > why left
                        guessedRead.replace(readLeftExtreme, commonInRead.length(), commonInNode);
                        middleFound = true;
                }
                if (expandFromRight) {
                        pair<int, int> bounds(readRightExtreme, readLength - readRightExtreme);
                        rightFound = expand(leftNode, erroneousRead, qualityProfile, guessedRead, bounds, true, status);
                }
                if (expandFromLeft) {
                        pair<int, int> bounds(0, readLeftExtreme);
                        leftFound = expand(leftNode, erroneousRead, qualityProfile, guessedRead, bounds, false, status);
                }
        } else {
                int dif=findDifference(commonInNode, erroneousRead, qualityProfile, 0);
                if (dif < maxError) {
                        guessedRead = commonInNode;
                        status = fullHealing;
                        return true;
                }
                else {
                        string newRead;
                        if (checkForIndels(commonInNode, erroneousRead, maxError, qualityProfile, newRead)) {
                                guessedRead = newRead;
                                status = fullHealing;
                                return true;
                        } else {
                                status = kmerfound;
                                return false;
                        }
                }
        }

        int dif = findDifference(guessedRead, erroneousRead, qualityProfile, 0);
        if (dif < maxError) {
                if ((expandFromLeft && !leftFound)
                                || (expandFromRight && !rightFound)
                                || !middleFound)
                        status = parHealing;
                return true;
        }
        else {
                string newRead;
                if (checkForIndels(guessedRead, erroneousRead, maxError, qualityProfile, newRead)) {
                        //TODO GILES > updated value of newRead is not used
                        if ((expandFromLeft && !leftFound)
                                        || (expandFromRight && !rightFound)
                                        || !middleFound)
                                status = parHealing;
                        return true;
                }
                //revert the guess to the state at the start of this method
                guessedRead = initialGuessedRead;
        }
        return false;
}

/**
 * picks closest result to erroneousread
 * @return true if bestMatch is changed
 */
bool ReadCorrection::findBestMatch(vector<string> const &results,
                string &erroneousRead,
                bool rightDir, string &bestMatch, int readLength) {
        bool find = false;
        if(!rightDir) {
                std::reverse(erroneousRead.begin(), erroneousRead.end());
        }
        double minSim = 20;
        NW_Alignment Nw;
        for(auto const &result : results) {
                string strItem = result;
                if (!rightDir) {
                        std::reverse(strItem.begin(), strItem.end());
                }
                double newSim = Nw.get_similarity_per(strItem, erroneousRead);
                if (newSim > minSim) {
                        bestMatch = result;
                        minSim = newSim;
                        find = true;
                        if (newSim == 100)
                                return find;
                }
        }
        return find;
}

/**
 * count number of differences between originalRead and guessedRead
 *      originalRead starting at pos startOfRead, guessedRead starting at pos 0
 * @return weighted sum of qualityvalues at differences
 */
int ReadCorrection::findDifference(string const &guessedRead,
                string const &originalRead,
                string const &qualityProfile,
                int startOfRead) {
        unsigned int i = startOfRead;
        unsigned int j = 0;
        unsigned int d = 0;
        unsigned int predif = 0;
        while(i < originalRead.length() && j < guessedRead.length()) {
                if (originalRead[i] != guessedRead[j]
                                && originalRead[i] != 'N'
                                && guessedRead[j] != 'N') {
                        if (i != 0)
                                d = d + ((int) qualityProfile[i]) / (i - predif);
                        else
                                d = d + ((int) qualityProfile[i]);
                        predif = i;
                }
                i++;
                j++;
        }
        return d;
}

/**
 * count number of differences between two strings
 * @return number of positions in which the strings differ
 */
int ReadCorrection::findDifference(string const &a, string const &b) {
        unsigned int d = 0;
        for(unsigned int i = 0; i < a.length() && i < b.length(); ++i) {
                if (a[i] != b[i] && a[i] != 'N' && b[i] != 'N') {
                        d++;
                }
        }
        int dif = a.length() > b.length() ? a.length() - b.length()
                : b.length() - a.length();
        return d + dif;

}

typedef std::pair<SSNode, string> nodePath;
typedef std::pair<nodePath, int> minHeapElement;

/**
 * create a list of all possible solutions in the graph
 * @return list of 
 */
vector<string> ReadCorrection::getAllSolutions(SSNode const &rootNode,
                string const &readPart,
                string const &qualityProfile,
                bool forward) {
        vector<string> results;
        string root;
        std::stack<nodePath> mystack;
        clock_t begin = clock();
        mystack.push(make_pair(rootNode, root));
        while (!mystack.empty()) {
                pair<SSNode, string> r = mystack.top();
                mystack.pop();
                SSNode leftNode = r.first;
                string currentPath = r.second;

                for (ArcIt it = (forward ? leftNode.rightBegin() : leftNode.leftBegin());
                                it !=  (forward ? leftNode.rightEnd() : leftNode.leftEnd());
                                ++it) {

                        SSNode rrNode = dbg.getSSNode(it->getNodeID());
                        if (rrNode.getNodeID() == -leftNode.getNodeID() || !rrNode.isValid())
                                continue;
                        string nodeContent = rrNode.getSequence();
                        string content;
                        string newPath;
                        if (forward) {
                                content = nodeContent.substr(kmerSize - 1, nodeContent.length());
                                newPath = currentPath + content;
                        } else {
                                content = nodeContent.substr(0, nodeContent.length() - (kmerSize - 1));
                                newPath = content + currentPath;
                        }

                        if (newPath.length() < readPart.length()) {
                                if (newPath.length() > kmerSize ) {
                                        double errorDif;
                                        double errorCommonThreshold = newPath.length() * .3 * 25;
                                        if (forward) {
                                                errorDif = findDifference(newPath, readPart, qualityProfile, 0);
                                        } else {
                                                string tempNewPath = newPath;
                                                string tempQualityProfile = qualityProfile;
                                                string tempReadPart = readPart;
                                                std::reverse(tempNewPath.begin(), tempNewPath.end());
                                                std::reverse(tempQualityProfile.begin(), tempQualityProfile.end());
                                                std::reverse(tempReadPart.begin(), tempReadPart.end());
                                                errorDif = findDifference(tempNewPath, tempReadPart, tempQualityProfile, 0);
                                        }
                                        if (errorDif < errorCommonThreshold) {
                                                mystack.push(make_pair(rrNode, newPath));
                                        }
                                } else {
                                        mystack.push(make_pair(rrNode, newPath));
                                }
                        } else {
                                string newStr = newPath.substr(0, readPart.length());
                                results.push_back(newStr);
                        }
                        if (results.size() > 200) //stop because it is taking too much effort
                                break;
                }
                clock_t endt = clock();
                double passedTime = (double) (endt - begin) / (double) (CLOCKS_PER_SEC * 60);
                if (passedTime > .5) { //stop because it is taking too much time
                        cout << "this read took more than half minute to process";
                        break;
                }
        }
        return results;
}






