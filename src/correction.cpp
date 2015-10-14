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
                                && getline(readsFile, readInfo.originalContent)) {
                        if (FASTQ == library->getFileType()) {
                                getline(readsFile,readInfo.orientation );
                                getline(readsFile,readInfo.qProfile );
                        } else {
                                readInfo.orientation='+';
                                int qvLength = readInfo.qProfile.length();
                                for (unsigned int i = 0; i < qvLength; i++) {
                                        readInfo.qProfile[i] = '#';
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
                if (correctRead(reads[i])) {
                        supportedReads++;
                }
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
                //write qProfile
                outFastq << readInfo.qProfile << endl;
        }
}

/**
 * Try to find a kmer from the read in the graph and correct the read
 * @return true if correction is found
 */
bool ReadCorrection::findKmer(readCorrectionStatus &status,
                string const &original, string &guess,
                string const &qProfile) {
        bool found = false;
        TString read = original;
        int startOfRead = 0;
        for (TStringIt it = read.begin(); it != read.end(); it++ ) {
                Kmer kmer = *it;
                if (!checkForAnswer(kmer, startOfRead, original, guess, qProfile, status)) {
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
                string const &original, string &guess,
                string const &qProfile) {
        bool found = false;
        int readLength = original.length();
        for (int kmerStart = 0; kmerStart + kmerSize < readLength
                        && kmerSize <= readLength; kmerStart += 5) {
                string tempstr = original.substr(kmerStart, kmerSize);
                if (!dbg.kmerExistsInGraph(tempstr)) {
                        if (recKmerCorrection(tempstr, qProfile, kmerStart, 1)) {
                                Kmer kmer = tempstr;
                                if (checkForAnswer(kmer, kmerStart, original, guess, qProfile, status)) {
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
bool ReadCorrection::correctionByKmer(readCorrectionStatus &status,
                string const &original, string &guess, string const &qProfile) {
        bool found = false;
        if (original.length() >= kmerSize) {
                found = findKmer(status, original, guess, qProfile);
        }
        if (status == kmerNotfound) {
                found = findSimilarKmer(status, original, guess, qProfile);
        }
        return found;
}

/**
 * Compute how far left we can shift the seed to the left
 */
int ReadCorrection::shiftSize(match_t const &m, string const &reference) {
        int shift = kmerSize - m.len < m.query ? kmerSize - m.len : m.query;
        int j = 0, k = 0;
        while (j < m.ref && reference[j - 1] != '>' && k < shift) {
                --j, ++k;
        }
        return k;
}

/**
 * check if the kmer starting at m.ref - shift is actually a kmer of nucleotides
 */
bool ReadCorrection::seedIsContainedInKmer(match_t const &m, string const &reference, int shift) {
        int j = m.ref - shift, k = 0;
        while (k < kmerSize
                        && (reference[j] == 'A'
                                || reference[j] == 'T'
                                || reference[j] == 'C'
                                || reference[j] == 'G'
                                || reference[j] == 'a'
                                || reference[j] == 't'
                                || reference[j] == 'c'
                                || reference[j] == 'g')) {
                ++j, ++k;
        }
        return k == kmerSize;
}
/**
 * Try to correct read using MEMs
 * @return true if correction is found
 */
bool ReadCorrection::correctionByMEM(vector<match_t> &matches, string const &reference,
                readCorrectionStatus &status, string const &original,
                string &guess, string const &qProfile) {
        for (int i = 0; i < matches.size(); i++) {
                match_t m = matches[i];
                if (m.len >= kmerSize) {
                        continue; //we already in previous step checked this kmer
                }
                int shift = shiftSize(m, reference);
                if (seedIsContainedInKmer(m, reference, shift)) {
                        string refstr = reference.substr(m.ref - shift, kmerSize);
                        string foundKmerStr = original.substr(m.query - shift, kmerSize);
                        float diff = findDifference(refstr, foundKmerStr);
                        int startOfRead = m.query - shift;
                        Kmer kmer = refstr;
                        if (kmerSize / (diff + 1) > 6 && checkForAnswer(kmer, startOfRead, original, guess, qProfile, status)) {
                                return true;
                        }
                }
        }
        return false;
}

/**
 * for every reference sequence, find MEMs and process them
 * @return true if correction is found
 */
bool ReadCorrection::correctionByMEM(readCorrectionStatus &status, string const &original, string &guess, string const &qProfile) {
        for (unsigned int r = 0; r < saVec.size(); ++r) {
                int min_len = kmerSize / 3;
                bool print = 0; //not sure what it prints if set to 1
                long memCounter = 0; //this does not really work
                string query = original;
                vector<match_t> matches; //will contain the matches
                saVec[r]->MEM(query, matches, min_len, print, memCounter, true, 1);
                sort(matches.begin(), matches.end(), sortMemResultbySize);
                if (correctionByMEM(matches, references[r], status, original, guess, qProfile)) {
                        return true;
                }
        }
        return false;
}

/**
 * correct a single read
 */
bool ReadCorrection::correctRead(readStructStr &readInfo) {
        unsigned int readLength=readInfo.originalContent.length();
        readCorrectionStatus status = kmerNotfound;
        string original =  readInfo.originalContent;
        string guess = original;
        string qProfile = readInfo.qProfile;
        //find in regular way
        bool found = correctionByKmer(status, original, guess, qProfile);
        //find with essaMEM
        if (status == kmerNotfound) {
                found = correctionByMEM(status, original, guess, qProfile);
        }
        //all attempt for error correction finished
        //now write the guessed Read into the file
        if (guess.length() > readLength) {
                guess = guess.substr(0, readLength);
        } else if (guess.length() < readLength) {
                cout << readInfo.intID << endl;
                guess = original;
        }
        if (guess.length() != readInfo.qProfile.length()) {
                if (guess.length() > readInfo.qProfile.length()) {
                        guess = guess.substr(0, readInfo.qProfile.length());
                }
                if (guess.length() < readInfo.qProfile.length()) {
                        readInfo.qProfile = readInfo.qProfile.substr(0, guess.length());
                }
        }
        //save the correction so we can write it to file later
        readInfo.corrctReadContent = guess;
        return found;
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
                string const &original,
                string &guess, string const &qProfile,
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
                        original, guess, qProfile, status);
}

/**
 *
 * @return
 */
bool ReadCorrection::recKmerCorrection(string &kmerStr,
                string const &qProfile,
                int kmerStart, int round) {
        if (round > 5) {
                return false;
        }
        string oriKmer = kmerStr;
        int pos = lowQualityPos(qProfile, kmerStart, oriKmer, round);
        char currentBase = kmerStr[pos - kmerStart];
        char bases[4] = {'A', 'C', 'G', 'T'};
        for (char const &x : bases) {
                if (currentBase != x) {
                        kmerStr[pos - kmerStart] = x;
                        if (dbg.kmerExistsInGraph(kmerStr)) {
                                return true;
                        }
                }
        }
        for (char const &x : bases) {
                if (currentBase != x) {
                        kmerStr[pos - kmerStart] = x;
                        if (recKmerCorrection(kmerStr, qProfile,
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
        for (int j = 0; j < round; ++j) {
                char lowValue = '~';
                int k = 0;
                for (unsigned int i = startOfRead;
                                i < (startOfRead + kmer.length()); ++i) {
                        if (lowValue > quality[i]
                                        && quality[i] >= preLow
                                        && preIndex != i ) {
                                lowIndex = i;
                                lowValue = quality[i];
                                if (kmer[k] == 'N') {
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
                } else {
                        if (cInRef == '-') {
                                //insertion occured delete char from read
                                //do nothing
                        } else {
                                if (cInRead == '-') {
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
                int const maxError, string const &qProfile,
                string &newRead) {
        NW_Alignment Nw;
        string refCopy1 = ref;

        double fullsim = Nw.alignment(refCopy1, refCopy1);
        double sim = Nw.alignment(refCopy1, query);
        if ((sim / fullsim) > .8 ) {
                newRead = applyINDchanges(refCopy1, query);
                double dif = findDifference(ref, newRead, qProfile, 0);
                if (dif < maxError) {
                        return true;
                }
        }
        return false;
}

bool ReadCorrection::expand(SSNode const &node, string const &original,
                string const &qProfile, string &guess,
                pair<int, int> bounds, bool forward,
                readCorrectionStatus &status) {
        bool found = false;
        int readLength = original.length() < qProfile.length()
                ? original.length()
                : qProfile.length();
        if (forward ? node.getNumRightArcs() : node.getNumLeftArcs() > 0) {
                string remainingInRead = original.substr(bounds.first, bounds.second);
                string remainingInQuality = qProfile.substr(bounds.first, bounds.second);
                vector<string> results = getAllSolutions(node, remainingInRead, remainingInQuality, forward);
                if (results.size() > 0) {
                        string bestMatch;
                        if (findBestMatch(results, remainingInRead, forward, bestMatch, readLength)) {
                                guess.replace(bounds.first, bestMatch.length(), bestMatch);
                                found = true;
                        }
                } else {
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
                string const &original,
                string &guess,
                string const &qProfile,
                readCorrectionStatus &status) {
        string initialguess = guess;
        int readLength = original.length() < qProfile.length() ? original.length() : qProfile.length();
        int nodeLength = nodeContent.length();
        bool leftFound = startOfNode >= startOfRead;
        bool rightFound = nodeLength - startOfNode >= readLength - startOfRead;
        int nodeRightExtreme = nodeLength - startOfNode >= readLength - startOfRead ? readLength - startOfRead + startOfNode : nodeLength;
        int nodeLeftExtreme = startOfNode > startOfRead ? startOfNode - startOfRead : 0;
        string commonInNode = nodeContent.substr(nodeLeftExtreme, nodeRightExtreme - nodeLeftExtreme );
        int maxError = readLength * avgQualityError * maxErrorRate;
        if (!leftFound || !rightFound) {
                //compute the common part in the read
                int readLeftExtreme = startOfRead > startOfNode ? startOfRead - startOfNode : 0;
                int readRightExtreme = nodeLength - startOfNode >= readLength - startOfRead ? readLength : nodeLength - startOfNode + startOfRead;
                string commonInRead = original.substr(readLeftExtreme, readRightExtreme - readLeftExtreme);
                //correct the left if needed
                if (!leftFound) {
                        pair<int, int> bounds(0, readLeftExtreme);
                        leftFound = expand(leftNode, original, qProfile, guess, bounds, false, status);
                }
                bool middleFound = false;
                int errorCommonThreshold = commonInRead.length()*avgQualityError * maxErrorRate;
                if (findDifference(commonInNode, original, qProfile, readLeftExtreme) < errorCommonThreshold) {
                        guess.replace(readLeftExtreme, commonInRead.length(), commonInNode);
                        middleFound = true;
                }
                if (!rightFound) {
                        pair<int, int> bounds(readRightExtreme, readLength - readRightExtreme);
                        rightFound = expand(leftNode, original, qProfile, guess, bounds, true, status);
                }
                int dif = findDifference(guess, original, qProfile, 0);
                if (dif < maxError) {
                        if (!leftFound || !rightFound || !middleFound) {
                                status = parHealing;
                        } else {
                                status = fullHealing;
                        }
                        return true;
                }
                string newRead;
                if (checkForIndels(guess, original, maxError, qProfile, newRead)) {
                        if (!leftFound || !rightFound || !middleFound) {
                                status = parHealing;
                        } else {
                                status = fullHealing;
                        }
                        guess = newRead;
                        return true;
                }
                //revert the guess to the state at the start of this method
                guess = initialguess;
                return false;
        }
        int dif = findDifference(commonInNode, original, qProfile, 0);
        if (dif < maxError) {
                guess = commonInNode;
                status = fullHealing;
                return true;
        }
        string newRead;
        if (checkForIndels(commonInNode, original, maxError, qProfile, newRead)) {
                guess = newRead;
                status = fullHealing;
                return true;
        }
        status = kmerfound;
        return false;
}


/**
 * picks closest result to original
 * @return true if bestMatch is changed
 */
bool ReadCorrection::findBestMatch(vector<string> const &results,
                string &original,
                bool rightDir, string &bestMatch, int readLength) {
        bool find = false;
        if (!rightDir) {
                std::reverse(original.begin(), original.end());
        }
        double minSim = 20;
        NW_Alignment Nw;
        for (auto const &result : results) {
                string strItem = result;
                if (!rightDir) {
                        std::reverse(strItem.begin(), strItem.end());
                }
                double newSim = Nw.get_similarity_per(strItem, original);
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
 * count number of differences between original and guess
 *      with original starting at pos start, guess starting at pos 0
 * @return weighted sum of qualityvalues at differences
 */
int ReadCorrection::findDifference(string const &guess, string const &original,
                string const &qProfile, int start) {
        unsigned int i = start;
        unsigned int j = 0;
        unsigned int d = 0;
        unsigned int predif = 0;
        while (i < original.length() && j < guess.length()) {
                if (original[i] != guess[j]
                                && original[i] != 'N'
                                && guess[j] != 'N') {
                        if (i != 0)
                                d = d + ((int) qProfile[i]) / (i - predif);
                        else
                                d = d + ((int) qProfile[i]);
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
        for (unsigned int i = 0; i < a.length() && i < b.length(); ++i) {
                if (a[i] != b[i] && a[i] != 'N' && b[i] != 'N') {
                        d++;
                }
        }
        int dif = a.length() > b.length() ? a.length() - b.length()
                : b.length() - a.length();
        return d + dif;

}

typedef pair<SSNode, string> nodePath;
typedef pair<nodePath, int> minHeapElement;

/**
 * create a list of all possible solutions in the graph
 * @return list of
 */
vector<string> ReadCorrection::getAllSolutions(SSNode const &rootNode,
                string const &readPart,
                string const &qProfile,
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
                                        double errorCommonThreshold = newPath.length() * maxErrorRate * avgQualityError;
                                        if (forward) {
                                                errorDif = findDifference(newPath, readPart, qProfile, 0);
                                        } else {
                                                string tempNewPath = newPath;
                                                string tempqProfile = qProfile;
                                                string tempReadPart = readPart;
                                                std::reverse(tempNewPath.begin(), tempNewPath.end());
                                                std::reverse(tempqProfile.begin(), tempqProfile.end());
                                                std::reverse(tempReadPart.begin(), tempReadPart.end());
                                                errorDif = findDifference(tempNewPath, tempReadPart, tempqProfile, 0);
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
                        if (results.size() > 200) {
                                //stop because it is taking too much effort
                                break;
                        }
                }
                clock_t endt = clock();
                double passedTime = (double) (endt - begin) / (double) (CLOCKS_PER_SEC * 60);
                if (passedTime > .1) { //stop because it is taking too much time
                        cout << "this read took more than .1 minute to process";
                        break;
                }
        }
        return results;
}
























