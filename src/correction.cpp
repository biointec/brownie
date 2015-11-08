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
        numOfReads(0),numOfAllReads(0),numOfSupportedReads(0),minSimPer(60),maxTimePerRead(.5),dijk(g)
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
                                size_t qvLength = readInfo.qProfile.length();
                                for (size_t i = 0; i < qvLength; i++) {
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
 * write the reads to the outputfile
 */
void ReadCorrection::writeOutputReads(vector<readStructStr> const &reads) {
        for (size_t i = 0; i < reads.size(); ++i) {
                readStructStr readInfo = reads[i];
                //write ID
                outFastq << readInfo.strID << endl;
                //write correction
                outFastq << readInfo.corrctReadContent << endl;
                //write orientation
                if (FASTQ == library->getFileType()) {
                        outFastq << readInfo.orientation << endl;
                        //write qProfile
                        outFastq << readInfo.qProfile << endl;
                }
        }
}

/**
 * Correct a batch of reads
 */
void ReadCorrection::correctReads(vector<readStructStr> &reads) {
        int supportedReads = 0;
        //#pragma omp parallel for reduction(+:supportedReads)
        for (size_t i = 0; i < reads.size(); ++i) {

                if (correctRead(reads[i]))
                        supportedReads++;
        }
}

bool ReadCorrection::correctRead(readStructStr &readInfo) {
/*
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
        checkReadSize(guess, readInfo);
        //now write the guessed Read into the file
        //save the correction so we can write it to file later
        readInfo.corrctReadContent = guess;
        return found;*/

        readCorrectionStatus status = kmerNotfound;
        string original =  readInfo.originalContent;
        string guess = original;
        string qProfile = readInfo.qProfile;

        bool found = false;
        if (original.length() >= kmerSize) {
                found = newCorrectRead(readInfo,status);
        }
        if (status == kmerNotfound) {
                found = findSimilarKmer(status, original, guess, qProfile);
        }
        if (status == kmerNotfound) {
                found = correctionByMEM(status, original, guess, qProfile);
        }
        return found;

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
                        if (status == anotherKmer|| status == kmerNotfound)
                                continue;
                        else
                                break;
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
        size_t readLength = original.length();
        for (size_t kmerStart = 0; kmerStart + kmerSize < readLength
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
        size_t j = m.ref - shift, k = 0;
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
        for (size_t i = 0; i < matches.size(); i++) {
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
                        if (kmerSize / (diff + 1) > 6 &&checkForAnswer(kmer, startOfRead, original, guess, qProfile, status)) {//
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
                << ((double)numOfSupportedReads /(double) numOfReads) * 100 << "%"
                << "\r";
        cout.flush();
}

/**
 * process one read file for error correction
 */
void ReadCorrection::CorrectErrorsInLibrary(ReadLibrary *input) {
        library = input;
        numOfAllReads = library->getNumReads();
        //opening file for input and output
        readsFile.open(library->getInputFilename().c_str(), ios::in);
        outFastq.open(library->getOutputFileName(), ios::out);
        clock_t begin=clock();
        while (readsFile.is_open()) {
                vector<readStructStr> reads;
                readInputReads(reads);
                clock_t start=clock();
                correctReads(reads);
                clock_t end=clock();
                cout<<"correction of these reads took "<< double(end-start)/(double) (CLOCKS_PER_SEC*60)<<" minutes"<<endl;
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
                        << ": " << input.getInputFilename() << ", type: "
                        << input.getFileType() << endl;
                CorrectErrorsInLibrary(&input);
        }
}
void ReadCorrection::checkReadSize(string &guess,readStructStr &readInfo){
        string original =  readInfo.originalContent;
        unsigned int readLength=readInfo.originalContent.length();
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
                for (size_t i = startOfRead;
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
bool ReadCorrection::expand(SSNode  &node, string const &original,
                string const &qProfile, string &guess,
                pair<int, int> bounds, bool forward,
                readCorrectionStatus &status) {
        bool found = false;
        if (forward ? node.getNumRightArcs() : node.getNumLeftArcs() > 0) {
                string remainingInRead = original.substr(bounds.first, bounds.second);
                vector<string> results = getAllSolutions(node, remainingInRead,  forward);
                if (results.size() > 0) {
                        string bestMatch;
                        if (findBestMatch(results, remainingInRead, forward, bestMatch)) {
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

bool ReadCorrection::recursiveCompare(SSNode &leftNode,
                string const &nodeContent,
                int startOfNode, int startOfRead,
                string const &original,
                string &guess,
                string const &qProfile,
                readCorrectionStatus &status) {
        string initialguess = guess;
        NW_Alignment Nw;
        int readLength = original.length() < qProfile.length() ? original.length() : qProfile.length();
        int nodeLength = nodeContent.length();
        bool leftFound = startOfNode >= startOfRead;
        bool rightFound = nodeLength - startOfNode >= readLength - startOfRead;
        int nodeRightExtreme = nodeLength - startOfNode >= readLength - startOfRead ? readLength - startOfRead + startOfNode : nodeLength;
        int nodeLeftExtreme = startOfNode > startOfRead ? startOfNode - startOfRead : 0;
        string commonInNode = nodeContent.substr(nodeLeftExtreme, nodeRightExtreme - nodeLeftExtreme );
        int readLeftExtreme = startOfRead > startOfNode ? startOfRead - startOfNode : 0;
        int readRightExtreme = nodeLength - startOfNode >= readLength - startOfRead ? readLength : nodeLength - startOfNode + startOfRead;
        string commonInRead = original.substr(readLeftExtreme, readRightExtreme - readLeftExtreme);
        guess.replace(readLeftExtreme, commonInRead.length(), commonInNode);
        if (!leftFound) {
                pair<int, int> bounds(0, readLeftExtreme);
                leftFound = expand(leftNode, original, qProfile, guess, bounds, false, status);
        }
        if (!rightFound) {
                pair<int, int> bounds(readRightExtreme, readLength - readRightExtreme);
                rightFound = expand(leftNode, original, qProfile, guess, bounds, true, status);
        }
        guess= checkForIndels(original, guess);
        if (!leftFound || !rightFound) {
                status = parHealing;
        } else {
                status = fullHealing;
        }
        if (Nw.get_similarity_perEnhanced(original , guess )>minSimPer){
                return true;
        }
        guess=initialguess;
        return false;
}


/**
 * picks closest result to original
 * @return true if bestMatch is changed
 */
bool ReadCorrection::findBestMatch(vector<string> const &results,
                string &original,
                bool rightDir, string &bestMatch) {
        bool find = false;
        NW_Alignment Nw;
        if (!rightDir) {
                std::reverse(original.begin(), original.end());
        }
        double minSim = minSimPer;
        for (auto const &result : results) {
                string strItem = result;
                if (!rightDir) {
                        std::reverse(strItem.begin(), strItem.end());
                }
                double newSim = Nw.get_similarity_perEnhanced(strItem, original);
                if (newSim > minSim) {
                        bestMatch = result;
                        minSim = newSim;
                        find = true;
                        if (newSim > 99)
                                return find;
                }
        }
        return find;
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
typedef pair<nodePath, double> HeapElement;

/**
 * create a list of all possible solutions in the graph
 * @return list of
 */


int ReadCorrection::findRecSolutionsForward(vector<string> &results, SSNode const rootNode,
                                      string const &readPart,
                                      string currentPath,clock_t& start
) {
        clock_t end=clock();
        double passedTime=double(end-start)/(double) (CLOCKS_PER_SEC*60);
        if(passedTime>maxTimePerRead)
                return 0;
        if (currentPath.length()>readPart.length()){
                string newStr = currentPath.substr(0, readPart.length());
                results.push_back(newStr);
        }
        else{
                NW_Alignment Nw;
                if (currentPath.length()<kmerSize|| Nw.get_similarity_perEnhanced(currentPath,readPart.substr(0,currentPath.length()) )>minSimPer){
                        for (ArcIt it=rootNode.rightBegin();it!=rootNode.rightEnd();it++){
                                SSNode rNode = dbg.getSSNode(it->getNodeID());
                                string nodeContent = rNode.getSequence();
                                string  path = currentPath + nodeContent.substr(kmerSize - 1, nodeContent.length());;
                                findRecSolutionsForward(results,rNode,readPart,path,start);
                        }

                }
        }
        return 1;
}


/**
 * search graph in backware direction to find the path mapping to the readpart, it initiate at rootNode and the lengh is the limition point
 * @return list of results which can be potentially mapped to the input readpart
 */

int ReadCorrection::findRecSolutionsBackward(vector<string> &results, SSNode const rootNode,
                                      string const &readPart,
                                      string currentPath,clock_t& start
) {
        clock_t end=clock();
        double passedTime=double(end-start)/(double) (CLOCKS_PER_SEC*60);
        if(passedTime>maxTimePerRead)
                return 0;
        if (currentPath.length()>readPart.length()){
                string newStr = currentPath.substr(0, readPart.length());
                results.push_back(newStr);
        }
        else{
                NW_Alignment Nw;
                std::reverse(currentPath.begin(), currentPath.end());
                if (currentPath.length()<kmerSize|| Nw.get_similarity_perEnhanced(currentPath,readPart.substr(0,currentPath.length()))>minSimPer){
                        std::reverse(currentPath.begin(), currentPath.end());
                        for (ArcIt it=rootNode.leftBegin();it!=rootNode.leftEnd();it++){
                                SSNode lNode = dbg.getSSNode(it->getNodeID());
                                string nodeContent = lNode.getSequence();
                                string  path =  nodeContent.substr(0, nodeContent.length() - (kmerSize - 1))+currentPath;
                                findRecSolutionsBackward(results,lNode,readPart,path,start);
                        }

                }
        }
        return 1;
}

struct comparator {
         bool operator()(const HeapElement& f, const HeapElement& s)
        {
                return s.second >= f.second;
        }
};
/**
 * get all solutions by traversing graph in forward and backward direction and save the result
 * @return list of results which can be potentially mapped to the input readpart
 */
vector<string> ReadCorrection::getAllSolutions(SSNode &rootNode,
                                               string const &readPart,
                                               bool forward)
{
        NW_Alignment Nw;
        vector<string> results;
        if (rootNode.isLoop())
                return results;
        priority_queue< HeapElement,vector<HeapElement> ,comparator> heap;
        clock_t start=clock();
        heap.push(make_pair( make_pair(rootNode, ""),0));
        while (!heap.empty()) {
                clock_t end=clock();
                double pass= double(end-start)/(double) (CLOCKS_PER_SEC*60);
                if (pass>maxTimePerRead){
                        if (results.size()==0)
                                rootNode.setLoop(true);
                        return results;
                }
                pair< pair<SSNode, string>, double> r =heap.top();
                heap.pop();;
                SSNode node = r.first.first;
                string currentPath = r.first.second;
                for (ArcIt it = (forward ? node.rightBegin() : node.leftBegin());
                     it !=  (forward ? node.rightEnd() : node.leftEnd());
                it++) {
                        SSNode rrNode = dbg.getSSNode(it->getNodeID());
                        if (rrNode.getNodeID() == -node.getNodeID() || !rrNode.isValid())
                                continue;
                        string nodeContent = rrNode.getSequence();
                        string newPath="";
                        if (forward)
                                newPath = currentPath + nodeContent.substr(kmerSize - 1, nodeContent.length());
                        else
                                newPath = nodeContent.substr(0, nodeContent.length() - (kmerSize - 1)) + currentPath;
                        if (!forward)
                                std::reverse(newPath.begin(), newPath.end());
                        if (newPath.length() < readPart.length()) {
                                double sim= Nw.get_similarity_perEnhanced(newPath,readPart.substr(0,newPath.length()));
                                if (newPath.length() < kmerSize||sim>minSimPer )
                                {
                                        if (!forward)
                                                std::reverse(newPath.begin(), newPath.end());
                                        heap.push(make_pair( make_pair(rrNode, newPath), sim));

                                }
                        } else {
                                string newStr = newPath.substr(0, readPart.length());
                                results.push_back(newStr);
                                if (findDifference(readPart, newStr)==0)
                                        return results;
                                if (results.size()>1000 )
                                        return results;
                        }
                }

        }
        return results;
}
/**
 * correct reads, searchng kmers in read to find hit then extend from right and left, then making bridges between unfound part of reads
 * go exhustivley at the begining and end of read to find the matched part in graph.
 *
 * @return true if it finds any kmer of read in graph.
 */

bool ReadCorrection::newCorrectRead(readStructStr &readInfo,readCorrectionStatus &status) {
        string read =  readInfo.originalContent;
        int startOfRead = 0;
        int readLength= readInfo.originalContent.length();
        string  guess = readInfo.originalContent;
        size_t numOfHit=0;
        int curReadLeftExtreme = 0,curReadRightExtreme=0,prevReadRightExtreme=0;
        int minLeftInRead=readLength;
        int maxRightInRead=0;
        SSNode currNode, prevNode,firstNode,lastNode;
        NW_Alignment Nw;

        while (startOfRead<=(readLength-kmerSize)) {
                Kmer kmer = read.substr(startOfRead, kmerSize);
                NodePosPair result = dbg.getNodePosPair(kmer);
                if (!result.isValid()){
                        startOfRead++;
                        continue;
                }
                numOfHit++;
                int startOfNode = result.getPosition();
                currNode = dbg.getSSNode(result.getNodeID());
                int nodeLength = currNode.getSequence().length();
                size_t nodeRightExtreme = nodeLength - startOfNode >= readLength - startOfRead ? readLength - startOfRead + startOfNode : nodeLength;
                size_t nodeLeftExtreme = startOfNode > startOfRead ? startOfNode - startOfRead : 0;
                curReadLeftExtreme = startOfRead > startOfNode ? startOfRead - startOfNode : 0;
                curReadRightExtreme = nodeLength - startOfNode >= readLength - startOfRead ? readLength : nodeLength - startOfNode + startOfRead;
                if(curReadLeftExtreme<prevReadRightExtreme){
                        nodeLeftExtreme=nodeLeftExtreme+prevReadRightExtreme-curReadLeftExtreme;
                        curReadLeftExtreme=prevReadRightExtreme;
                }
                string commonInNode = currNode.getSequence().substr(nodeLeftExtreme, nodeRightExtreme - nodeLeftExtreme );
                bool newBiggerMap= curReadRightExtreme-curReadLeftExtreme>maxRightInRead-minLeftInRead && curReadLeftExtreme<minLeftInRead;
                if (curReadLeftExtreme-prevReadRightExtreme>1 &&prevNode.isValid()){
                        vector<string> bridges;
                        string lostPart=read.substr(prevReadRightExtreme, curReadLeftExtreme-prevReadRightExtreme);
                        findBridge(bridges,prevNode,currNode,lostPart,"");
                        if (bridges.size() > 0) {
                                string bestMatch="";
                                if (findBestMatch(bridges,lostPart, true, bestMatch)) {
                                        guess.replace(prevReadRightExtreme, bestMatch.length(), bestMatch);
                                }
                        }
                         prevReadRightExtreme=curReadLeftExtreme;

                }
                if (((curReadLeftExtreme==prevReadRightExtreme )||(newBiggerMap))
                        &&(Nw.get_similarity_perEnhanced(read.substr(curReadLeftExtreme, curReadRightExtreme - curReadLeftExtreme),commonInNode )>minSimPer))
                {
                        if (newBiggerMap){
                                guess= readInfo.originalContent;
                                minLeftInRead=readLength;//to make sure it will be decreased later
                                maxRightInRead=0;
                        }
                        guess.replace(curReadLeftExtreme, curReadRightExtreme - curReadLeftExtreme, commonInNode);
                        startOfRead=curReadRightExtreme+1;
                        prevReadRightExtreme=curReadRightExtreme;
                        minLeftInRead=minLeftInRead>curReadLeftExtreme?curReadLeftExtreme:minLeftInRead;
                        firstNode=minLeftInRead==curReadLeftExtreme ?currNode:firstNode;
                        maxRightInRead=curReadRightExtreme;
                        prevNode=currNode;
                }
                else
                        startOfRead++;
        }
        if(numOfHit==0){
                status=ReadCorrection::kmerNotfound;
                return false;
        }
        if (firstNode.isValid()>0 && minLeftInRead>0)
        {
                vector<string> leftResult= getAllSolutions(firstNode,read.substr(0, minLeftInRead),false);
                if (leftResult.size() > 0) {
                        string bestMatch,leftPart=read.substr(0, minLeftInRead);
                        if (findBestMatch(leftResult,leftPart, false, bestMatch)) {
                                guess.replace(0, bestMatch.length(), bestMatch);
                        }
                }
        }
        if(currNode.isValid() && maxRightInRead<readLength){
                vector<string> rightResult= getAllSolutions(currNode,read.substr(maxRightInRead, readLength-maxRightInRead+1),true);
                if (rightResult.size() > 0) {
                        string bestMatch, rightPart=read.substr(maxRightInRead, readLength-maxRightInRead+1);
                        if (findBestMatch(rightResult, rightPart, true, bestMatch)) {
                                guess.replace(maxRightInRead, bestMatch.length(), bestMatch);
                        }
                }

        }

        guess= checkForIndels(read, guess);
        if (Nw.get_similarity_perEnhanced(read , guess )>minSimPer){
                readInfo.corrctReadContent=guess;
                status=ReadCorrection::fullHealing;
                return true;
        }
        return false;
}



void ReadCorrection::findBridge(vector<string> &results  ,SSNode startNode,SSNode &endNode,string& readPart,string currentPath){

        if (currentPath.length()>readPart.length()){
                for (ArcIt it=startNode.rightBegin();it!=startNode.rightEnd();it++){
                        SSNode rNode = dbg.getSSNode(it->getNodeID());
                        if (rNode.getNodeID()==endNode.getNodeID()){
                                string newStr = currentPath.substr(0, readPart.length());
                                results.push_back(newStr);
                                break;
                        }
                }
        }
        else{
                NW_Alignment Nw;
                if (currentPath.length()<kmerSize|| Nw.get_similarity_perEnhanced(currentPath,readPart.substr(0,currentPath.length()) )>minSimPer){
                        for (ArcIt it=startNode.rightBegin();it!=startNode.rightEnd();it++){
                                SSNode rNode = dbg.getSSNode(it->getNodeID());
                                string nodeContent = rNode.getSequence();
                                string  path = currentPath + nodeContent.substr(kmerSize - 1, nodeContent.length());;
                                //double length= dijk.shortestPath(rNode,endNode);
                                //if (readPart.length()-currentPath.length()>=length)
                                findBridge(results,rNode,endNode,readPart,path);
                        }

                }
        }
}

/**this procedure compares number of mismaches between two input string before and after alignment
 *      if number of mismaches decrease after alignment then it return the new aligned string
 *
 * @return the modified read if taking into account indels decrease number of errors, otherwise the initial input read
 */
string  ReadCorrection::checkForIndels(string const &ref, string query){
        NW_Alignment Nw;
        string refCopy=ref;
        string queryCopy=query;
        size_t initialDiff= findDifference(ref,query);
        Nw.enhancedAlignment(refCopy,queryCopy);
        size_t allignedDiff= findDifference(refCopy,queryCopy);
        if (allignedDiff<initialDiff){
                return( applyINDchanges(refCopy, queryCopy));
        }
        return query;
}
