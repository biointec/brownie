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


ReadCorrection::ReadCorrection(DBGraph &g, Settings &s) :
        dbg(g), library(NULL), settings(s), kmerSize(s.getK())
{
        cout << endl << "welcome to error correction part" << endl;
        //*******************************************************
        cout << "populateTable" << endl;
        dbg.populateTable();
        //*******************************************************
        cout << "making reference for essaMEM" << endl;
        references = makeRefForessaMEM(dbg);
        cout << "using essaMEM" << endl;
        std::string meta = "noise";
        for (unsigned int i = 0; i < references.size(); i++) {
                sparseSA* sa1 = init_essaMEM(references[i], meta);
                saVec.push_back(sa1);
        }
}

void ReadCorrection::readInputReads(vector<readStructStr> &reads, double &numOfReads) {
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
                                for (unsigned int i = 0; i < readInfo.qualityProfile.length(); i++) {
                                       readInfo.qualityProfile[i]='#';
                                }
                        }
                        numOfReads++;
                        readInfo.strID.erase(0, 1);
                        readInfo.strID = "@" + readInfo.strID;
                        readInfo.intID=numOfReads;
                        reads.push_back(readInfo);
                } else {
                        readsFile.close();
                }
        }
}

void ReadCorrection::correctRead(readStructStr &readInfo, int &numOfSupportedReads) {
        unsigned int readLength=readInfo.erroneousReadContent.length();
        bool found=false;
        if (readLength<50) {
                int stop=0;
                stop++;
        }
        readCorrectionStatus status=kmerNotfound;
        SSNode leftNode;
        string erroneousRead=readInfo.erroneousReadContent;
        string correctRead=readInfo.corrctReadContent;
        string qualityProfile=readInfo.qualityProfile;
        int kmerStart=0;
        TString read =erroneousRead;
        string guessedRead=erroneousRead;
        //find in regular way
        int startOfRead=0;
        if (erroneousRead.length()>=kmerSize) {
                for ( TStringIt it = read.begin(); it != read.end(); it++ ) {
                        Kmer kmer = *it;
                        if (!checkForAnswer( kmer,  startOfRead, correctRead,  erroneousRead,guessedRead , qualityProfile, status ) ) {
                                startOfRead++;
                                if (status==anotherKmer||status==kmerNotfound) {
                                        continue;
                                }
                        } else {
                                found=true;
                                break;
                        }
                }
        }
        if (status==kmerNotfound) {
                kmerStart=0;
                while (kmerStart+kmerSize<readLength && erroneousRead.length()>=kmerSize) {
                        string tempstr=erroneousRead.substr(kmerStart,kmerSize);
                        if (!dbg.kmerExistsInGraph(tempstr)) {
                                if(recKmerCorrection(tempstr, qualityProfile, kmerStart, 1)) {
                                        Kmer kmer = tempstr;
                                        if (checkForAnswer( kmer,  kmerStart, correctRead,  erroneousRead,guessedRead , qualityProfile,status )) {
                                                found=true;
                                                break;
                                        }
                                }
                        }
                        kmerStart=kmerStart+5;
                }
        }
        //find with essaMEM
        if (status==kmerNotfound) {
                unsigned int r=0;
                while(r<saVec.size()) {
                        int min_len =kmerSize/3;
                        bool print = 0; //not sure what it prints if set to 1
                        long memCounter = 0; //this does not really work
                        std::string query = erroneousRead;
                        vector<match_t> matches; //will contain the matches
                        saVec[r]->MEM(query, matches, min_len, print, memCounter, true, 1);
                        sort ( matches.begin(), matches.end(), sortMemResultbySize );
                        if (matches.size()>0) {
                                match_t bestMatch;
                                string bestRefstr="";
                                int bestShiftLeft=0;
                                unsigned int cutLength=0;
                                for (unsigned int i=0; i<matches.size(); i++) {
                                        match_t m = matches[i];
                                        unsigned int shiftLeft=kmerSize-m.len<m.query?  kmerSize-m.len :m.query;
                                        int j=m.ref;
                                        unsigned int k=0;
                                        int lastValue=0;
                                        if (m.len==kmerSize) {
                                                continue;//we already in previous step checked this kmer
                                        }
                                        while(j>0&& references[r][j-1]!='>' && k<shiftLeft) {
                                                lastValue=k;
                                                k++;
                                                j--;
                                        }
                                        shiftLeft=k;
                                        lastValue=0;
                                        j=m.ref-shiftLeft;
                                        k=0;
                                        while(references[r][j]!='<' && k<=kmerSize && (references[r][j]=='A'||references[r][j]=='T'||references[r][j]=='C'
                                                ||references[r][j]=='G'||references[r][j]=='a'||references[r][j]=='t'||references[r][j]=='c'||references[r][j]=='g')) {
                                                lastValue=k;
                                                j++;
                                                k++;
                                        }
                                        cutLength=lastValue;
                                        string refstr=references[r].substr(m.ref-shiftLeft, cutLength);
                                        string foundKmerStr=erroneousRead.substr(m.query-shiftLeft,cutLength);
                                        float d=findDifference(refstr, foundKmerStr);
                                        if (refstr.length()/(d+1)>6&&  cutLength==kmerSize && foundKmerStr.length()==kmerSize && refstr.length()==kmerSize) {
                                                bestMatch=m;
                                                bestRefstr=refstr;
                                                bestShiftLeft=shiftLeft;
                                                Kmer kmer = bestRefstr;
                                                startOfRead=bestMatch.query-bestShiftLeft;
                                                if (checkForAnswer( kmer,  startOfRead, correctRead,  erroneousRead,guessedRead , qualityProfile ,status)) {
                                                        found=true;
                                                        break;
                                                }
                                        }
                                }
                                if (found) {
                                        break;
                                }
                        }
                        if (found) {
                                break;
                        }
                        r++;
                }
        }
        //all attempt for error correction finished, now write the guessed Read into the file
        if (guessedRead.length()!=readInfo.qualityProfile.length()) {
                if (guessedRead.length()>readInfo.qualityProfile.length()) {
                        string str =guessedRead.substr(0, readInfo.qualityProfile.length());
                guessedRead=str;
                }
                if (guessedRead.length()<readInfo.qualityProfile.length()) {
                        string str =readInfo.qualityProfile.substr(0, guessedRead.length());
                        readInfo.qualityProfile=str;
                }
        }
        if (found) {
        numOfSupportedReads += 1;
        }
        if(guessedRead.length()>readLength) {
                string refstr=guessedRead.substr(0, readLength);
                guessedRead=refstr;
        }
        if (guessedRead.length()<readLength) {
                cout <<readInfo.intID<<endl;
                guessedRead=readInfo.erroneousReadContent;
        }
        //save the correction so we can write it to file later
        readInfo.corrctReadContent = guessedRead;
}


void ReadCorrection::CorrectErrorsInLibrary(ReadLibrary *input) {
        library = input;
        size_t numOfAllreads = library->getNumOfReads();
        double numOfReads = 0;
        int numOfSupportedReads = 0;
        //opening file for input and output
        readsFile.open(library->getFilename().c_str(), ios::in);
        outFastq.open(library->getOutputFileName(), ios::out);
        
        clock_t begin=clock();


        while (readsFile.is_open()) {
                vector<readStructStr> reads;
                readInputReads(reads, numOfReads);
                #pragma omp parallel for reduction(+:numOfSupportedReads)
                for (int i = 0; i < reads.size(); ++i) {
                        correctRead(reads[i], numOfSupportedReads);
                }
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
                
                clock_t endt=clock();
                double passedTime=double(endt-begin)/(double) (CLOCKS_PER_SEC*60);
                double progress = ((double)numOfReads/(double)numOfAllreads)*100;
                double remainingTime=((100-progress)*passedTime)/progress;
                cout << "Processing read number " << numOfReads << "/" << numOfAllreads
                         << " which  means: " << progress
                         << "%  progress. Approximate remaining time is " << remainingTime
                         << " minutes. Error correction is less than "
                         << (numOfSupportedReads / numOfReads) * 100 << "%" << "\r";
                cout.flush();
        }
        outFastq.close();
        cout << endl
             << "Number of reads which are found in graph: " << numOfSupportedReads
             << " out of " << numOfReads
             << " which is " << (numOfSupportedReads / numOfReads) * 100
             << "%" << endl;
}

void ReadCorrection::errorCorrection(LibraryContainer &libraries) {

        for (size_t i = 0; i <libraries.getSize(); i++) {
                ReadLibrary &input =libraries.getInput(i);

                cout << "Processing file " << i+1 << "/" << libraries.getSize() << ": "
                << input.getFilename() << ", type: " << input.getFileType()
                << endl;
                CorrectErrorsInLibrary(&input);
        }
}
bool ReadCorrection::checkForAnswer( Kmer kmer, int  startOfRead, string & correctRead, string & erroneousRead,string & guessedRead , string &qualityProfile ,readCorrectionStatus &status) {
        NodePosPair result = dbg.getNodePosPair(kmer);
        if (!result.isValid()) {
                return false;
        }
        if (!dbg.kmerExistsInGraph(kmer)) {
                return false;
        }
        int startOfNode=result.getPosition();
        SSNode leftNode=dbg.getSSNode(result.getNodeID());
        string nodeContent=leftNode.getSequence();
        status=kmerfound;
        if (recursiveCompare(leftNode,nodeContent, startOfNode, startOfRead,correctRead, erroneousRead, guessedRead,qualityProfile, status )) {
                return true;
        } else {
                return false;
        }
}

bool ReadCorrection::recKmerCorrection(string &kmerStr,const string & qualityProfile ,int kmerStart ,int round) {

        if(round>5) {
                return false;
        }
        string oriKmer=kmerStr;

        int pos=lowQualityPos(qualityProfile,kmerStart,oriKmer,round);
        char currentBase=kmerStr[pos-kmerStart] ;
        if (currentBase!='A') {
                kmerStr[pos-kmerStart]='A';
                

                if (dbg.kmerExistsInGraph(kmerStr)) {
                        return true;
                }
        }
        if (currentBase!='T') {
                kmerStr[pos-kmerStart]='T';
                if (dbg.kmerExistsInGraph(kmerStr)) {
                        return true;
                }
        }

        if(currentBase!='C') {
                kmerStr[pos-kmerStart]='C';
                if (dbg.kmerExistsInGraph(kmerStr)) {
                        return true;
                }
        }

        if(currentBase!='G') {
                kmerStr[pos-kmerStart]='G';
                if (dbg.kmerExistsInGraph(kmerStr)) {
                        return true;
                }
        }
        if (currentBase!='A') {
                kmerStr[pos-kmerStart]='A';
                return recKmerCorrection(kmerStr,qualityProfile,kmerStart,round+1);
        }
        if (currentBase!='T') {
                kmerStr[pos-kmerStart]='T';
                return recKmerCorrection(kmerStr,qualityProfile,kmerStart,round+1);
        }

        if(currentBase!='C') {
                kmerStr[pos-kmerStart]='C';
                return recKmerCorrection(kmerStr,qualityProfile,kmerStart,round+1);
        }

        if(currentBase!='G') {
                kmerStr[pos-kmerStart]='G';
                return recKmerCorrection(kmerStr,qualityProfile,kmerStart,round+1);
        }
        return true;
}
int ReadCorrection::lowQualityPos(string quality,int startOfRead,string kmer, int round) {
        int i=startOfRead;
        char lowValue='~';
        int lowIndex=i;
        int preLow=0;
        int preIndex=-1;
        int j=1;
        while(j<=round) {
                unsigned int i=startOfRead;
                lowValue='~';
                int k=0;
                while(i<(startOfRead+kmer.length())) {
                        if (lowValue>quality[i]&&quality[i]>=preLow && preIndex!=i ) {
                                lowIndex=i;
                                lowValue=quality[i];
                                if(kmer[k]=='N') {
                                        quality[i]='!';
                                }
                        }
                        i++;
                        k++;
                }
                preLow=lowValue;
                preIndex=lowIndex;
                j++;
        }
        return lowIndex;
}
string ReadCorrection::applyINDchanges(string reference, string read) {
        unsigned int  i=0;
        unsigned int j=0;
        string newRead="";
        while(i<reference.length()&&j<read.length()) {
                char cInRef=reference[i];
                char cInRead=read[j];
                if (cInRead!='-'&& cInRef!='-') {
                        newRead=newRead+cInRead;
                }
                else {
                        if(cInRef=='-') {
                                //insertion occured delete char from read
                                //do nothing
                        } else {
                                if(cInRead=='-') {//deletion occured inser char to read
                                        newRead=newRead+cInRef;
                                }
                        }
                }
                i++;
                j++;
        }
        return newRead;
}
bool ReadCorrection::checkForIndels(string ref, string query, const int maxError, string& qualityProfile, string& newRead ) {
        NW_Alignment Nw;
        string refCopy1=ref;

        double fullsim=Nw.alignment(refCopy1,refCopy1);
        double sim=Nw.alignment(refCopy1,query);
        if ((sim/fullsim)>.8 ) {
                newRead= applyINDchanges(refCopy1,query);
                double dif=findDifference( ref,newRead, qualityProfile,0);
                if (dif<maxError) {
                        return true;
                }
        }
        return false;
}
bool ReadCorrection::recursiveCompare(SSNode leftNode,string nodeContent,int  startOfNode,int startOfRead, string & correctRead, string & erroneousRead,string & guessedRead , string &qualityProfile,readCorrectionStatus &status )
{

        NW_Alignment Nw;
        string initialGuessedRead=guessedRead;

        int readLenth=erroneousRead.length()<qualityProfile.length()? erroneousRead.length():qualityProfile.length();
        int nodeLenght=nodeContent.length();

        int nodeRightExtreme =nodeLenght-startOfNode>=readLenth-startOfRead ?readLenth-startOfRead+startOfNode:nodeLenght;
        int nodeLeftExtreme =startOfNode-startOfRead>=0 ?startOfNode- startOfRead:0;

        int readRightExtreme=nodeLenght-startOfNode>=readLenth-startOfRead ?readLenth:nodeLenght-startOfNode+startOfRead;
        int readLeftExtream=startOfNode-startOfRead>=0? 0:startOfRead-startOfNode;

        bool expandFromRight=nodeLenght-startOfNode>=readLenth-startOfRead? false:true;
        bool expandFromLeft=startOfNode-startOfRead>=0 ?false:true;
        bool rightFound=false;
        bool leftFound=false;
        bool midleFound=false;
        string commonInRead="";
        commonInRead=erroneousRead.substr(readLeftExtream,readRightExtreme-readLeftExtream);
        string commonInNode="";
        commonInNode=nodeContent.substr(nodeLeftExtreme,nodeRightExtreme-nodeLeftExtreme );

        string leftRemainingInRead="";
        leftRemainingInRead=erroneousRead.substr(0,readLeftExtream);
        string rightRemainigInRead= "";
        rightRemainigInRead=erroneousRead.substr(readRightExtreme,readLenth-readRightExtreme);

        string leftRemainingInQuality=qualityProfile.substr(0,readLeftExtream);
        string rightRemainingInQuality=qualityProfile.substr(readRightExtreme,readLenth-readRightExtreme);

        int errorCommonThreshold=sqrt( commonInRead.length())*33;
        if ((expandFromLeft || expandFromRight) &&findDifference(commonInNode,erroneousRead, qualityProfile, readLeftExtream)<errorCommonThreshold) {
                guessedRead.replace(readLeftExtream,commonInRead.length(),commonInNode);
                midleFound=true;
        }

        int maxError=dbg.getReadLength()*33*.2;
        if (!expandFromLeft && !expandFromRight) {
                int dif=findDifference( commonInNode,erroneousRead, qualityProfile,0);
                if (dif<maxError) {
                        guessedRead=commonInNode;
                        status=fullHealing;
                        return true;
                }
                else {
                        string newRead="";
                        if (checkForIndels(commonInNode, erroneousRead,maxError,qualityProfile,newRead) ) {
                                guessedRead=newRead;
                                status=fullHealing;
                                return true;
                        } else {
                                status=kmerfound;
                                return false;
                        }
                }
        }

        string bestrightMatch="";
        std::vector<std::string> rightResults;
        if (expandFromRight) {
                if (leftNode.getNumRightArcs()>0) {
                        getAllRightSolutions(leftNode,rightRemainigInRead,rightRemainingInQuality, rightRemainigInRead.length(),rightResults );
                        if (rightResults.size()>0) {
                                if (findBestMatch(rightResults,rightRemainigInRead ,rightRemainingInQuality,true ,bestrightMatch, leftNode, readLenth)) {
                                        guessedRead.replace(readRightExtreme,bestrightMatch.length(),bestrightMatch);
                                        rightFound=true;
                                }
                        }
                        else {
                                status=graphIsMissing;
                        }
                } else {
                        status=graphIsMissing;
                }

        }

        string bestLeftMatch="";
        std::vector<std::string> leftResults;
        if(expandFromLeft) {
                if (leftNode.getNumLeftArcs()>0) {
                        getAllLeftSolutions(leftNode,leftRemainingInRead,leftRemainingInQuality ,leftRemainingInRead.length(),leftResults);
                        if (leftResults.size()>0) {
                                if( findBestMatch(leftResults,leftRemainingInRead ,leftRemainingInQuality,false ,bestLeftMatch, leftNode, readLenth)) {
                                        guessedRead.replace(0,bestLeftMatch.length(),bestLeftMatch);
                                        leftFound=true;
                                }
                        } else {
                                status=graphIsMissing;
                        }

                } else {
                        status=graphIsMissing;
                }
        }

        int dif=findDifference(guessedRead,erroneousRead, qualityProfile, 0);

        if (dif<maxError) {
                if ((expandFromLeft&&!leftFound)||(expandFromRight&&!rightFound)||(!midleFound))
                        status=parHealing;
                return true;
        }
        else {
                string newRead="";
                if (checkForIndels(guessedRead,erroneousRead,maxError,qualityProfile,newRead)) {
                        if ((expandFromLeft&&!leftFound)||(expandFromRight&&!rightFound)||!midleFound)
                                status=parHealing;
                        return true;
                }
                guessedRead=initialGuessedRead;
        }
        return false;


}
bool ReadCorrection::findBestMatch(vector<string>& results, string erroneousRead, string qualityProfile,bool rightDir , string& bestrightMatch, SSNode &leftNode,  int readLength) {
        std::string first= results[0];



        bool find=false;
        if(rightDir==false) {
                std::reverse(qualityProfile.begin(), qualityProfile.end());
                std::reverse(erroneousRead.begin(), erroneousRead.end());

        }
        double minSim=20;
        NW_Alignment Nw;
        for(std::vector<string>::iterator it =results.begin()  ; it != results.end(); ++it) {

                string strItem=*it;
                if (rightDir==false) {
                        std::reverse(strItem.begin(), strItem.end());
                }
                double newSim=Nw.get_similarity_per(strItem,erroneousRead);

                if (newSim>minSim) {
                        bestrightMatch=*it;
                        minSim=newSim;
                        find= true;
                        if (newSim==100)
                                return find;
                }
        }
        return find;
}

int ReadCorrection::findDifference(string guessedRead, string originalRead, string &qualityProfile, int startOfRead) {
        unsigned int i=startOfRead;
        unsigned int j=0;
        unsigned int d=0;
        unsigned int predif=0;
        while(i<originalRead.length()&&j<guessedRead.length()) {
                if (originalRead[i]!=guessedRead[j] &&originalRead[i]!='N'&& guessedRead[j]!='N') {
                        if (i!=0)
                                d=d+ ((int)qualityProfile[i])/(i-predif);
                        else
                                d=d+ ((int)qualityProfile[i]);
                        predif=i;
                }
                i++;
                j++;
        }
        return d;

}
int ReadCorrection::findDifference(string a, string b) {
        unsigned int i=0;
        unsigned int j=0;
        unsigned int d=0;
        while( i<a.length()&&j<b.length()) {
                if (a[i]!=b[j] && a[i]!='N' && b[i]!='N') {
                        d++;
                }
                i++;
                j++;
        }
        int dif=a.length()>b.length()?a.length()-b.length():b.length()-a.length();
        d=d+dif;
        return d;

}
typedef std::pair < SSNode,string> nodePath;
typedef std::pair<nodePath, int> minHeapElement;
void ReadCorrection::getAllRightSolutions(SSNode rootNode, string readPart,string qualityProfile ,unsigned int depth, std::vector<std::string> & results) {


        string root="";
        nodePath p = std::make_pair(rootNode ,root  );
        std::stack<nodePath> mystack;
        clock_t begin=clock();

        mystack.push(make_pair(rootNode,root));

        while (!mystack.empty()) {

                pair<SSNode, string> r=mystack.top();
                mystack.pop();
                clock_t endt=clock();
                double passedTime=double(endt-begin)/(double) (CLOCKS_PER_SEC*60);
                if (passedTime>.5) {
                        cout<<"this read took more than half minute to process";
                        break;
                }
                SSNode leftNode= r.first;
                string currentPath=r.second;
                for ( ArcIt it = leftNode.rightBegin(); it != leftNode.rightEnd(); it++ ) {
                        SSNode rrNode = dbg.getSSNode ( it->getNodeID() );
                        if (rrNode.getNodeID()==-leftNode.getNodeID()|| !rrNode.isValid())
                                continue;
                        string ndoeContent=rrNode.getSequence();
                        string content=ndoeContent.substr(kmerSize-1, ndoeContent.length());
                        string newPath=currentPath+content;
                        if (newPath.length() <depth) {
                                double errorCommonThreshold=newPath.length()*.3*25;
                                double errorDif=findDifference(newPath,readPart, qualityProfile, 0);
                                if (newPath.length()>kmerSize ) {
                                        if (errorDif<errorCommonThreshold) {
                                                nodePath newPair = std::make_pair( rrNode,newPath  );
                                                mystack.push(newPair);
                                        }
                                } else
                                        mystack.push(make_pair(rrNode,newPath));

                        }
                        else {
                                string newStr=newPath.substr(0,depth);
                                results.push_back(newStr);
                        }

                        if (results.size()>200)
                                break;
                }
        }

}
void ReadCorrection::getAllLeftSolutions(SSNode rootNode,string readPart,string qualityProfile , unsigned int depth, std::vector<std::string> & results) {

        unsigned int kmerSize=settings.getK();
        typedef std::pair < SSNode,string> nodePath;
        std::stack<nodePath> mystack;
        string root="";
        nodePath p = std::make_pair(rootNode ,root  );
        mystack.push(p);
        clock_t begin=clock();
        while (!mystack.empty()) {
                pair<SSNode, string> r=mystack.top();
                mystack.pop();
                clock_t endt=clock();
                double passedTime=double(endt-begin)/(double) (CLOCKS_PER_SEC*60);
                if (passedTime>.5) {
                        cout<<"this read took more than half minute to process";
                        break;
                }
                SSNode leftNode=r.first;
                string currentPath=r.second;
                for ( ArcIt it = leftNode.leftBegin(); it != leftNode.leftEnd(); it++ ) {
                        SSNode llNode = dbg.getSSNode ( it->getNodeID() );
                        if (llNode.getNodeID()==-leftNode.getNodeID()|| !llNode.isValid())
                                continue;
                        string ndoeContent=llNode.getSequence();
                        string content=ndoeContent.substr(0, ndoeContent.length()-kmerSize+1);
                        string newPath=content+currentPath;
                        if (newPath.length() <depth) {
                                string tempNewPath=newPath;
                                string tempQualityProfile=qualityProfile;
                                string tempReadPart=readPart;
                                double errorDif= findDifference(tempNewPath,tempReadPart, tempQualityProfile, 0);
                                double errorCommonThreshold=newPath.length()*.2*25;
                                if (newPath.length()>kmerSize ) {

                                        std::reverse(tempNewPath.begin(), tempNewPath.end());
                                        std::reverse(tempQualityProfile.begin(), tempQualityProfile.end());
                                        std::reverse(tempReadPart.begin(), tempReadPart.end());
                                        if (errorDif<errorCommonThreshold) {
                                                nodePath newPair = std::make_pair( llNode,newPath  );
                                                mystack.push(newPair);

                                        } else {
                                                int stop=0;
                                                stop++;
                                        }
                                } else {
                                        nodePath newPair = std::make_pair( llNode,newPath  );
                                        mystack.push(newPair);
                                }
                        } else {
                                string newStr=newPath.substr(newPath.length()-depth,depth);
                                results.push_back(newStr);
                        }
                }
                if (results.size()>200)
                        break;
        }
}


