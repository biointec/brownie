#include "alignment.h"
#include "graph.h"
#include "kmernode.h"
#include "settings.h"
#include <list>
#include <queue>
#include <string.h>
#include <stack>
#include "library.h"



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
bool sortMemResultbySize ( match_t left, match_t right )
{
        return left.len>right.len;
}
void DBGraph::CorrectErrorsInLibrary(ReadLibrary &library) {

        cout <<endl<<"welcome to error correction part"<<endl;
        //*******************************************************
        FileType type =library.getFileType();
        string readFileName=library.getFilename();
        size_t numOfAllreads= library.getNumOfReads();
        unsigned int kmerSize=settings.getK();
        //*******************************************************
        table = new KmerNodeTable ( settings, numNodes );
        cout<<"populateTable"<<endl;
        table->populateTable ( nodes );
        //*******************************************************
        cout<<"making reference for essaMEM"<<endl;
        vector<string> references =makeRefForessaMEM();
        cout<<"using essaMEM"<<endl;
        std::string meta = "noise";
        vector<sparseSA*> saVec;
        for (unsigned int i=0; i<references.size(); i++) {
                sparseSA* sa1 = init_essaMEM(references[i], meta);
                saVec.push_back(sa1);
        }
        //opening file for input and output
        ifstream readsFile;
        readsFile.open(readFileName.c_str(),ios::in);
        ofstream outFastq;

        outFastq.open(library.getOutputFileName(), ios::out);

        readStructStr readInfo;
        string nodeContent;
        string guessedRead;
        string correctRead;
        string qualityProfile;
        string erroneousRead;
        string  strID;



        double numOfReads=0;
        int numOfSupportedReads=0;


        unsigned int readLength=0;
        clock_t begin=clock();

        while(getline(readsFile,readInfo.strID )&&getline(readsFile,readInfo.erroneousReadContent )) {


                if (FASTQ==type)
                {
                        getline(readsFile,readInfo.orientation );
                        getline(readsFile,readInfo.qualityProfile );
                } else {
                        readInfo.orientation='+';
                        //readInfo.qualityProfile=readInfo.erroneousReadContent;
                        for (unsigned int i=0; i<readInfo.qualityProfile.length(); i++) {
                                readInfo.qualityProfile[i]='#';
                        }
                }
                numOfReads++;

                if ( (int)numOfReads % OUTPUT_FREQUENCY == 0 ) {
                        clock_t endt=clock();
                        double passedTime=double(endt-begin)/(double) (CLOCKS_PER_SEC*60);
                        double progress = ((double)numOfReads/(double)numOfAllreads)*100;
                        double remainingTime=((100-progress)*passedTime)/progress;
                        cout << "Processing read number " << numOfReads<<"/"<<numOfAllreads
                        <<" which  means: " << progress<<"%  progress. Approximate remaining time is " <<remainingTime <<" minutes. Error correction is less than " <<(numOfSupportedReads/numOfReads)*100 <<"%"<<"\r";
                        cout.flush();
                }


                readInfo.intID=numOfReads;
                readLength=readInfo.erroneousReadContent.length();
                bool found=false;

                if (readLength<50)
                {
                        int stop=0;
                        stop++;
                }
                readCorrectionStatus status=kmerNotfound;
                //fullHealing,
                //parHealing,
                //kmerNotfound,
                //graphIsMissing

                SSNode leftNode;
                erroneousRead=readInfo.erroneousReadContent;
                correctRead=readInfo.corrctReadContent;
                qualityProfile=readInfo.qualityProfile;
                strID=readInfo.strID;


                int kmerStart=0;

                TString read =erroneousRead;
                guessedRead=erroneousRead;

                //find in regular way



                int startOfRead=0;
                if (erroneousRead.length()>=kmerSize) {
                        for ( TStringIt it = read.begin(); it != read.end(); it++ ) {
                                Kmer kmer = *it;
                                if (!cheakForAnswer( kmer,  startOfRead, correctRead,  erroneousRead,guessedRead , qualityProfile, status ) ) {
                                        startOfRead++;
                                        if (status==anotherKmer||status==kmerNotfound)
                                                continue;
                                }
                                else {
                                        found=true;
                                        break;
                                }
                        }
                }
                if (status==kmerNotfound) {
                        kmerStart=0;
                        while (kmerStart+kmerSize<readLength && erroneousRead.length()>=kmerSize) {
                                string tempstr=erroneousRead.substr(kmerStart,kmerSize);
                                NodePosPair result=table->find(tempstr);
                                if(!result.isValid())
                                        if( recKmerCorrection(tempstr,qualityProfile,kmerStart ,1 )) {
                                                Kmer kmer = tempstr;
                                                if (cheakForAnswer( kmer,  kmerStart, correctRead,  erroneousRead,guessedRead , qualityProfile,status )) {
                                                        found=true;
                                                        break;
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
                                                if (m.len==kmerSize)
                                                        continue;//we already in previous step checked this kmer

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
                                                                if (cheakForAnswer( kmer,  startOfRead, correctRead,  erroneousRead,guessedRead , qualityProfile ,status)) {
                                                                        found=true;
                                                                        break;
                                                                }
                                                        }
                                        }
                                        if (found)
                                                break;
                                }
                                if (found)
                                        break;
                                r++;
                        }
                }



                //all attempt for error correction finished, now write the guessed Read into the file

                if (guessedRead.length()!=readInfo.qualityProfile.length()) {
                        // cout <<readInfo.intID<<endl;
                        if (guessedRead.length()>readInfo.qualityProfile.length()) {
                                string str =guessedRead.substr(0, readInfo.qualityProfile.length());
                                guessedRead=str;
                        }
                        if (guessedRead.length()<readInfo.qualityProfile.length()) {
                                string str =readInfo.qualityProfile.substr(0, guessedRead.length());
                                readInfo.qualityProfile=str;
                        }

                }
                if (found)
                        numOfSupportedReads++;
                string fastqStrID=readInfo.strID;
                fastqStrID.erase(0, 1);
                fastqStrID="@"+fastqStrID;
                if(guessedRead.length()>readLength) {
                        //cout <<readInfo.intID<<endl;
                        string refstr=guessedRead.substr(0, readLength);
                        guessedRead=refstr;
                }
                if (guessedRead.length()<readLength) {
                        cout <<readInfo.intID<<endl;
                        guessedRead=readInfo.erroneousReadContent;
                }
                outFastq<<fastqStrID<<endl<<guessedRead<<endl<<readInfo.orientation<<endl<<readInfo.qualityProfile<<endl;
        }
        cout<<endl <<"Number of reads which are found in graph: "<< numOfSupportedReads<<" out of "<< numOfReads <<" which is "<< (numOfSupportedReads/numOfReads)*100 <<"%"<<endl;
        readsFile.close();
        outFastq.close();
}
void DBGraph::errorCorrection(LibraryContainer &libraries) {

        for (size_t i = 0; i <libraries.getSize(); i++) {
                ReadLibrary &input =libraries.getInput(i);

                cout << "Processing file " << i+1 << "/" << libraries.getSize() << ": "
                << input.getFilename() << ", type: " << input.getFileType()
                << endl;
                CorrectErrorsInLibrary(input);
        }
}
bool DBGraph::cheakForAnswer( Kmer kmer, int  startOfRead, string & correctRead, string & erroneousRead,string & guessedRead , string &qualityProfile ,readCorrectionStatus &status) {



        NodePosPair result = table->find ( kmer );
        if ( !result.isValid() )
                return false;

        int startOfNode=result.getPosition();
        SSNode leftNode=getSSNode(result.getNodeID());
        string nodeContent=leftNode.getSequence();
        status=kmerfound;
        if (recursiveCompare(leftNode,nodeContent, startOfNode, startOfRead,correctRead, erroneousRead, guessedRead,qualityProfile, status )) {

                return true;
        } else {

                return false;
        }

}
bool sortMemResultbyRefPos ( match_t left, match_t right )
{
        return (left.ref<right.ref) ;
}
bool sortMemResultbyQueryPos ( match_t left, match_t right )
{
        return (left.query<right.query) ;
}
sparseSA*  DBGraph::init_essaMEM(string &ref, std::string meta) {

        std::vector<std::string> refdescr;
        refdescr.push_back(meta);
        std::vector<long> startpos;
        startpos.push_back(0); //only one reference
        bool printSubstring = false;
        bool printRevCompForw = false;
        sparseSA * sa;
        sa = new sparseSA(ref,				//reference
                          refdescr,			//
                          startpos,			//
                          false,			//to use or not to use 4 column format
                          1,				//sparseness
                          false,			//suffixlinks
                          true,				//child arrays
                          1,				//skip parameter
                          printSubstring,	//
                          printRevCompForw	//
        );
        sa->construct();
        return sa;
}
vector<string> DBGraph::makeRefForessaMEM() {

        vector<string> list;
        string content1="";
        unsigned int maxSize=5000000;
        int i=-numNodes;
        while(i<numNodes) {
                if (i==0)
                        i++;
                SSNode n = getSSNode(i);
                if(!n.isValid()) {
                        i++;
                        continue;
                }
                stringstream convert;
                convert<<i;
                if (content1.size()+n.getSequence().size()<maxSize)
                        content1=content1+"<"+convert.str()+ ">"+ n.getSequence();
                else {
                        cout<<"split to new string"<<endl;
                        string newContent=content1;
                        content1="";
                        list.push_back(newContent);
                }
                i++;
        }
        list.push_back(content1);
        return list;

}
string DBGraph::findCorrectKmerWithessaMEM(TString &read, string & erroneousRead, string reference) {
        std::string ref = "";//"CATGGACTGACGTGCTTCTACTACATCATGCGACTTACTAC";

        std::string meta = "noise";
        ref=reference;
        int min_len = 20;
        //init the enhanced suffix array
        //mem finding
        int r=0;
        string guessedRead=read.getSequence();


        for ( TStringIt it = read.begin(); it != read.end(); it++ ) {
                int kmerSize=settings.getK();
                Kmer kmer = *it;
                string kmerstr=kmer.str();
                std::string query = kmer.str();
                std::string RCqueery=kmer.getReverseComplement().str();
                vector<match_t> matches; //will contain the matches
                bool print = 0; //not sure what it prints if set to 1
                long memCounter = 0; //this does not really work
                //sa.MEM(query, matches, min_len, print, memCounter, true, 1);
                sparseSA *sa = init_essaMEM(ref, meta);
                sa->MEM(query, matches, min_len, print, memCounter, true, 1);

                if (matches.size()>0) {
                        for (unsigned int i=0; i<matches.size(); i++) {
                                match_t m = matches[i];
                                if(m.len>=27 && m.query==0 ) {
                                        string refstr=reference.substr(m.ref, kmerSize);
                                        int d=findDifference(refstr, kmer.str());
                                        if (d<2) {
                                                guessedRead.replace(r,kmerSize,refstr);
                                                erroneousRead=guessedRead;
                                                return guessedRead;
                                        }
                                }
                                if (m.query==0) {
                                        int stop=0;
                                        stop++;
                                }
                        }
                }
                //delete sa;
                r++;
                if (r>10) {
                        break;
                }
        }
        return "";
}
bool DBGraph::recKmerCorrection(string &kmerStr,const string & qualityProfile ,int kmerStart ,int round) {


        if(round>5) {
                return false;
        }
        string oriKmer=kmerStr;

        int pos=lowQualityPos(qualityProfile,kmerStart,oriKmer,round);
        char currentBase=kmerStr[pos-kmerStart] ;
        if (currentBase!='A') {
                kmerStr[pos-kmerStart]='A';
                NodePosPair result=table->find(kmerStr);
                if (result.isValid())
                        return true;
        }
        if (currentBase!='T') {
                kmerStr[pos-kmerStart]='T';
                NodePosPair result=table->find(kmerStr);
                if (result.isValid())
                        return true;
        }

        if(currentBase!='C') {
                kmerStr[pos-kmerStart]='C';
                NodePosPair result=table->find(kmerStr);
                if (result.isValid())
                        return true;
        }

        if(currentBase!='G') {
                kmerStr[pos-kmerStart]='G';
                NodePosPair result=table->find(kmerStr);
                if (result.isValid())
                        return true;
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
int DBGraph::lowQualityPos(string quality,int startOfRead,string kmer, int round) {
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
string DBGraph::applyINDchanges(string reference, string read) {
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
bool DBGraph::checkForIndels(string ref, string query,const int maxError,string& qualityProfile, string& newRead ) {
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
bool DBGraph::recursiveCompare(SSNode leftNode,string nodeContent,int  startOfNode,int startOfRead, string & correctRead, string & erroneousRead,string & guessedRead , string &qualityProfile,readCorrectionStatus &status )
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

        int maxError=readLength*33*.2;
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
                                status=DBGraph::kmerfound;
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
                        //writeLocalCytoscapeGraph(leftNode.getNodeID(),leftNode.getNodeID(),150);
                        //SSNode newLeftNode=  getSSNode(leftNode.getNodeID()*-1);
                        //getAllRightSolutions2(newLeftNode,leftRemainingInRead,leftRemainingInQuality,leftRemainingInRead.length(),leftResults);
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
bool DBGraph::findBestMatch(vector<string>& results, string erroneousRead, string qualityProfile,bool rightDir , string& bestrightMatch, SSNode &leftNode,  int readLength) {
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
                //double newErrorDist=findDifference(strItem,erroneousRead, qualityProfile, 0);
                double newSim=Nw.get_similarity_per(strItem,erroneousRead);

                if (newSim>minSim) {
                        bestrightMatch=*it;
                        //cout<<bestrightMatch<<endl;
                        minSim=newSim;
                        find= true;
                        if (newSim==100)
                                return find;
                }
        }
        return find;
}

int DBGraph::findDifference(string guessedRead, string originalRead, string &qualityProfile, int startOfRead) {
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
int DBGraph::findDifference(string a, string b) {
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
        return d;//+dif;

}
typedef std::pair < SSNode,string> nodePath;
typedef std::pair<nodePath, int> minHeapElement;
void DBGraph::getAllRightSolutions(SSNode rootNode, string readPart,string qualityProfile ,unsigned int depth, std::vector<std::string> & results) {


        unsigned int kmerSize=settings.getK();

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
                SSNode leftNode= r.first; // currentRoot.first;
                string currentPath=r.second;
                for ( ArcIt it = leftNode.rightBegin(); it != leftNode.rightEnd(); it++ ) {
                        SSNode rrNode = getSSNode ( it->getNodeID() );
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
void DBGraph:: getAllLeftSolutions(SSNode rootNode,string readPart,string qualityProfile , unsigned int depth, std::vector<std::string> & results) {

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
                        SSNode llNode = getSSNode ( it->getNodeID() );
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


