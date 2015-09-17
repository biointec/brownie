#ifndef ESSAMEM_H
#define ESSAMEM_H

bool sortMemResultbyRefPos ( match_t left, match_t right )
{
        return (left.ref<right.ref) ;
}

bool sortMemResultbyQueryPos ( match_t left, match_t right )
{
        return (left.query<right.query) ;
}
bool sortMemResultbySize ( match_t left, match_t right )
{
        return left.len>right.len;
}

sparseSA*  init_essaMEM(string &ref, std::string meta) {

        std::vector<std::string> refdescr;
        refdescr.push_back(meta);
        std::vector<long> startpos;
        startpos.push_back(0); //only one reference
        bool printSubstring = false;
        bool printRevCompForw = false;
        sparseSA * sa;
        sa = new sparseSA(ref,                  //reference
                          refdescr,             //
                          startpos,             //
                          false,                //4 column format or not
                          1,                    //sparseness
                          false,                //suffixlinks
                          true,                 //child arrays
                          1,                    //skip parameter
                          printSubstring,       //
                          printRevCompForw      //
        );
        sa->construct();
        return sa;
}

vector<string> makeRefForessaMEM(DBGraph &dbg) {
        vector<string> list;
        string content1="";
        unsigned int maxSize=5000000;
        int i = -dbg.getNumNodes();
        while(i < dbg.getNumNodes()) {
                if (i==0)
                        i++;
                SSNode n = dbg.getSSNode(i);
                if(!n.isValid()) {
                        i++;
                        continue;
                }
                stringstream convert;
                convert<<i;
                if (content1.size() + n.getSequence().size() < maxSize)
                        content1 = content1 + "<" + convert.str() +  ">"
                                   + n.getSequence();
                else {
                        cout << "split to new string" << endl;
                        string newContent = content1;
                        content1 = "";
                        list.push_back(newContent);
                }
                i++;
        }
        list.push_back(content1);
        return list;
}

/*
string findCorrectKmerWithessaMEM(TString &read, string & erroneousRead, string reference) {
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
                std::string RCqueery = kmer.getReverseComplement().str();
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
*/
#endif
