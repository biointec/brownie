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
#endif
