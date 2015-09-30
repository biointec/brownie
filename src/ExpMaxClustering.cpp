
#include "ExpMaxClustering.h"

using namespace std;


float ExpMaxClustering::calculateMean(string clusterFileName){
        ifstream  frequencyStream;
        frequencyStream.open(clusterFileName.c_str());
        string line;
        double sumOfAll=0;
        double totalNum=0;
        float avg=0;
        double num=0;
        double index=0;

        //cout<<clusterFileName<<endl;
        while (!frequencyStream.eof())
        {
                frequencyStream>>line;
                vector<std::string> tokens= util.digestLine(line);
                index =util.strToDouble(tokens[0]);
                num= util.strToDouble(tokens[1]);
                //cout<<num<<endl;
                sumOfAll=sumOfAll+index*num;
                totalNum=totalNum+num;
                if (totalNum!=0)
                avg=(avg*(totalNum-num)+index*num)/(totalNum);

        }
        frequencyStream.close();
        return avg;
}
bool ExpMaxClustering::isErroneous(double index){
        float erronousProb=gsl_ran_poisson_pdf(index,curErronousClusterMean);
        float correctProb=gsl_ran_poisson_pdf(index,curCorrectClusterMean);
        if (erronousProb/correctProb>1)
                return true;
        return false;

}
bool ExpMaxClustering::classifier(){
        double  numOfcorrect = 0,  numOfIncorrect = 0, sumOfBoth=0;
        size_t index;
        ifstream  frequencyStream;
        ofstream correctNodeStream;
        ofstream erronousNodeStream;
        correctNodeStream.open(correctNodesFileName.c_str(), ios::out);
        erronousNodeStream.open(erronousNodesFileName.c_str(), ios::out);
        frequencyStream.open(frequencyFileName.c_str());
        correctNodeStream.setf(ios_base::fixed);
        erronousNodeStream.setf(ios_base::fixed);
        if(!frequencyStream.good()){
                std::cerr << "Error opening file "<< frequencyFileName << std::endl;
                return false;
        }
        string line;
        frequencyStream>>line;
        bool firstHasItem=false;
        bool secondHasItem=false;
        while (!frequencyStream.eof())
        {
                frequencyStream>>line;
                vector<string> tokens= util.digestLine(line);
                index =util.strToDouble(tokens[0]);
                numOfcorrect= util.strToDouble(tokens[2]);
                numOfIncorrect= util.strToDouble( tokens[3]);
                //cout<<fixed<<setprecision(2)<<numOfIncorrect<<endl;
                sumOfBoth=numOfcorrect+numOfIncorrect;
                if (isErroneous(index)){
                        erronousNodeStream<<fixed<<setprecision(2)<<index<<"," <<util.ConvertToString(sumOfBoth)<<endl;
                        firstHasItem=true;
                }
                else{
                        correctNodeStream<<fixed<<setprecision(2)<<index<<","<<fixed <<util.ConvertToString(sumOfBoth)<<endl;
                        secondHasItem=true;
                }

        }
        if (!firstHasItem ||!secondHasItem)
                this->numOfClusters=1;
        correctNodeStream.close();
        erronousNodeStream.close();
        frequencyStream.close();
        return true;
}
void ExpMaxClustering::updateParameter(){
        perCorrectClusterMean=curCorrectClusterMean;
        perErronousClusterMean=curErronousClusterMean;
        curErronousClusterMean=calculateMean(erronousNodesFileName);
        curCorrectClusterMean=calculateMean(correctNodesFileName);
}


void ExpMaxClustering::findIntersectionPoint(){
        double e=2.718281;
        double c=curErronousClusterMean/curCorrectClusterMean;
        intersectionPoint=(curErronousClusterMean-curCorrectClusterMean)* (log(e)/log(c));
        cout<<"The intersection of these curves is :"<<intersectionPoint<<endl;
}
void ExpMaxClustering::findIntersectionPoint(double curErronousClusterMean , double curCorrectClusterMean){
        double e=2.718281;
        double c=curErronousClusterMean/curCorrectClusterMean;
        intersectionPoint=(curErronousClusterMean-curCorrectClusterMean)* (log(e)/log(c));
        cout<<"The intersection of these curves is :"<<intersectionPoint<<endl;
}
void ExpMaxClustering::doClassification(){
        size_t i=0;
        while(fabs( perErronousClusterMean-curErronousClusterMean)>divergenceThreshold ||fabs( perCorrectClusterMean-curCorrectClusterMean)>divergenceThreshold){
                //cout<<"correctMean :      "<<curCorrectClusterMean<<endl;
                //cout<<"erronousMean:      "<<curErronousClusterMean<<endl;
                //cout<<"next round------------------------"<<endl;
                classifier();
                updateParameter();
                i++;
                if (i>20)
                        break;
        }
        cout<<"Mean of nodes in correct nodes cluster   :       "<<curCorrectClusterMean<< endl;
        cout<<"Mean of nodes in erronous nodes cluster  :       "<<curErronousClusterMean<<endl;
        findIntersectionPoint();
}
