#include <google/sparse_hash_set>
#include <iostream>
#include <fstream>
#include <tr1/functional>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <boost/lexical_cast.hpp>
#include <iomanip>
using namespace std;

class commonUtil{

        char delimiter ;
public:
        commonUtil(char seperator){
                delimiter=seperator;
        }
        commonUtil(){
                delimiter=',';
        }
        vector<string> digestLine(string line) {
                vector<string> elems;
                string token="";
                size_t i=0;
                while(i<line.length()){
                        if (line[i]==delimiter){
                                elems.push_back(token);
                                token="";
                        }
                        else
                                token=token+line[i];
                        i++;
                }
                elems.push_back(token);
                return elems;
        }
        string ConvertToString (float number){
                std::ostringstream buff;
                buff<< fixed <<setprecision(2)<<number;
                return buff.str();
        }
        double strToDouble(string str){
                double value;
                try
                {
                        value = boost::lexical_cast<double>(str);
                }
                catch (boost::bad_lexical_cast const&)
                {
                        value = 0;
                }
                return value;
        }
        float strToFloat(string str){
                return (atof(str.c_str()));
        }
};



class ExpMaxClustering{


private:
        float perErronousClusterMean;
        float perCorrectClusterMean;

        string erronousNodesFileName;
        string correctNodesFileName;
        std::string frequencyFileName;
        float divergenceThreshold;
        commonUtil util;


        void initialization(string FileName, string firstClusterFileName, string secondClusterFileName,float divergenceValue,float initialMeanOfFirstCluster,float initialMeanOfSecondCluster  ){
                divergenceThreshold=divergenceValue;
                perCorrectClusterMean=initialMeanOfSecondCluster;
                perErronousClusterMean=initialMeanOfFirstCluster;
                curErronousClusterMean=perErronousClusterMean+divergenceThreshold*2;
                curCorrectClusterMean=perCorrectClusterMean+divergenceThreshold*2;
                frequencyFileName=FileName;
                erronousNodesFileName=firstClusterFileName;//"erronousCluster.dat";
                correctNodesFileName=secondClusterFileName;//"correctCluster.dat";
                numOfClusters=2;
        }
        float calculateMean(string clusterFileName);
        bool isErroneous(double index);
        bool classifier();
        void updateParameter();
        void findIntersectionPoint();

public:
        double intersectionPoint;
        size_t numOfClusters;
        float curErronousClusterMean;
        float curCorrectClusterMean;
        ExpMaxClustering(string FileName, string firstClusterFileName, string secondClusterFileName,float divergenceValue,float initialMeanOfFirstCluster,float initialMeanOfSecondCluster){
                initialization(FileName, firstClusterFileName,secondClusterFileName,divergenceValue,initialMeanOfFirstCluster, initialMeanOfSecondCluster );
        }

        ExpMaxClustering(string fileName){
                initialization(fileName, "erronousNodes.dat","correctNodes.dat",.01,1, 50 );
        }
        ExpMaxClustering(){
                //initialization("scov_001.dat", "erronousNodes.dat","correctNodes.dat",.01,1, 50 );
        }
        void doClassification();
        void findIntersectionPoint(double curErronousClusterMean , double curCorrectClusterMean);
};
