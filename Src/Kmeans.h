/*
 *  Adyton: A Network Simulator for Opportunistic Networks
 *  Copyright (C) 2015  Nikolaos Papanikos, Dimitrios-Georgios Akestoridis,
 *  and Evangelos Papapetrou
 *
 *  This file is part of Adyton.
 *
 *  Adyton is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Adyton is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Adyton.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Written by Evangelos Papapetrou.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <float.h>
#include <vector>

using namespace std;

//---------------------------------------------Added by giorg---------------------------------------------//
struct RandomValues{
    int id;
    std::vector<double> values;
};
//--------------------------------------------------------------------------------------//

class Kmeans
{
private:

    int nodeID; //id of node maintaining the structure
    int DEST; //id of the destination for which structure is maintained

    int N; //maximum number of elements
    int D; //dimension
    int numberofrepetitions; //the number of repetitions for minimizing the effect of random initialization
    int MaxNumofClusters; //Maximum number of clusters
    int BestNumofClusters; //the best number of clusters




    //-------------------Added for supporting normalization----------------------------------//
    bool normalizedvalues;
    double averageofutils;
    double stdevofutils;
    //--------------------------------------------------------------------------------------//
    double n; //weight value for LVQ
    int period_size;
    bool LVQ_on;
    bool weightedKMeans;

    double **Utilities; // Utilities[N][D]: array used for storing the elements
    double **NonNormUtilities; // NonNormUtilities[N][D]: array used for storing the elements in non-normalized mode
    double *weights; //weights[N]: aaray used for storing the weight for each element - used in weightedKMeans
    double **FinalCenter; //FinalCenter[MaxNumofClusters][D]: array used for storing the coordinates of the cluster centers (note: MaxNumofClusters->maximum number of clusters)
    int *FinalCluster; //FinalCluster[N]: array used for storing the index of the cluster in which an element belongs to


    //---------------------------------------------Added by giorg---------------------------------------------//
    double avr;
    vector<RandomValues> dataArray;
    double *oneHopUtilities;
    double *twoHopUtilities;
    int oneHopListSize;
    int twoHopListSize;
    int oneHopListCnt;
    int twoHopListCnt;
    bool trainingFlag;
    //--------------------------------------------------------------------------------------//


    int elementcnt; //current number of stored elements
    int elementLVQcnt; //current number of elements considered in LVQ

    bool training; //indicates whether we are still in the training phase

    void KMeansForDifferentClusterNumbers();
    double KMeansAlgo(double **tempCenters, int *tempClusters, int numofk);
    bool InitializeCenters(double **tempCenters, int numofk);
    double CalculateSilhouetteScore(double **ttCenters, int *ttClusters, int kvalue);

public:
    Kmeans(int NID,int Destination);
    //-------------------Added for supporting normalization----------------------------------//
    Kmeans(int ID, int Destination, bool enablenorm, bool lvqon, int psize, bool wKmeans);
    //--------------------------------------------------------------------------------------//
    ~Kmeans();
    int getBestNumofClusters();
    int getDest();
    int getNumberOfUtilities();
    bool isInTraining();
    double** getFinalCenter();

    void UpdatewithNewElement(double elementValue);
    void UpdatewithNewElement(double *elementValue);
    void setElementOfUtilities(double elementValue);
    void setElementOfUtilities(double *elementValue);
    bool checkUtilitiesAreAllTheSame();
    int returnClusterRank(double UtilityValue);
    int returnClusterRank(double *UtilityValue);
    void  LVQ(double UtilityValue);
    void  LVQ(double *UtilityValue);



    //---------------------------------------------Added by giorg---------------------------------------------//

    void CustomSetElementOfUtilities(double elementValue,double average,double stdDev,int sample);
    void CustomUpdatewithNewElement(double elementValue, double average, double stdDev, int sample);
    bool CustomCheckUtilitiesAreAllTheSame() ;

    void CustomKMeansForDifferentClusterNumbers();
    double CustomKMeansAlgo(double **tempCenters, int *tempClusters, int numofk);
    bool CustomInitializeCenters(double **tempCenters, int numofk);
    double CustomCalculateSilhouetteScore(double **ttCenters, int *ttClusters, int kvalue);
    double average();

    double stdev();

    int getElementcnt() ;

    void generateRandomValuesFromStats(double average, double standardDeviation, int sampleNumber);
    void lvqGenerateRandomValuesFromStats(double average, double standardDeviation, int sampleNumber);
    //---------------------------------------------------------------------------------------------------------//
};
