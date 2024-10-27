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

#ifndef KMEANS_H
#define KMEANS_H

#include <random>
#include "Kmeans.h"

#endif

Kmeans::Kmeans(int ID, int Destination) {

    nodeID = ID;
    DEST = Destination;
    N = 50;
    D = 1;
    numberofrepetitions = 20;
    MaxNumofClusters = 4;
    BestNumofClusters = 0;
    n = 0.05;
    period_size = 0;
    LVQ_on = false;
    weightedKMeans = false;
    //-------------------Added for supporting normalization----------------------------------//
    normalizedvalues = false;
    averageofutils = -1.0;
    stdevofutils = -1.0;
    //--------------------------------------------------------------------------------------//

    Utilities = (double **) malloc(N * sizeof(double *));
    NonNormUtilities = (double **) malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        Utilities[i] = (double *) malloc(D * sizeof(double));
        NonNormUtilities[i] = (double *) malloc(D * sizeof(double));
    }
    weights = (double *) malloc(N * sizeof(double));
    FinalCenter = (double **) malloc(MaxNumofClusters * sizeof(double *));
    for (int i = 0; i < MaxNumofClusters; i++) {
        FinalCenter[i] = (double *) malloc(D * sizeof(double));
    }
    FinalCluster = (int *) malloc(N * sizeof(int));

    for (int i = 0; i < D; i++) {
        for (int j = 0; j < N; j++) {
            Utilities[j][i] = -100.0;
            NonNormUtilities[j][i] = -100.0;
            FinalCluster[j] = -1;
            weights[j] = 1.0;
        }
        for (int j = 0; j < MaxNumofClusters; j++) {
            FinalCenter[j][i] = -100.0;
        }
    }

    elementcnt = 0;
    elementLVQcnt = 0;
    training = true;
}

//-------------------Added for supporting normalization----------------------------------//
Kmeans::Kmeans(int ID, int Destination, bool enablenorm, bool lvqon, int psize, bool wKmeans) {

    nodeID = ID;
    DEST = Destination;
    D = 1;
    //---------------------------------------------Added by giorg---------------------------------------------//
    //to N einai iso me to athroisma(oneHopListSize + twoHopListSize).
    N = 200;
    oneHopListSize = 200;
    twoHopListSize = 0;
    oneHopListCnt = 0;
    twoHopListCnt = 0;
    trainingFlag = false;
    //-------------------------------------------------------------------------------------------------------//

    numberofrepetitions = 20;
    MaxNumofClusters = 4;
    BestNumofClusters = 0;
    normalizedvalues = enablenorm;
    averageofutils = -1.0;
    stdevofutils = -1.0;
    n = 0.05;
    period_size = psize;
    LVQ_on = lvqon;
    weightedKMeans = wKmeans;

    if (enablenorm && (D > 1)) {
        printf("Error: Currently we support normalization only for one dimension.\n");
        exit(1);
    }
    if ((period_size != 0) && LVQ_on) {
        printf("Error: Currently we do not support LQV and periodic update at the same time.\n");
        exit(1);
    }
    if (weightedKMeans && LVQ_on) {
        printf("Error: LVQ and Weighted KMeans can not be implemented concurrently.\n");
        exit(1);
    }

    Utilities = (double **) malloc(N * sizeof(double *));
    NonNormUtilities = (double **) malloc(N * sizeof(double *));

    //---------------------------------------------Added by giorg---------------------------------------------//
    //arxikopoihsh ton 2 pinakon (oneHopUtilities , twoHopUtilities) poy tha apothikeoun ta utilities

    oneHopUtilities = (double *) malloc(oneHopListSize * sizeof(double *));
    if (twoHopListSize != 0) {
        twoHopUtilities = (double *) malloc(twoHopListSize * sizeof(double *));
    }


    //-------------------------------------------------------------------------------------------------------//

    weights = (double *) malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        Utilities[i] = (double *) malloc(D * sizeof(double));
        NonNormUtilities[i] = (double *) malloc(D * sizeof(double));
    }
    FinalCenter = (double **) malloc(MaxNumofClusters * sizeof(double *));
    for (int i = 0; i < MaxNumofClusters; i++) {
        FinalCenter[i] = (double *) malloc(D * sizeof(double));
    }
    FinalCluster = (int *) malloc(N * sizeof(int));


    //---------------------------------------------Added by giorg---------------------------------------------//
    //arxikopoihsh ton 2 pinakon (oneHopUtilities , twoHopUtilities)   me tis times -100
    for (int i = 0; i < oneHopListSize; i++) {
        oneHopUtilities[i] = -100;
    }

    for (int i = 0; i < twoHopListSize; i++) {
        twoHopUtilities[i] = -100;
    }
    //-------------------------------------------------------------------------------------------------------//

    for (int i = 0; i < D; i++) {
        for (int j = 0; j < N; j++) {
            Utilities[j][i] = -100.0;
            NonNormUtilities[j][i] = -100.0;
            FinalCluster[j] = -1;
            weights[j] = 1.0;
        }
        for (int j = 0; j < MaxNumofClusters; j++) {
            FinalCenter[j][i] = -100.0;
        }
    }
    elementcnt = 0;
    elementLVQcnt = 0;
    training = true;
}
//--------------------------------------------------------------------------------------//

Kmeans::~Kmeans() {

/*	for (int i = 0; i < N; ++i) {
		free(Utilities[i]);
		free(NonNormUtilities[i]);
	}
	free(Utilities);
	free(NonNormUtilities);
	free(weights);
	for (int i = 0; i < MaxNumofClusters; ++i) {
		free(FinalCenter[i]);
	}
	free(FinalCenter);
	free(FinalCluster);
*/
    return;
}

int Kmeans::getBestNumofClusters() { return BestNumofClusters; }

int Kmeans::getDest() { return DEST; }

int Kmeans::getNumberOfUtilities() { return elementcnt + elementLVQcnt; }

bool Kmeans::isInTraining() { return training; }

double **Kmeans::getFinalCenter() { return FinalCenter; }


void Kmeans::UpdatewithNewElement(double elementValue) {
    if (LVQ_on) LVQ(elementValue);
    if (period_size != 0) setElementOfUtilities(elementValue);
}

void Kmeans::UpdatewithNewElement(double *elementValue) {
//To be implemented. Update with an element when dimension >1.
}

void Kmeans::setElementOfUtilities(double elementValue) {

    int current_index = elementcnt % N;

    int weightindex;

    if ((training) || (period_size != 0)) {

        NonNormUtilities[current_index][0] = elementValue;

        //weights[current_index]//Set here the weight for this element
        if (weightedKMeans) {
            for (int i = 0; i < N; i++) {
                weightindex = (current_index - i < 0) ? (N + current_index - i) : (current_index - i);
                weights[i] = exp(-(weightindex + 1) / 400.0);
            }
        }

        elementcnt++;

    }

    int kmeans_trigger = elementcnt % N;
    if ((period_size != 0) && (elementcnt > N)) kmeans_trigger = (elementcnt - N) % period_size;

    if (kmeans_trigger == 0) {
        if (checkUtilitiesAreAllTheSame()) {
            elementcnt = elementcnt - 1;


            NonNormUtilities[current_index][0] = -100.0;


            weights[current_index] = 1.0;

        } else {
            training = false;
            KMeansForDifferentClusterNumbers();
        }
    }

}
//---------------------------------------------------------------Added by giorg---------------------------------------------------------------------------------//
//Add two custom methods (CustomSetElementOfUtilities , CustomUpdatewithNewElement) which take 3 data.

void Kmeans::CustomSetElementOfUtilities(double elementValue, double average, double stdDev, int sample) {


    //---------------------------------------------Added by giorg---------------------------------------------//
    //ipologizo thn thesi poy tha balo to neo elementValue sthn lista oneHopUtilities
    int oneHopCurrentIndex = oneHopListCnt % oneHopListSize;
    //--------------------------------------------------------------------------------------------------------//

    int weightindex;

    if ((training) || (period_size != 0)) {

        //---------------------------------------------Added by giorg---------------------------------------------//
        oneHopUtilities[oneHopCurrentIndex] = elementValue;
        //--------------------------------------------------------------------------------------------------------//


        //weights[current_index]//Set here the weight for this element
        if (weightedKMeans) {
            for (int i = 0; i < oneHopListSize; i++) {
                weightindex = (oneHopCurrentIndex - i < 0) ? (oneHopListSize + oneHopCurrentIndex - i) : (
                        oneHopCurrentIndex - i);
                weights[i] = exp(-(weightindex + 1) / 400.0);
            }
        }

        //---------------------------------------------Added by giorg---------------------------------------------//
        oneHopListCnt++;

        //parago ta synthetic data moy kai ta bazo sthn twoHopUtilities lista mou
        if (twoHopListSize != 0) {

            generateRandomValuesFromStats(average, stdDev, sample);
        }

        //--------------------------------------------------------------------------------------------------------//
    }

    int kmeans_trigger = oneHopListCnt % oneHopListSize;
    if ((period_size != 0) && (oneHopListCnt > oneHopListSize))
        kmeans_trigger = (oneHopListCnt - oneHopListSize) % period_size;


    if (kmeans_trigger == 0) {
        if (CustomCheckUtilitiesAreAllTheSame()) {
            oneHopUtilities[oneHopCurrentIndex] = -100;
            oneHopListCnt--;
            elementcnt--;
            weights[oneHopCurrentIndex] = 1.0;

        } else {
            training = false;


            //---------------------------------------------Added by giorg---------------------------------------------//
            // ginontai merge ta 2 list mou(oneHopUtilities , twoHopUtilities) sthn NonNormUtilities.
            for (int i = 0; i < oneHopListSize; i++) {
                NonNormUtilities[i][0] = oneHopUtilities[i];
            }

            if (twoHopListCnt > twoHopListSize) {
                for (int i = 0; i < (twoHopListSize); i++) {
                    NonNormUtilities[i + oneHopListSize][0] = twoHopUtilities[i];
                }

            } else if (twoHopListSize != 0) {
                for (int i = 0; i < (twoHopListCnt); i++) {
                    NonNormUtilities[i + oneHopListSize][0] = twoHopUtilities[i];
                }

            }

            if (oneHopListCnt + twoHopListCnt < N) {
                elementcnt = oneHopListCnt + twoHopListCnt;

            } else {
                elementcnt = N;
            }

            //--------------------------------------------------------------------------------------------------------//
            CustomKMeansForDifferentClusterNumbers();

        }


    }
}


void Kmeans::CustomUpdatewithNewElement(double elementValue, double average, double stdDev, int sample) {

    // klhsh ths lvq  gia to oneHopUtility alla kai gia ta synthetic data poy paragontai me thn xrhsh ton data(average ,stdDev ,sample)
    if (LVQ_on) {
        int oneHopCurrentIndex = oneHopListCnt % oneHopListSize;
        //ananeono to oneHopUtilities etsi oste na paragontai ananaiomena data (average ,stdDev ,sample) akoma kai meta thn training periodo.

        oneHopUtilities[oneHopCurrentIndex] = elementValue;
        oneHopListCnt++;
        LVQ(elementValue);

        if (twoHopListSize != 0) {
            lvqGenerateRandomValuesFromStats(average, stdDev, sample);
        }

    }
    if (period_size != 0) CustomSetElementOfUtilities(elementValue, average, stdDev, sample);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------//

void Kmeans::setElementOfUtilities(double *elementValue) {

//To be implemented. Insert an element when dimension >1.

}

bool Kmeans::checkUtilitiesAreAllTheSame() {

    for (int i = 0; i < N; i++) {
        for (int k = i + 1; k < N; k++) {
            for (int j = 0; j < D; j++) {
                if (NonNormUtilities[i][j] != NonNormUtilities[k][j]) return false;
            }
        }
    }
    return true;
}

//-------------------------------------------------------------Added by giorg-----------------------------------------------------------------------------------//

//metatroph ths checkUtilitiesAreAllTheSame etsi oste na elegxei thn oneHopUtilities lista.

bool Kmeans::CustomCheckUtilitiesAreAllTheSame() {

    for (int i = 0; i < oneHopListSize; i++) {
        for (int k = i + 1; k < oneHopListSize; k++) {

            if (oneHopUtilities[i] != oneHopUtilities[k]) return false;
        }
    }
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------//



int Kmeans::returnClusterRank(double UtilityValue) {

    //Returns the rank of the cluster the 'UtilityValue' belongs to. Valid rank values: 1-BestNumofClusters. The rank of the top-level cluster is 1.
    double min = 999999999.0;
    double distance = 0.0;
    int mincluster = 0;
    int rank = 1;

    //-------------------Added for supporting normalization----------------------------------//
    if (normalizedvalues) {
        UtilityValue = (UtilityValue - averageofutils) / stdevofutils;
    }
    //--------------------------------------------------------------------------------------//

    for (int i = 0; i < this->BestNumofClusters; i++) {
        distance = sqrt(pow(FinalCenter[i][0] - UtilityValue, 2));
        if (min - distance > DBL_EPSILON) {
            min = distance;
            mincluster = i;
        }
    }

    for (int i = 0; i < this->BestNumofClusters; i++) {
        if ((FinalCenter[i][0] - FinalCenter[mincluster][0] > DBL_EPSILON) && (mincluster != i)) {
            rank++;
        }
    }

    return rank;
}


int Kmeans::returnClusterRank(double *UtilityValue) {

    //Returns the rank of the cluster the 'UtilityValue' belongs to. Valid rank values: 1-BestNumofClusters. The rank of the top-level cluster is 1.
    //To be implemented. Return the rank when dimension >1.
    return 0;
}


void Kmeans::LVQ(double UtilityValue) {


    double min = 999999999.0;
    double distance = 0.0;
    int mincluster = 0;
    //-------------------Added for supporting normalization----------------------------------//
    if (normalizedvalues) {
        for (int i = 0; i < this->BestNumofClusters; i++) {
            FinalCenter[i][0] = FinalCenter[i][0] * stdevofutils + averageofutils;
        }
        double currentM2 = stdevofutils * (elementcnt + elementLVQcnt);
        elementLVQcnt++;
        double delta = UtilityValue - averageofutils;
        averageofutils = averageofutils + delta / (elementcnt + elementLVQcnt);
        currentM2 = currentM2 + delta * (UtilityValue - averageofutils);
        stdevofutils = currentM2 / (elementcnt + elementLVQcnt);
        UtilityValue = (UtilityValue - averageofutils) / stdevofutils;
        for (int i = 0; i < this->BestNumofClusters; i++) {
            FinalCenter[i][0] = (FinalCenter[i][0] - averageofutils) / stdevofutils;
        }
    }
    //--------------------------------------------------------------------------------------//



    for (int i = 0; i < this->BestNumofClusters; i++) {
        distance = sqrt(pow(FinalCenter[i][0] - UtilityValue, 2));
        if (min - distance > DBL_EPSILON) {
            min = distance;
            mincluster = i;
        }
    }

    FinalCenter[mincluster][0] += n * (UtilityValue - FinalCenter[mincluster][0]);

}

void Kmeans::LVQ(double *UtilityValue) {

    //To be implemented.
    //Update using LVQ when dimension > 1.

}

void Kmeans::KMeansForDifferentClusterNumbers() {

    double **Center;
    double **tCenter;
    int *tCluster;
    int *Cluster;

    double error;
    double minimumerror = 999999999.9;
    double silhouette;
    double max_silhouette = -999999999.9;

    //-------------------Added for supporting normalization----------------------------------//
    if (normalizedvalues) {
        double sum = 0.0;
        for (int i = 0; i < N; i++) {
            sum = sum + NonNormUtilities[i][0];
        }
        averageofutils = (sum / (double) (N));
        sum = 0.0;
        for (int i = 0; i < N; i++) {
            sum = sum + pow(NonNormUtilities[i][0] - averageofutils, 2);
        }
        stdevofutils = sqrt(sum / (double) (N));

        for (int i = 0; i < N; i++) {
            Utilities[i][0] = (NonNormUtilities[i][0] - averageofutils) / stdevofutils;
        }
    } else {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < D; j++) {
                Utilities[i][j] = NonNormUtilities[i][j];
            }
        }
    }
    //--------------------------------------------------------------------------------------//

    tCenter = (double **) malloc(MaxNumofClusters * sizeof(double *));
    Center = (double **) malloc(MaxNumofClusters * sizeof(double *));
    for (int i = 0; i < MaxNumofClusters; i++) {
        tCenter[i] = (double *) malloc(D * sizeof(double));
        Center[i] = (double *) malloc(D * sizeof(double));
    }
    Cluster = (int *) malloc(N * sizeof(int));
    tCluster = (int *) malloc(N * sizeof(int));

    //Repeat calculation of centers (for different maximum number of clusters) and choose the BestNumofClusters based on Silhouette score
    for (int k = 2; k <= MaxNumofClusters; k++) {

        for (int i = 0; i < MaxNumofClusters; i++) {
            for (int j = 0; j < D; j++) {
                Center[i][j] = -100.0;
            }
        }
        for (int i = 0; i < N; i++) {
            Cluster[i] = -1;
        }
        minimumerror = 999999999.9;

        //Repeat calculation of centers (for various initializations) and choose the one with the minimum error
        for (int rep = 0; rep < numberofrepetitions; rep++) {
            error = KMeansAlgo(tCenter, tCluster, k);
            if (minimumerror - error > DBL_EPSILON) {
                minimumerror = error;
                for (int i = 0; i < k; i++) {
                    for (int j = 0; j < D; j++) {
                        Center[i][j] = tCenter[i][j];
                    }
                }
                for (int i = 0; i < N; i++) {
                    Cluster[i] = tCluster[i];
                }
            }
        }

        if (minimumerror < 999999999.9) {
            //Choose best k based on Silhouette score
            silhouette = CalculateSilhouetteScore(Center, Cluster, k);
            if (silhouette - max_silhouette > DBL_EPSILON) {
                max_silhouette = silhouette;
                for (int i = 0; i < k; i++) {
                    FinalCluster[i] = Cluster[i];
                    for (int j = 0; j < D; j++) {
                        FinalCenter[i][j] = Center[i][j];
                    }
                }
                BestNumofClusters = k;
            }
        }

    }

    if (BestNumofClusters == 0) {
        printf("This should never have happened. At least two clusters should always exist.\n");
        exit(1);
    }

    for (int i = 0; i < MaxNumofClusters; ++i) {
        free(Center[i]);
        free(tCenter[i]);
    }
    free(Center);
    free(tCenter);
    free(Cluster);
    free(tCluster);

}


double Kmeans::KMeansAlgo(double **tempCenters, int *tempClusters, int numofk) {
    int kmeansiterations;
    bool flagchanges;
    double distance;
    double minimumdistance;
    int minimumcluster;
    double error;
    double errorhelper;
    bool successful_initialization = false;

    for (int i = 0; i < MaxNumofClusters; i++) {
        for (int j = 0; j < D; j++) {
            tempCenters[i][j] = -100.0;
        }
    }
    for (int i = 0; i < N; i++) {
        tempClusters[i] = -1;
    }

    successful_initialization = InitializeCenters(tempCenters, numofk);

    if (successful_initialization) {
        //Calculate distances from centers
        kmeansiterations = 0;
        flagchanges = true;
        while (flagchanges) {

            flagchanges = false;

            minimumcluster = -1;
            for (int i = 0; i < N; i++) {
                minimumdistance = 999999999;
                for (int j = 0; j < numofk; j++) {
                    distance = 0;
                    for (int k = 0; k < D; k++) {
                        //weights do not participate in assigning points to clusters
                        distance = distance + pow(tempCenters[j][k] - Utilities[i][k], 2);
                    }
                    distance = sqrt(distance);
                    if (minimumdistance - distance > DBL_EPSILON) {
                        minimumdistance = distance;
                        minimumcluster = j;
                    }
                }
                if (tempClusters[i] != minimumcluster) flagchanges = true;
                tempClusters[i] = minimumcluster;
            }

            //Calculate new centers
            double sumofelements[numofk][D];
            double countofelements[numofk];

            for (int i = 0; i < numofk; i++) {
                countofelements[i] = 0;
                for (int j = 0; j < D; j++) {
                    sumofelements[i][j] = 0;
                }
            }

            for (int i = 0; i < N; i++) {
                //countofelements[tempClusters[i]]++;
                countofelements[tempClusters[i]] = countofelements[tempClusters[i]] + weights[i];
                for (int j = 0; j < D; j++) {
                    sumofelements[tempClusters[i]][j] =
                            sumofelements[tempClusters[i]][j] + weights[i] * Utilities[i][j];
                }
            }

            for (int i = 0; i < numofk; i++) {
                for (int j = 0; j < D; j++) {
                    if (countofelements[i] != 0) {
                        tempCenters[i][j] = sumofelements[i][j] / countofelements[i];
                    } else {
                        tempCenters[i][j] = -100.0;
                    }
                }
            }

            //Calculate error
            error = 0.0;
            for (int i = 0; i < N; i++) {
                errorhelper = 0.0;
                for (int j = 0; j < D; j++) {
                    errorhelper = errorhelper + weights[i] * pow(tempCenters[tempClusters[i]][j] - Utilities[i][j], 2);
                }
                error = error + errorhelper;
            }
            error = sqrt(error);

            kmeansiterations++;
        }
        return error;
    } else {
        return 999999999.91;
    }

}


bool Kmeans::InitializeCenters(double **tempCenters, int numofk) {
    //Random Initialization
    int seed = rand();
    srand(seed);
    int rnd;

    int slidingrnd;
    bool center_is_different;
    int compared_dimensions;

    for (int i = 0; i < numofk; i++) {

        rnd = rand() % N;
        slidingrnd = rnd;

        do {
            center_is_different = true;
            for (int x = 0; x < i; x++) {
                compared_dimensions = 0;
                for (int y = 0; y < D; y++) {
                    if (abs(Utilities[slidingrnd][y] - tempCenters[x][y]) <= 3.0 * DBL_EPSILON) compared_dimensions++;
                }
                if (compared_dimensions == D) {
                    center_is_different = false;
                    slidingrnd = (slidingrnd + 1) % N;
                    break;
                }
            }
        } while ((slidingrnd != rnd) && (!center_is_different));

        if (center_is_different) {
            for (int j = 0; j < D; j++) {
                tempCenters[i][j] = Utilities[slidingrnd][j];
            }
        } else {
            return false;
        }

    }

    return true;

}

double Kmeans::CalculateSilhouetteScore(double **ttCenters, int *ttClusters, int kvalue) {

    double silh;
    double si;
    double ai;
    double bi;
    int currentCluster;
    int lowAvgDisimPosition;
    double minDisim;
    double max;

    double distsample;

    double *SimDis = (double *) malloc(kvalue * sizeof(double));
    int *SimDisCounter = (int *) malloc(kvalue * sizeof(int));


    si = 0.0;
    silh = 0.0;
    ai = 0.0;
    bi = 0.0;
    currentCluster = 0;
    max = 0.0;

    for (int i = 0; i < N; i++) {
        currentCluster = ttClusters[i];

        for (int a = 0; a < kvalue; a++) {
            SimDis[a] = 0.0;
            SimDisCounter[a] = 0;
        }

        for (int l = 0; l < N; l++) {
            distsample = 0.0;
            for (int j = 0; j < D; j++) {
                distsample += pow(Utilities[i][j] - Utilities[l][j], 2);
            }

            distsample = sqrt(distsample);

            if (currentCluster == ttClusters[l]) {
                SimDis[currentCluster] = SimDis[currentCluster] + distsample;
                SimDisCounter[currentCluster]++;
            } else {
                SimDisCounter[ttClusters[l]]++;
                SimDis[ttClusters[l]] = SimDis[ttClusters[l]] + distsample;
            }
        }

        for (int c = 0; c < kvalue; c++) {
            SimDis[c] = sqrt(SimDis[c]) / SimDisCounter[c];
        }

        ai = SimDis[currentCluster];

        lowAvgDisimPosition = -1;
        minDisim = 999999999.9;
        for (int b = 0; b < kvalue; b++) {
            if (b != currentCluster) {
                if (SimDis[b] < minDisim) {
                    lowAvgDisimPosition = b;
                    minDisim = SimDis[b];
                }
            }
        }
        bi = SimDis[lowAvgDisimPosition];
        max = (ai < bi) ? bi : ai;
        si = (bi - ai) / max;
        silh += si;
    }

    free(SimDis);
    free(SimDisCounter);

    return silh / N;
}

//---------------------------------------------Added by giorg---------------------------------------------//


//Metatroph ton 4 methodon etsi oste na diatrexoyn mono tis theseis poy exoyn dextei timi apo thn lista NonNormUtilities


void Kmeans::CustomKMeansForDifferentClusterNumbers() {

    double **Center;
    double **tCenter;
    int *tCluster;
    int *Cluster;

    double error;
    double minimumerror = 999999999.9;
    double silhouette;
    double max_silhouette = -999999999.9;

    //-------------------Added for supporting normalization----------------------------------//
    if (normalizedvalues) {
        double sum = 0.0;
        for (int i = 0; i < N; i++) {
            sum = sum + NonNormUtilities[i][0];
        }
        averageofutils = (sum / (double) (N));
        sum = 0.0;
        for (int i = 0; i < N; i++) {
            sum = sum + pow(NonNormUtilities[i][0] - averageofutils, 2);
        }
        stdevofutils = sqrt(sum / (double) (N));

        for (int i = 0; i < N; i++) {
            Utilities[i][0] = (NonNormUtilities[i][0] - averageofutils) / stdevofutils;
        }
    } else {
        for (int i = 0; i < elementcnt; i++) {
            for (int j = 0; j < D; j++) {
                Utilities[i][j] = NonNormUtilities[i][j];
            }
        }
    }
    //--------------------------------------------------------------------------------------//

    tCenter = (double **) malloc(MaxNumofClusters * sizeof(double *));
    Center = (double **) malloc(MaxNumofClusters * sizeof(double *));
    for (int i = 0; i < MaxNumofClusters; i++) {
        tCenter[i] = (double *) malloc(D * sizeof(double));
        Center[i] = (double *) malloc(D * sizeof(double));
    }
    Cluster = (int *) malloc(elementcnt * sizeof(int));
    tCluster = (int *) malloc(elementcnt * sizeof(int));

    //Repeat calculation of centers (for different maximum number of clusters) and choose the BestNumofClusters based on Silhouette score
    for (int k = 2; k <= MaxNumofClusters; k++) {

        for (int i = 0; i < MaxNumofClusters; i++) {
            for (int j = 0; j < D; j++) {
                Center[i][j] = -100.0;
            }
        }
        for (int i = 0; i < elementcnt; i++) {
            Cluster[i] = -1;
        }
        minimumerror = 999999999.9;

        //Repeat calculation of centers (for various initializations) and choose the one with the minimum error
        for (int rep = 0; rep < numberofrepetitions; rep++) {
            error = CustomKMeansAlgo(tCenter, tCluster, k);
            if (minimumerror - error > DBL_EPSILON) {
                minimumerror = error;
                for (int i = 0; i < k; i++) {
                    for (int j = 0; j < D; j++) {
                        Center[i][j] = tCenter[i][j];
                    }
                }
                for (int i = 0; i < elementcnt; i++) {
                    Cluster[i] = tCluster[i];
                }
            }
        }

        if (minimumerror < 999999999.9) {
            //Choose best k based on Silhouette score
            silhouette = CustomCalculateSilhouetteScore(Center, Cluster, k);
            if (silhouette - max_silhouette > DBL_EPSILON) {
                max_silhouette = silhouette;
                for (int i = 0; i < k; i++) {
                    FinalCluster[i] = Cluster[i];
                    for (int j = 0; j < D; j++) {
                        FinalCenter[i][j] = Center[i][j];
                    }
                }
                BestNumofClusters = k;
            }
        }

    }

    if (BestNumofClusters == 0) {
        printf("This should never have happened. At least two clusters should always exist.\n");
        exit(1);
    }

    for (int i = 0; i < MaxNumofClusters; ++i) {
        free(Center[i]);
        free(tCenter[i]);
    }
    free(Center);
    free(tCenter);
    free(Cluster);
    free(tCluster);

}


double Kmeans::CustomKMeansAlgo(double **tempCenters, int *tempClusters, int numofk) {
    int kmeansiterations;
    bool flagchanges;
    double distance;
    double minimumdistance;
    int minimumcluster;
    double error;
    double errorhelper;
    bool successful_initialization = false;

    for (int i = 0; i < MaxNumofClusters; i++) {
        for (int j = 0; j < D; j++) {
            tempCenters[i][j] = -100.0;
        }
    }
    for (int i = 0; i < elementcnt; i++) {
        tempClusters[i] = -1;
    }

    successful_initialization = CustomInitializeCenters(tempCenters, numofk);

    if (successful_initialization) {
        //Calculate distances from centers
        kmeansiterations = 0;
        flagchanges = true;
        while (flagchanges) {

            flagchanges = false;

            minimumcluster = -1;
            for (int i = 0; i < elementcnt; i++) {
                minimumdistance = 999999999;
                for (int j = 0; j < numofk; j++) {
                    distance = 0;
                    for (int k = 0; k < D; k++) {
                        //weights do not participate in assigning points to clusters
                        distance = distance + pow(tempCenters[j][k] - Utilities[i][k], 2);
                    }
                    distance = sqrt(distance);
                    if (minimumdistance - distance > DBL_EPSILON) {
                        minimumdistance = distance;
                        minimumcluster = j;
                    }
                }
                if (tempClusters[i] != minimumcluster) flagchanges = true;
                tempClusters[i] = minimumcluster;
            }

            //Calculate new centers
            double sumofelements[numofk][D];
            double countofelements[numofk];

            for (int i = 0; i < numofk; i++) {
                countofelements[i] = 0;
                for (int j = 0; j < D; j++) {
                    sumofelements[i][j] = 0;
                }
            }

            for (int i = 0; i < elementcnt; i++) {
                //countofelements[tempClusters[i]]++;
                countofelements[tempClusters[i]] = countofelements[tempClusters[i]] + weights[i];
                for (int j = 0; j < D; j++) {
                    sumofelements[tempClusters[i]][j] =
                            sumofelements[tempClusters[i]][j] + weights[i] * Utilities[i][j];
                }
            }

            for (int i = 0; i < numofk; i++) {
                for (int j = 0; j < D; j++) {
                    if (countofelements[i] != 0) {
                        tempCenters[i][j] = sumofelements[i][j] / countofelements[i];
                    } else {
                        tempCenters[i][j] = -100.0;
                    }
                }
            }

            //Calculate error
            error = 0.0;
            for (int i = 0; i < elementcnt; i++) {
                errorhelper = 0.0;
                for (int j = 0; j < D; j++) {
                    errorhelper = errorhelper + weights[i] * pow(tempCenters[tempClusters[i]][j] - Utilities[i][j], 2);
                }
                error = error + errorhelper;
            }
            error = sqrt(error);

            kmeansiterations++;
        }
        return error;
    } else {
        return 999999999.91;
    }

}

bool Kmeans::CustomInitializeCenters(double **tempCenters, int numofk) {
    //Random Initialization
    int seed = rand();
    srand(seed);
    int rnd;

    int slidingrnd;
    bool center_is_different;
    int compared_dimensions;

    for (int i = 0; i < numofk; i++) {

        rnd = rand() % elementcnt;
        slidingrnd = rnd;

        do {
            center_is_different = true;
            for (int x = 0; x < i; x++) {
                compared_dimensions = 0;
                for (int y = 0; y < D; y++) {
                    if (abs(Utilities[slidingrnd][y] - tempCenters[x][y]) <= 3.0 * DBL_EPSILON) compared_dimensions++;
                }
                if (compared_dimensions == D) {
                    center_is_different = false;
                    slidingrnd = (slidingrnd + 1) % elementcnt;
                    break;
                }
            }
        } while ((slidingrnd != rnd) && (!center_is_different));

        if (center_is_different) {
            for (int j = 0; j < D; j++) {
                tempCenters[i][j] = Utilities[slidingrnd][j];
            }
        } else {
            return false;
        }

    }

    return true;

}


double Kmeans::CustomCalculateSilhouetteScore(double **ttCenters, int *ttClusters, int kvalue) {

    double silh;
    double si;
    double ai;
    double bi;
    int currentCluster;
    int lowAvgDisimPosition;
    double minDisim;
    double max;

    double distsample;

    double *SimDis = (double *) malloc(kvalue * sizeof(double));
    int *SimDisCounter = (int *) malloc(kvalue * sizeof(int));


    si = 0.0;
    silh = 0.0;
    ai = 0.0;
    bi = 0.0;
    currentCluster = 0;
    max = 0.0;

    for (int i = 0; i < elementcnt; i++) {
        currentCluster = ttClusters[i];

        for (int a = 0; a < kvalue; a++) {
            SimDis[a] = 0.0;
            SimDisCounter[a] = 0;
        }

        for (int l = 0; l < elementcnt; l++) {
            distsample = 0.0;
            for (int j = 0; j < D; j++) {
                distsample += pow(Utilities[i][j] - Utilities[l][j], 2);
            }

            distsample = sqrt(distsample);

            if (currentCluster == ttClusters[l]) {
                SimDis[currentCluster] = SimDis[currentCluster] + distsample;
                SimDisCounter[currentCluster]++;
            } else {
                SimDisCounter[ttClusters[l]]++;
                SimDis[ttClusters[l]] = SimDis[ttClusters[l]] + distsample;
            }
        }

        for (int c = 0; c < kvalue; c++) {
            SimDis[c] = sqrt(SimDis[c]) / SimDisCounter[c];
        }

        ai = SimDis[currentCluster];

        lowAvgDisimPosition = -1;
        minDisim = 999999999.9;
        for (int b = 0; b < kvalue; b++) {
            if (b != currentCluster) {
                if (SimDis[b] < minDisim) {
                    lowAvgDisimPosition = b;
                    minDisim = SimDis[b];
                }
            }
        }
        bi = SimDis[lowAvgDisimPosition];
        max = (ai < bi) ? bi : ai;
        si = (bi - ai) / max;
        silh += si;
    }

    free(SimDis);
    free(SimDisCounter);

    return silh / elementcnt;
}





//Dhmiourgo 3 methodous oi opoies bgazoyn ta 3 data mas me thn xrhsh ths oneHopUtilities lista mas.

double Kmeans::average() {
    double sum = 0.0;
    int total = getElementcnt();

    for (int i = 0; i < total; i++) {
        sum += oneHopUtilities[i];
    }

    avr = sum / total;

    return avr;
}

double Kmeans::stdev() {

    double squaredDifferencesSum = 0.0;
    double difference = 0.0;

    int total = getElementcnt();

    for (int i = 0; i < total; i++) {
        difference = oneHopUtilities[i] - avr;
        squaredDifferencesSum += difference * difference;

    }

    // Calculate the variance (average of squared differences)
    double variance = squaredDifferencesSum / total;

    // Calculate the standard deviation
    double standardDeviation = sqrt(variance);

    return standardDeviation;
}

int Kmeans::getElementcnt() {
    if (oneHopListCnt <= oneHopListSize) {
        return oneHopListCnt;
    } else {
        return oneHopListSize;
    }
}


//Dhmiourgo methodo h opoia me thn xrhsh ton data mas(avr , stdDev , sampleNumber) paragei ta synthetic data mas ta opoia apothikeuontai sthn lista twoHopUtilities.

void Kmeans::generateRandomValuesFromStats(double average, double standardDeviation, int sampleNumber) {
    double syntheticData;
    int current_index;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(average, standardDeviation);


    for (int i = 0; i < sampleNumber; i++) {
        current_index = twoHopListCnt % twoHopListSize;

        do {
            syntheticData = distribution(gen);
        } while (syntheticData < 0);

        twoHopUtilities[current_index] = syntheticData;
        twoHopListCnt++;

    }

}


//methodos me thn opoia kano klhsh ths methodou  lvq gia kathe synthetic data poy dhmioyrgo.
void Kmeans::lvqGenerateRandomValuesFromStats(double average, double standardDeviation, int sampleNumber) {

    double syntheticData;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(average, standardDeviation);

    for (int i = 0; i < sampleNumber; i++) {

        do {
            syntheticData = distribution(gen);
        } while (syntheticData < 0);

        LVQ(syntheticData);

    }
}

//-----------------------------------------------------------------------------------------------------//