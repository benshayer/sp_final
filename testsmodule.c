#include "spkmeansmodule.h"
int testWAM()
{
    double** matrix;
    double observ1[3] = {1,2,3};
    double observ2[3] = {2,3,4};
    double observ3[3] = {4,5,6};
    double* observations[3] = {observ1, observ2, observ3};
    matrix = CreateWeightedAdjacencyMatrix(observations,3,3);
    return 0;
}
int testLaplacian()
{
    
    double WAMLine0[3] = {0, 4, 3}, WAMLine1[3] = {4,0,5}, WAMLine2[3] = {3,5,0};
    double DDMLine0[3] = {1/pow(5,0.5),0,0}, DDMLine1[3] = {0,1/pow(6,0.5),0}, DDMLine2[3] = {0,0,1/pow(8,0.5)}; 
    double *matrixWAM[3] = {WAMLine0,WAMLine1,WAMLine2}, *matrixDDM[3] = {DDMLine0, DDMLine1, DDMLine2}, **Laplacian;
    Laplacian = ComputeNormalizedGraphLaplacian(matrixWAM,matrixDDM,3);
}
