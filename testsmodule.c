#include "spkmeansmodule.c"
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
/*
int main(int argc, char *argv[])
{
    double LINE0[4] ={4,-30,60,-35}, LINE1[4] = {-30,300,-675,420}, LINE2[4] = {60,-675,1620,-1050}, LINE3[4] = {-35,420,-1050,700};
    double LINE0_[4] ={1,0.5302358,0.7561642, 0.3645064}, LINE1_[4] = { 0.5302358, 1.0000000,0.3779162, 0.4705346}, LINE2_[4] = { 0.7561642,0.3779162,1.0000000, 0.4844589}, LINE3_[4] = { 0.3645064, 0.4705346,0.4844589, 1.0000000};
    double LINE01[5] ={1,-1,-2,1,1}, LINE11[5] = {-1,0,1,3,2}, LINE21[5] = {-2,1,3,1,1}, LINE31[5] = {1,3,1,4,0}, LINE41[5] = {1,2,1,0,5};
    double *S[4] ={LINE0,LINE1,LINE2,LINE3};
    double *S2[4] = {LINE0_,LINE1_,LINE2_,LINE3_};
    double *S3[5] = {LINE01,LINE11,LINE21,LINE31,LINE41};
    flowJacobiAlgo(S3,5);
}
*/