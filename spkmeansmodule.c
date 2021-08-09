#include <math.h>
#include <stdlib.h>

double calcWeight(double* observ1, double* observ2, int dim)
{
    double sum = 0, currDiff, norm,weight;
    int i;
    for (i=0;i<dim;i++)
    {
        currDiff = (observ1[i]-observ2[i])*(observ1[i]-observ2[i]);
        sum = sum + currDiff;
    }
    norm = pow(sum,0.5);
    norm = -norm/2;
    return exp(norm);
}

double** CreateWeightedAdjacencyMatrix(double** observations, int dim, int n)
{
    double** wam;
    double currWeight;
    int i,j;
    wam = calloc(n,sizeof(int*));
    for(i=0;i<n;i++)
    {
        wam[i] = calloc(n,sizeof(int));
    }
    for(i=0;i<n;i++)
    {
        for(j=i+1;j<n;j++)
        {
            currWeight = calcWeight(observations[i],observations[j],dim);
            wam[i][j] = currWeight;
            wam[j][i] = currWeight;
        }
    }
    return wam;
}

int main(int argc, char *argv[])
{
    double** matrix;
    double observ1[3] = {1,2,3};
    double observ2[3] = {2,3,4};
    double observ3[3] = {4,5,6};
    double* observations[3] = {observ1, observ2, observ3};
    matrix = CreateWeightedAdjacencyMatrix(observations,3,3);

    return 0;
}

