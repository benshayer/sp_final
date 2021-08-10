#include <math.h>
#include <stdlib.h>
#include "testsmodule.c"

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
void MatrixMultiply_helper(double** matrixA, double** matrixB, double** result, int n)
{
    int i,j,k;
    double currItem;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            currItem = 0;
            for(k=0;k<n;k++)
            {
                currItem+=matrixA[i][k] * matrixB[k][j];
            }
            result[i][j] = currItem;
        }
    }
}
void LaplacianNorm_helper(double** Lnorm, int n)
{
    int i=0,j=0;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if (i==j)
            {
                Lnorm[i][i] = 1 - Lnorm[i][i];
            }
            else
            {
                Lnorm[i][j] = -Lnorm[i][j];
            }
        }
    }
}
double** ComputeNormalizedGraphLaplacian(double** wam, double** ddm_square, int n)
{
    double** Lnorm, **midMatrix;
    int i,j;
    Lnorm = calloc(n,sizeof(int*));
    midMatrix = calloc(n,sizeof(int*));
    for(i=0;i<n;i++)
    {
        Lnorm[i] = calloc(n,sizeof(int));
        midMatrix[i] = calloc(n,sizeof(int));
    }
    MatrixMultiply_helper(ddm_square,wam,midMatrix,n);
    MatrixMultiply_helper(midMatrix,ddm_square,Lnorm,n);
    LaplacianNorm_helper(Lnorm, n);
    return Lnorm;

}

int main(int argc, char *argv[])
{
    testWAM();
    testLaplacian();
    return 0;
}

