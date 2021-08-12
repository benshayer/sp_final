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
void findMaxOffDiagonal(double** matrixA, int* k, int* l, int n)
{  
    double max = 0;
    int i=0, j=0;
    for (i=0;i<n;i++)
    {
        for (j=i+1; j<n;j++)
        {
            if (matrixA[i][j]>max)
            {
                max = matrixA[i][j];
                *k = i;
                *l = j;
            }
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

void calculateRotateValues(double** matrixA, double* c, double* s, int i, int j)
{
    double tetha, signTetha, absTetha, t;
    tetha = (matrixA[j][j]-matrixA[i][i])/2*matrixA[i][j];
    (tetha<0) ? (signTetha = -1) : (signTetha=1);
    (tetha<0) ? (absTetha = -tetha) : (absTetha = tetha);
    t = signTetha / (absTetha + (pow(tetha,2)+1));
    *c = 1/(pow((pow(t,2)+1),0.5));
    *s = t*(*c);
}
void transform(double** matrixA, double** matrixAt, int n)
{
    int i=0, j=0;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            matrixAt[i][j] = matrixA[j][i];
        }
    }
}
void RotateMatrix_helper(double** matrixA, double** matrixAtag,double c, double s,int i, int j, int n)
{
    int r;
    for(int r=0; r<n;r++)
    {
        if(r!=i && r!=j)
        {
            matrixAtag[r][i] = c*matrixA[r][i] - s*matrixA[r][j];
            matrixAtag[i][r] = matrixAtag[r][i];
            matrixAtag[r][j] = c*matrixA[r][j] + s* matrixA[r][i];
            matrixAtag[j][r] = matrixAtag[r][j];
        }
    }

    matrixAtag[i][i] = pow(c,2)*matrixA[i][i] + pow(s,2)*matrixA[j][j] - 2 * s * c * matrixA[i][j];
    matrixAtag[j][j] = pow(s,2) *matrixA[i][i] + pow(c,2) * matrixA[j][j] + 2 * s * c * matrixA[i][j];
    matrixAtag[i][j] = 0;
    matrixAtag[j][i] = 0;
}
void freeMatrix(double** matrix, int n)
{
    int i;
    for (i=0;i<n;i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}
void calcOFFMatrix(double** matrixA,double* offA, int n)
{
    double offMatrix = 0;
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=i+1;j<n;j++)
        {
            offMatrix += pow(2*matrixA[i][j],2);
        }
    }
    *offA = offMatrix;
}

void JacobiAlgorithm(double** matrixA, int n)
{

    double **matrixAtag, **matrixP, **matrixV, c, s,offA , offAtag, epsilon = 0.001;
    int maxElementOffDiagonalI, maxElementOffDiagonalJ,i,j;
    do{
        matrixAtag = callock(n, sizeof(int*));
        assert(matrixAtag!=NULL);
        matrixP = callock(n,sizeof(int*));
        assert(matrixP!=NULL);
        for(i=0;i<n;i++)
        {
            matrixAtag[i] = calloc(n, sizeof(int));
            assert(matrixAtag[i]!=NULL);
            matrixP[i] = calloc(n,sizeof(int));
            assert(matrixP[i]!=NULL);
            matrixP[i][i] = 1;    
        }
        findMaxOffDiagonal(matrixA,&maxElementOffDiagonalI,&maxElementOffDiagonalJ,n);
        calculateRotateValues(matrixA,&c, &s,i,j);
        matrixP[i][i] = c;
        matrixP[j][j] = c;
        matrixP[i][j] = s;
        matrixP[j][i] = -s;
        RotateMatrix_helper(matrixA, matrixAtag, c, s,maxElementOffDiagonalI,maxElementOffDiagonalJ, n);
        calcOFFMatrix(matrixA, &offA, n);
        calcOFFMatrix(matrixAtag, &offAtag, n);
        freeMatrix(matrixA,n);
        matrixA = matrixAtag;

    } while(1);
    
}

int main(int argc, char *argv[])
{
    testWAM();
    testLaplacian();
    return 0;
}

