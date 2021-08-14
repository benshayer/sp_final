#include <math.h>
#include <stdlib.h>
#include "testsmodule.c"
#include <assert.h>
#include <stdio.h>

/* test */

typedef struct EignValue{
    double value;
    int index;
}EignValue;

int cmpfunc (const void * a, const void * b) {
   EignValue firstE = *(EignValue*)a;
   EignValue secondE = *(EignValue*)b;
   if (firstE.value > secondE.value){
       return 1;
   }
   if (secondE.value > firstE.value) {
       return -1;
   }
   else {
       if (firstE.index <= secondE.index) {
           return -1;
       }
       else {
           return 1;
       }
   }
}

int * initDataPoints(char* filename, double ***data_vectors)
{
    FILE *ifp = NULL;
    ifp = fopen(filename,"r");
    double *vector;
    double current_value;
    int d = 2;
    int n = 10;
    int i = 0;
    int j = 0;
    int *value;
    char c;
    vector = (double *)calloc(d, sizeof(double));
    assert(vector != NULL);
    while (fscanf(ifp,"%lf%c", &current_value, &c) == 2)
    {
        if (i == (d - 1))
        {
            d *= 2;
            vector = (double *)realloc(vector, d * sizeof(double));
            assert(vector != NULL);
        }
        vector[i] = current_value;
        i++;
        if (c == '\n')
        {
            vector = (double *)realloc(vector, (i) * sizeof(double));
            assert(vector != NULL);
            if (j == (n - 1))
            {
                n *= 10;
                *data_vectors = (double **)realloc(*data_vectors, n * sizeof(double *));
                assert(*data_vectors != NULL);
            }

            (*data_vectors)[j] = vector;
            j++;
            vector = (double *)calloc((i), sizeof(double));
            assert(vector != NULL);
            d = i;
            i = 0;
        } 
    }
    free(vector);
    *data_vectors = (double **)realloc(*data_vectors, (j) * sizeof(double *));
    assert(*data_vectors != NULL);
    value = (int *)calloc(2, sizeof(int));
    assert(value != NULL);
    value[0] = j;
    value[1] = d;
    return value;
}

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
    wam = calloc(n,sizeof(double*));
    for(i=0;i<n;i++)
    {
        wam[i] = calloc(n,sizeof(double));
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

double** DiagonalDegreeMatrix(double** matrix, int n)
{
    double** ddm;
    double sumOfRow;
    int i,j;
    ddm = calloc(n,sizeof(double*));
    for(i=0;i<n;i++)
    {
        ddm[i] = calloc(n,sizeof(double));
    }
    for (i=0;i<n;i++){
        sumOfRow = 0;
        for (j=0;j<n;j++){
            sumOfRow += matrix[i][j];
        }
        ddm[i][i] = 1/sqrt(sumOfRow);
    }
    return ddm;
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
    double max = 0, absCurrElem = 0;
    int i=0, j=0;
    for (i=0;i<n;i++)
    {
        for (j=i+1; j<n;j++)
        {
            absCurrElem = pow(pow(matrixA[i][j],2),0.5);
            if (absCurrElem>max)
            {
                max = absCurrElem;
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
    Lnorm = calloc(n,sizeof(double*));
    midMatrix = calloc(n,sizeof(double*));
    for(i=0;i<n;i++)
    {
        Lnorm[i] = calloc(n,sizeof(double));
        midMatrix[i] = calloc(n,sizeof(double));
    }
    MatrixMultiply_helper(ddm_square,wam,midMatrix,n);
    MatrixMultiply_helper(midMatrix,ddm_square,Lnorm,n);
    LaplacianNorm_helper(Lnorm, n);
    return Lnorm;

}

void calculateRotateValues(double** matrixA, double* c, double* s, int i, int j)
{
    double tetha, signTetha, absTetha, t;
    tetha = (matrixA[j][j]-matrixA[i][i])/(2*matrixA[i][j]);
    (tetha<0) ? (signTetha = -1) : (signTetha=1);
    (tetha<0) ? (absTetha = -tetha) : (absTetha = tetha);
    t = signTetha / (absTetha + pow((pow(tetha,2)+1),0.5));
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
void RotateMatrix_helper(double** matrixA,double c, double s,int i, int j, int n)
{   
    double *Icol, *Jcol;
    Icol = calloc(n,sizeof(double));
    Jcol = calloc(n,sizeof(double));
    int r;
    for (r=0;r<n;r++)
    {
        Icol[r] = matrixA[r][i];
        Jcol[r] = matrixA[r][j];
    }
    for(int r=0; r<n;r++)
    {
        if(r!=i && r!=j)
        {
            matrixA[r][i] = c*Icol[r] - s*Jcol[r];
            matrixA[i][r] = matrixA[r][i];
            matrixA[r][j] = c*Jcol[r] + s* Icol[r];
            matrixA[j][r] = matrixA[r][j];
        }
    }

    matrixA[i][i] = pow(c,2)*Icol[i] + pow(s,2)*Jcol[j] - 2 * s * c * Jcol[i];
    matrixA[j][j] = pow(s,2) *Icol[i] + pow(c,2) * Jcol[j] + 2 * s * c * Jcol[i];
    matrixA[i][j] = 0;
    matrixA[j][i] = 0;
    free(Icol);
    free(Jcol);
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
            offMatrix += 2 * pow(matrixA[i][j],2);
        }
    }
    *offA = offMatrix;
}

void copyMatrix(double** matrixSource, double** matrixDest, int n)
{
    int i=0,j=0;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            matrixDest[i][j] = matrixSource[i][j];
        }
    }
}

void JacobiAlgorithm(double** matrixA,double** matrixV, int n)
{

    double **matrixAtag, **matrixP, **matrixNewV, c, s,offA , offAtag, epsilon = 0.001;
    int maxElementOffDiagonalI, maxElementOffDiagonalJ,i,j, stopCondition=1;
    EignValue *eignValues;
    for (i=0; i<n;i++)
    {
        matrixV[i][i] = 1; //Init V to be identity matrix
        for(j=0;j<n;j++)
        {
            if(i!=j)
            {
                matrixV[i][j] =0;
            }
        }
    }
    do{
        calcOFFMatrix(matrixA, &offA, n);
        matrixNewV = calloc(n,sizeof(double*));
        assert(matrixNewV!=NULL);
        matrixP = calloc(n,sizeof(double*));
        assert(matrixP!=NULL);
        for(i=0;i<n;i++)
        {
            matrixNewV[i] = calloc(n, sizeof(double));
            assert(matrixNewV[i]!=NULL);
            matrixP[i] = calloc(n,sizeof(double));
            assert(matrixP[i]!=NULL);
            matrixP[i][i] = 1;
            for(j=0;j<n;j++)
            {
                if(i!=j)
                {
                    matrixP[i][j] =0;
                }
            }    
        }
        findMaxOffDiagonal(matrixA,&maxElementOffDiagonalI,&maxElementOffDiagonalJ,n);
        calculateRotateValues(matrixA,&c, &s,maxElementOffDiagonalI,maxElementOffDiagonalJ);
        matrixP[maxElementOffDiagonalI][maxElementOffDiagonalI] = c;
        matrixP[maxElementOffDiagonalJ][maxElementOffDiagonalJ] = c;
        matrixP[maxElementOffDiagonalI][maxElementOffDiagonalJ] = s;
        matrixP[maxElementOffDiagonalJ][maxElementOffDiagonalI] = -s; //Fill relevant values of P
        RotateMatrix_helper(matrixA,c, s,maxElementOffDiagonalI,maxElementOffDiagonalJ, n); //Calculate A' matrix
        MatrixMultiply_helper(matrixV,matrixP,matrixNewV,n); //Get the current matrix V
        calcOFFMatrix(matrixA, &offAtag, n);
        copyMatrix(matrixNewV,matrixV,n);
        if ((offA-offAtag)<epsilon)
        {
            stopCondition = 0;
        }
        freeMatrix(matrixP,n);
        freeMatrix(matrixNewV,n);
        if (offAtag==0) {break;}
    } while(stopCondition);
}

int TheEigengapHeuristic(double* eigenvalues, int len) {
    double deltaI = 0;
    double currMax = 0;
    int position=0;
    int i,j;
    for(i=1; i<=(len/2);i++){
        deltaI = fabs(eigenvalues[i-1]-eigenvalues[i]);
        if (deltaI > currMax){
            currMax = deltaI;
            position = i-1;
        }
    }
    return position;
}


void getMatrixSortedEignVectors(double** matrixA, double** matrixV, double** matrixU, int n, int k)
{
    EignValue* eignValues;
    int i,j, indexElem;
    eignValues = calloc(n,sizeof(EignValue));
    for (i=0; i<n; i++){
        eignValues[i].value = matrixA[i][i];
        eignValues[i].index = i;
    }
    qsort(eignValues,n, sizeof(EignValue),cmpfunc);
    for(j=0;j<k;j++)
    {
        indexElem = eignValues[j].index;
        for(i=0;i<n;i++)
        {
            matrixU[i][j] = matrixA[i][indexElem];
        }
    }
}
void normalizedMatrixUtoMatrixT(double** matrixU,int n, int k)
{
    int i,j;
    double sumRow;
    for(i=0;i<n;i++)
    {   
        sumRow=0;
        for(j=0;j<k;j++)
        {
            sumRow+=pow(matrixU[i][j],2);
        }
        for(j=0;j<k;j++)
        {
            matrixU[i][j] = matrixU[i][j]/pow(sumRow,0.5);
        }
    }
}

int main(int argc, char *argv[])
{   
    int* values;
    int n;
    double** data_vectors;
    n=10;
    data_vectors = (double **)calloc(n, sizeof(double *));
    initDataPoints("/Users/giladf/Desktop/Semester B/sp_final/tests/input_1.txt",&data_vectors);

    return 0;
}

