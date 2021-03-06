#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "spkmeans.h"
#include "kmeans.c"


typedef struct EignValue{
    double value;
    int index;
}EignValue;


int convertStringIntoGoalEnum(char* UserGoal)
{
    if (!strcmp(UserGoal,"spk"))
    {
        return 0;
    } else if (!strcmp(UserGoal,"wam")){
        return 1;
    } else if (!strcmp(UserGoal, "ddg")){
        return 2;
    } else if (!strcmp(UserGoal,"lnorm")){
        return 3;
    } else if (!strcmp(UserGoal,"jacobi")){
        return 4;
    } else 
    {
        return 5;
    };
}

int cmpfunc (const void * a, const void * b) { /*Comperator for qsort of EignValues struct method*/
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

int doubleCmpr (const void *a, const void * b) /*Comperator for qsort of EignValues from type double method*/
{ 
    double first = *(double*)a;
    double second = *(double*)b;
    if (first > second)
    {
        return 1;
    }
    else {
        if (first <second) 
        {
            return -1;
        }
        else {
            return 0;
        }
    }
}

int * initDataPoints(char* filename, double ***data_vectors)
{
    FILE *ifp = NULL;
    double *vector;
    double current_value;
    int *value;
    char c;
    int d = 2;
    int n = 10;
    int i = 0;
    int j = 0;
    ifp = fopen(filename,"r");
    vector = (double *)calloc(d, sizeof(double));
    assert(vector != NULL && "An Error Has Occured");
    while (fscanf(ifp,"%lf%c", &current_value, &c) >= 1)
    {
        if (i == (d - 1))
        {
            d *= 2;
            vector = (double *)realloc(vector, d * sizeof(double));
            assert(vector != NULL && "An Error Has Occured");
        }
        vector[i] = current_value;
        i++;
        if (c == '\n' || c=='\r')
        {
            vector = (double *)realloc(vector, (i) * sizeof(double));
            assert(vector != NULL && "An Error Has Occured");
            if (j == (n - 1))
            {
                n *= 10;
                *data_vectors = (double **)realloc(*data_vectors, n * sizeof(double *));
                assert(*data_vectors != NULL && "An Error Has Occured");
            }

            (*data_vectors)[j] = vector;
            j++;
            vector = (double *)calloc((i), sizeof(double));
            assert(vector != NULL && "An Error Has Occured");
            d = i;
            i = 0;
        } 
    }
    if (j == (n - 1))
            {
                n *= 10;
                *data_vectors = (double **)realloc(*data_vectors, n * sizeof(double *));
                assert(*data_vectors != NULL && "An Error Has Occured");
            }
    if (c != '\n' && c!='\r'){
        (*data_vectors)[j] = vector;
        j++;
        d=i;
    }
    *data_vectors = (double **)realloc(*data_vectors, (j) * sizeof(double *));
    assert(*data_vectors != NULL && "An Error Has Occured");
    value = (int *)calloc(2, sizeof(int));
    assert(value != NULL && "An Error Has Occured");
    value[0] = j;
    value[1] = d;
    fclose(ifp);
    return value;
}

double calcWeight(double* observ1, double* observ2, int dim)
{
    double sum = 0, currDiff=0, norm=0;
    int i=0;
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
    assert(wam!=NULL && "An Error Has Occured");
    for(i=0;i<n;i++)
    {
        wam[i] = calloc(n,sizeof(double));
        assert( wam[i]!=NULL && "An Error Has Occured");
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

double** DiagonalDegreeMatrix(double** matrix, int n, int sqrtMode)
{
    double** ddm;
    double sumOfRow;
    int i,j;
    ddm = calloc(n,sizeof(double*));
    assert(ddm!=NULL && "An Error Has Occured");
    for(i=0;i<n;i++)
    {
        ddm[i] = calloc(n,sizeof(double));
        assert(ddm[i]!=NULL && "An Error Has Occured");
    }
    for (i=0;i<n;i++){
        sumOfRow = 0;
        for (j=0;j<n;j++){
            sumOfRow += matrix[i][j];
        }
        if (sqrtMode){
            ddm[i][i] = 1/sqrt(sumOfRow);
        }
        else{
            ddm[i][i] =(sumOfRow);
        }
        
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
    int i=0;
    Lnorm = calloc(n,sizeof(double*));
    assert(Lnorm!=NULL && "An Error Has Occured");
    midMatrix = calloc(n,sizeof(double*));
    assert(midMatrix!=NULL && "An Error Has Occured");
    for(i=0;i<n;i++)
    {
        Lnorm[i] = calloc(n,sizeof(double));
        assert(Lnorm[i]!=NULL && "An Error Has Occured");
        midMatrix[i] = calloc(n,sizeof(double));
        assert(midMatrix[i]!=NULL && "An Error Has Occured");
    }
    MatrixMultiply_helper(ddm_square,wam,midMatrix,n);
    MatrixMultiply_helper(midMatrix,ddm_square,Lnorm,n);
    LaplacianNorm_helper(Lnorm, n);
    freeMatrix(midMatrix,n);
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
    int r;
    Icol = calloc(n,sizeof(double));
    assert(Icol!=NULL && "An Error Has Occured");
    Jcol = calloc(n,sizeof(double));
    assert(Jcol!=NULL && "An Error Has Occured");
    for (r=0;r<n;r++)
    {
        Icol[r] = matrixA[r][i];
        Jcol[r] = matrixA[r][j];
    }
    for(r=0; r<n;r++)
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

void UpdateMatrixV(double **matrixV, double c, double s, int i, int j, int n)
{
    int r=0;
    double *Icol, *Jcol;
    Icol = calloc(n,sizeof(double));
    assert(Icol!=NULL && "An Error Has Occured");
    Jcol = calloc(n,sizeof(double));
    assert(Jcol!=NULL && "An Error Has Occured");
    for (r=0;r<n;r++)
    {
        Icol[r] = matrixV[r][i];
        Jcol[r] = matrixV[r][j];
    }
    for (r=0;r<n;r++)
    {
        matrixV[r][i] = c*Icol[r] - s*Jcol[r];
        matrixV[r][j] = s*Icol[r] + c*Jcol[r];
    }
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


void getEignValues(double** matrixA,double* EignValues, int n)
{
    int i;
    for(i=0;i<n;i++)
    {
        EignValues[i] = matrixA[i][i];
    }
}

void JacobiAlgorithm(double** matrixA,double** matrixV, int n)
{

    double c, s,offA , offAtag, epsilon = pow(10,-15);
    int maxElementOffDiagonalI, maxElementOffDiagonalJ,i,j, stopCondition=1, countIter=0;
    for (i=0; i<n;i++)
    {
        matrixV[i][i] = 1; /*Init V to be identity matrix*/
        for(j=0;j<n;j++)
        {
            if(i!=j)
            {
                matrixV[i][j] =0;
            }
        }
    }
    do{
        countIter++;
        calcOFFMatrix(matrixA, &offA, n);
        findMaxOffDiagonal(matrixA,&maxElementOffDiagonalI,&maxElementOffDiagonalJ,n);
        calculateRotateValues(matrixA,&c, &s,maxElementOffDiagonalI,maxElementOffDiagonalJ);
        RotateMatrix_helper(matrixA,c, s,maxElementOffDiagonalI,maxElementOffDiagonalJ, n); /*Calculate A' matrix*/
        UpdateMatrixV(matrixV,c,s,maxElementOffDiagonalI,maxElementOffDiagonalJ,n); /*Update MatrixV - The matrix with EignVectors*/
        calcOFFMatrix(matrixA, &offAtag, n);
        if ((offA-offAtag)<epsilon)
        {
            stopCondition = 0;
        }
        if (offAtag==0) {break;}
    } while(stopCondition && countIter<100);
}

int TheEigengapHeuristic(double* eigenvalues, int len) {
    double deltaI = 0;
    double currMax = 0;
    int position=0;
    int i;
    qsort(eigenvalues,len, sizeof(double),doubleCmpr);
    for(i=1; i<=(len/2);i++){
        deltaI = eigenvalues[i]-eigenvalues[i-1];
        if (deltaI > currMax){
            currMax = deltaI;
            position = i;
        }
    }
    return position;
}


void getMatrixSortedEignVectors(double** matrixA, double** matrixV, double** matrixU, int n, int k)
{
    EignValue* eignValues;
    int i,j, indexElem;
    eignValues = calloc(n,sizeof(EignValue));
    assert(eignValues!=NULL && "An Error Has Occured");
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
            matrixU[i][j] = matrixV[i][indexElem];
        }
    }
    free(eignValues);
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
        if (sumRow!=0)
        {
            for(j=0;j<k;j++)
                {
                    matrixU[i][j] = matrixU[i][j]/pow(sumRow,0.5);
                }
        }
    }

}

double** getNewDataPointsDimK(double** observations, int n, int dim, int* k)
{
    int i;
    double **weightedAdjMatrix, **Lnorm, **ddMatrix, **EignVectorsMatrix,**matrixNewPointsToKmeans, *EignValues;
    weightedAdjMatrix = CreateWeightedAdjacencyMatrix(observations,dim,n);
    ddMatrix = DiagonalDegreeMatrix(weightedAdjMatrix,n,1);
    Lnorm = ComputeNormalizedGraphLaplacian(weightedAdjMatrix,ddMatrix,n);
    EignVectorsMatrix = (double**)calloc(n,sizeof(double*));
    assert(EignVectorsMatrix!=NULL && "An Error Has Occured");
    matrixNewPointsToKmeans = (double**)calloc(n,sizeof(double*));
    assert(matrixNewPointsToKmeans!=NULL && "An Error Has Occured");
    for(i=0;i<n;i++)
    {
        EignVectorsMatrix[i] = calloc(n,sizeof(double));
        assert(EignVectorsMatrix[i]!=NULL && "An Error Has Occured");
        matrixNewPointsToKmeans[i] = calloc(n,sizeof(double));
        assert(matrixNewPointsToKmeans[i]!=NULL && "An Error Has Occured");
    }
    JacobiAlgorithm(Lnorm,EignVectorsMatrix,n);
    if(*k==0)
    {
        EignValues = calloc(n,sizeof(double));
        getEignValues(Lnorm,EignValues,n);
        *k = TheEigengapHeuristic(EignValues,n);
        free(EignValues);
    }
    getMatrixSortedEignVectors(Lnorm,EignVectorsMatrix,matrixNewPointsToKmeans,n,*k);
    normalizedMatrixUtoMatrixT(matrixNewPointsToKmeans,n,*k);
    freeMatrix(weightedAdjMatrix,n);
    freeMatrix(Lnorm,n);
    freeMatrix(ddMatrix,n);
    freeMatrix(EignVectorsMatrix,n);
    return matrixNewPointsToKmeans;
}

void getFirstKCentroids(double** dataPoints, double** centroidsToFill, int k)
{
    int i,j;
    for(i=0;i<k;i++)
    {
        for(j=0;j<k;j++)
        {
            centroidsToFill[i][j] = dataPoints[i][j];
        }
    }
}

int checkNegativeZero(double value)
{
    if (value>-0.00005 && value<0){ return 1;}
    return 0;
}
void printMatrix(double** matrix, int a, int b) {
    int i,j;
    for (i = 0; i < a; i++)
    {
        for (j = 0; j <  b - 1; j++)
        {
            if (checkNegativeZero(matrix[i][j])){
                printf("0.0000,");
            } else{
                printf("%.4f,", matrix[i][j]);
            }
            
        }
        if (checkNegativeZero(matrix[i][b-1])){
            printf("0.0000");
        } else{
            printf("%.4f", matrix[i][b-1]);
        }
        putchar('\n');
    }
}


void flowSPKforC(double** observations, int n, int dim, int k,int max_iter)
{
    double** centroids, **newDataPoints;
    int i;
    newDataPoints = getNewDataPointsDimK(observations,n,dim, &k); /*Get the new Data Points from K dimension into newDataPoints*/
    centroids = (double**)calloc(k,sizeof(double*));
    assert(centroids!=NULL && "An Error Has Occured");
    for(i=0;i<k;i++)
    {
        centroids[i] = calloc(k,sizeof(double));
        assert(centroids[i]!=NULL && "An Error Has Occured");
    }
    getFirstKCentroids(newDataPoints,centroids,k); /*Start KMeans Algorithm*/
    calculate_kmeans(newDataPoints,centroids,n,k,k,max_iter);
    printMatrix(centroids,k,k);
    freeMatrix(newDataPoints,n);
    freeMatrix(centroids,k);
}

void flowJacobiAlgo(double** matrix,int n)
{   
    double** matrixEignVectors;
    int i,j;
    matrixEignVectors = (double**)calloc(n,sizeof(double*));
    assert(matrixEignVectors!=NULL && "An Error Has Occured");
    for(i=0;i<n;i++)
    {
        matrixEignVectors[i] = calloc(n,sizeof(double));
        assert(matrixEignVectors[i]!=NULL && "An Error Has Occured");
    }
    JacobiAlgorithm(matrix,matrixEignVectors,n); /*Run Jacobi Algorithm to get EignValues and EignVectors*/
    for(i=0;i<n;i++)
    {
        if (i<(n-1)){
            printf("%0.4f,",matrix[i][i]);
        } else{
            printf("%0.4f\n",matrix[i][i]);
        }
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if (j<(n-1))
            {
            if (checkNegativeZero(matrixEignVectors[j][i])){
                printf("0.0000,");
            } else{
                printf("%.4f,", matrixEignVectors[j][i]);
            }
            }
            else
            {
                if (checkNegativeZero(matrixEignVectors[n-1][i])){
                    printf("0.0000");
                } else{
                    printf("%.4f\n", matrixEignVectors[n-1][i]);
                    }
            }
        }
    }
}

void flowWam(double** data_vectors,int d, int n) {
    double** WeightedAdjMatrix;
    WeightedAdjMatrix = CreateWeightedAdjacencyMatrix(data_vectors, d, n);
    printMatrix(WeightedAdjMatrix,n,n);
    freeMatrix(WeightedAdjMatrix,n);
}

void flowDdg(double** data_vectors,int d, int n) {
    double** WeightedAdjMatrix, **DDMatrix;
    WeightedAdjMatrix = CreateWeightedAdjacencyMatrix(data_vectors, d, n);
    DDMatrix = DiagonalDegreeMatrix(WeightedAdjMatrix,n,0);
    printMatrix(DDMatrix,n,n);
    freeMatrix(WeightedAdjMatrix,n);
    freeMatrix(DDMatrix,n);
}

void flowLnorm(double** data_vectors,int d, int n) {
    double** WeightedAdjMatrix, **DDMatrix, **lNormMatrix;
    WeightedAdjMatrix = CreateWeightedAdjacencyMatrix(data_vectors, d, n);
    DDMatrix = DiagonalDegreeMatrix(WeightedAdjMatrix,n,1);
    lNormMatrix = ComputeNormalizedGraphLaplacian(WeightedAdjMatrix,DDMatrix,n);
    printMatrix(lNormMatrix,n,n);
    freeMatrix(WeightedAdjMatrix,n);
    freeMatrix(DDMatrix,n);
    freeMatrix(lNormMatrix,n);
}

int main(int argc, char *argv[])
{   
    int* values;
    int n,k,d;
    char* flow, *nameOfFile;
    double** data_vectors;
    n = 10;
    data_vectors = (double **)calloc(n, sizeof(double *));
    assert(data_vectors!=NULL && "An Error Has Occured");
    k = atoi(argv[1]);
    if (argc!=4)
    {
        printf("Invalid Input!");
        exit(0);
    }
    flow = argv[2];
    nameOfFile = argv[3];
    values = initDataPoints(nameOfFile, &data_vectors);
    n = values[0];
    d = values[1];
    free(values);
    switch(convertStringIntoGoalEnum(flow))
    {
        case spk:
            if (k>=n || k<0)
                {
                printf("Invalid Input!");
                freeMatrix(data_vectors,n);
                exit(0);
                }
            flowSPKforC(data_vectors,n,d,k,300);
            break;
        case wam:
            flowWam(data_vectors, d, n);
            break;
        case ddg:
            flowDdg(data_vectors,d,n);
            break;
        case lnorm:
            flowLnorm(data_vectors,d,n);
            break;
        case jacobi:
            flowJacobiAlgo(data_vectors,n);
            break;
        default:
            printf("Invalid Input!");
            freeMatrix(data_vectors,n);
            exit(0);
    }
    freeMatrix(data_vectors,n);
  return 0;
}

