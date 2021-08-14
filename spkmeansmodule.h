double calcWeight(double* observ1, double* observ2, int dim);
double** ComputeNormalizedGraphLaplacian(double** wam, double** ddm_square, int n);
double** CreateWeightedAdjacencyMatrix(double** observations, int dim, int n);
void flowSPKforC(double** observations, int n, int dim, int k,int max_iter);
void getMatrixSortedEignVectors(double** matrixA, double** matrixV, double** matrixU, int n, int k);
double** getNewDataPointsDimK(double** observations, int n, int dim, int* k);
int* initDataPoints(char* filename, double ***data_vectors);
void getMatrixSortedEignVectors(double** matrixA, double** matrixV, double** matrixU, int n, int k);
void JacobiAlgorithm(double** matrixA,double** matrixV, int n);
void LaplacianNorm_helper(double** Lnorm, int n);
void MatrixMultiply_helper(double** matrixA, double** matrixB, double** result, int n);
int TheEigengapHeuristic(double* eignvalues, int len);
double** DiagonalDegreeMatrix(double** matrix, int n)
