#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void calculate_kmeans(double** data_vectors, double** centroids, int n, int d, int k,int max_iter);
double diffFromCentroid(double *curr_vector, double *centroid, int d);
int findMatchCluster(double *curr_vector, double **centroids, int d, int k);
void calculateCentroids(double ***new_centroids, int *sizes, int d, int k);