#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "kmeans.h"

void calculateCentroids(double ***new_centroids, int *sizes, int d, int k)
{
    int i;
    int j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            if (sizes[i]!=0)
            {
                            (*new_centroids)[i][j] = (*new_centroids)[i][j] / sizes[i];
            }
        }
    }
}

void calculate_kmeans(double** data_vectors, double** centroids, int n, int d, int k,int max_iter){
    int count_iter;
    double **new_centroids;
    int *sizes;
    int flag_centroids_different;
    int i;
    int j;
    int m;
    flag_centroids_different = 0;
    count_iter = 0;
    while (count_iter < max_iter && !flag_centroids_different)
    {
        new_centroids = (double **)calloc(k, sizeof(double *));
        assert(new_centroids!=NULL);
        sizes = (int *)calloc(k, sizeof(int));
        for (j = 0; j < k; j++)
        {
            new_centroids[j] = (double *)calloc(d, sizeof(double));
            assert(new_centroids[j]!=NULL);
        }
        for (i = 0; i < n; i++)
        {
            int match_cluster = findMatchCluster(data_vectors[i], centroids, d, k);
            sizes[match_cluster] += 1;
            for (m = 0; m < d; m++)
            {
                new_centroids[match_cluster][m] += data_vectors[i][m];
            }
        }
        calculateCentroids(&new_centroids, sizes, d, k);
        flag_centroids_different = 1;
        for (i = 0; i < k; i++)
        {
            for (j = 0; j < d; j++)
            {
                if (new_centroids[i][j] != centroids[i][j])
                {
                    flag_centroids_different = 0;
                    break;
                }
            }
            if (!flag_centroids_different)
            {
                break;
            }
        }
        for (i = 0; i < k; i++)
        {
            for (j = 0; j < d; j++)
            {
                centroids[i][j] = new_centroids[i][j];
            }
        }
        count_iter++;
        for (i = 0; i < k; i++)
        {
            free(new_centroids[i]);
        }
        free(new_centroids);
        free(sizes);
    }
}


double diffFromCentroid(double *curr_vector, double *centroid, int d)
{
    int i;
    double distance = 0;
    for (i = 0; i < d; i++)
    {
        distance += (curr_vector[i] - centroid[i]) * (curr_vector[i] - centroid[i]);
    }
    return distance;
}

int findMatchCluster(double *curr_vector, double **centroids, int d, int k)
{
    double min_distance = 1000000;
    double curr_distance;
    int min_index = -1;
    int i;
    for (i = 0; i < k; i++)
    {
        curr_distance = diffFromCentroid(curr_vector, centroids[i], d);
        if (curr_distance < min_distance)
        {
            min_distance = curr_distance;
            min_index = i;
        }
    }
    return min_index;
}

