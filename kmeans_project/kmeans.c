#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


int *initDataPoints(double ***data_vectors);
double **initCentroids(double **data_vectors, int k, int d);
double diffFromCentroid(double *curr_vector, double *centroid, int d);
int findMatchCluster(double *curr_vector, double **centroids, int d, int k);
void calculateCentroids(double ***new_centroids, int *sizes, int d, int k);
void checkArguments(int k, int max_iter, int amount_vectors);
int *initDataPoints(double ***data_vectors)
{
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
    while (scanf("%lf%c", &current_value, &c) == 2)
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
            /*
            for (int k=0; k<i;k++)
                printf("%lf,",vector[k]);
            */
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
    /*
    for (int m=0; m<j;m++){
        for (int n=0; n<d;n++){
            printf("%Lf,",(*data_vectors)[m][n]);
        }
        putchar('\n');
    }
    */
    value = (int *)calloc(2, sizeof(int));
    assert(value != NULL);
    value[0] = j;
    value[1] = d;
    return value;
}

double **initCentroids(double **data_vectors, int k, int d)
{
    double **centroids = (double **)calloc(k, sizeof(double *));
    int j;
    int i;
    for (j = 0; j < k; j++)
    {
        centroids[j] = (double *)calloc(d, sizeof(double));
    }
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            centroids[i][j] = data_vectors[i][j];
        }
    }
    return centroids;
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

void calculateCentroids(double ***new_centroids, int *sizes, int d, int k)
{
    int i;
    int j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            (*new_centroids)[i][j] = (*new_centroids)[i][j] / sizes[i];
        }
    }
}

void checkArguments(int k, int max_iter, int amount_vectors)
{
    if (k <= 0 || max_iter < 0 || k>=amount_vectors)
    {
        printf("Arguments should be in the correct format!");
        exit(0);
    }
}

int main(int argc, char *argv[])
{
    int max_iter;
    int k;
    int n;
    double** data_vectors;
    int d;
    int count_iter;
    int amount_vector;
    int* values;
    double** centroids;
    double **new_centroids;
    int *sizes;
    int flag_centroids_different;
    int i;
    int j;
    int m;
    max_iter = 200;
    if (argc == 3)
    {
        max_iter = atoi(argv[2]);
    }
    k = atoi(argv[1]);
    n = 10;
    data_vectors = (double **)calloc(n, sizeof(double *));
    values = initDataPoints(&data_vectors);
    amount_vector = values[0];
    checkArguments(k, max_iter,amount_vector);
    d = values[1];
    count_iter = 0;
    free(values);
    centroids = initCentroids(data_vectors, k, d);
    flag_centroids_different = 0;

    while (count_iter < max_iter && !flag_centroids_different)
    {
        new_centroids = (double **)calloc(k, sizeof(double *));
        sizes = (int *)calloc(k, sizeof(int));
        for (j = 0; j < k; j++)
        {
            new_centroids[j] = (double *)calloc(d, sizeof(double));
        }
        for (i = 0; i < amount_vector; i++)
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
    }
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d - 1; j++)
        {
            printf("%0.4f,", centroids[i][j]);
        }
        printf("%0.4f", centroids[i][d - 1]);
        putchar('\n');
    }
    free(centroids);
    return 0;
}