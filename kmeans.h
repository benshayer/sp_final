#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


static PyObject* fit_capi(PyObject *self, PyObject *args);
static void convert_data_to_c(PyObject* data_points_py, double** data_points, int n, int d);
static void convert_point_to_c(PyObject* point_py, double* point, int d);
static void calculate_kmeans(double** data_points, double** centroids, int n, int d, int k,int max_iter);
static double diffFromCentroid(double *curr_vector, double *centroid, int d);
static int findMatchCluster(double *curr_vector, double **centroids, int d, int k);
static void calculateCentroids(double ***new_centroids, int *sizes, int d, int k);