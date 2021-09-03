#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spkmeans.h"
#include "kmeans.h"

static PyObject* kmeans_pp_capi(PyObject *self, PyObject *args);
static PyObject* getSPKDataPoints_capi(PyObject *self, PyObject *args);
static PyObject* jacobi_capi(PyObject *self, PyObject *args);
static PyObject* wam_capi(PyObject *self, PyObject *args);
static PyObject* ddg_capi(PyObject *self, PyObject *args);
static PyObject* lnorm_capi(PyObject *self, PyObject *args);
static void convert_data_to_c(PyObject* data_points_py, double** data_points, int n, int d);
static void convert_point_to_c(PyObject* point_py, double* point, int d);

static PyObject* kmeans_pp_capi(PyObject *self, PyObject *args)
{
    PyObject* data_points_py;
    PyObject* centroids_py;
    int n;
    int d;
    int k;
    int i;
    int j;
    int max_iter;
    double** data_points;
    double** centroids;
    PyObject * item;
    PyObject* final_centroids;
    PyObject* current_centroid;
    if(!PyArg_ParseTuple(args,"iiiiOO",&n,&d,&k,&max_iter,&data_points_py,&centroids_py)){
        return NULL;
    }
    final_centroids = PyList_New(k);
    data_points = (double **)calloc(n, sizeof(double*));
    assert(data_points!=NULL);
    centroids = (double **)calloc(k,sizeof(double*));
    assert(centroids!=NULL);
    convert_data_to_c(data_points_py,data_points,n,d);
    convert_data_to_c(centroids_py,centroids,k,d);
    calculate_kmeans(data_points,centroids,n,d,k,max_iter);
    for (i=0;i<k;i++){
        current_centroid = PyList_New(d);
        for (j=0;j<d;j++){
            item = PyFloat_FromDouble(centroids[i][j]);
            PyList_SetItem(current_centroid,j,item);
        }
        PyList_SetItem(final_centroids,i,current_centroid);
        free(centroids[i]);
    }
    for (i=0;i<n;i++){
        free(data_points[i]);
    }
    free(centroids);
    free(data_points);
    return final_centroids;

}

static PyObject* getSPKDataPoints_capi(PyObject *self, PyObject *args)
{
    PyObject* data_points_py;
    PyObject* PyNewDataPoints;
    PyObject *item;
    PyObject* currentLine;
    int n,d,k,i,j;
    double** data_points, **newDataPoints;
    if(!PyArg_ParseTuple(args,"iiiO",&n,&d,&k,&data_points_py)){
        return NULL;
    }
    PyNewDataPoints = PyList_New(n);
    data_points = (double **)calloc(n, sizeof(double*));
    assert(data_points!=NULL);
    convert_data_to_c(data_points_py,data_points,n,d);
    newDataPoints = getNewDataPointsDimK(data_points,n,d,&k);
    for (i=0;i<n;i++){
        currentLine = PyList_New(k);
        for (j=0;j<k;j++){
            item = PyFloat_FromDouble(newDataPoints[i][j]);
            PyList_SetItem(currentLine,j,item);
        }
        PyList_SetItem(PyNewDataPoints,i,currentLine);
        free(newDataPoints[i]);
    }
    for (i=0;i<n;i++){
        free(data_points[i]);
    }
    free(newDataPoints);
    free(data_points);
    return PyNewDataPoints;

}

static PyObject* jacobi_capi(PyObject *self, PyObject *args)
{
    PyObject* data_points_py;
    int n;
    double** data_points;
    if(!PyArg_ParseTuple(args,"iO",&n,&data_points_py)){
        return NULL;
    }
    data_points = (double **)calloc(n, sizeof(double*));
    assert(data_points!=NULL);
    convert_data_to_c(data_points_py,data_points,n,n);
    flowJacobiAlgo(data_points,n);
    freeMatrix(data_points,n);
    
    Py_RETURN_NONE;
}


static PyObject* ddg_capi(PyObject *self, PyObject *args)
{

    PyObject* data_points_py;
    int n, d;
    double** data_points;
    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&data_points_py)){
        return NULL;
    }
    data_points = (double **)calloc(n, sizeof(double*));
    assert(data_points!=NULL);
    convert_data_to_c(data_points_py,data_points,n,d);
    flowDdg(data_points,d,n);
    freeMatrix(data_points,n);
    Py_RETURN_NONE;
}

static PyObject* lnorm_capi(PyObject *self, PyObject *args)
{
    PyObject* data_points_py;
    int n, d;
    double** data_points;
    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&data_points_py)){
        return NULL;
    }
    data_points = (double **)calloc(n, sizeof(double*));
    assert(data_points!=NULL);
    convert_data_to_c(data_points_py,data_points,n,d);
    flowLnorm(data_points,d,n);
    freeMatrix(data_points,n);
    Py_RETURN_NONE;
}

static PyObject* wam_capi(PyObject *self, PyObject *args)
{
    PyObject* data_points_py;
    int n, d;
    double** data_points;
    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&data_points_py)){
        return NULL;
    }
    data_points = (double **)calloc(n, sizeof(double*));
    assert(data_points!=NULL);
    convert_data_to_c(data_points_py,data_points,n,d);
    flowWam(data_points,d,n);
    Py_RETURN_NONE;
}


static void convert_data_to_c(PyObject* data_points_py, double** data_points, int n, int d){
    int i;
    int j;
    double current_item;
    PyObject* point_py;
    double* point;
    for (i=0;i<n;i++){
        point_py = PyList_GET_ITEM(data_points_py,i);
        assert(point!=NULL);
        point = (double*)calloc(d,sizeof(double));
        data_points[i] = (double *)calloc(d,sizeof(double));
        convert_point_to_c(point_py, point, d);
        for(j=0;j<d;j++){
            current_item = point[j];
            data_points[i][j] = current_item;
        }
        free(point);
        
    }

}

static void convert_point_to_c(PyObject* point_py, double* point, int d){
    int i;
    double temp;
    PyObject* item;
    for (i=0;i<d;i++){
        item = PyList_GetItem(point_py,i);
        temp = PyFloat_AsDouble(item);
        point[i] = temp;
        
    }
}

static PyMethodDef _capiMethods[] = {
    {"kmeans_pp", (PyCFunction) kmeans_pp_capi, METH_VARARGS, PyDoc_STR("Get datapoints and centroids and calc kmeans")},
    {"getnew_datapoints", (PyCFunction) getSPKDataPoints_capi, METH_VARARGS, PyDoc_STR("A function to run Jacobi algorithm to find Eignvalues and Eignvectors")},
    {"jacobi", (PyCFunction) jacobi_capi, METH_VARARGS, PyDoc_STR("A function to run Jacobi algorithm to find Eignvalues and Eignvectors")},
    {"wam", (PyCFunction) wam_capi, METH_VARARGS, PyDoc_STR("A function to find WAM for set of data points")},
    {"ddg", (PyCFunction) ddg_capi, METH_VARARGS, PyDoc_STR("A function to find DDG for set of data points")},
    {"lnorm", (PyCFunction) lnorm_capi, METH_VARARGS, PyDoc_STR("A function to find LNORM for set of data points")},
    {NULL,NULL,0,NULL}
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans",
    NULL,
    -1,
    _capiMethods
};

PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if(!m){
        return NULL;
    }
    return m;
}