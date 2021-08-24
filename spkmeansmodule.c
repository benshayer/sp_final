#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spkmeans.c"

static void jacobi_capi(PyObject *self, PyObject *args);
static void wam_capi(PyObject *self, PyObject *args);
static void ddg_capi(PyObject *self, PyObject *args);
static void lnorm_capi(PyObject *self, PyObject *args);
static void convert_data_to_c(PyObject* data_points_py, double** data_points, int n, int d);
static void convert_point_to_c(PyObject* point_py, double* point, int d);

static void jacobi_capi(PyObject *self, PyObject *args)
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
}


static void ddg_capi(PyObject *self, PyObject *args)
{
    PyObject* data_points_py;
    int n, d;
    double** data_points;
    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&data_points_py)){
        return NULL;
    }
    data_points = (double **)calloc(n, sizeof(double*));
    assert(data_points!=NULL);
    convert_data_to_c(data_points_py,data_points,n,n);
    flowDdg(data_points,d,n);
}

static void lnorm_capi(PyObject *self, PyObject *args)
{
    PyObject* data_points_py;
    int n, d;
    double** data_points;
    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&data_points_py)){
        return NULL;
    }
    data_points = (double **)calloc(n, sizeof(double*));
    assert(data_points!=NULL);
    convert_data_to_c(data_points_py,data_points,n,n);
    flowLnorm(data_points,d,n);
}

static void wam_capi(PyObject *self, PyObject *args)
{
    PyObject* data_points_py;
    int n, d;
    double** data_points;
    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&data_points_py)){
        return NULL;
    }
    data_points = (double **)calloc(n, sizeof(double*));
    assert(data_points!=NULL);
    convert_data_to_c(data_points_py,data_points,n,n);
    flowWam(data_points,d,n);
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
    {"jacobi", (PyCFunction) jacobi_capi, METH_VARARGS, PyDoc_STR("A function to run Jacobi algorithm to find Eignvalues and Eignvectors")},
    {"wam", (PyCFunction) wam_capi, METH_VARARGS, PyDoc_STR("A function to run Jacobi algorithm to find Eignvalues and Eignvectors")},
    {"ddg", (PyCFunction) ddg_capi, METH_VARARGS, PyDoc_STR("A function to run Jacobi algorithm to find Eignvalues and Eignvectors")},
    {"lnorm", (PyCFunction) lnorm_capi, METH_VARARGS, PyDoc_STR("A function to run Jacobi algorithm to find Eignvalues and Eignvectors")},
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