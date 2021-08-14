#include "kmeans.h"

static PyObject* fit_capi(PyObject *self, PyObject *args)
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


static void calculate_kmeans(double** data_vectors, double** centroids, int n, int d, int k,int max_iter){
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


static double diffFromCentroid(double *curr_vector, double *centroid, int d)
{
    int i;
    double distance = 0;
    for (i = 0; i < d; i++)
    {
        distance += (curr_vector[i] - centroid[i]) * (curr_vector[i] - centroid[i]);
    }
    return distance;
}

static int findMatchCluster(double *curr_vector, double **centroids, int d, int k)
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

static void calculateCentroids(double ***new_centroids, int *sizes, int d, int k)
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

static PyMethodDef _capiMethods[] = {
    {"fit", (PyCFunction) fit_capi, METH_VARARGS, PyDoc_STR("A function to calculate kmeans of given data points and centroids")},
    {NULL,NULL,0,NULL}
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    _capiMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if(!m){
        return NULL;
    }
    return m;
}