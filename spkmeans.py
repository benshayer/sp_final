import numpy as np
import csv
import spkmeans
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("args", nargs="+")
parser = parser.parse_args()
args = parser.args
if len(args) == 3:
    k,  goal, file_name = int(args[0]), str(args[1]), str(args[2])
else:
    print("Invalid Input!")
    exit()


def calculate_points_distance(row_i, row_j):
    distance = np.sum([(x - y) ** 2 for (x, y) in zip(row_i, row_j)])
    return distance

def initial_kmeans_pp_centroids(k, full_data_points):
    n = len(full_data_points)
    if k >= n:
        print("Invalid Input!")
        exit()
    initial_centroids = calculate_initial_centroids(k, full_data_points)
    centroids_points = get_centroids_points(initial_centroids, full_data_points)
    return [centroids_points, full_data_points, initial_centroids]

def calculate_initial_centroids(k, data_points):
    np.random.seed(0)
    centroids = []
    n = len(data_points)
    first_centroid = np.random.choice(n)
    centroids.append(first_centroid)
    probabilities = np.zeros((n), dtype=np.float64)
    min_distances = np.full(n, np.inf)
    total_distances = 0
    curr_distance = 0
    z = 0
    while z != k - 1:
        new_centroid = data_points.loc[centroids[z]]
        for curr_row in data_points.itertuples():
            index, curr_vector = int(curr_row.Index), curr_row[1:]
            curr_distance = calculate_points_distance(curr_vector, new_centroid)
            if curr_distance < min_distances[index]:
                min_distances[index] = curr_distance
        total_distances = min_distances.sum()
        for i in range(n):
            probabilities[i] = min_distances[i] / total_distances
        next_centroid = np.random.choice(n, size=None, p=probabilities)
        centroids.append(next_centroid)
        z += 1

    return centroids

def get_centroids_points(centroids_indexes, data_points):
    centroids_points = data_points.loc[centroids_indexes]
    return centroids_points


def print_centroids(centroids):
    k = len(centroids)
    d = len(centroids[0])
    for i in range(k):
        for j in range(d - 1):
            print(f'{format(centroids[i][j],".4f")}',end=",")
        print(f'{format(centroids[i][d-1], ".4f")}')


def initDataPointsPython(filename):
     dataVectors = []
     file = open(filename, 'r')
     for row in file:
         rowWithOutComma = row.split(',')
         newRow=[]
         for item in rowWithOutComma:
             newRow.append(float(item))
         dataVectors.append(newRow)
     return dataVectors

if __name__ == '__main__':
    dataVectors = initDataPointsPython(file_name)
    n = len(dataVectors)
    d = len(dataVectors[0])
    max_iter=300
    if goal=="spk":
        if (k<0 or k>=n):
            print("Invalid Input!")
            exit()
        newDataPoints = spkmeans.getnew_datapoints(n,d,k,dataVectors)
        newDataPoints_df = pd.DataFrame(newDataPoints)
        k=len(newDataPoints[0])
        centroids_indexes = calculate_initial_centroids(k, newDataPoints_df)
        centroids_points = get_centroids_points(centroids_indexes, newDataPoints_df)
        centroids_list =centroids_points.values.tolist()
        for j in range(k - 1):
            print(centroids_indexes[j], end=',')
        print(centroids_indexes[k - 1])
        final_centroids = spkmeans.kmeans_pp(n, k, k, max_iter, newDataPoints, centroids_list)
    elif goal=="wam":
        spkmeans.wam(n,d,dataVectors)
    elif goal=="ddg":
        spkmeans.ddg(n,d,dataVectors)
    elif goal == "lnorm":
        spkmeans.lnorm(n,d,dataVectors)
    elif goal=="jacobi":
        spkmeans.jacobi(n,dataVectors)



