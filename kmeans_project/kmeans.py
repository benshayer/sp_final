import numpy as np
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument("args", nargs="+")
parser = parser.parse_args()
args = parser.args
if len(args) == 2:
    k, max_iter = int(args[0]), int(args[1])
else:
    k = int(args[0])
    max_iter = 200


def init_data_points():
    data_lines = []
    data_list = []
    while True:
        try:
            data_lines.append(input())
        except EOFError:
            break
    for line in data_lines:
        vector = []
        for value in line.split(','):
            if value:
                vector.append(float(value))
        if vector:
            data_list.append(vector)
    return data_list


def init_centroids(k, data_points_list):
    centroids = []
    for i in range(k):
        centroids.append(data_points_list[i])
    return centroids


def diff_from_centroid(x, centroid):
    distance = 0
    d = len(x)
    for i in range(d):
        distance += (x[i] - centroid[i]) ** 2
    return distance


def mean_of_cluster(cluster):
    d = len(cluster[0])
    centroid = [0 for i in range(d)]
    for i in range(d):
        for vector in cluster:
            try:
                centroid[i] += vector[i]
            except:
                print(i, vector)
                raise Exception

        centroid[i] = centroid[i] / len(cluster)
    return centroid


def calculate_k_means(k, max_iter):
    if (k <= 0 or max_iter < 0):
        print("The arguments should be in correct format!")
        exit()
    if not max_iter:
        max_iter = 200
    data_points_list = init_data_points()
    centroids = init_centroids(k, data_points_list)
    prev_centroids = init_centroids(k, data_points_list)
    clusters = [[] for i in range(k)]
    count_iter = 0
    while (count_iter < max_iter):
        for current_vector in data_points_list:
            distances_from_centroids = [diff_from_centroid(current_vector, centroids[i]) for i in range(k)]
            cluster_position = distances_from_centroids.index(min(distances_from_centroids))
            clusters[cluster_position].append(current_vector)
        for i in range(k):
            centroids[i] = mean_of_cluster(clusters[i])
        clusters = [[] for i in range(k)]
        count_iter += 1
        if prev_centroids == centroids:
            break
        prev_centroids = list.copy(centroids)
    for i in range(k):
        d = len(centroids[i])
        for j in range(d - 1):
            print('{:.4f}'.format(centroids[i][j]), end=',')
        print('{:.4f}'.format(centroids[i][d - 1]))


if __name__ == '__main__':
    calculate_k_means(k, max_iter)
