import pandas as pd
import numpy as np
import argparse
import mykmeanssp
import time

parser = argparse.ArgumentParser()
parser.add_argument("args", nargs="+")
parser = parser.parse_args()
args = parser.args
if len(args) == 4:
    k, max_iter, file_name_1, file_name_2 = int(args[0]), int(args[1]), str(args[2]), str(args[3])
else:
    k, file_name_1, file_name_2 = int(args[0]), str(args[1]), str(args[2])
    max_iter = 300


def convert_txt_df(txt_path):
    data_points = pd.read_csv(txt_path, header=None, index_col=0)
    return data_points


def get_full_data_points(file_path_1, file_path_2):
    data_points_1 = convert_txt_df(file_path_1)
    data_points_2 = convert_txt_df(file_path_2)
    full_data_points = pd.merge(data_points_1, data_points_2, left_index=True, right_index=True)
    return full_data_points


def calculate_points_distance(row_i, row_j):
    distance = np.sum([(x - y) ** 2 for (x, y) in zip(row_i, row_j)])
    return distance


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
        start_function_time = time.time()
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


def initial_kmeans_pp_centroids(k, max_iter, path_1, path_2):
    full_data_points = get_full_data_points(path_1, path_2)
    initial_centroids = calculate_initial_centroids(k, full_data_points)
    centroids_points = get_centroids_points(initial_centroids, full_data_points)
    return [centroids_points, full_data_points, initial_centroids]


def print_centroids(centroids):
    k = len(centroids)
    d = len(centroids[0])
    for i in range(k):
        for j in range(d - 1):
            print('{:.4f}'.format(centroids[i][j]), end=',')
        print('{:.4f}'.format(centroids[i][d - 1]))


if __name__ == '__main__':
    centroids, full_data_points, centroids_indexes = initial_kmeans_pp_centroids(k, max_iter, file_name_1, file_name_2)
    centroids_list = centroids.values.tolist()
    full_data_points_list = full_data_points.values.tolist()
    n = len(full_data_points_list)
    d = len(full_data_points_list[0])
    for j in range(k - 1):
        print(centroids_indexes[j], end=',')
    print(centroids_indexes[k - 1])
    final_centroids = mykmeanssp.fit(n, d, k, max_iter, full_data_points_list, centroids_list)
    print_centroids(final_centroids)
