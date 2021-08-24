import numpy as np
import csv
import spkmeans
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("args", nargs="+")
parser = parser.parse_args()
args = parser.args
if len(args) == 3:
    k,  goal, file_name = int(args[0]), str(args[1]), str(args[2])
else:
    print("Error!")
    exit()


# WAM = np.array([[0, 4, 3], [4, 0, 5], [3, 5, 0]])
# DDM = np.array([[1 / 5 ** 0.5, 0, 0], [0, 1 / 6 ** 0.5, 0], [0, 0, 1 / 8 ** 0.5]])
# mid = (np.matmul(DDM, WAM))
# x = np.matmul(mid,DDM)
# id = np.identity(3)
# l_norm = id-x
# for i in range(2):
#     for j in range(2):
#         print(x[i][j], end='')
#     print('\n')


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
    spkmeans.jacobi(n,dataVectors)



