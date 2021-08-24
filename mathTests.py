# import numpy as np
# import csv


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


# def initDataPointsPython(filename):
#     dataVectors = []
#     file = open("g.txt", 'r')
#     for row in file:
#         rowWithOutComma = row.split(',')
#         newRow=[]
#         for item in rowWithOutComma:
#             newRow.append(int(item))
#         dataVectors.append(newRow)
#     return dataVectors




