import math

x = 3.1415929336


def truncate(f, n):
    return math.floor(f * 10 ** n) / 10 ** n


print(truncate(x, 5))
