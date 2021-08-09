from setuptools import setup, Extension


setup(name='mykmeanssp',
      version='1.0',
      description='My implemantation to calculate kmeans from given Data Points and Centroids. By Ben and Gilad',
      ext_modules=[Extension('mykmeanssp', sources=['kmeans.c'])])
