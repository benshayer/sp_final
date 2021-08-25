from setuptools import setup, Extension


setup(name='spkmeans',
      version='1.0',
      description='My implemantation to calculate kmeans from given Data Points and Centroids. By Ben and Gilad',
      ext_modules=[Extension('spkmeans', sources=['spkmeansmodule.c','spkmeans.c'])])
