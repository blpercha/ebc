Ensemble Biclustering for Classification
==============

A python implementation of the Ensemble Biclustering for Classification (EBC) algorithm. EBC is a co-clustering algorithm that allows you to perform co-clustering on very large sparse N-dimensional matrix. For details and examples of using EBC please reference [this paper](http://www.ncbi.nlm.nih.gov/pubmed/26219079).

##Files:

- `ebc.py`: the implementation of EBC algorithm.
- `matrix.py`: a N-dimensional sparse matrix implementation using python dict that allows basic get/set/sum operations.
- `resources/`: this directory contains many example data that could be used to run EBC on.

##Usage:

Using ebc is easy. First you need a sparse matrix constructed as a SparseMatrix instance defined in `matrix.py`. You can use `matrix.py`'s built in `read_data()` method to construct the sparse matrix.

Once you have the sparse matrix and import the ebc module, you can simply do the following:

    ebc = EBC(sparse_matrix, n_clusters=[30, 125], max_iterations=10)
    cXY, objective, iter = ebc.run()

The returned `cXY` contains the clustering assignments along each axis in a list of lists, `objective` contains the final objective function value, `iter` contains the iterations that it ran for to converge.

##Dependency:

To run this implementation of EBC algorithm, you need to have `numpy` installed.

##Efficiency:

Need to be updated.