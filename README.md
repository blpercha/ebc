Ensemble Biclustering for Classification
==============

A python implementation of the Ensemble Biclustering for Classification (EBC) algorithm. Although having "biclustering" in its name, EBC is a co-clustering algorithm that allows you to perform co-clustering on very large sparse N-dimensional matrix. For details and examples of using EBC please reference [this paper](http://www.ncbi.nlm.nih.gov/pubmed/26219079).

## Files

- `ebc.py`: an implementation of EBC algorithm using sparse matrix.
- `ebc2d.py`: a vectorized implementation of EBC algorithm using (2D) numpy dense array.
- `matrix.py`: a N-dimensional sparse matrix implementation using python dict that allows basic get/set/sum operations.
- `ebc-sample-matrices/`: this directory contains many example data that could be used to run EBC on.

## Usage

#### General Usage

Using ebc is easy. First you need a sparse matrix constructed as a SparseMatrix instance defined in `matrix.py`. You can use `matrix.py`'s built in `read_data()` method to construct the sparse matrix.

Once you have the sparse matrix and import the ebc module, you can simply do the following:

    ebc = EBC(sparse_matrix, n_clusters=[30, 125], max_iterations=10)
    cXY, objective, iter = ebc.run()

The returned `cXY` contains the clustering assignments along each axis in a list of lists, `objective` contains the final objective function value, `iter` contains the iterations that it ran for to converge.

#### Efficiency Considerations

In short, `EBC` is built for highly sparse multi-dimensional matrix, while `EBC2D` is a vectorized version of `EBC`, and is built for less sparse 2D matrix. The running time of `EBC` increases linearly with the number of non-zero elements, while the running time of `EBC2D` will only increase when matrix size increases.

- If your input matrix is highly sparse or you need support for N-dimensional matrix with N > 2, you should use `EBC` class in `ebc.py`.
- If your input matrix is not very sparse (e.g. 5% density), and can fit into memory with `numpy`, you can choose to use `EBC2D` class in `ebc2d.py`. Note that `EBC2D` only supports 2D matrix, since large multi-dimensional dense matrix can hardly fit into memory.

To give you more sense of the efficiency, the `EBC` implementation runs for ~50 seconds on our 14052 x 7272 example matrix (which is >99% sparse) on a MacBook, and in this case it is faster than `EBC2D`. However, in our testing, for a 5000 x 5000 matrix of 95% sparsity, `EBC2D` can get a `x5` speedup over `EBC`. And this speedup grows as the sparsity decreases.

#### Dependencies

To run this implementation of EBC algorithm, you need to have `numpy` installed.

## References

- Percha, Bethany, and Russ B. Altman. "Learning the Structure of Biomedical Relationships from Unstructured Text." PLoS computational biology 11.7 (2015).
- Dhillon, Inderjit S., Subramanyam Mallela, and Dharmendra S. Modha. "Information-theoretic co-clustering." Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining. ACM, 2003.

## Questions?

This code was written by Beth Percha and Yuhao Zhang. We welcome any questions or comments, and would appreciate it if you would let us know if you make any substantial modifications or improvements. You can reach us at blpercha@stanford.edu.
