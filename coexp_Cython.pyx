# coexp_Cython.pyx
# cython: language_level=3
# cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False


import numpy as np
cimport numpy as np
from scipy.stats import kendalltau
from libc.math cimport sqrt
from scipy.sparse import issparse
from .external.corals import CorALS
np.import_array()
from numpy.random import default_rng


# ----------- Cython In-Block Shuffle -----------

# Shuffles a row of matrix
# This code is used to validate the significance of correlations  probabilities. 
# The actual coorilation is compared to a randomized correlation to make sure that the correlation given is not merely gotten by chance



def _indices_block_shuffle(int n, int block_size, object rng):
  #cdef is a cython optimization and helps with integrating the underlying c, defining types.
   cdef int i
   blocks = np.arange(0, n, block_size)
   rng.shuffle(blocks)
   indices = []
   for i in blocks:
       for j in range(block_size):
           if i + j < n:
               indices.append(i + j)
   return np.array(indices, dtype=np.int32)


def inblock_shuffle(np.ndarray[np.float64_t, ndim=2] a, int block_size=20, int axis=0, object seed=365):
# cython definition of types
   cdef int n = a.shape[axis]
   cdef object rng = default_rng(seed)


   if block_size > n:
       raise ValueError("block_size must be less than or equal to a.shape[axis]")
   elif block_size < 1:
       raise ValueError("block_size must be greater than 0")
   elif block_size == 1:
       return a
   elif block_size == n:
       return rng.permutation(a, axis=axis)
   else:
       idx = _indices_block_shuffle(n, block_size, rng)
       return np.take(a, idx, axis=axis)


def cy_correlation_matrix(
    str method,
    np.ndarray[np.float64_t, ndim=2] x,
    np.ndarray[np.float64_t, ndim=2] y = None,
):
   cdef int i, j, n_x, n_y
   cdef np.ndarray[np.float64_t, ndim=2] corrs


   if issparse(x):
       x = x.toarray()
   if y is not None and issparse(y):
       y = y.toarray()
   elif y is None:
       y = x


   n_x = x.shape[1]
   n_y = y.shape[1]

# Pearson and Spearman are two types of correlation mathamatical formulas and methods of calculating correlations. 
# This code could be expanded to other correlation methods in the future (expandable)
   if method == "pearson" or method == "spearman":
        # CorALS is a software package that works to create large-scale correlation networks when working with high-dimensional data
       corrs = CorALS(X=x, Y=y, correlation_type=method)
       return corrs

    # utilizes numpy for speed optmization ( numpy works underlyingly in C )
   corrs = np.empty((n_x, n_y), dtype=np.float64)
   # Cython decleration which improves speed
   cdef np.ndarray[np.float64_t, ndim=1] col_x
   cdef np.ndarray[np.float64_t, ndim=1] col_y


   for i in range(n_x):
       col_x = x[:, i]
       for j in range(n_y):
           col_y = y[:, j]
        # utilizes Scipy for fast implementation 
           corrs[i, j] = kendalltau(col_x, col_y)[0]


   if y is x:
       corrs = np.triu(corrs) + np.triu(corrs, k=1).T
       np.fill_diagonal(corrs, 1)


   return corrs


# ----------- Squareform Helper (upper triangle) ----------- 
# Optimization to only create the lowe diagonal of the correlation matrix so that there are no repitions when creating coorrelations. 
def squareform_cy(np.ndarray[np.float64_t, ndim=2] mat):
   cdef int i, j, idx = 0, n = mat.shape[0]
   cdef int num = (n * (n - 1)) // 2
   cdef np.ndarray[np.float64_t, ndim=1] out = np.empty(num, dtype=np.float64)
   for i in range(n):
       for j in range(i + 1, n):
           out[idx] = mat[i, j]
           idx += 1
   return out

# The full loop utilizing all the functions above 
# ----------- Main Permutation Loop -----------
def cython_loop(
   np.ndarray[np.float64_t, ndim=1] streaming_mean,
   np.ndarray[np.float64_t, ndim=1] streaming_m2,
   int num_permutations,
   int block_size,
   np.ndarray[np.float64_t, ndim=2] X,
   object inblock_permutation,
   object rng


):
    # Cython type definitions
   cdef int streaming_count = 0
   cdef int i, j, length = streaming_mean.shape[0]
   cdef double delta_val
   cdef object x_shuff, corr_mat
   cdef np.ndarray[np.float64_t, ndim=1] corrs

    # create a correlation matrix made up from the randomly shuffled matrix of Genes and Genes 
    # This would mean that the code would try to find the coorilation between 
   for i in range(num_permutations):
       x_shuff = inblock_shuffle(X, block_size=block_size, seed=rng)
       corr_mat = cy_correlation_matrix("pearson", X, x_shuff)

       corrs = squareform_cy(corr_mat)
       streaming_count += 1
       for j in range(length):
           delta_val = corrs[j] - streaming_mean[j]
           streaming_mean[j] += delta_val / streaming_count
           streaming_m2[j] += delta_val * (corrs[j] - streaming_mean[j])


       corrs = squareform_cy(corr_mat.T)
       streaming_count += 1
       for j in range(length):
           delta_val = corrs[j] - streaming_mean[j]
           streaming_mean[j] += delta_val / streaming_count
           streaming_m2[j] += delta_val * (corrs[j] - streaming_mean[j])

   return streaming_count
