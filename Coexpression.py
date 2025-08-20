# TODO: wrap some low-level code in here
# and keep the functionality of this function
def gene_coexpressions(
    X: np.ndarray,
    correlation_method: str = settings.correlation_method,
    num_permutations: int = settings.n_permutations,
    block_size: int = settings.permutations_block_size,
    permutation_method: str = settings.permutation_method,
    sort_by: np.ndarray = None,
    seed: AnyRandom = 603550,
) -> tuple[np.ndarray, np.ndarray]:
    

    """
    Calculate gene coexpressions with permutation testing.

    Parameters:
    -----------
    X : np.ndarray
        Gene expression matrix with shape (cells, genes).
    correlation_method : str, optional
        Method to compute correlation (i.e 'pearson' or 'spearman').
        Default is from settings.
    num_permutations : int, optional
        Number of permutations for the permutation test.
        Each permutation results in two correlation matrices.
        Default is from settings.
    block_size : int, optional
        Block size for in-block shuffling during permutation.
        Default is from settings.
    sort_by : np.ndarray, optional
        Array to sort cells by, to reduce sequencing depth biases.
        Default is None.
    random_state : int, optional
        Seed for random number generator.
        Default is 0.

    Returns:
    --------
    Tuple containing:
    - corrs: 1D array (squareformed) of true correlations.
    - pvals: 1D array (squareformed) of p-values for the correlations.
    """ 
    
    rng = get_rng(seed)
    X = X.copy()
    m, n = X.shape  # m = cells, n = genes
    # For spearman correlation, rank X and than use pearson
    if correlation_method == "spearman":
        X = rankdata(X, axis=0)
        correlation_method = "pearson"
    # NOTE for Ariel: you don't have to implement "derangement" for now.
    if permutation_method == "derangement":
        inblock_permutation = utils.inblock_derangement
    elif permutation_method == "shuffle":
        inblock_permutation = utils.inblock_shuffle
    # Sort cells by n_counts
    # to reduce sequencing depth biases
    if sort_by is not None:
        X = X[np.argsort(sort_by)]
    # Init streaming variables, 1d squareformed arrays
    corrs_len = (n * n - n) // 2
    

    
    # Python Implementation StartsHere 
    streaming_count = 0
    streaming_mean  = np.zeros(corrs_len, dtype=np.float64)
    streaming_m2    = np.zeros(corrs_len, dtype=np.float64)
    for _ in range(num_permutations):
        x_shuff = inblock_permutation(X, block_size, seed=rng)
        corr_mat = correlation_matrix(correlation_method, x=X, y=x_shuff)
        # ask if the squareform is worth it. 
        corrs = squareform(corr_mat, checks=False)# check that upper is diff than lower diagonal. (false it isnt the same as the lower) 1D array
        streaming_count += 1
        delta = corrs - streaming_mean
        streaming_mean += delta / streaming_count
        streaming_m2 += delta * (corrs - streaming_mean)
        # the transpose of the correalation matrix
        # is another set of permutations
        corrs = squareform(corr_mat.T, checks=False)
        streaming_count += 1
        delta = corrs - streaming_mean
        streaming_mean += delta / streaming_count
        streaming_m2 += delta * (corrs - streaming_mean)
    streaming_std = np.sqrt(streaming_m2 / (streaming_count - 1))# Numpy knows that streaming_count is a number and
    # True corrs (1d squareformed array)
    corrs = squareform(
        correlation_matrix(method=correlation_method, x=X), checks=False
    )
    # Measure distances from mean
    distance = np.abs(corrs - streaming_mean)


 
# Cython-optimized Code Starts Here 
#     X = np.ascontiguousarray(X, dtype=np.float64)
#     streaming_mean = np.zeros(corrs_len, dtype=np.float64)
#     streaming_m2 = np.zeros(corrs_len, dtype=np.float64)
    
#     streaming_count = coexp_Cython.cython_loop(
#     streaming_mean,         # <-- 1
#     streaming_m2,           # <-- 2
#     num_permutations,       # <-- 3
#     block_size,             # <-- 4
#     X,                      # <-- 5
#     inblock_permutation,    # <-- 6
#     rng                     # <-- 7
# )
#     streaming_std = np.sqrt(streaming_m2 / (streaming_count - 1)) #Numpy knows that streaming_count is a number and
#     # True corrs (1d squareformed array)
#   # calling the cython code should generate a matrix of the random correlations between one gene and another randomly selected gene. 
#    # The below line of code should create the actual correlation between genes found in the same cell 
#     corrs = squareform(
#         correlation_matrix(method=correlation_method, x=X), checks=False
#     )
#    # Measure distances from mean
#     distance = np.abs(corrs - streaming_mean)


    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(streaming_std == 0, 10, distance / streaming_std)
        w = np.where(distance == 0, 0, w)

    pvals = 2 * (1 - norm.cdf(w))
    
    # corrs: return the true correlation found from computing the correlation between genes found in the same cell. 

    # pvals: return the "p" value, or computed significance of the true correlation found. 
    return corrs, pvals
