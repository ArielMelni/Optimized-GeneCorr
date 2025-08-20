# Optimized-GeneCorr
Optimized gene correlation matrix generation from count data  
*Summer Research @ Bar-Ilan Bioinformatics Lab*

---

## 🔬 Bioinformatics Lab Focus

**Overall Objective**  
Analyze gene coexpression networks (GCNs) in elderly cells compared to younger cells, with the broader aim of applying these methodologies to cancer research. Similar to aging cells, cancer cells exhibit difficulties in regulating cell division and are increasingly prone to genetic mutations.  

**Methodology**  
1. Develop a framework to build and analyze gene coexpression networks in scRNA-seq.  
2. Generate a cell-type-specific GCN atlas across young and aged tissues.  
3. Compare GCNs and their dynamics between young and aged tissues.  

**Bioinformatics Significance**  
Accurately building and analyzing coexpression networks is essential for identifying key genes, understanding their interactions, and drawing meaningful biological conclusions.  

**Network Analysis Metrics**  
To compare networks, we examine structural and functional properties such as:  
- Graph density  
- Node expression  
- Edge coexpression  
- LCC size and density  
- Effective number of components  
- Degree assortativity  
- Percolation threshold  
- Node degree distribution  
- Node harmonic centrality  

---

## 💻 Code Description

### **Coexpression.py**

**Parameters**  
- **X** (*np.ndarray*)  
  Gene expression matrix with shape `(cells, genes)`.  

- **correlation_method** (*str*, optional)  
  Method to compute correlation (e.g., `'pearson'` or `'spearman'`).  
  Default is from settings.  

- **num_permutations** (*int*, optional)  
  Number of permutations for the permutation test.  
  Each permutation results in two correlation matrices.  
  Default is from settings.  

- **block_size** (*int*, optional)  
  Block size for in-block shuffling during permutation.  
  Default is from settings.  

- **sort_by** (*np.ndarray*, optional)  
  Array to sort cells by, to reduce sequencing depth biases.  
  Default is `None`.  

- **random_state** (*int*, optional)  
  Seed for random number generator.  
  Default is `0`.  

**Outputs**  
1. **Coexpression Matrix** – the probability that two genes are coexpressed in the same cell.  
2. **P-values Matrix** – the statistical significance of the coexpression matrix probabilities.  

**Generating and Utilizing P-values**  
- P-values are computed through permutation-based statistical testing. This involves shuffling the gene matrix and recalculating coexpression probabilities for randomly selected genes.  
- The probability of two random genes being coexpressed is compared against the observed coexpression probability of two genes found in the same cell.  
- **P-value generation enables networks to capture statistical significance by correcting for random coexpression, rather than relying solely on raw coexpression probabilities.**

**Additional Features**  
- Includes code to run **either Cython-based loops** or **pure Python loops**, with both implementations available and commented for flexibility.  

---

### **Cython.py**  
**Key Features That Make the Cython Code Faster Than Python**

| Feature | Benefit |
|---------|---------|
| `cdef` static typing | Converts Python objects/loops to C-level operations |
| Typed NumPy arrays | Fast, contiguous memory access |
| Cython loops vs Python loops | Avoids Python interpreter overhead |
| Streaming mean/variance | Saves memory and avoids storing all permutation matrices |
| Optimized libraries (`CorALS`) | Fast correlation calculation for high-dimensional data |


---

## ⏱ Runtime Comparison & Optimization

The performance improvement of the Cython implementation over the pure Python code was evaluated using Python’s built-in `timeit` module. Each implementation was executed repeatedly, and the overall execution time was recorded to ensure that the observed speedup is consistent and not due to random variation.

### **Results**

Below are screenshots of the runtime tests, with Cython showing an overall 30 second speed optimization:

![Python Runtime](path/to/python_runtime.png)  
*Python implementation runtime across multiple runs.*

![Cython Runtime](path/to/cython_runtime.png)  
*Cython implementation runtime across multiple runs.*


---




