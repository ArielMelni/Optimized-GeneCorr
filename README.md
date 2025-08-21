# Optimized-GeneCorr
Optimized gene correlation matrix generation from count data  
*Summer Research @ Bar-Ilan Bioinformatics Lab*

---

## 🔬 Bioinformatics Lab Focus

**Abstract**  

Quantifying heterogeneity, the diversity of gene expression patterns across cells, is a critical tool in bioinformatics because it provides insights into the mechanisms underlying genetic change. For instance, older cells show elevated heterogeneity, revealing a declined performance of their cellular processes such as DNA replication and repair. Since cancerous cells exhibit similarly Impaired capacity, this study utilizes aging as a model system for cancer research. Methodologically, this study develops a framework to build single-cell gene coexpression networks (GCNs) from scRNA-seq data. In these networks, genes are represented as nodes and edges between nodes indicate the likelihood of coexpression. Comparing the divergence of network characteristics from the GCNs of young and old cells provide insights into the biological systems that drive genetic alterations during aging. Applying this methodology to analyze single-cell coexpression networks in cancerous cells represents a promising avenue for future research in bioinformatics and cancer biology.
**The Problem** 
- The code that creates GCNs is done in Python creating overhead that slows down the creation of GCN's Significantly.

**My Objective** 
- Develop a solution to Python's impaired capacity and optimize the current speed of creating a GCN. 

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

### **coexp_Cython.pyx**  
**Key Features That Make the Cython Code Faster Than Python**

| Feature | Benefit |
|---------|---------|
| `cdef` static typing | Converts Python objects/loops to C-level operations |
| Typed NumPy arrays | Fast, contiguous memory access |
| Cython loops vs Python loops | Avoids Python interpreter overhead |
| Streaming mean/variance | Saves memory and avoids storing all permutation matrices |
| Optimized libraries (`CorALS`) | Fast correlation calculation for high-dimensional data |

- In order to run the code, first compilation is required, which allows for much of the speed optimizations. 

---

## ⏱ Runtime Comparison & Optimization

The performance improvement of the Cython implementation over the pure Python code was evaluated using Python’s built-in `timeit` module. Each implementation was executed repeatedly, and the overall execution time was recorded to ensure that the observed speedup is consistent and not due to random variation.

### **Results**

Cython saw a 34 Second Speed optimization: 

*Python implementation runtime across multiple runs: **1 min and 57 seconds** *

![Python Runtime](https://github.com/ArielMelni/Optimized-GeneCorr/blob/main/8C38DD30-33A5-4112-BC58-02A0257EE8B9.jpeg )  

*Cython implementation runtime across multiple runs. **1 min and 23 seconds** * 

![Cython Runtime](https://github.com/ArielMelni/Optimized-GeneCorr/blob/main/44D4B4DF-6195-4957-BF11-4B4944283085.jpeg)  



---




