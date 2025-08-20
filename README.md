# Optimized-GeneCorr
Optimized gene correlation matrix generation from count data  
*Summer Research @ Bar-Ilan Bioinformatics Lab*

---

## Bioinformatics Lab Focus

**Overall Objective:**  
Analyze gene coexpression networks (GCNs) in elderly cells compared to younger cells, with the broader aim of applying these methodologies to cancer research. Similar to aging cells, cancer cells exhibit difficulties in regulating cell division and are increasingly prone to genetic mutations.  

**Methodology:**  
1. Develop a framework to build and analyze gene co-expression networks in scRNA-seq.  
2. Generate a cell-type-specific GCN atlas across young and aged tissues.  
3. Compare GCNs and their dynamics from young to aged tissues.  

**Bioinformatics significance:**  
Accurately building and analyzing coexpression networks is essential for identifying key genes, understanding their interactions, and drawing meaningful biological conclusions. 

**Network Analysis:**  
By constructing GCNs from gene expression data of both age groups, we can investigate structural and functional differences between the networks.  

**Features of Networks for Analysis:**  
- Graph density  
- Node expression  
- Edge co-expression  
- LCC size + density  
- Effective number of components  
- Degree assortativity  
- Percolation threshold  
- Node degree  
- Node harmonic

---
## Code Description 

**Coexpression.py**

*Outputs* 
1. Coexpression Matrix: a representation of the probability of two genes coexpressed in the same cell.
2. A Corresponding P-vals matrix: an assessment of the statistical significance of the Coexpression Matrix probabilities

*Importance of Generating and utilizing P-Values* 

( This is imperetive in ensurinh that coexpression probabilities are understood, and larger networks are also understood with the proper praportions in mind and understood as they are related to the pvals(
- Pvals are created through statistical testings, including shuffling the matrix and understanding the coexpression probability for genes that are randommly selected and finding the statistical difference from the mean of the found probabilities to see to what exrtent these probabilities found should be undersood as significant.
- Includes both starteer codes to run it in python and in cython.
- Directions to run either one are commented. 
  
   
**Cython.py** 
- utilizs underlying C programming, and combines it with python to effective combine the larger python code that it is written in with the other coe. so that this couldse results in speed. 

**origional python code** 
- includes the full python for loop to create the random coexpression matrix.
- This loop includes much python overhead and slows the code down incredibly. 



