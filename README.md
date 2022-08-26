# SiFS
SiFS is an algorithm for the identification of disease-associated genes. It integrates the information obtained from protein-protein interaction network (PPI) and gene expression profiles. SiFS curates a set of potential disease-causing genes from microarray data as disease genes by maximizing both significance and functional similarity of the selected gene subset. 

# Using SiFS on Sample Data
You can test the source using the sample data set as follows:
1. gcc SiFS.c -lm -o sifs.out
2. ./sifs.out Data/GeneExpr.txt Data/PPI_Network.txt Data/Subset.txt <#genes to be selected> <weight parameter> <Output File>

# Using SiFS on Your Own Data
The gene expression and PPI network must depict the expression profiles and interaction between of genes that lie within their intersection set. The gene expression data must have the following representation:
  ROW_1   #samples  #Genes  #Classes
  ROW_2   HGNC_ID_1  EXPR_1 EXPR_2 ... EXPR_P
  ROW_3   HGNC_ID_2  EXPR_1 EXPR_2 ... EXPR_P
  .
  .
  .
  ROW_N   HGNC_ID_2  EXPR_1 EXPR_2 ... EXPR_P
The sample expression file can also be referred in order to understand the representation format.

The PPI network file must have a numeric representation as follows:
  ROW1    #Genes  #Interactions
  ROW2       1            2         Weight of interaction between HGNC_ID_I and HGNC_ID_2 (if weight > 0)
  ROW3       1            3         Weight of interaction between HGNC_ID_1 and HGNC_ID_3
  .
  .
  .
  
So, each row depicts the interaction partners and the weight of their interaction. The sample PPI Network file can also be referred in order to understand the representation format.

The SiFS can be used to curate genes from a subset of genes of the gene expression data or from the entire set. If using SiFS to select genes from a subset of genes, then they must be provided in Subset.txt or path to relevant file must be provided. On the other hand, to use SiFS on the entire data set, Subset.txt must have the complete set of genes or a path to the file must be given.
  
The weight parameter controls the the relative importance of significance and functional similarity of the candidate gene with respect to the already-selected genes. Its values must lie within the range [0,1].
