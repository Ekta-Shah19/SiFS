# SiFS
SiFS is an algorithm for the identification of disease-associated genes. It integrates the information obtained from protein-protein interaction network (PPI) and gene expression profiles. SiFS curates a set of potential disease-causing genes from microarray data as disease genes by maximizing both significance and functional similarity of the selected gene subset. 

# Using SiFS on Sample Data
You can test the source using the sample data set as follows:<br />
1. gcc&nbsp;SiFS.c&nbsp;-lm&nbsp;-o&nbsp;sifs.out<br />
2. ./sifs.out&nbsp;Data/GeneExpr.txt&nbsp;Data/PPI_Network.txt&nbsp;Data/Subset.txt&nbsp;#Genes_to_be_selected &nbsp;weight_parameter&nbsp;Output_File<br />

# Using SiFS on Your Own Data
The gene expression and PPI network must depict the expression profiles and interaction between of genes that lie within their intersection set. The gene expression data must have the following representation:<br />
  ROW_1 &emsp;&emsp;  #samples &emsp; #Genes &emsp; #Classes<br />
  ROW_2  &emsp;&emsp; HGNC_ID_1 &nbsp; EXPR_1 &nbsp; EXPR_2 &nbsp;...&nbsp; EXPR_P<br />
  ROW_3 &emsp;&emsp;  HGNC_ID_2 &nbsp; EXPR_1 &nbsp; EXPR_2 &nbsp;...&nbsp; EXPR_P<br />
  .<br />
  .<br />
  .<br />
  ROW_N &emsp;&emsp;  HGNC_ID_2 &nbsp; EXPR_1 &nbsp; EXPR_2 &nbsp;...&nbsp; EXPR_P<br />
The sample expression file can also be referred in order to understand the representation format.<br />

The PPI network file must have a numeric representation as follows:<br />
  ROW1&emsp;&emsp;#Genes&emsp;#Interactions<br />
  ROW2&emsp;&emsp;1&emsp;2&emsp;Weight of interaction between HGNC_ID_I and HGNC_ID_2 (if weight > 0)<br />
  ROW3&emsp;&emsp;1&emsp;3&emsp;Weight of interaction between HGNC_ID_1 and HGNC_ID_3<br />
  .<br />
  .<br />
  .<br />
  <br />
So, each row depicts the interaction partners and the weight of their interaction. The sample PPI Network file can also be referred in order to understand the representation format.<br />

The SiFS can be used to curate genes from a subset of genes of the gene expression data or from the entire set. If using SiFS to select genes from a subset of genes, then they must be provided in Subset.txt or path to relevant file must be provided. On the other hand, to use SiFS on the entire data set, Subset.txt must have the complete set of genes or a path to the file must be given.<br />
  
The weight parameter controls the the relative importance of significance and functional similarity of the candidate gene with respect to the already-selected genes. Its values must lie within the range [0,1].<br />
