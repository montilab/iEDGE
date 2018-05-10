# iEDGE
**I**ntegrative analysis of (**E**pi-)**D**NA and **G**ene **E**xpression Data

iEDGE is a R-package for performing integrative analysis of epi-DNA and gene expression data. This analysis was based on [Monti, Chapuy, et. al., Cancer Cell (2012).](http://www.ncbi.nlm.nih.gov/pubmed/22975378), which describes an approach to integrate copy number alteration to gene expression dataset in Diffuse Large B cell lymphoma samples and extends the methodology to study in detail the effect of p53 and cell cycle deregulation in DLBCL. This package provides the generalized pipeline to perform integrative analysis of epi-genomic data (methylation, copy number alteration, mutation, microRNA, etc) and gene-expression data arising from overlapping samples. The package consistent of 3 major modules:
1. cis/trans differential expression based analysis of epi-DNA and gene expression data
2. pathway enrichment analysis of significant cis and trans genes correlated with the epi-DNA alterations
3. pruning of cis/trans epi-DNA to gene expression connections using a statistical mediation model

Reporting options are available in table (.txt) format as well as in html format, supported by interactive plots and tables

# Installing iEDGE

```

library(devtools)

#contact maintainer for PAT
install_github("montilab/iEDGE",auth_token = PAT)

library(iEDGE)

```

