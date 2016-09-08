# iEDGE
**I**ntegrative analysis of (**E**pi-)**D**NA and **G**ene **E**xpression Data

iEDGE is a R-package for performing integrative analysis of epi-DNA and gene expression data. This analysis was based on [Monti, Chapuy, et. al., Cancer Cell (2012).](http://www.ncbi.nlm.nih.gov/pubmed/22975378), which describes an approach to integrate copy number alteration to gene expression dataset in Diffuse Large B cell lymphoma samples and extends the methodology to study in detail the effect of p53 and cell cycle deregulation in DLBCL. This package provides the generalized pipeline to perform integrative analysis of epi-genomic data (methylation, copy number alteration, mutation, microRNA, etc) and gene-expression data arising from overlapping samples. The package consistent of 3 major modules:
1. cis/trans differential expression based analysis of epi-DNA and gene expression data
2. pathway enrichment analysis of significant cis/trans correlations
3. pruning of cis/trans epi-DNA to gene expression connections using conditional independence criteria

Reporting options are available in table (.txt) format as well as in html format, supported by interactive plots and tables

# Installing iEDGE

(1) Please install the following Bioconductor packages first before installing iEDGE:

```
source("http://bioconductor.org/biocLite.R")
biocLite(c("Biobase", "limma", "org.Hs.eg.db", "sva"))

```

(2) Installing iEDGE

```

library(devtools)

PAT <- "04fe676593e46b6bda5a5d09431156e8a500349a"

install_github("montilab/iEDGE",auth_token = PAT)

library(iEDGE)

```

