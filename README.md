# iEDGE
**I**ntegrative analysis of (**E**pi-)**D**NA and **G**ene **E**xpression Data

iEDGE is an R-package for performing integrative analysis of (epi-)DNA
and gene expression data that builds upon, and significantly expands,
the methodology first described in [Monti, Chapuy, et. al., Cancer
Cell (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22975378). This
package provides the generalized pipeline to perform integrative
analysis of epi-genomic data (methylation, copy number alteration,
mutation, microRNA, etc.) and gene-expression data from
paired samples, and consists of 3 major modules:
1. cis/trans differential expression analysis of epi-DNA and gene expression data
2. pathway enrichment analysis of significant cis and trans genes associated with the epi-DNA alterations
3. pruning of cis/trans epi-DNA to gene expression connections using a statistical mediation model

Reporting options are available in table (.txt) format as well as in html format, supported by interactive plots and tables

# Installing iEDGE

```

library(devtools)

install_github("montilab/iEDGE")

library(iEDGE)

```

