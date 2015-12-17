
if (FALSE){
brca<-readRDS("../../datasets/TCGA_brca/TCGA_brca_eset_hasannot.RDS")


fdat<-fData(brca)


genes<-as.character(fdat$accession)
genes.er<-c("ESR1") #6q25.1
genes.pr<-c("PGR")
genes.her2<-c("ERBB2") 

ind<-which(genes %in% c(genes.er, genes.pr, genes.her2))

brca.sub<-brca[rev(match(c(genes.er, genes.pr, genes.her2), genes)),]

require(CBMRtools)

exprs(brca.sub)<-log2(exprs(brca.sub)+1)
p1<-heatmap.ggplot2(eSet= brca.sub, 
	col.clust = TRUE, 
	row.clust = FALSE,
	col.clust.hc = NA, row.clust.hc = NA,
	col.lab = c("er", "pr", "her2", "subtype"), 
	row.lab = c("accession"),
	heatmap.y.text = TRUE, heatmap.x.text = FALSE,
	heatmap.colorlegend.name = "RNASeq_expression",
	title.text = "TCGA BRCA log2 RNA-seq expression, z-score row normalized",
	col.legend.name = c("er", "pr", "her2", "subtype"),
	row.legend.name = c("accession"),
	row.scaling = "none",
	z.norm = FALSE,
	cuttree.col = "", cuttree.row = "",
	verbose = FALSE, show = FALSE)

p1




}









