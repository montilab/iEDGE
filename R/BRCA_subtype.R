###call BRCA subtypes from clinical file

if(FALSE){


clinical<-read.table("/Users/amyli/Desktop/git_projects/datasets/TCGA_brca/clinical_dir/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt",
 sep = "\t", header = TRUE, fill = TRUE)

mat<-clinical[-(1:2), -(1)]
colnames(mat)<-colnames(clinical)[-(1)]
rownames(mat)<-mat[,1]
mat<-data.frame(mat)

er <- as.character(mat$er_status_by_ihc)
er[!(er %in% c("Positive", "Negative"))]<- "Unknown"
pr <- as.character(mat$pr_status_by_ihc)
pr[!(pr %in% c("Positive", "Negative"))]<- "Unknown"

her2_ihc <- as.character(mat$her2_status_by_ihc)
her2_fish <- as.character(mat$her2_fish_status)
her2<-sapply(1:length(her2_ihc), 
	function(i) if (her2_ihc[i] == "Positive" | her2_fish[i] == "Positive"){
			return("Positive")
		} else if (her2_ihc[i] == "Negative" | her2_fish[i] == "Negative"){
			return("Negative")	
		} else{
			return("Unknown")
		}
	)

subtype<-sapply(1:length(er), 
	function(i) {
		er<-er[i]
		pr<-pr[i]
		her2<-her2[i]
		if( (er == "Positive" | pr == "Positive") & her2 == "Negative"){
			return("LuminalA")
		} else if ( (er == "Positive" | pr == "Positive") & her2 == "Positive"){
			return("LuminalB")
		} else if ( er == "Negative" & pr == "Negative" & her2 == "Negative"){
			return("TripleNegative")
		} else if (er == "Negative" & pr == "Negative" & her2 == "Positive"){
			return("Her2")
		} else {
			return("Unknown")
		}
	})

patient.id<-mat$bcr_patient_barcode


data.dir<-"/Users/amyli/Desktop/git_projects/datasets/TCGA_brca"
ge<-readRDS(paste(data.dir, "BRCA_rnaseqv2_rsem.RDS", sep = "/"))


unique.id<-as.character(colnames(ge))
patient.id<-as.character(sapply(unique.id,
	function(x) paste(strsplit(x, split = "-")[[1]][1:3], collapse = "-")))
sample.id<-as.character(sapply(unique.id,
	function(x) strsplit(x, split = "-")[[1]][4]))
tumor.status<-as.character(sapply(sample.id, 
	function(x){
		if (grepl("01|06", x)) return("Tumor")
		else if(grepl("11", x)) return("Normal")
		else return("Unknown")
		}))






plate<-as.character(sapply(unique.id, 
	function(x) strsplit(x, split = "-")[[1]][6]))

batch.dat<-read.table("../../datasets/TCGA_brca/BatchData.tsv", sep = "\t", fill = TRUE, header = TRUE, colClasses="character")
head(batch.dat)

batch.df<-data.frame(unique.id=batch.dat$Sample, batch.id = batch.dat$BatchId)

require(plyr)

df.ge<-data.frame(unique.id = unique.id, 
	patient.id = patient.id, 
	sample.id  = sample.id,
	tumor.status = tumor.status, 
	plate = plate)
df.ge<-join(df.ge, batch.df, by = "unique.id")

df.clinical<-data.frame(patient.id = mat$bcr_patient_barcode, er = er, pr = pr, her2 = her2, subtype = subtype)
df.clinical<-cbind(df.clinical, mat)



require(plyr)
df.clinical.join<-join(df.ge, df.clinical, by = "patient.id")

rownames(df.clinical.join)<-colnames(ge)

library(biomaRt)

df.genes<-data.frame(accession = rownames(ge))
hg19<-useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bmResult <-getBM(attributes = c('hgnc_symbol', 'wikigene_description'),
              filters = 'hgnc_symbol', values = df.genes$accession, mart = hg19)

colnames(bmResult)<-c("accession", "description")
df.genes.join<-join(df.genes, bmResult, by = "accession", 
	match = "first")
rownames(df.genes.join)<-df.genes.join$accession

brca.eset<-to_eset(assay.data = ge, pdat = df.clinical.join, fdat = df.genes.join)


saveRDS(brca.eset, file = paste(data.dir, "TCGA_brca_eset_full.RDS", sep = "/"))

#filter out samples without clinical annotation
brca.eset<-brca.eset[, !is.na(brca.eset$subtype)]
saveRDS(brca.eset, file = paste(data.dir, "TCGA_brca_eset_hasannot.RDS", sep = "/"))


brca.eset.tn<-brca.eset[, brca.eset$subtype == "TripleNegative"]

brca.eset.luminalA<-brca.eset[, brca.eset$subtype == "LuminalA"]
brca.eset.luminalB<-brca.eset[, brca.eset$subtype == "LuminalB"]
brca.eset.her2<-brca.eset[, brca.eset$subtype == "Her2"]

saveRDS(brca.eset.tn, file = paste(data.dir, "TCGA_brca_eset_tn.RDS", sep = "/"))
saveRDS(brca.eset.luminalA, file = paste(data.dir, "TCGA_brca_eset_luminalA.RDS", sep = "/"))
saveRDS(brca.eset.luminalB, file = paste(data.dir, "TCGA_brca_eset_luminalB.RDS", sep = "/"))
saveRDS(brca.eset.her2, file = paste(data.dir, "TCGA_brca_eset_her2.RDS", sep = "/"))


##save tumor samples only
brca.eset<-brca.eset[, brca.eset$tumor.status == "Tumor"]
brca.eset.tn<-brca.eset[, brca.eset$subtype == "TripleNegative"]
brca.eset.luminalA<-brca.eset[, brca.eset$subtype == "LuminalA"]
brca.eset.luminalB<-brca.eset[, brca.eset$subtype == "LuminalB"]
brca.eset.her2<-brca.eset[, brca.eset$subtype == "Her2"]

saveRDS(brca.eset.tn, file = paste(data.dir, "TCGA_brca_eset_tn_tumor.RDS", sep = "/"))
saveRDS(brca.eset.luminalA, file = paste(data.dir, "TCGA_brca_eset_luminalA_tumor.RDS", sep = "/"))
saveRDS(brca.eset.luminalB, file = paste(data.dir, "TCGA_brca_eset_luminalB_tumor.RDS", sep = "/"))
saveRDS(brca.eset.her2, file = paste(data.dir, "TCGA_brca_eset_her2_tumor.RDS", sep = "/"))



brca.eset.mad<-variationFilter(dat = brca.eset, score = "mad", 
	transform = "log2", ngenes = 500)


require(CBMRtools)

p1<-heatmap.ggplot2(eSet= brca.eset.mad, 
	col.clust = TRUE, 
	row.clust = TRUE,
	col.clust.hc = NA, row.clust.hc = NA,
	col.lab = c("plate", "er", "pr", "her2", "subtype", "tumor.status"), row.lab = "",
	heatmap.y.text = FALSE, heatmap.x.text = FALSE,
	heatmap.colorlegend.name = "RNASeq_expression",
	title.text = "TCGA BRCA log2 RNA-seq expression, z-score row normalized",
	col.legend.name = c("er", "pr", "her2", "subtype"),
	row.legend.name = "",
	row.scaling = "z-score.capped",
	z.norm = FALSE,
	cuttree.col = 4, cuttree.row = 3,
	verbose = TRUE, show = FALSE)


 	meta.c.lab<-levels(unique(p1$meta.c$id))
  	meta.c.color.string<-c("yellow", "khaki3", "gold", "chocolate", 
  		"darkred", "cyan", "white", 
  		"red", "blue", "green", "pink")
    meta.c.color<-as.character(sapply(meta.c.color.string, to.hex))


p2<-heatmap.ggplot2(eSet= brca.eset.mad, 
	col.legend.brewer = meta.c.color,
	col.clust = TRUE, 
	row.clust = TRUE,
	col.clust.hc = NA, row.clust.hc = NA,
	col.lab = c("er", "pr", "her2", "subtype"), row.lab = "",
	heatmap.y.text = FALSE, heatmap.x.text = FALSE,
	heatmap.colorlegend.name = "RNASeq_expression",
	title.text = "TCGA BRCA log2 RNA-seq expression, z-score row normalized",
	col.legend.name = c("er", "pr", "her2", "subtype"),
	row.legend.name = "",
	row.scaling = "z-score.capped",
	z.norm = FALSE,
	cuttree.col = 4, cuttree.row = 3,
	verbose = TRUE, show = FALSE)
  
ggsave(p2$heatmap, file = "BRCA_rnaseq_topmad500.pdf")



##perform paired edgeR
require(edgeR)
brca<-readRDS("../../datasets/TCGA_brca/TCGA_brca_eset_full.RDS")
brca.haspair<-brca[, brca$patient.id[which(brca$tumor.status == "Normal")]]

Subject<-factor(brca.haspair$patient.id)
Treat<-factor(brca.haspair$tumor.status)
design<-model.matrix(~Subject+Treat)

brca.haspair.mad1k<-variationFilter(dat = brca.haspair, score = "mad", 
	transform = "none", ngenes = 1000)

y <- DGEList(counts=exprs(brca.haspair.mad1k))
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)







}

