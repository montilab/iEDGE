library(iEDGE)

args<-commandArgs(TRUE)

#f_in is one of RDS file created by 1.preprocess_all_TCGA.R
f_in<-args[1]
#f_out is the path to the output directory
f_out<-args[2]

cat(paste("f_in: ", f_in, "\n", sep = ""))
cat(paste("f_out: ", f_out, "\n", sep = ""))

header<-paste(gsub(".RDS", "", basename(f_in)), sep = "")
cat(paste("header: ", header, "\n", sep = ""))

dat<-readRDS(f_in)

#change to actual path
#download gmt files from MSigDB
gs.dir<-"path_to_geneset_gmt_files"
gs.names<-c("h.all.v5.0.symbols.gmt","c2.cp.v5.0.symbols.gmt", "c3.tft.v5.0.symbols.gmt")
gs.names<-paste(gs.dir, gs.names, sep = "/")

res<-run_iEDGE(dat, header, f_out, 
	gs.file = gs.names, 
	gepid = "gene_symbol", cnid = "Unique.Name", 
	cndir = "alteration_direction", fdr.cis.cutoff = 0.25, 
	fdr.trans.cutoff = 0.01, fc.cis = 1.2, fc.trans = 1.5, min.drawsize = 3, onesided.cis = TRUE, 
	onesided.trans = FALSE, uptest = "Amplification", downtest = "Deletion", 
	min.group = 2,
	prune.col = "fdr", prune.thres = 0.05, hyperthres = 0.25,
	cis.boxplot = TRUE, trans.boxplot = FALSE, bipartite = TRUE, enrichment = TRUE, html = TRUE)


