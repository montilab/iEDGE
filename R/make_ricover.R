#make_ricover.R
library(Biobase)
library(biomaRt)
f.dir.in<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/wes_7_16_15"

eset.focal<-readRDS(paste(f.dir.in, "/", "ricover_wes_gistic_all_focal.RDS", sep = ""))
eset.focal.binary<-readRDS(paste(f.dir.in, "/", "ricover_wes_gistic_all_focal_binary.RDS", sep = ""))
eset.arm<-readRDS(paste(f.dir.in, "/", "ricover_wes_gistic_all_arm.RDS", sep = ""))
eset.arm.binary<-readRDS(paste(f.dir.in, "/", "ricover_wes_gistic_all_arm_binary.RDS", sep = ""))
eset.or<-readRDS(paste(f.dir.in, "/", "ricover_wes_gistic_all_or.RDS", sep = ""))
eset.or.binary<-readRDS(paste(f.dir.in, "/", "ricover_wes_gistic_all_or_binary.RDS", sep = ""))

eset.gep<-readRDS(paste(f.dir.in, "/", "ricover_gep.RDS", sep = ""))
eset.gep.log<-readRDS(paste(f.dir.in, "/", "ricover_gep_log.RDS", sep = ""))
eset.gep.log.corrected<-readRDS(paste(f.dir.in, "/", "ricover_gep_log_correct.RDS", sep = ""))

eset.cisgenes<-readRDS(paste(f.dir.in, "/", "ricover_wes_gistic_cisgenes.RDS", sep = ""))

##biomart for genes on each arm

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
 
# Only use standard human chromosomes
normal.chroms <- c(1:22, "X", "Y", "M")
 
# Filter on HGNC symbol and chromosome, retrieve genomic location and band
my.symbols <- fData(eset.gep)$accession
 
my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                    filters = c("hgnc_symbol", "chromosome_name"),
                    values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),
                                  mart = ensembl)

my.regions$arm<-paste(my.regions$chromosome_name,  sapply(my.regions$band, function(x) substr(x, 1, 1)), sep ="")

get_genes_in_arm<-function(x,y){
	#x<-"5p"
	return(y$hgnc_symbol[which(y$arm == x)])
}

get_genes_in_arm_in_focal<-function(x, eset.focal, eset.cisgenes){
	cytoband.focal<-fData(eset.focal)$Descriptor
	ind.focal<-grep(x, cytoband.focal)
	genes.focal<-eset.cisgenes$genes[ind.focal]
	return(unique(do.call(c,genes.focal)))
}

eset.cisgenes.arm<-list()
eset.cisgenes.arm$meta<-fData(eset.arm)
eset.cisgenes.arm$genes<-lapply(1:nrow(eset.cisgenes.arm$meta), 
	function(i){
		get_genes_in_arm(x=eset.cisgenes.arm$meta$Descriptor[i], y = my.regions)
		})


##in arm also in focal

dir.focal<-sapply(fData(eset.focal)$Unique.Name, function(x){
	if (grepl("Amp", x)) return("Amplification") else return("Deletion")
	})
arm.focal<-sapply(fData(eset.focal)$Descriptor, function(x){
	if(grepl("p", x))
		paste(strsplit(x, split = "p")[[1]][1], "p", sep = "")
	else 
		paste(strsplit(x, split = "q")[[1]][1], "q", sep = "")
	
	})
both.focal<-paste(dir.focal, arm.focal, sep = "_")

dir.arm<-sapply(fData(eset.arm)$Unique.Name, function(x){
	if (grepl("Amp", x)) return("Amplification") else return("Deletion")
	})
arm.arm<-sapply(fData(eset.arm)$Descriptor, function(x){
	if(grepl("p", x))
		paste(strsplit(x, split = "p")[[1]][1], "p", sep = "")
	else 
		paste(strsplit(x, split = "q")[[1]][1], "q", sep = "")
	
	})
both.arm<-paste(dir.arm, arm.arm, sep = "_")

match.focal<-match(both.focal, both.arm)
ind.focal<-which(!is.na(match.focal))

eset.cisgenes.arm.infocal<-list()
eset.cisgenes.arm.infocal$meta<-fData(eset.focal)[ind.focal,]
eset.cisgenes.arm.infocal$genes<-eset.cisgenes$genes[ind.focal]

ind.arm.infocal<-match.focal[!is.na(match.focal)]

samples.gep<-colnames(eset.gep)
samples.cn<-colnames(eset.focal)

samples.map<-read.table("/Users/amyli/Desktop/git_projects/datasets/ricover2015/IDmapping.txt",
	sep = "\t", header = TRUE)

for(i in colnames(samples.map)){
	samples.map[,i]<-gsub("-", ".", samples.map[,i])
}

samples.map$wes_id<-match(samples.map$WES.ID, samples.cn)
samples.map$gep_id<-match(samples.map$caseID, samples.gep)

samples.map.sub<-subset(samples.map, !is.na(wes_id) & !is.na(gep_id))

eset.arm.infocal<-eset.arm[ind.arm.infocal,]
fData(eset.arm.infocal)<-fData(eset.focal)[ind.focal,]

eset.arm.infocal.binary<-eset.arm.binary[ind.arm.infocal,]
fData(eset.arm.infocal.binary)<-fData(eset.focal)[ind.focal,]


ricover_wes<-list(
	cn.focal = eset.focal[, samples.map.sub$wes_id],
	cn.focal.binary = eset.focal.binary[, samples.map.sub$wes_id],
	cn.arm = eset.arm[, samples.map.sub$wes_id],
	cn.arm.binary = eset.arm.binary[, samples.map.sub$wes_id],
	cn.focal.or.arm = eset.or[, samples.map.sub$wes_id],
	cn.focal.or.arm.binary = eset.or.binary[, samples.map.sub$wes_id],
	cn.arm.infocal = eset.arm.infocal[, samples.map.sub$wes_id],
	cn.arm.infocal.binary = eset.arm.infocal.binary[, samples.map.sub$wes_id],
	gep = eset.gep[, samples.map.sub$gep_id],
	gep.log = eset.gep.log[, samples.map.sub$gep_id],
	gep.log.corrected = eset.gep.log.corrected[, samples.map.sub$gep_id],
	gep.cisgenes = eset.cisgenes,
	gep.cisgenes.arm = eset.cisgenes.arm,
	gep.cisgenes.arm.infocal = eset.cisgenes.arm.infocal
	)

saveRDS(ricover_wes, file = paste(f.dir.in,"/", "ricover_wes_all.RDS", sep = ""))

