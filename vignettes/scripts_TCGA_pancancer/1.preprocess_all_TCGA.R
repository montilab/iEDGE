library(Biobase)
library(iEDGE)


#preprocess data 

#change to actual path
gepdir<-"path_to_rnaseq_esets_rds_files"
gepfns<-list.files(gepdir, pattern = ".rds$")

cancertypes<-unique(sapply(gepfns, function(i) strsplit(i, split = "_")[[1]][1]))
gepdate<-"2015_02_04"
geplist<-lapply(cancertypes, function(i){
	gep<-readRDS(paste(gepdir, "/", i, "_", gepdate, "_ES.rds", sep = ""))
	gep<-gep[!is.na(fData(gep)$gene_symbol),]
	gep<- gep[, pData(gep)[, "tissue_type"] %in% "Tumor"]
	gep<-gep[which( apply(exprs(gep), 1, var) != 0),]
	exprs(gep)<-log2(exprs(gep)+1)
	return(gep)
	})

cndir<-lapply(cancertypes, function(i){
	#change to actual path
	fn<-paste("path_to_firehose_tcga_copynumber_data",
		i,"/analyses__*/",
		i, "/*/gdac.broadinstitute.org_", i, 
		"-TP.CopyNumber_Gistic2.Level_4.*00.0.0", sep = "")
	#find directories matching pattern
	res<-Sys.glob(file.path(pattern =fn))
	if(length(res) == 0) return(NA)
	#if one than more, return last, should be sorted by date
	if(length(res) > 1) return(res[length(res)]) 
	else return(res)
	})
cancertypeswcn<-cancertypes[!is.na(cndir)]
cndirf<-cndir[!is.na(cndir)]
gepf<-geplist[!is.na(cndir)]

cnf<-lapply(1:length(cndirf), function(i){
	print(i)
	gistic_in<-cndirf[[i]]
	my.symbols<-fData(gepf[[i]])[, "gene_symbol"]
	gisticdata<-make_GISTIC2(gistic_in = gistic_in, 
		all_genes = my.symbols, 
		amp_thres_arm = 0.1, del_thres_arm = -0.1, 
		amp_qvalue_arm = 0.25, del_qvalue_arm = 0.25, 
		qvalue_focal = 0.25,
		all_lesions = paste(gistic_in, "/", "all_lesions.conf_99.txt", sep = ""), 
		f.amp = paste(gistic_in,"/", "amp_genes.conf_99.txt", sep = ""), 
		f.del = paste(gistic_in, "/", "del_genes.conf_99.txt", sep =""),
		broad_thres = TRUE
		)
	return(gisticdata)
	})

#run local only
cnffocal<-lapply(cnf, function(i) return(i$focal))

OUTDIR<-paste("../data/combined")
dir.create(OUTDIR, recursive = TRUE)

#save combined
lapply(1:length(gepf), function(i){
	gep<-gepf[[i]]
	gepid<-as.character(pData(gep)[, "bcr_sample_barcode"])
	colnames(gep)<-gepid

	cn<-cnffocal[[i]]
	cnid<-colnames(cn$cn)
	cnshort<-gsub("\\.", "-", substr(cnid, 1, 15))
	mapping<-data.frame(CNID = cnid, GEPID = gepid[match(cnshort, gepid)])
	mapping<-mapping[apply(mapping, 1, function(i) !any(is.na(i))), ]

	header <- paste("TCGA_", cancertypeswcn[i], sep ="")
	dat<-iEDGE_combine(alt = cn, gep = gep, mapping = mapping,
		mapping.cn = "CNID", mapping.gep = "GEPID", uppercase = TRUE)

	OUTDIR<-paste("../data/combined")
	dir.create(OUTDIR, recursive = TRUE)
	OUTFILE<-paste(OUTDIR, "/", header, ".RDS", sep = "")
	saveRDS(dat, OUTFILE)
	})

