require(CBMRtools)
require(Biobase)

read_GISTIC1<-function(all_lesions){

	x<-read.table(all_lesions,sep = "\t", header =T)
	x.meta.ind<-1:which(colnames(x) == "Amplitude.Threshold")
	##remove actual copy number rows
	x<-x[x$Amplitude.Threshold != "Actual Log2 Ratio Given",]
	#make eset
	fdat<-x[, x.meta.ind]
	for(i in colnames(fdat)){
		fdat[, i]<-gsub(" ", "", fdat[,i])
	}
	dat<-x[, -x.meta.ind]
	pdat<-data.frame(sample_id = colnames(dat))
	eset<-to.eSet(mat = dat, pdat = pdat, fdat = fdat)
	return(eset)
}

focal_or_arm<-function(eset.focal, eset.arm){
	eset.or<-eset.focal
	for(ind.focal in 1:nrow(eset.focal)){
		name.focal<-fData(eset.focal)$Unique.Name[ind.focal]
		name.arm<-paste(name.focal, "-", "ContributingBroadEvents", sep = "")
		ind.arm<-which(fData(eset.arm)$Unique.Name == name.arm)
		if(length(ind.arm)!=0){
			##has matching arm level copy number alteration
			x.both<-rbind(exprs(eset.focal[ind.focal,]), exprs(eset.arm[ind.arm,]))
			x.or<-apply(x.both, 2, max) #take max of copy number status
			exprs(eset.or)[ind.focal,]<-x.or
			fData(eset.or)$Unique.Name[ind.focal]<-paste(fData(eset.or)$Unique.Name[ind.focal],
				"_or_arm", sep = "")
		}
	}
	return(eset.or)
}


get_cis_genes<-function(f.genes, eset.focal, direction){

	##get amp genes 
	res<-read.table(f.genes, sep = "\t", header = FALSE, fill = T)
	res<- res[, apply(res, 2, function(x) !all(is.na(x)))]
	last.meta.ind<-which(res[,1] == "genes in wide peak")-1
	
	res.meta<-t(res[1:last.meta.ind,])
	colnames(res.meta)<-res.meta[1,]
	res.meta<-data.frame(res.meta[-1,])
	for(i in colnames(res.meta)){
		res.meta[,i]<-gsub(" ", "", res.meta[,i])
	}
	res.meta<-data.frame(direction = direction, res.meta)
	res.genes<-res[-(1:last.meta.ind), -1]
	res.genes<-sapply(1:ncol(res.genes),
		function(x) setdiff(unique(as.character(res.genes[, x])), ""))

	res.genes<-lapply(res.genes, function(x) gsub("\\[|\\]", "", x))
	
	return(list(meta = res.meta, genes = res.genes))
}


if(TRUE){

	f.dir.in<-"/Users/amyli/Desktop/git_projects/datasets/lymphoma2010/GISTIC_20110620.4amy"
	f.dir.out<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/cancercell2012"
	all_lesions<-paste(f.dir.in, "/", "all_lesions_file.txt", sep = "")
	eset<-read_GISTIC1(all_lesions = all_lesions)

	#separate into focal, arm or either
	eset.focal<-eset[!grepl("ContributingBroadEvents", fData(eset)$Unique.Name),]
	eset.arm<-eset[grepl("ContributingBroadEvents", fData(eset)$Unique.Name),]
	eset.or<-focal_or_arm(eset.focal, eset.arm)

	#save gistic files as R objects (RDS)
	saveRDS(eset.focal, file = paste(f.dir.out, "/", "cancercell2012_gistic_all_focal.RDS", sep = ""))
	saveRDS(eset.arm, file = paste(f.dir.out, "/", "cancercell2012_gistic_all_arm.RDS", sep = ""))
	saveRDS(eset.or, file = paste(f.dir.out, "/", "cancercell2012_gistic_all_or.RDS", sep = ""))

	#binary phenotype

	mat<-exprs(eset)
	mat[mat == 2]<-1
	exprs(eset)<-mat
	eset.focal<-eset[!grepl("ContributingBroadEvents", fData(eset)$Unique.Name),]
	eset.arm<-eset[grepl("ContributingBroadEvents", fData(eset)$Unique.Name),]
	eset.or<-focal_or_arm(eset.focal, eset.arm)

	#save gistic files as R objects (RDS)
	saveRDS(eset.focal, file = paste(f.dir.out, "/", "cancercell2012_gistic_all_focal_binary.RDS", sep = ""))
	saveRDS(eset.arm, file = paste(f.dir.out, "/", "cancercell2012_gistic_all_arm_binary.RDS", sep = ""))
	saveRDS(eset.or, file = paste(f.dir.out, "/", "cancercell2012_gistic_all_or_binary.RDS", sep = ""))

	#read cis genes

	eset.focal<-readRDS(paste(f.dir.out, "/", "cancercell2012_gistic_all_focal_binary.RDS", sep = ""))

	f.amp<-paste(f.dir.in, "/", "Amp_genes.txt", sep = "")
	f.del<-paste(f.dir.in, "/", "Del_genes.txt", sep = "")
	genes.amp<-get_cis_genes(f.genes = f.amp, eset.focal, direction = "Amplification")
	genes.del<-get_cis_genes(f.genes = f.del, eset.focal, direction = "Deletion")
	
	genes.both<-list(meta = rbind(genes.amp[["meta"]], genes.del[["meta"]]),
		genes = c(genes.amp[["genes"]], genes.del[["genes"]]))

	sample_meta<-fData(eset.focal)
	sample_meta$wide.peak.boundaries<-as.character(sapply(sample_meta$Wide.Peak.Limits, 
		function(x) strsplit(x, split = "\\(")[[1]][1]))
	gene_meta<-genes.both[['meta']]
	idx<-match(sample_meta$wide.peak.boundaries, gene_meta$wide.peak.boundaries)
	genes.both[['meta']]<-genes.both[['meta']][idx,]
	genes.both[['genes']]<-genes.both[['genes']][idx]
	
	saveRDS(genes.both, file = paste(f.dir.out, "/", "cancercell2012_gistic_cisgenes.RDS", sep = ""))



}