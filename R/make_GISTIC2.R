#require(Biobase)

#' @import Biobase
to.eSet<-function(mat, pdat, fdat){
	#require(Biobase)
	mat<-as.matrix(mat)

	#checking data type and dimensions
	if (!is.data.frame(pdat))
		stop("pdat must be a data frame")
	
	if (!is.data.frame(fdat))
		stop("fdat must be a data frame")
	
	if ( nrow(fdat) != nrow(mat))
		stop("nrow(fdat) must equal nrow(mat)")
	
	if( nrow(pdat) != ncol(mat))
		stop("nrow(pdat) must equal ncol(mat)")

	if (!all(rownames(fdat) == rownames(mat))){
		warning("fdat rownames and mat rownames do not match, setting fdat rownames to mat rownames")
		rownames(fdat) <- rownames(mat)
	}

	if (!all(rownames(pdat) == colnames(mat))){
		warning("pdat rownames and mat colnames do not match, setting fdat rownames to mat rownames")
		rownames(pdat) <- colnames(mat)
	}

	fMetaData<-data.frame(labelDescription = colnames(fdat), row.names = colnames(fdat))
	featureData<-new("AnnotatedDataFrame", data= fdat, varMetadata=fMetaData) 

	pMetaData<-data.frame(labelDescription = colnames(pdat), row.names = colnames(pdat))
	phenoData<-new("AnnotatedDataFrame", data= pdat, varMetadata=pMetaData) 
	
	eSet<-ExpressionSet(assayData=mat, featureData = featureData, phenoData = phenoData,  annotation = "")
	return(eSet)
}

read_GISTIC2_focal<-function(all_lesions){
	x<-read.table(all_lesions,sep = "\t", header =T)
	x<-x[, !apply(x, 2, function(i){all(is.na(i))})]
	x.meta.ind<-1:which(colnames(x) == "Amplitude.Threshold")
	x<-x[x$Amplitude.Threshold != "Actual Copy Change Given",]
	fdat<-x[, x.meta.ind]
	for(i in colnames(fdat)){
		fdat[, i]<-gsub(" ", "", fdat[,i])
	}
	dat<-x[, -x.meta.ind]
	dat[dat == 2]<-1
	pdat<-data.frame(sample_id = colnames(dat))
	eset<-to.eSet(mat = dat, pdat = pdat, fdat = fdat)
	return(eset)
}

#' @import Biobase
read_GISTIC2_focal_with_thres<-function(all_lesions, amp_thres =0.3, del_thres = -0.3){
	x<-read.table(all_lesions,sep = "\t", header =T)
	x<-x[, !apply(x, 2, function(i){all(is.na(i))})]
	x.meta.ind<-1:which(colnames(x) == "Amplitude.Threshold")
	x<-x[x$Amplitude.Threshold == "Actual Copy Change Given",]
	fdat<-x[, x.meta.ind]
	for(i in colnames(fdat)){
		fdat[, i]<-gsub(" ", "", fdat[,i])
	}
	dat<-x[, -x.meta.ind]
	pdat<-data.frame(sample_id = colnames(dat))
	eset<-to.eSet(mat = dat, pdat = pdat, fdat = fdat)

	fData(eset)$Unique.Name<-gsub("-CNvalues", "", fData(eset)$Unique.Name)
	dirs<-grep("Amp", fData(eset)$Unique.Name)

	ncols<-ncol(eset)
	mat<-t(sapply(1:nrow(eset), function(i){
		is_amp<-grepl("Amp", fData(eset)$Unique.Name[i])
		x<-as.numeric(exprs(eset)[i,, drop = FALSE])
		res<-rep(0, ncols)

		if(is_amp)
		res[x>amp_thres]<-1
		else
		res[x<del_thres]<-1
		return(res)
	}))
	exprs(eset)<-mat
	return(eset)	
}

read_GISTIC2_broad<-function(broad_lesions){
	x<-read.table(broad_lesions,sep = "\t", header =T)
	x<-x[, !apply(x, 2, function(i){all(is.na(i))})]
	x.meta.ind<-which(colnames(x) == "Chromosome.Arm")
	fdat<-data.frame(Descriptor = x[, x.meta.ind])
	for(i in colnames(fdat)){
		fdat[, i]<-gsub(" ", "", fdat[,i])
	}
	dat<-x[, -x.meta.ind]
	pdat<-data.frame(sample_id = colnames(dat))
	eset<-to.eSet(mat = dat, pdat = pdat, fdat = fdat)

	return(eset)
}


#' @import Biobase
read_GISTIC2_broad_with_thres<-function(broad_lesions, 
	broad_lesions_threshold,
	amp_qvalue = 0.25, del_qvalue = 0.25,
	amp_thres =0.3, del_thres = -0.3){
	x<-read.table(broad_lesions,sep = "\t", header =T)
	x<-x[, !apply(x, 2, function(i){all(is.na(i))})]
	x.meta.ind<-which(colnames(x) == "Chromosome.Arm")
	fdat<-data.frame(Descriptor = x[, x.meta.ind])
	for(i in colnames(fdat)){
		fdat[, i]<-gsub(" ", "", fdat[,i])
	}
	dat<-x[, -x.meta.ind]
	rownames(dat)<-paste("R", 1:nrow(dat), sep = "")
	pdat<-data.frame(sample_id = colnames(dat))
	eset<-to.eSet(mat = dat, pdat = pdat, fdat = fdat)

	#fData(eset)$Unique.Name<-gsub("-CNvalues", "", fData(eset)$Unique.Name)
	#dirs<-grep("Amp", fData(eset)$Unique.Name)
	if(!is.na(broad_lesions_threshold))
	broad_thres<-read.csv(broad_lesions_threshold, sep = "\t", header = T)
	
	eset.amp<-eset
	eset.del<-eset

	if(!is.na(broad_lesions_threshold)){
		broad_thres<-read.csv(broad_lesions_threshold, sep = "\t", header = T)
		eset.amp.ind<-which(fData(eset.amp)$Descriptor %in% broad_thres$Arm[broad_thres$Amp.q.value < amp_qvalue])
		eset.amp<-eset.amp[eset.amp.ind,]
		eset.del.ind<-which(fData(eset.del)$Descriptor %in% broad_thres$Arm[broad_thres$Del.q.value < del_qvalue])
		eset.del<-eset.del[eset.del.ind,]
	}

	if(nrow(eset.amp)>0)
	fData(eset.amp)$Unique.Name<-paste("AmplificationArm", 1:nrow(eset.amp), sep = "")
	else fData(eset.amp)$Unique.Name<-character(0)
	if(nrow(eset.del)>0)
	fData(eset.del)$Unique.Name<-paste("DeletionArm", 1:nrow(eset.del), sep = "")
	else fData(eset.del)$Unique.Name<-character(0)

	mat<-rbind(exprs(eset.amp), exprs(eset.del))
	fdat<-rbind(fData(eset.amp), fData(eset.del))
	pdat<-pData(eset.amp)
	rownames(mat)<-paste("R", 1:nrow(mat), sep = "")
	eset<-to.eSet(mat = mat, pdat = pdat, fdat = fdat)

	ncols<-ncol(eset)
	mat<-t(sapply(1:nrow(eset), function(i){
		is_amp<-grepl("Amp", fData(eset)$Unique.Name[i])
		x<-as.numeric(exprs(eset)[i,, drop = FALSE])
		res<-rep(0, ncols)

		if(is_amp)
		res[x>amp_thres]<-1
		else
		res[x<del_thres]<-1
		return(res)
	}))
	exprs(eset)<-mat

	return(eset)
}

get_chr<-function(x){
	arm<-"p"
	x.split<-strsplit(x, split = arm)[[1]]

	if (length(x.split) == 1){
		arm<-"q"
		x.split<-strsplit(x, split = arm)[[1]]
	}
	return(paste(x.split[1], arm, sep = ""))
}

focal_or_arm<-function(eset.focal, eset.arm){
	eset.or<-eset.focal
	vec<-as.character(fData(eset.or)$Unique.Name)
	vecstr<-gsub("+[0-9]", "", vec)
	vecnum<-gsub("+[a-z]|+[A-Z]", "", vec)

	for(ind.focal in 1:nrow(eset.focal)){
		name.focal<-fData(eset.focal)$Unique.Name[ind.focal]
		desc.focal<-fData(eset.focal)$Descriptor[ind.focal]
		desc.arm<-get_chr(desc.focal)
		if(grepl("Amplification", fData(eset.focal)$Unique.Name[ind.focal]))
			ind.arm<-which(fData(eset.arm)$Descriptor == desc.arm & grepl("Amplification", fData(eset.arm)$Unique.Name))	
		else
			ind.arm<-which(fData(eset.arm)$Descriptor == desc.arm & grepl("Deletion", fData(eset.arm)$Unique.Name))
		if(length(ind.arm)!=0){
			##has matching arm level copy number alteration
			x.both<-rbind(exprs(eset.focal[ind.focal,]), exprs(eset.arm[ind.arm,]))
			x.or<-apply(x.both, 2, max) #take max of copy number status
			exprs(eset.or)[ind.focal,]<-x.or
			fData(eset.or)$Unique.Name[ind.focal]<-paste(vecstr[ind.focal],
				"orarm", vecnum[ind.focal],sep = "")
		}
	}
	return(eset.or)
}

get_cis_genes<-function(f.genes,direction){
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

	#remove characters after "|": alternative ensembl name
	res.genes<-lapply(res.genes, function(x) gsub("\\|.*", "",  x))
	res.genes<-lapply(res.genes, function(x) gsub(" ", "",  x))
	res.genes<-lapply(res.genes, function(x) toupper(x))

	#res.genes<-lapply(res.genes, function(x) gsub("\\[|\\]", "", x))
	return(list(meta = res.meta, genes = res.genes))
}


#get_genes_in_arm<-function(x,y){
#	return(y$hgnc_symbol[which(y$arm == x)])
#}

#' @import org.Hs.eg.db AnnotationDbi
#' @export
get_genes_in_arm<-function(x){
	#library("org.Hs.eg.db")
	db<-as.list( revmap(org.Hs.eg.db::org.Hs.egMAP) )

	res<-lapply(x, function(i){	
		inds<-grep(paste("^", i, sep = ""), names(db))
		ids<-as.character(unlist(db[inds]))
		symbols<-select(org.Hs.eg.db,
	       keys = ids,
	       columns=c("ENSEMBL","ENTREZID","SYMBOL","GENENAME"),
	       keytype="ENTREZID")
		return(unique(symbols[, "SYMBOL"]))
		})

	return(res)
}

add_direction<-function(cn, remove.cols = TRUE){
	cnid<-"Unique.Name"
	#make a vector of directions of alteration and add to fData(cn)
	cndirvector<-as.character(sapply(fData(cn)[, cnid], function(x){
		if (grepl("Amp", x)) return("Amplification")
		else return("Deletion")
		}))
	fData(cn)<-data.frame(fData(cn), alteration_direction = cndirvector)
	#display only three columns from fData(cn) in limma tables
	if(remove.cols){
		fData(cn)<-fData(cn)[, c("Unique.Name", "Descriptor", "alteration_direction")]
	}
	return(cn)
}



#wrapper for make_GISTIC2
#' @import org.Hs.eg.db
#' @export
make_GISTIC2_with_thres_focal_and_arm<-function(gistic_in, all_genes, 
	amp_thres_arm, del_thres_arm, 
	amp_qvalue_arm, del_qvalue_arm, 
	qvalue_focal,
	all_lesions = NA, f.amp = NA, f.del = NA, 
	broad_thres = TRUE){

	cat("Reading GISTIC2 data..\n")

	if(is.na(all_lesions))
	all_lesions<-paste(gistic_in, "/", "all_lesions.conf_99.txt", sep = "")

	cat("Reading focal alterations...\n")
	eset.focal<-read_GISTIC2_focal(all_lesions = all_lesions)
	eset.focal<-eset.focal[as.numeric(fData(eset.focal)$q.values)<qvalue_focal,]

	cat("Reading arm alterations...\n")
	broad_lesions<-paste(gistic_in, "/", "broad_values_by_arm.txt", sep = "")
	if(broad_thres){
		broad_lesions_threshold<-paste(gistic_in, "/", "broad_significance_results.txt", sep = "")
	} else {
		broad_lesions_threshold<-NA
	}	
	eset.arm<-read_GISTIC2_broad_with_thres(broad_lesions = broad_lesions, 
	broad_lesions_threshold = broad_lesions_threshold,
	amp_qvalue = amp_qvalue_arm, del_qvalue = del_qvalue_arm, 
	amp_thres = amp_thres_arm, del_thres = del_thres_arm)

	eset.or<-focal_or_arm(eset.focal, eset.arm)

	##get cis genes
	cat("Fetching cis genes in focal..\n")

	if(is.na(f.amp))
	f.amp<-paste(gistic_in, "/", "amp_genes.conf_99.txt", sep = "")

	if(is.na(f.del))
	f.del<-paste(gistic_in, "/", "del_genes.conf_99.txt", sep = "")

	genes.amp<-get_cis_genes(f.genes = f.amp, direction = "Amplification")
	genes.del<-get_cis_genes(f.genes = f.del, direction = "Deletion")
	
	genes.both<-list(meta = rbind(genes.amp[["meta"]], genes.del[["meta"]]),
		genes = c(genes.amp[["genes"]], genes.del[["genes"]]))

	sample_meta<-fData(eset.focal)
	sample_meta$wide.peak.boundaries<-as.character(sapply(sample_meta$Wide.Peak.Limits, 
		function(x) strsplit(x, split = "\\(")[[1]][1]))
	gene_meta<-genes.both[['meta']]
	idx<-match(sample_meta$wide.peak.boundaries, gene_meta$wide.peak.boundaries)
	cisgenes<-genes.both[['genes']][idx]

	cisgenes.arm<-get_genes_in_arm(fData(eset.arm)[, "Descriptor"])

	eset.focal<-add_direction(eset.focal)
	eset.arm<-add_direction(eset.arm)
	eset.or<-add_direction(eset.or)

	res<-list(focal = list(cn=eset.focal, cis=cisgenes), 
		arm = list(cn=eset.arm, cis=cisgenes.arm),
		or = list(cn=eset.or, cis=cisgenes))
	return(res)
}

