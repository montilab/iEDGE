require(Biobase)
require(biomaRt)

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
#	eset.amp.ind<-which(fData(eset.amp)$Descriptor %in% broad_thres$Arm[broad_thres$Amp.q.value < amp_qvalue])
#	eset.amp<-eset.amp[eset.amp.ind,]
	
	eset.del<-eset
	#eset.del.ind<-which(fData(eset.del)$Descriptor %in% broad_thres$Arm[broad_thres$Del.q.value < del_qvalue])
	#eset.del<-eset.del[eset.del.ind,]
	
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


get_genes_in_arm<-function(x,y){
	return(y$hgnc_symbol[which(y$arm == x)])
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
make_GISTIC2<-function(f.dir.in, binarize = TRUE, my.symbols){

	cat("Reading GISTIC2 data..\n")
	all_lesions<-paste(f.dir.in, "/", "all_lesions.conf_99.txt", sep = "")
	broad_lesions<-paste(f.dir.in, "/", "broad_values_by_arm.txt", sep = "")
	eset.focal<-read_GISTIC2_focal(all_lesions = all_lesions)
	eset.arm<-read_GISTIC2_broad(broad_lesions = broad_lesions)

	fdat<-fData(eset.focal)
	amp.thres<-fdat$Amplitude.Threshold[grepl("Amplification", fdat$Unique.Name)][1]
	amp.thres<-sapply(strsplit(amp.thres, split = ";")[[1]], 
		function(x) strsplit(x, split = ":")[[1]][2])[2]
	amp.breaks<-c(-100, as.numeric(unlist(strsplit(gsub("<t<", ",", amp.thres), ","))), 100)
	eset.arm.amp<-eset.arm
	arm.amp.mat<-t(apply(exprs(eset.arm.amp), 1, 
		function(x) 
			as.numeric(as.character(cut(x = x, breaks = amp.breaks, labels = c(0, 1, 2))))))
	exprs(eset.arm.amp)<-arm.amp.mat

	del.thres<-fdat$Amplitude.Threshold[grepl("Deletion", fdat$Unique.Name)][1]
	del.thres<-sapply(strsplit(del.thres, split = ";")[[1]], 
		function(x) strsplit(x, split = ":")[[1]][2])[2]
	del.breaks<-rev(c(100, as.numeric(unlist(strsplit(gsub(">t>", ",", del.thres), ","))), -100))
	eset.arm.del<-eset.arm
	arm.del.mat<-t(apply(exprs(eset.arm.del), 1, 
		function(x) 
			as.numeric(as.character(cut(x = x, breaks = del.breaks, labels = c(2, 1, 0))))))
	exprs(eset.arm.del)<-arm.del.mat

	amp.arm.uniquename<-paste("AmplificationArm", 1:nrow( fData(eset.arm.amp)), sep = "")
	del.arm.uniquename<-paste("DeletionArm", 1:nrow( fData(eset.arm.del)), sep = "")
	
	fdat.arm<-rbind(cbind(Unique.Name = amp.arm.uniquename, fData(eset.arm.amp)),
		cbind(Unique.Name = del.arm.uniquename, fData(eset.arm.del)))
	exprs.arm<-rbind(exprs(eset.arm.amp), exprs(eset.arm.del))
	rownames(exprs.arm)<-rownames(fdat.arm)
	colnames(exprs.arm)<-pData(eset.arm)[,1]
	eset.arm<-to.eSet(mat = exprs.arm, pdat = pData(eset.arm), fdat = fdat.arm)

	eset.or<-focal_or_arm(eset.focal, eset.arm)

	#binary phenotype
	if(binarize == TRUE){
		cat("Binarizing alterations..\n")
		mat<-exprs(eset.focal)
		mat[mat == 2]<-1
		exprs(eset.focal)<-mat
		mat<-exprs(eset.arm)
		mat[mat == 2]<-1
		exprs(eset.arm)<-mat
		mat<-exprs(eset.or)
		mat[mat == 2]<-1
		exprs(eset.or)<-mat
	}

	##get cis genes
	cat("Fetching cis genes in focal..\n")
	f.amp<-paste(f.dir.in, "/", "amp_genes.conf_99.txt", sep = "")
	f.del<-paste(f.dir.in, "/", "del_genes.conf_99.txt", sep = "")
	genes.amp<-get_cis_genes(f.genes = f.amp, direction = "Amplification")
	genes.del<-get_cis_genes(f.genes = f.del, direction = "Deletion")
	
	genes.both<-list(meta = rbind(genes.amp[["meta"]], genes.del[["meta"]]),
		genes = c(genes.amp[["genes"]], genes.del[["genes"]]))

	#cisgenes<-c(genes.amp[["genes"]], genes.del[["genes"]])

	sample_meta<-fData(eset.focal)
	sample_meta$wide.peak.boundaries<-as.character(sapply(sample_meta$Wide.Peak.Limits, 
		function(x) strsplit(x, split = "\\(")[[1]][1]))
	gene_meta<-genes.both[['meta']]
	idx<-match(sample_meta$wide.peak.boundaries, gene_meta$wide.peak.boundaries)
	#cisgenes[['meta']]<-cisgenes[['meta']][idx,]
	cisgenes<-genes.both[['genes']][idx]


	##get cis genes in arm
	# Filter on HGNC symbol and chromosome, retrieve genomic location and band
	cat("Fetching cis genes in arm..\n")
		#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	ensembl<- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

	# Only use standard human chromosomes
	normal.chroms <- c(1:22, "X", "Y", "M")
	my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                    filters = c("hgnc_symbol", "chromosome_name"),
                    values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),
                                  mart = ensembl)
	my.regions$arm<-paste(my.regions$chromosome_name,  
		sapply(my.regions$band, function(x) substr(x, 1, 1)), sep ="")


	cisgenes.arm<-lapply(1:nrow(fData(eset.arm)), 
		function(i){
			get_genes_in_arm(x=fData(eset.arm)$Descriptor[i], y = my.regions)
		})
	##get cis genes in arm and in focal
	cat("Fetching cis genes in arm and in focal..\n")
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
	cisgenes.arm.infocal<-cisgenes[ind.focal]

	eset.focal<-add_direction(eset.focal)
	eset.arm<-add_direction(eset.arm)
	eset.or<-add_direction(eset.or)

	res<-list(focal = list(cn=eset.focal, cis=cisgenes), 
		arm = list(cn=eset.arm, cis=cisgenes.arm),
		or = list(cn=eset.or, cis=cisgenes.arm.infocal))
	return(res)
}



#wrapper for make_GISTIC2
make_GISTIC2_with_thres_arm<-function(gistic_in, all_genes, 
	amp_thres_arm, del_thres_arm, 
	amp_qvalue_arm, del_qvalue_arm, all_lesions = NA, f.amp = NA, f.del = NA, 
	broad_thres = TRUE){

	cat("Reading GISTIC2 data..\n")

	if(is.na(all_lesions))
	all_lesions<-paste(gistic_in, "/", "all_lesions.conf_99.txt", sep = "")

	cat("Reading focal alterations...\n")
	eset.focal<-read_GISTIC2_focal(all_lesions = all_lesions)

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

	#cisgenes<-c(genes.amp[["genes"]], genes.del[["genes"]])

	sample_meta<-fData(eset.focal)
	sample_meta$wide.peak.boundaries<-as.character(sapply(sample_meta$Wide.Peak.Limits, 
		function(x) strsplit(x, split = "\\(")[[1]][1]))
	gene_meta<-genes.both[['meta']]
	idx<-match(sample_meta$wide.peak.boundaries, gene_meta$wide.peak.boundaries)
	cisgenes<-genes.both[['genes']][idx]

	##get cis genes in arm
	# Filter on HGNC symbol and chromosome, retrieve genomic location and band
	cat("Fetching cis genes in arm..\n")
		#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	ensembl<- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

	# Only use standard human chromosomes
	normal.chroms <- c(1:22, "X", "Y", "M")
	my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                    filters = c("hgnc_symbol", "chromosome_name"),
                    values = list(hgnc_symbol=all_genes, chromosome_name=normal.chroms),
                                  mart = ensembl)
	my.regions$arm<-paste(my.regions$chromosome_name,  
		sapply(my.regions$band, function(x) substr(x, 1, 1)), sep ="")


	cisgenes.arm<-lapply(1:nrow(fData(eset.arm)), 
		function(i){
			get_genes_in_arm(x=fData(eset.arm)$Descriptor[i], y = my.regions)
		})
	
	eset.focal<-add_direction(eset.focal)
	eset.arm<-add_direction(eset.arm)
	eset.or<-add_direction(eset.or)

	res<-list(focal = list(cn=eset.focal, cis=cisgenes), 
		arm = list(cn=eset.arm, cis=cisgenes.arm),
		or = list(cn=eset.or, cis=cisgenes))
	return(res)
}


#wrapper for make_GISTIC2
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

	#cisgenes<-c(genes.amp[["genes"]], genes.del[["genes"]])

	sample_meta<-fData(eset.focal)
	sample_meta$wide.peak.boundaries<-as.character(sapply(sample_meta$Wide.Peak.Limits, 
		function(x) strsplit(x, split = "\\(")[[1]][1]))
	gene_meta<-genes.both[['meta']]
	idx<-match(sample_meta$wide.peak.boundaries, gene_meta$wide.peak.boundaries)
	cisgenes<-genes.both[['genes']][idx]

	##get cis genes in arm
	# Filter on HGNC symbol and chromosome, retrieve genomic location and band
	cat("Fetching cis genes in arm..\n")
		#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	ensembl<- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

	# Only use standard human chromosomes
	normal.chroms <- c(1:22, "X", "Y", "M")
	my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                    filters = c("hgnc_symbol", "chromosome_name"),
                    values = list(hgnc_symbol=all_genes, chromosome_name=normal.chroms),
                                  mart = ensembl)
	my.regions$arm<-paste(my.regions$chromosome_name,  
		sapply(my.regions$band, function(x) substr(x, 1, 1)), sep ="")


	cisgenes.arm<-lapply(1:nrow(fData(eset.arm)), 
		function(i){
			get_genes_in_arm(x=fData(eset.arm)$Descriptor[i], y = my.regions)
		})
	
	eset.focal<-add_direction(eset.focal)
	eset.arm<-add_direction(eset.arm)
	eset.or<-add_direction(eset.or)

	res<-list(focal = list(cn=eset.focal, cis=cisgenes), 
		arm = list(cn=eset.arm, cis=cisgenes.arm),
		or = list(cn=eset.or, cis=cisgenes))
	return(res)
}





#wrapper for make_GISTIC2
make_GISTIC2_with_thres<-function(gistic_in, all_genes, 
	amp_thres_focal, amp_thres_arm, del_thres_focal, del_thres_arm, 
	amp_qvalue_arm, del_qvalue_arm, all_lesions = NA, f.amp = NA, f.del = NA){

	cat("Reading GISTIC2 data..\n")

	if(is.na(all_lesions))
	all_lesions<-paste(gistic_in, "/", "all_lesions.conf_99.txt", sep = "")

	eset.focal<-read_GISTIC2_focal_with_thres(all_lesions, 
		amp_thres = amp_thres_focal, del_thres = del_thres_focal)

	broad_lesions<-paste(gistic_in, "/", "broad_values_by_arm.txt", sep = "")
	broad_lesions_threshold<-paste(gistic_in, "/", "broad_significance_results.txt", sep = "")

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

	#cisgenes<-c(genes.amp[["genes"]], genes.del[["genes"]])

	sample_meta<-fData(eset.focal)
	sample_meta$wide.peak.boundaries<-as.character(sapply(sample_meta$Wide.Peak.Limits, 
		function(x) strsplit(x, split = "\\(")[[1]][1]))
	gene_meta<-genes.both[['meta']]
	idx<-match(sample_meta$wide.peak.boundaries, gene_meta$wide.peak.boundaries)
	#cisgenes[['meta']]<-cisgenes[['meta']][idx,]
	cisgenes<-genes.both[['genes']][idx]

	##get cis genes in arm
	# Filter on HGNC symbol and chromosome, retrieve genomic location and band
	cat("Fetching cis genes in arm..\n")
		#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	ensembl<- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

	# Only use standard human chromosomes
	normal.chroms <- c(1:22, "X", "Y", "M")
	my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                    filters = c("hgnc_symbol", "chromosome_name"),
                    values = list(hgnc_symbol=all_genes, chromosome_name=normal.chroms),
                                  mart = ensembl)
	my.regions$arm<-paste(my.regions$chromosome_name,  
		sapply(my.regions$band, function(x) substr(x, 1, 1)), sep ="")


	cisgenes.arm<-lapply(1:nrow(fData(eset.arm)), 
		function(i){
			get_genes_in_arm(x=fData(eset.arm)$Descriptor[i], y = my.regions)
		})
	
	eset.focal<-add_direction(eset.focal)
	eset.arm<-add_direction(eset.arm)
	eset.or<-add_direction(eset.or)

	res<-list(focal = list(cn=eset.focal, cis=cisgenes), 
		arm = list(cn=eset.arm, cis=cisgenes.arm),
		or = list(cn=eset.or, cis=cisgenes))
	return(res)
}



#wrapper for make_GISTIC2
make_GISTIC2_with_thres_nocisgenes<-function(gistic_in, all_genes, 
	amp_thres_focal, amp_thres_arm, del_thres_focal, del_thres_arm, 
	amp_qvalue_arm, del_qvalue_arm, all_lesions = NA, f.amp = NA, f.del = NA){

	cat("Reading GISTIC2 data..\n")

	if(is.na(all_lesions))
	all_lesions<-paste(gistic_in, "/", "all_lesions.conf_99.txt", sep = "")

	eset.focal<-read_GISTIC2_focal_with_thres(all_lesions, 
		amp_thres = amp_thres_focal, del_thres = del_thres_focal)

	broad_lesions<-paste(gistic_in, "/", "broad_values_by_arm.txt", sep = "")
	broad_lesions_threshold<-paste(gistic_in, "/", "broad_significance_results.txt", sep = "")

	eset.arm<-read_GISTIC2_broad_with_thres(broad_lesions = broad_lesions, 
	broad_lesions_threshold = broad_lesions_threshold,
	amp_qvalue = amp_qvalue_arm, del_qvalue = del_qvalue_arm, 
	amp_thres = amp_thres_arm, del_thres = del_thres_arm)

	eset.or<-focal_or_arm(eset.focal, eset.arm)

	eset.focal<-add_direction(eset.focal)
	eset.arm<-add_direction(eset.arm)
	eset.or<-add_direction(eset.or)

	res<-list(focal = list(cn=eset.focal), 
		arm = list(cn=eset.arm),
		or = list(cn=eset.or))
	return(res)
}



###get GISTIC2 data from TCGA brca dataset

make_TCGA_brca_gistic<-function(f.header, f.dir.in, f.dir.out) {

	all_lesions<-paste(f.dir.in, "/", "all_lesions.conf_99.txt", sep = "")
	broad_lesions<-paste(f.dir.in, "/", "broad_values_by_arm.txt", sep = "")
	eset.focal<-read_GISTIC2_focal(all_lesions = all_lesions)
	eset.arm<-read_GISTIC2_broad(broad_lesions = broad_lesions)

	fdat<-fData(eset.focal)
	amp.thres<-fdat$Amplitude.Threshold[grepl("Amplification", fdat$Unique.Name)][1]
	amp.thres<-sapply(strsplit(amp.thres, split = ";")[[1]], 
		function(x) strsplit(x, split = ":")[[1]][2])[2]
	amp.breaks<-c(-100, as.numeric(unlist(strsplit(gsub("<t<", ",", amp.thres), ","))), 100)
	eset.arm.amp<-eset.arm
	arm.amp.mat<-t(apply(exprs(eset.arm.amp), 1, 
		function(x) 
			as.numeric(as.character(cut(x = x, breaks = amp.breaks, labels = c(0, 1, 2))))))
	exprs(eset.arm.amp)<-arm.amp.mat

	del.thres<-fdat$Amplitude.Threshold[grepl("Deletion", fdat$Unique.Name)][1]
	del.thres<-sapply(strsplit(del.thres, split = ";")[[1]], 
		function(x) strsplit(x, split = ":")[[1]][2])[2]
	del.breaks<-rev(c(100, as.numeric(unlist(strsplit(gsub(">t>", ",", del.thres), ","))), -100))
	eset.arm.del<-eset.arm
	arm.del.mat<-t(apply(exprs(eset.arm.del), 1, 
		function(x) 
			as.numeric(as.character(cut(x = x, breaks = del.breaks, labels = c(2, 1, 0))))))
	exprs(eset.arm.del)<-arm.del.mat

	amp.arm.uniquename<-paste("AmplificationPeakArm", 1:nrow( fData(eset.arm.amp)), sep = "")
	del.arm.uniquename<-paste("DeletionPeakArm", 1:nrow( fData(eset.arm.del)), sep = "")
	fdat.arm<-rbind(cbind(Unique.Name = amp.arm.uniquename, fData(eset.arm.amp)),
		cbind(Unique.Name = del.arm.uniquename, fData(eset.arm.del)))
	exprs.arm<-rbind(exprs(eset.arm.amp), exprs(eset.arm.del))
	rownames(exprs.arm)<-rownames(fdat.arm)
	colnames(exprs.arm)<-pData(eset.arm)[,1]
	eset.arm<-to.eSet(mat = exprs.arm, pdat = pData(eset.arm), fdat = fdat.arm)

	eset.or<-focal_or_arm(eset.focal, eset.arm)
	#save gistic files as R objects (RDS)
	saveRDS(eset.focal, file = paste(f.dir.out, "/", f.header, "_gistic_all_focal.RDS", sep = ""))
	saveRDS(eset.arm, file = paste(f.dir.out, "/", f.header,  "_gistic_all_arm.RDS", sep = ""))
	saveRDS(eset.or, file = paste(f.dir.out, "/", f.header, "_gistic_all_or.RDS", sep = ""))

	#binary phenotype

	mat<-exprs(eset.focal)
	mat[mat == 2]<-1
	exprs(eset.focal)<-mat

	mat<-exprs(eset.arm)
	mat[mat == 2]<-1
	exprs(eset.arm)<-mat

	mat<-exprs(eset.or)
	mat[mat == 2]<-1
	exprs(eset.or)<-mat

	#save gistic files as R objects (RDS)
	saveRDS(eset.focal, file = paste(f.dir.out, "/", f.header, "_gistic_all_focal_binary.RDS", sep = ""))
	saveRDS(eset.arm, file = paste(f.dir.out, "/", f.header, "_gistic_all_arm_binary.RDS", sep = ""))
	saveRDS(eset.or, file = paste(f.dir.out, "/", f.header, "_gistic_all_or_binary.RDS", sep = ""))

	##get cis genes
	eset.focal<-readRDS(paste(f.dir.out, "/", f.header, "_gistic_all_focal_binary.RDS", sep = ""))

	f.amp<-paste(f.dir.in, "/", "amp_genes.conf_99.txt", sep = "")
	f.del<-paste(f.dir.in, "/", "del_genes.conf_99.txt", sep = "")
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
	saveRDS(genes.both, file = paste(f.dir.out, "/", f.header, "_gistic_cisgenes.RDS", sep = ""))
}