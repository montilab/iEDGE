to_eset<-function(assay.data, pdat, fdat){
	require(Biobase)
	exprs<-as.matrix(assay.data)
	pMetadata<-data.frame(labelDescription = colnames(pdat), row.names= colnames(pdat))
	phenoData<-new("AnnotatedDataFrame", data = pdat, varMetadata = pMetadata)
	fMetaData<-data.frame(labelDescription = colnames(fdat), row.names = colnames(fdat))
	featureData<-new("AnnotatedDataFrame", data= fdat, varMetadata=fMetaData)  
	annotation <-""
	eSet<-ExpressionSet(assayData=exprs, phenoData = phenoData,  annotation = annotation, 
		featureData = featureData)
	return(eSet)
}



#' \code{cls2df} reads from .cls file into data frame
#' @param f .cls file name
cls2df<-function(f){
	con <- file(f) 
	open(con)
	att.name <- list()
	att.value <- list()

	first.line<-readLines(con, n = 1, warn = FALSE) #skip
	current.line <- 1
	while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
 		if (grepl("^#", line)){
 			curr.name<-strsplit(line, split = " ")[[1]]
 			curr.name<-paste(curr.name[2:length(curr.name)], collapse= "_vs_")
 			att.name[[current.line]]<-curr.name
 		} else{
 			curr.value<-strsplit(line, split = " ")[[1]]
 			att.value[[current.line]]<-curr.value
 			current.line <- current.line + 1
 		}
	} 
	close(con)

	df<-do.call(cbind, att.value)
	df<-data.frame(df)
	names(df)<-unlist(att.name)
	return(df)
}





#' \code{read_cis_genes} reads cis genes from GISTIC files, helper function for read_gistic
#' @param file.in.peak is the file containing in_peak genes
#' @param file.in.region is the file containing in_region genes
#' @param cn.dir must be one of "amp" - Amplification or "del" - Deletion
read_cis_genes<-function(file.in.peak, #file cis genes in peak
	file.in.region, #file containing cis genes in region
	cn.dir #must be one of "amp" or "del"
	){
	#file.in<-"~/Desktop/monti_lab/year2/DataIntegration/ricover/data/ricover/datasets/gistic2_analysis_allelic_capseg/job.50618380/del_genes.conf_99.sorted_in_region.txt"
	res.peak<-read.table(file.in.peak, sep = "\t", header = FALSE, fill = T)
	res.peak<- res.peak[, apply(res.peak, 2, function(x) !all(is.na(x)))]
	last.meta.ind<-which(res.peak[,1] == "genes in wide peak")-1
	
	res.meta.peak<-t(res.peak[1:last.meta.ind,])
	colnames(res.meta.peak)<-res.meta.peak[1,]
	res.meta.peak<-data.frame(res.meta.peak[-1,])
	res.genes.peak<-res.peak[-(1:last.meta.ind), -1]
	res.genes.peak<-sapply(1:ncol(res.genes.peak),
		function(x) setdiff(unique(as.character(res.genes.peak[, x])), ""))

	#Unique.Name column is found GISTIC2 file formats only, add to GISTIC1 file if needed
	if(cn.dir == "amp"){
		dir.header <-"AmplificationPeak"
	} else if (cn.dir == "del"){
		dir.header <- "DeletionPeak"
	} else{
		stop("in read_cis_genes, cn.dir must be one of \"amp\" or \"del\"")
	}
	if(!("Unique.Name" %in% colnames(res.meta.peak))){
		res.meta.peak<-cbind(Unique.Name=paste(dir.header, 1:nrow(res.meta.peak), sep = ""), 
			res.meta.peak)
	}

	if (is.na(file.in.region)){ #peak files provided only
		warning("only peak file is provided, region file missing")
		res.list<-list(genesets.peak = res.genes.peak, 
			meta.peak = res.meta.peak)
		return(res.list)
	}
	else {
		res.region<-read.table(file.in.region, sep = "\t", header = FALSE, fill = T)
		last.meta.ind<-which(res.region[,1] == "genes in wide peak")-1
		res.meta.region<-t(res.region[1:last.meta.ind,])
		colnames(res.meta.region)<-res.meta.region[1,]
		res.meta.region<-data.frame(res.meta.region[-1,])

		res.genes.region<-res.region[-(1:last.meta.ind), -1]
		res.genes.region<-sapply(1:ncol(res.genes.region),
			function(x) setdiff(unique(as.character(res.genes.region[, x])), ""))

	#	if (!all(as.character(res.meta.peak[, "cytoband"]) == as.character(res.meta.region[, "cytoband"]))){
	#		stop("cytobands in peak and region are inconsistent!")
	#	}

		if(!("Unique.Name" %in% colnames(res.meta.region))){
			res.meta.region<-cbind(Unique.Name=paste(dir.header, 1:nrow(res.meta.region), sep = ""), 
				res.meta.region)
		}
		res.list<-list(genesets.peak = res.genes.peak, 
			genesets.region = res.genes.region, 
			meta.peak = res.meta.peak, 
			meta.region = res.meta.region)
		return(res.list)
	}
}


#' \code{read_gistic} reads from GISTIC files
#' @param lesions.dir is the file containing all_lesions
#' @param genes.peak.dir is the file containing in_peak_conf## genes
#' @param genes.region.dir is the file containing in_region_conf## genes
#' @param match.by c(first, second), first = unique id column name in lesions.dir (default = Descriptor), second = unique id column name in genes.peak.dir (default = cytoband)
#' @param cn.dir must be one of "amp" - Amplification or "del" - Deletion
#' @export
read_gistic<-function(lesions.dir, genes.peak.dir, genes.region.dir,
	#binary = T, 
	match.by = c("Descriptor", "cytoband"), 
	#match.by[1] = unique id in lesions, 
	#match.by[2] = unique id in genes$meta.peak
	cn.dir #amp or del
	){

	lesions.all<-read.table(lesions.dir,sep = "\t", header =T)
	lesions.all<-lesions.all[, apply(lesions.all, 2, function(x) !all(is.na(x)))]
	#if (binary == TRUE)
	lesions<-lesions.all[!grepl("Actual",lesions.all$Amplitude.Threshold), ]
	#else
	lesions.cont<-lesions.all[grepl("Actual",lesions.all$Amplitude.Threshold), ]
	
	lesions$Descriptor<-gsub(" ","", lesions$Descriptor)
	lesions$Unique.Name<-gsub(" ", "", lesions$Unique.Name)

	m<-nrow(lesions)
	n<-ncol(lesions)
	last_desc_ind<-which(colnames(lesions) == "Amplitude.Threshold")

	res.genes<-read_cis_genes(file.in.peak = genes.peak.dir, 
		file.in.region = genes.region.dir, 
		cn.dir = cn.dir)

	lesions.rows<-match(res.genes$meta.peak[, match.by[2]],lesions[, match.by[1]])
	lesions.mat<-lesions[lesions.rows, (last_desc_ind + 1):n]
	lesions.mat.cont<-lesions.cont[lesions.rows, (last_desc_ind + 1):n]

	lesions.meta<-lesions[lesions.rows, 1:last_desc_ind]

	return(list(mat = lesions.mat, mat.cont = lesions.mat.cont, 
		meta = lesions.meta, cisgenes = res.genes))
}



#' \code{read_res} converts .res data into internal gene expression data format
#' @param f gene expression .res file
#' @param row.ind whether first row has numeric row indices, default = FALSE
read_res<-function(f, row.ind = FALSE){
	x<-read.table(f, sep = "\t", header = F, skip = 3, quote = "")
	con<-file(f)
	x.header<-strsplit( readLines(con, n=1), split = "\t")[[1]][-(1:2)]
	close(con)

	##next line added if only files has first column as numeric row indices
	if (row.ind == TRUE) x <-x[,-1]
	description<-x[,1]
	accession<-toupper(x[,2])
	x<-x[, -(1:2)]
	odd_col<-sapply(1:dim(x)[2], function(x) x %% 2 == 1)
	x<-x[, odd_col]
	rownames(x)<-accession
	colnames(x)<-x.header[odd_col]
	return(list(mat = x, description = description, accession = accession))
}

read_eset<-function(f){
	require(Biobase)
	x<-readRDS(f)
	return(list(mat = exprs(x), description = fData(x)$description, accession = fData(x)$accession))

}

#note: ge.confounder is only included in annotation => expression values are uncorrected
#confounder correction can be done in the differential expression step

#' \code{match_samples} match samples found in ge and gistic objects
#' @param ge ge object
#' @param gistic gistic object
#' @param fmap.file file containing mapping from ge to gistic samples
#' @param fmap.gistic.col column name in fmap.file containing gistic ids
#' @param fmap.ge.col column name in fmap.file containing ge ids
#' @param ge.confounder sample confounder indicator file name, default = NA, no confounder

match_samples<-function(gistic, #gistic object 
	ge, #ge object
	fmap.file, #fmap file name 
	fmap.gistic.col = "snp.array.ID", 
	fmap.ge.col = "mrna.array.ID", 
	ge.confounder.file = NA){

	x<-colnames(gistic$mat)
	y<-colnames(ge$mat)
	
	fmap<-read.table(fmap.file, header = TRUE, sep = "\t")

	fmap<-fmap[fmap[, fmap.gistic.col] %in% x,]
	fmap<-fmap[fmap[, fmap.ge.col] %in% y,]
	
	x.match<-match(fmap[, fmap.gistic.col], x)

	#hotfix for names with "-" in sample map file
	#x.match<-match(gsub("-", ".", fmap[, fmap.gistic.col]), x)

	if (any(is.na(x.match)) | any(duplicated(x.match)))
		stop("duplicated or missing sample mapping")

	gistic$mat<-gistic$mat[, x.match]
	gistic$mat.cont<-gistic$mat.cont[, x.match]
	y.match<-match(fmap[, fmap.ge.col], y)

	if (any(is.na(y.match)) | any(duplicated(y.match)))
		stop("duplicated or missing sample mapping")
	ge$mat<-ge$mat[,y.match]

	if (!is.na(ge.confounder.file)){
		ge.confounder<-cls2df(ge.confounder.file)
		ge$confounder <-ge.confounder[y.match, , drop = FALSE]
		rownames(ge$confounder)<-1:length(y.match)
	}
	return(list(gistic = gistic, ge = ge))
}



#note: ge.confounder is only included in annotation => expression values are uncorrected
#confounder correction can be done in the differential expression step

#' \code{process_data} process gistic and ge data from text files
#' @param amp.peak.file GISTIC amp genes in peak file
#' @param amp.region.file GISTIC amp genes in region file
#' @param del.peak.file GISTIC del genes in peak file
#' @param del.region.file GISTIC del genes in region file
#' @param lesions.file GISTIC all_lesions file
#' @param ge.file Gene Expression res file
#' @param fmap.file file mapping tabular file from GISTIC to GE sample names
#' @param ge.confounder.file confounder indicator cls file, must match order in fmap.file
#' @param ... see addition arguments passed to \code{match_samples, read_res, read_cis_genes, read_gistic}
#' @return list containing :
#' \item{amp}{amplification gistic and ge}
#' \item{del}{deletion gistic and ge}
#' @export
process_data<-function(
	amp.peak.file,
	amp.region.file = NA,
	del.peak.file,
	del.region.file = NA,
	lesions.file,
	ge.file,
	fmap.file,
	ge.confounder.file = NA,
	match.by = c("Descriptor", "cytoband"),
	fmap.gistic.col = "snp.array.ID",
	fmap.ge.col = "mrna.array.ID", 
	row.ind = FALSE, 
	ge.format = "res" #or eset
	){

	cat("reading gistic: amplifications\n")
	gistic.amp<-read_gistic(lesions.dir = lesions.file, 
		genes.peak.dir = amp.peak.file, 
		genes.region.dir = amp.region.file,
		match.by =  match.by, 
		cn.dir = "amp")
	cat("reading gistic: deletions\n")
	gistic.del<-read_gistic(lesions.dir = lesions.file, 
		genes.peak.dir = del.peak.file, 
		genes.region.dir = del.region.file, 
		match.by =  match.by, 
		cn.dir = "del")
	cat("reading gene expression\n")

	if (ge.format == "res"){
		ge<-read_res(f = ge.file, row.ind = row.ind) #this data is unlogged!
	} else if (ge.format == "eset"){
		ge<-read_eset(f = ge.file)
	} else {
		stop ("invalid ge.format: must be one of res or eset")
	}
	cat("matching ge expression and gistic samples:amplification\n")
	gistic_ge_amp<-match_samples(gistic = gistic.amp, 
		ge = ge, 
		fmap.file = fmap.file, 
		fmap.gistic.col = fmap.gistic.col, 
		fmap.ge.col = fmap.ge.col, 
		ge.confounder = ge.confounder.file)
	cat("matching ge expression and gistic samples:deletion\n")
	gistic_ge_del<-match_samples(gistic = gistic.del, 
		ge = ge, 
		fmap.file = fmap.file, 
		fmap.gistic.col = fmap.gistic.col, 
		fmap.ge.col = fmap.ge.col, 
		ge.confounder = ge.confounder.file)
	cat("data processing done\n")
	return(list(amp =gistic_ge_amp,
		del = gistic_ge_del))
}














