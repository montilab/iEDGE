require(CBMRtools)
require(Biobase)

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


#' \code{read_res} converts .res data into internal gene expression data format
#' @param f gene expression .res file
#' @param row.ind whether first row has numeric row indices, default = FALSE
make_gep<-function(f, row.ind = FALSE, f.confounder = NA){
	x<-read.table(f, sep = "\t", header = F, skip = 3, quote = "")
	con<-file(f)
	x.header<-strsplit( readLines(con, n=1), split = "\t")[[1]][-(1:2)]
	close(con)

	#next line added if only files has first column as numeric row indices
	if (row.ind == TRUE) x <-x[,-1]
	description<-x[,1]
	accession<-toupper(x[,2])
	x<-x[, -(1:2)]
	odd_col<-sapply(1:dim(x)[2], function(x) x %% 2 == 1)
	x<-x[, odd_col]
	rownames(x)<-accession
	colnames(x)<-x.header[odd_col]
	
	fdat<-data.frame(description = description, accession = accession)
	if (!is.na(f.confounder)){
		batch<-cls2df(f = f.confounder)
		pdat<-data.frame(sample_id = colnames(x), batch = as.factor(batch[,1]))
	} else {
		pdat<-data.frame(sample_id = colnames(x))
	}
	#pdat<-cbind(pdat, batch)
	eset<-to.eSet(mat = x, pdat = pdat, fdat = fdat)
	return(eset)
}






lm_correct<-function(y,x){
	res<-lm(y ~ x)
	return(res$residuals + res$coefficients[1])
}
correct_batch<-function(eset, batch.name){
	mat<-exprs(eset)
	batch<-pData(eset)[, batch.name]
	mat.correct<-t(apply(mat, 1, function(i) lm_correct(y = i, x = batch)))
	exprs(eset)<-mat.correct
	return(eset)
}

if(TRUE){

	#cancercell 2012
	f.dir.in<-"/Users/amyli/Desktop/git_projects/datasets/lymphoma2010"
	f.file<-"lymphoma2010_wold.entrez.mas.gnorm.GS.res"
	f.confounder<-"lymphoma2010_wold.03vs10.cls"
	f.dir.out<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/cancercell2012"
	
	eset<-make_gep(f = paste(f.dir.in, "/", f.file, sep = ""), row.ind = FALSE, 
		f.confounder = paste(f.dir.in, "/", f.confounder, sep = ""))
	saveRDS(eset, file = paste(f.dir.out, "/", "cancercell2012_gep.RDS", sep = ""))

	eset.log<-eset
	mat<-exprs(eset.log)
	mat[mat<0]<-0
	mat<-log2(mat + 1)
	exprs(eset.log)<-mat

	saveRDS(eset.log, file = paste(f.dir.out, "/", "cancercell2012_gep_log.RDS", sep = ""))

	eset.log.correct<-correct_batch(eset = eset.log, batch.name = "batch")
	saveRDS(eset.log.correct, file = paste(f.dir.out, "/", "cancercell2012_gep_log_correct.RDS", sep = ""))

	#ricover
	f.dir.in<-"/Users/amyli/Desktop/git_projects/datasets/ricover2015"
	f.file<-"ricover2014.rma.hgu133plus2hsensgcdf.res"
	#f.confounder<-paste(f.dir.in,"/",  "ricover2014.batches.cls", sep = "")
	f.dir.out<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/wes_7_16_15"
	eset<-make_gep(f = paste(f.dir.in, "/", f.file, sep = ""), row.ind = TRUE,
		f.confounder = NA)


	pData(eset)$batch<-as.numeric(grepl("c_D",colnames(eset)))
	#a$batch<-cls2df(f = f.confounder)[,1]
	
	saveRDS(eset, file = paste(f.dir.out, "/", "ricover_gep.RDS", sep = ""))

	eset.log<-eset
	mat<-exprs(eset.log)
	mat[mat<0]<-0
	mat<-log2(mat + 1)
	exprs(eset.log)<-mat

	saveRDS(eset.log, file = paste(f.dir.out, "/", "ricover_gep_log.RDS", sep = ""))

	eset.log.correct<-correct_batch(eset = eset.log, batch.name = "batch")
	saveRDS(eset.log.correct, file = paste(f.dir.out, "/", "ricover_gep_log_correct.RDS", sep = ""))

}