#' \code{base.summary} adds additional columns to DE summary table from gene expression data
#' @param ge.mat gene expression matrix
#' @param groups sample indicator for groups
#' @param groups.label c(a, b) a = label for positive samples, b = label for negative samples
#' @param logged TRUE if gene expression is log2-transformed, base.summary provides statistics on unlogged data

#additional columns for DE summary table
base.summary<-function(ge.mat, groups, 
	groups.label = c("case", "control"), #label for positive, followed by label for negative samples
	logged = TRUE #are expression values logged?
	){
	if (logged) ge.mat<-2^(ge.mat)
	gene_id<-rownames(ge.mat)
	case<-groups.label[1]
	control<-groups.label[2]
	ge.mat.case<-ge.mat[, groups == case]
	ge.mat.control<-ge.mat[, groups == control]
	mean.case<-apply(ge.mat.case, 1, mean)
	mean.control<-apply(ge.mat.control, 1, mean)
	sd.case<-apply(ge.mat.case, 1, sd)
	sd.control<-apply(ge.mat.control, 1, sd)
	length.case<-ncol(ge.mat.case)
	length.control<-ncol(ge.mat.control)
	class<-sapply(mean.case > mean.control, 
		function(x){if(x == TRUE){case} else{control}})

	FC.actual<-mean.case/mean.control
	logFC.actual<-log2(FC.actual)

	res<-data.frame(gene_id, 
		mean.control, mean.case, 
		sd.control, sd.case, 
		FC.actual, logFC.actual,
		length.control,length.case, class)
	return(res)
}

#' \code{base.summary.wrapper} wrapper for \code{base.summary}
#' @param ge internal ge object
#' @param gistic internal gistic object
#' @param logged TRUE if gene expression is log2-transformed, base.summary.wrapper provides statistics on unlogged data
#' @export
base.summary.wrapper<-function(ge, gistic, logged = TRUE){
	ge.mat<-ge$mat
	gistic.mat<-gistic$mat
	gistic.mat[gistic.mat ==2]<-1
	n_alt<-length(gistic$meta$Descriptor)
	res<-list()
	for (i in 1:n_alt){
		alteration_id<-gistic$meta$Unique.Name[i]
		groups <-as.factor(as.numeric(gistic.mat[i, ]))
		res[[i]]<-base.summary(ge.mat, groups, groups.label = c(1, 0), logged = logged)
		res[[i]]<-cbind(res[[i]], alteration_id)
	}
	return(do.call(rbind, res))
}


#' \code{do_test_limma} performs limma differential expression test
#' @import limma
#' @param mat gene expression matrix
#' @param design design matrix
#' @param alternative must be one of "greater", "less" or "equal", P.Value and adj.P.Val turned to one-sided results if specified
#' @param logged TRUE if gene expression is log2-transformed, limma model is performed on log2 expression values
#' @return limma topTable data frame
#' @export
do_test_limma<-function(mat, design, alternative = "equal", logged = TRUE){
	require(limma)
	#require(miRtest)
	if (!logged){ 
		mat.min<-min(mat)
		if (mat.min <= 0) mat<- mat - mat.min + 1
		mat<-log2(mat)
	}
	fit <- lmFit(mat, design)
	astr<-paste("case", "-", "control", sep = "")
		prestr="makeContrasts("
		poststr=",levels=design)"
		commandstr=paste(prestr,astr,poststr,sep="")
	contrast.matrix<-eval(parse(text = commandstr))
	fit2 <- contrasts.fit(fit, contrast.matrix) #coef case - coef control
	fit2 <- eBayes(fit2)

	fit2.top<-topTable(fit2, coef = 1, adjust.method = "fdr", number = dim(fit2)[1], sort.by = "none")
	n<-length(fit2$t)

	#report one-sided p values
	if (alternative == "less" | alternative == "greater"){
		p.val.onesided <- sapply(1:n, 
			function(i){
				if (alternative == "greater")
					lower.tail<-fit2$t[i] < 0
				else 
					lower.tail <-fit2$t[i] >= 0
				return(pt(q = abs(fit2$t[i]), df = fit2$df.total[i], lower.tail = lower.tail))
			})
		p.val.onesided.adjust<-p.adjust(p.val.onesided, method = "fdr")
		fit2.top$P.Value <- p.val.onesided 
		fit2.top$adj.P.Val <-p.val.onesided.adjust
	}
	return(fit2.top)
}

#' \code{do_gistic2ge} wrapper for performing differential expression test
#' @import plyr
#' @param ge internal ge object
#' @param gistic internal gistic object
#' @param case values to be considered altered samples in gistic lesions file, default c(1,2)
#' @param control values to be considered unaltered samples in gstic lesions file, default c(0)
#' @param interaction_type "cis" or "trans" type of analysis to performed
#' @param test "limma" or "ttest" DE test to be used, default "limma"
#' @param alternative "equal", "greater", or "less" is DE test one-sided? 
#'  rule of thumb for trans alternative = "greater", for cis, alterative = "greater" for amp and "less" for del
#' @param logged TRUE if gene expression is log2-transformed, DE is performed on log2 expression values
#' @return data frame with DE results
#' @export
do_gistic2ge<-function(ge, gistic, 
	case=c(1, 2), 
	control=(0), 
	#alteration_type = "amp", #amp or del
	interaction_type = "cis", 
	test = "limma", #ttest or limma 
	alternative = "equal", #equal, greater, or less
	boundary = "peak", #or "region" or "peak_region"
	logged = TRUE #is ge logged?
	){

	require(plyr)
	
	res<-list()

	#get ge data
	ge.mat<-ge$mat
	ge.mat<-as.matrix(ge.mat)
	ge.mat.col<-colnames(ge.mat)	
	#get and reformat gistic data
	gistic.mat<-gistic$mat
	gistic.sub<-gistic.mat
	gistic.sub.col<-colnames(gistic.sub)
	for(i in case) gistic.sub[gistic.sub == i]<-TRUE
	for(i in control) gistic.sub[gistic.sub == i]<-FALSE
	gistic.sub[gistic.sub == TRUE]<-1
	gistic.sub[gistic.sub == FALSE]<-0

	if (boundary == "region"){
		gistic.cisgenes.peak<-gistic$cisgenes[["genesets.peak"]]
		gistic.cisgenes.region<-gistic$cisgenes[["genesets.region"]]
		gistic.cisgenes<-sapply(1:length(gistic.cisgenes.peak), 
			function(i) {union(gistic.cisgenes.peak[[i]], gistic.cisgenes.region[[i]])})
	}
	else {
		gistic.cisgenes<-gistic$cisgenes[[paste("genesets", ".", boundary, sep = "")]]
	}
	alteration_ids<-gistic$meta$Unique.Name
	alteration_desc<-gistic$meta$Descriptor
	
	all_genes<-rownames(ge.mat)
	cat(paste("total of ", length(alteration_ids), " alterations: \n\n", sep = ""))

	#i: for each alteration
	for (i in 1:dim(gistic.mat)[1]){ 
		cat(paste(i, "\n", sep = ""))
		if (interaction_type == "cis")
			genes.keep<-intersect(all_genes, unique(gistic.cisgenes[[i]]))
		else 
			genes.keep<-setdiff(all_genes, unique(gistic.cisgenes[[i]]))
		
		if (length(genes.keep)<1) res[[i]]<-NA

		else {
			ge.sub<-subset(ge.mat, rownames(ge.mat) %in% genes.keep)
			treatment<-factor(as.numeric(gistic.sub[i,]))
			levels(treatment)<-c("control", "case")
			design.treatment<-data.frame(treatment = treatment)

			if (!("confounder" %in% names(ge))){ # no batch confounder
				design.formula<-paste("~0",  paste(colnames(design.treatment), 
					collapse = " + "), 
					sep = " + ")

				design<-model.matrix(data = design.treatment, object = as.formula(design.formula))
				rownames(design)<-gistic.sub.col
				colnames(design)[1:2]<-c("control", "case") #0 always first column if levels are c(0,1)

			} else {
				design.confounder<-ge$confounder
				n.confounders<-ncol(design.confounder)
				colnames(design.confounder)<-paste("confounder_", 1:n.confounders, sep = "")
				design.combined<-cbind(design.treatment, design.confounder)
				design.formula<-paste("~0",  paste(colnames(design.treatment),collapse = " + "), 
					paste(colnames(design.confounder),collapse = " + "),
					sep = " + ")

				design<-model.matrix(data = design.combined, object = as.formula(design.formula))
				colnames(design)<-c("control", "case", colnames(design.confounder))
			}

			if (any(apply(design[, c(1,2)], 2, sum) < 3)){
				cat(paste("less than three samples in case or control for ", alteration_ids[i], "\n\n", sep = ""))
				res[[i]]<-NA
			}
			else {
				alteration_id<-alteration_ids[i]
				alteration_descriptor<-alteration_desc[i]
				gene_id<-rownames(ge.sub)
				if (test == "limma"){
					res[[i]]<-do_test_limma(mat = ge.sub, design = design, alternative = alternative, logged = logged)
					if (is.data.frame(res[[i]]))
						res[[i]]<-cbind(alteration_id, alteration_descriptor, interaction_type, gene_id, res[[i]])					
				} else if (test == "ttest"){
					res[[i]]<-do_test_ttest(mat = ge.sub, design = design, alternative = alternative, logged = logged)
					if (is.data.frame(res[[i]]))
						res[[i]]<-cbind(alteration_id, alteration_descriptor, interaction_type, gene_id, res[[i]])
				} else {
					res[[i]]<-NA
					cat ("invalid test\n")
				}
			}
		}
	}
	res.df<-do.call(rbind, res[!is.na(res)])
	#res.df$adj.P.Val.all<-p.adjust(res.df$P.Value, method = "fdr")
	return(res.df)
}

#' \code{do_gistic2ge_all} wrapper for \code{do_gistic2ge} for performing multiple differential expression tests
#' @param gistic_ge_amp internal gistic_ge object for amplifications
#' @param logged TRUE if gene expression is log2-transformed, DE is performed on log2 expression values
#' @return list of data frames with DE results
#' @export
do_gistic2ge_all<-function(gistic_ge_amp, 
	gistic_ge_del, 
	logged = FALSE){

	gistic2ge.map<-list()
	gistic2ge.map[["ge_amp"]]<-gistic_ge_amp$ge
	gistic2ge.map[["ge_del"]]<-gistic_ge_del$ge
	gistic2ge.map[["gistic_amp"]]<-gistic_ge_amp$gistic
	gistic2ge.map[["gistic_del"]]<-gistic_ge_del$gistic
	gistic2ge.map[["summary_amp"]]<-base.summary.wrapper(ge = gistic2ge.map[["ge_amp"]], gistic = gistic2ge.map[["gistic_amp"]], logged = logged)
	gistic2ge.map[["summary_del"]]<-base.summary.wrapper(ge = gistic2ge.map[["ge_del"]], gistic = gistic2ge.map[["gistic_del"]], logged = logged)

	gistic2ge<-list()
	#gistic2ge.sig<-list()

	gistic2ge_new<-list()
	for (boundary in c("region", "peak")){	  
	   	for (j in c("cis", "trans")){
	   		alt_dir<-c("amp", "del")
	   		for (i in alt_dir){

	        if (j == "trans"){
	          alternative <- "equal"
	        } else {
	          if (i == "amp") alternative <- "greater"
	          else alternative <- "less"
	        }

	        id<-paste(boundary, i,j, sep = "_")
	        print(id)

	        gistic2ge[[id]]<-do_gistic2ge(ge = gistic2ge.map[[paste("ge_", i, sep = "")]], 
	          gistic = gistic2ge.map[[paste("gistic_", i, sep = "")]], 
	          case=c(1, 2), 
	          control=(0), 
	          interaction_type = j, 
	          test = "limma", #t or limma  
	          alternative = alternative,           
	          boundary = boundary, 
	          logged = logged)

	        summary.name<-paste("summary_", i, sep ="")
	        df.inner<-join(gistic2ge[[id]], gistic2ge.map[[summary.name]], by = c("alteration_id", "gene_id"))
	        genes.ind<-match( df.inner$gene_id, gistic2ge.map[[paste("ge_", i, sep = "")]]$accession)
	        df.inner$gene_description <-gistic2ge.map[[paste("ge_", i, sep = "")]]$description[genes.ind]
	        gistic2ge[[id]]<-df.inner
	    }	    

	    df<-rbind(gistic2ge[[paste(boundary, alt_dir[1], j, sep = "_")]],
	    		gistic2ge[[paste(boundary, alt_dir[2], j, sep = "_")]])
	    df$adj.P.Val.all<-p.adjust(df$P.Value, method = "fdr")
	   	df$FC<-2^(df$logFC) #fold change
	    df$FC.lm<-df$FC
        df$logFC.lm<-df$logFC
        col.order<-c("alteration_id", "alteration_descriptor", "interaction_type","gene_id",
	        	"adj.P.Val", "adj.P.Val.all", "P.Value", "class", "FC.lm",  "logFC.lm", "t", 
	        	"mean.case", "mean.control", "sd.case", "sd.control", "FC.actual", "logFC.actual", 
	        	"length.case", "length.control", "gene_description")

	  	df<-df[, col.order]
	    gistic2ge_new[[paste(boundary, j, sep = "_")]]<-df
   
	  }
	}
	return(gistic2ge_new)
}

subset_gistic2ge<-function(x, #list of data frames
	cis.fdr, trans.fdr, sig.FC){
	res<-list()
	for (boundary in c("region", "peak")){	  
	   	for (j in c("cis", "trans")){
	   		if (j == "cis") sig.fdr <- cis.fdr
	    	else sig.fdr <- trans.fdr
	    	df<-x[[paste(boundary, j, sep = "_")]]
	    	if (is.na(sig.FC)) 
	    	df.sig<-subset(df, adj.P.Val < sig.fdr)
	    	else 
	    	df.sig<-subset(df, adj.P.Val < sig.fdr & (logFC.lm > log2(sig.FC) | logFC.lm < -log2(sig.FC)))
	    	res[[paste(boundary, j, sep = "_")]]<-df.sig

	   	}
	}
	return(res)

}