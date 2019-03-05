
#' Construct iEDGE input object
#' @param cn matrix of epi-DNA alterations (m rows of alterations by n samples)
#' @param gep matrix of gene expression (g rows of genes by n samples)
#' @param cisgenes list of cis genes corresponding to each row in cn
#' @param transgenes optional list of trans genes corresponding to each row in cn
#' @param cn.fdat data frame of row annotations for cn, must have column of row identifiers, and optionally column of direction identifiers if one-sided DE tests are applied
#' @param gep.fdat data frame of row annotaions for gep, must have column of row identifiers
#' @param cn.pdat optional data frame of column annotations for cn
#' @param gep.pdat optional data frame of column annotations for gep
#' @export 
construct_iEDGE<-function(cn, gep, cisgenes, transgenes = NA, cn.fdat, gep.fdat, cn.pdat = NA, gep.pdat = NA){
	if(suppressWarnings(is.na(cn.pdat)))
		cn.pdat<-data.frame(colid = colnames(cn))
	if(suppressWarnings(is.na(gep.pdat)))
		gep.pdat<-data.frame(colid = colnames(gep))
	
	cn<-to.eSet(mat = cn, pdat = cn.pdat, fdat = cn.fdat)
	gep<-to.eSet(mat = gep, pdat = gep.pdat, fdat = gep.fdat)

	if(suppressWarnings(is.na(transgenes)))
		dat<-list(cn = cn, gep = gep, cisgenes = cisgenes)
	else 
		dat<-list(cn = cn, gep = gep, cisgenes = cisgenes, transgenes = transgenes)

	return(dat)
}


#' Main wrapper for iEDGE, assumes input data is constructed by construct_iEDGE
#' @import parallel
#' @param   dat the input iEDGE object, see construct_iEDGE for specification
#' @param	header header string for result file names
#' @param	outdir output directory
#' @param	gs.file default = NA, vector of characters giving full path of gmt files for enrichment analysis
#' @param	gepid default = "SYMBOL", colname in fData(dat$gep) indicating unique gene ids to display
#' @param	cnid default = "Unique.Name", colname in fData(dat$cn) indicating unique alterations to display
#' @param	cndesc default = "Descriptor", column in fData(dat$cn), can be NA if not available
#' @param	cndir default = "alteration_direction", column in fData(dat$cn) indicating directions of one sided differential expression, use NA if onesided.cis = FALSE and onesided.trans = FALSE
#' @param	fdr.cis.cutoff default = 0.25, fdr cis cutoff
#' @param	fdr.trans.cutoff default = 0.05, fdr trans cutoff
#' @param	fc.cis default = NA, fold change cutoff cis
#' @param	fc.trans default = NA, fold change cutoff trans
#' @param	min.drawsize default = 3, min drawsize for geneset enrichment
#' @param	onesided.cis default = TRUE, one sided cis differential expression
#' @param	onesided.trans default = FALSE, one sided trans differential expression
#' @param	uptest default = "Amplification", alterations for which upregulation in alteration is desired, value must be on fData(dat$cn)[, cndir], must specify if onesided.cis =TRUE
#' @param	downtest = "Deletion", alterations for which downregulation in alteration is desired
#' @param	min.group default = 2, mininum sample size in each group 
#' @param	prune.col default = "pvalue", column name to indicate significance of sobel test
#' @param	prune.thres default = 0.05, significance cutoff of sobel test
#' @param	hyperthres default = 0.25, threshold for hyperEnrichment (fdr cutoff)
#' @param	cis.boxplot default = TRUE, display boxplot for cis genes
#' @param	trans.boxplot default = TRUE, display boxplot for trans genes
#' @param 	enrich.heatmap = TRUE, display heatmap for pathway enrichment results
#' @param	bipartite default = TRUE, do bipartite graph/pruning
#' @param	html default = TRUE, do html report
#' @param	jsdir default = NA, default directory of iEDGE js files
#' @param	numcores default = NA (non-parallel), if specified uses mclapply with specified mc.cores
#' @param	cache default = list(DE = NULL, prunning = NULL, ui = NULL), optional, cached result of previous run_iEDGE 
#' @param includeheatmap = TRUE, optional, include heatmap for pathway enrichment in UI, must be enrichment=TRUE for this to work
#' @export
run_iEDGE<-function(dat, #iEDGE object
	header, #header string for result file names
	outdir, #directory for results
	gs.file = NA, #vector of characters giving full path of gmt files for enrichment analysis
	gepid = "SYMBOL", #colname in fData(dat$gep) indicating unique gene ids to display
	cnid = "Unique.Name", #colname in fData(dat$cn) indicating unique alterations to display
	cndesc = "Descriptor", #column in fData(dat$cn), can be NA if not available
	cndir = "alteration_direction", #column in fData(dat$cn) indicating directions of one sided differential expression
	#use NA if onesided.cis = FALSE and onesided.trans = FALSE
	fdr.cis.cutoff = 0.25, #fdr cis cutoff
	fdr.trans.cutoff = 0.05, #fdr trans cutoff
	fc.cis = NA, #fold change cutoff cis
	fc.trans = NA, #fold change cutoff trans
	min.drawsize = 3, #min drawsize for geneset enrichment
	onesided.cis = TRUE, #is cis DE one sided, default TRUE
	onesided.trans = FALSE, #is trans DE one sided, default FALSE
	uptest = "Amplification", #alterations for which upregulation in alteration is desired
		#value must be on fData(dat$cn)[, cndir]
		#must specify if onesided.cis =TRUE
	downtest = "Deletion", #alterations for which downregulation in alteration is desired
		#must specify if onesided.cis =TRUE
		#value must be on fData(dat$cn)[, cndir]
	min.group = 2, #mutinfo.seed = 7, mutinfo.nsamples = 500, mutinfo.bins = 5,  
	prune.col = "pvalue", #column name to indicate significance of sobel test
	prune.thres = 0.05, #significance cutoff of sobel test
	hyperthres = 0.25, #threshold for hyperEnrichment (fdr cutoff)
	cis.boxplot = TRUE, #display boxplot for cis genes
	trans.boxplot = TRUE, #display boxplot for trans genes
	#enrich.heatmap = TRUE, #display heatmap for pathway enrichment results
	enrichment = TRUE, #do pathway enrichment
	bipartite = TRUE, #do bipartite graph/pruning
	html = TRUE, #do html report
	jsdir = NA, #default directory of iEDGE js files
	numcores = NA, #default NA non-parallel, if specified uses mclapply with specified mc.cores
	cache = list(DE = NULL, DE.enrich = NULL, prunning = NULL, pruning.enrich = NULL, ui = NULL), #optional, cached result of previous run_iEDGE 
	includeheatmap = TRUE
	){

	res<-NULL
	de.enrich<-NULL
	pruning<-NULL
	ui<-NULL
	if(!enrichment) includeheatmap = FALSE

	do.DE<-TRUE
	do.de.enrich<-TRUE
	do.pruning<-TRUE
	do.pruning.enrich<-TRUE
	do.ui<-TRUE



	if(!is.null(cache[["de"]])) {
		de<-cache[["de"]]
		do.DE<-FALSE
	}

	if(!is.null(cache[["de.enrich"]])){
		de.enrich<-cache[["de.enrich"]]
		do.de.enrich<-FALSE
	}

	if(!is.null(cache[["pruning"]])){
		pruning<-cache[["pruning"]]
		do.pruning<-FALSE
	}

	if(!is.null(cache[["pruning.enrich"]])){
		pruning.enrich<-cache[["pruning.enrich"]]
		do.pruning.enrich<-FALSE
	}

	#if(!is.null(cache[["ui"]])){
	#	ui<-cache[["ui"]]
	#	do.ui<-FALSE
	#}

	do.gs<-TRUE


	if(is.na(gs.file[1])){ 
		do.gs<-FALSE
	} else {
		if(any(!file.exists(gs.file))){
			stop("gs.file must be vector of gmt file names, full paths only")
		}
	}

	cn<-dat$cn

	#check duplicated altids
	if(any(duplicated(fData(dat$cn)[, cnid]))){
		stop("Error: duplicated entries in fData(dat$cn)[, cnid], make unique...")
	}

	gep<-dat$gep
	cisgenes<-dat$cisgenes
	
	if("transgenes" %in% names(dat))
		transgenes<-dat$transgenes
	else transgenes <- NA

	suppressWarnings(dir.create(outdir, recursive = TRUE))
	base_dir<-paste(outdir,"/", header, sep = "")

	suppressWarnings(dir.create(base_dir, recursive = TRUE))
	de_dir<-paste(base_dir, "/tables", sep = "")

	if(do.DE){
		cat(paste("Running iEDGE for data set: ", header, "\n",sep = ""))

		cat("Making Differential Expression tables...\n")
		de<-iEDGE_DE(cn, gep, cisgenes, transgenes,
			header,
			gepid, cnid,  	
			f.dir.out = de_dir, 
			fdr.cis.cutoff = fdr.cis.cutoff, fdr.trans.cutoff = fdr.trans.cutoff, 
			min.group = min.group,
			min.drawsize = min.drawsize,  
			cndir = cndir,
			cis.onesided = onesided.cis, 
			trans.onesided = onesided.trans,
			fc.cis = fc.cis,
			fc.trans = fc.trans,
			uptest = uptest,
			downtest = downtest, 
			numcores = numcores
			)
	}
	if(do.de.enrich){
		if(enrichment == TRUE){
			print("enrichment ")

			gs<-NA
			if(do.gs){
				cat("Reading genesets..\n")
				gs<-lapply(gs.file, function(i){
					res<-read_gmt(i)$genesets
					res<-sapply(res, function(x){
						a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
						return(unique(do.call(c,a)))
					})
					cat(paste("Geneset loaded:", i, "\n", sep = ""))
					return(res)
					})
				names(gs)<-gsub(".gmt", "",basename(gs.file))


			}

			de.enrich<-do_de_enrichment(res = de, header, gep, gepid, cnid, gs, f.dir.out = de_dir)
		}
	}

	de.cis.sig<-de[["cis"]][["sig"]]
	de.cis.full<-de[["cis"]][["full"]]
	de.trans.sig<-de[["trans"]][["sig"]]

	pruning_dir<-paste(base_dir, "/pruning", sep = "")

	if(do.pruning){
		if(bipartite == TRUE){	
			pruning<-prune(f_cis_tab =  de.cis.sig, 
				f_trans_tab = de.trans.sig, 
				cn = cn, 
				gep = gep,
				alteration_id = cnid,
				gene_id = gepid,
				pruning_dir = pruning_dir, 
				gs = gs,
				prunecol = prune.col, prunethres = prune.thres, 
				min.drawsize = min.drawsize, 
				hypercol = "fdr", 
				hyperthres = hyperthres)
		} 
	}

	if(do.ui){
		if(html == TRUE){
			html_dir<-paste(base_dir, "/html", sep = "")
			if(is.na(jsdir))
				jsdir<-file.path(path.package("iEDGE"), "javascript")
	
			ui<-iEDGE_UI(cistab = de.cis.sig, cisfulltab = de.cis.full,
				transtab = de.trans.sig, cn = cn, gep = gep, cisgenes = cisgenes,
				outdir = html_dir, jsdir = jsdir, 
				pruning = pruning, pruningjsdir = paste(pruning_dir, "/js", sep = ""),
				altid = cnid, altdesc = cndesc, geneid = gepid, 
				cis.boxplot = cis.boxplot, trans.boxplot = trans.boxplot, bipartite = bipartite,
				heatmap = includeheatmap, heatmapdatadir = paste0(base_dir, "/tables"),
				header = header, gs.names = gs.file)
		} 
	}

	return(list(de = de, de.enrich = de.enrich, pruning = pruning, ui = ui))
}

#' Performs sobel test of mediation
#' @param x binary numeric vector of epi-DNA states
#' @param y matrix of cis gene expression ncol(x) must be equal to length(x)
#' @param z matrix of trans gene expression ncol(z) must be equal to length(x)
#' @param y.names names of cis genes in y
#' @param z.names names of trans genes in z
#' @export
calc_sobel<-function(x,y,z, y.names, z.names){
	ny<-nrow(y)
	nz<-nrow(z)
	res<-NULL
	if(ny>1 & nz>1){
		res<-calc_sobel_mat(x,y,z)
	} else if (ny>1 & nz == 1){
		res<-calc_sobel_y_z1(x,y,z)
	} else if (ny == 1 & nz > 1){
		res<-calc_sobel_y1_z(x,y,z)
	} else if (ny == 1 & nz == 1){
		res<-calc_sobel_y1_z1(x,y,z)
	}
	res<-data.frame(cis = y.names[res[, "yind"]], trans = z.names[res[, "zind"]], res)
	return(res)
}

# lapply wrapper preserving the names
lapplyn<-function(dat, fn){
	datn<-names(dat)
	res<-lapply(dat, fn)
	names(res)<-datn
	return(res)
}

#' @import Biobase
run_limma_accessory<-function(eset, design, log_offset = 1, log_base = 2){
	eset.sub.cond1<-eset[,design[,1] == 1]
	eset.sub.cond2<-eset[,design[,2] == 1]

	exprs.cond1<-exprs(eset.sub.cond1)
	exprs.cond2<-exprs(eset.sub.cond2)

	mean.cond1<-apply(exprs.cond1, 1, function(i) mean(i, na.rm = TRUE))
	sd.cond1<-apply(exprs.cond1, 1, function(i) sd(i, na.rm = TRUE))
	n.cond1<-apply(exprs.cond1, 1, function(i) sum(!is.na(i)))

	mean.cond2<-apply(exprs.cond2, 1, function(i) mean(i, na.rm = TRUE))
	sd.cond2<-apply(exprs.cond2, 1, function(i) sd(i, na.rm = TRUE))
	n.cond2<-apply(exprs.cond2, 1, function(i) sum(!is.na(i)))

	eset.unlog<-eset
	exprs(eset.unlog)<-log_base^exprs(eset.unlog) - log_offset

	eset.sub.cond1<-eset.unlog[,design[,1] == 1]
	eset.sub.cond2<-eset.unlog[,design[,2] == 1]

	exprs.cond1.unlog<-exprs(eset.sub.cond1)
	exprs.cond2.unlog<-exprs(eset.sub.cond2)

	mean.cond1.unlog<-apply(exprs.cond1.unlog,1, function(i) mean(i, na.rm = TRUE) )
	sd.cond1.unlog<-apply(exprs.cond1.unlog,1, function(i) sd(i, na.rm = TRUE))
	
	mean.cond2.unlog<-apply(exprs.cond2.unlog,1, function(i) mean(i, na.rm = TRUE) )
	sd.cond2.unlog<-apply(exprs.cond2.unlog,1, function(i) sd(i, na.rm = TRUE))

	high.class<-mean.cond2 >= mean.cond1
	high.class[high.class == TRUE]<- as.character(colnames(design)[2])
	high.class[high.class == FALSE]<- as.character(colnames(design)[1])
	res<-data.frame(mean0.log = mean.cond1, mean1.log = mean.cond2, 
		sd0.log = sd.cond1, sd1.log = sd.cond2,
		n0 = n.cond1, n1 = n.cond2,
		mean0.unlog = mean.cond1.unlog, mean1.unlog = mean.cond2.unlog,
		sd0.unlog = sd.cond1.unlog, sd1.unlog = sd.cond2.unlog,
		high.class = high.class)
	res<-cbind(res, fData(eset))
	return(res)
}

#' @import limma Biobase plyr
run_limma<-function(eset, design, 
	 onesided = TRUE,
	 cndir = "alteration_direction", 
	 uptest = "Amplification", 
	 downtest = "Deletion"){

	fit <- lmFit(eset, design)
	command_str<-paste("makeContrasts(",
	"(", colnames(design)[2] , "-", colnames(design)[1], ")", 
	 ",levels = design)", sep = "")
	contrast.matrix<-eval(parse(text =command_str)) 
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	
	#bug in number(length(fit2)) due to limma new version? revised to nrow(fit2)
	#fit2.table<-topTable(fit2, coef=1,adjust.method="BH", number =length(fit2), sort.by = "none")
	
	fit2.table<-topTable(fit2, coef=1,adjust.method="BH", number =nrow(fit2), sort.by = "none")
	

	if (onesided){

		if(!(cndir %in% colnames(fData(eset)))){
			stop("cndir must be in colnames(fData(eset))")
		}
		if(any(!(fData(eset)[, cndir] %in% c(uptest, downtest)))){
			stop("cndir must contain only labels from uptest and downtest")
		}
		p.val.onesided <- sapply(1:nrow(fit2.table), 
			function(i){
				if (fit2.table[i, cndir] == uptest){
				
					lower.tail<-fit2$t[i] < 0
				}
				else if(fit2.table[i, cndir] == downtest){
					
					lower.tail <-fit2$t[i] >= 0

				}
				return(pt(q = abs(fit2$t[i]), df = fit2$df.total[i], lower.tail = lower.tail))
			})
		p.val.onesided.adjust<-p.adjust(p.val.onesided, method = "fdr")
		fit2.table$P.Value <- p.val.onesided 
		fit2.table$adj.P.Val <-p.val.onesided.adjust
	}
	fit2.accessary<-run_limma_accessory(eset, design)
	fit2.table<-join(x=fit2.table, y=fit2.accessary, 
		by = colnames(fData(eset)), type = "left", match = "first")
	fit2.table$fold.change<-2^fit2.table$logFC
	fit2.table<-fit2.table[order(fit2.table$adj.P.Val, decreasing = FALSE),]
	return(fit2.table)
}

#' @import Biobase
iEDGE_DE_inner<-function(gep, #eset containing log2 gene expression
 cn, #eset containing copy number by sample, must be sorted with respect to columns in gep
 cisgenes, #list of cisgenes in gep with respect to cn ordered by rows in cn
 transgenes = NA, #optionally list of transgenes with respect ot cn, if not specified will consider all genes
 			 #with expeption of cisgenes
 gepid, #column name in fData of gep identifying the genes in cisgenes
 cnid, #column name in fData of cn identifying the alterations in names(cn)
 cndir, #column name in fData of cn identifying direction of test 
 #(Amplification = high expression in cn = 1 low in 0, Deletion = low in 1 high in 0)
 interaction_type = "cis", #cis or trans
 fdr.cutoff = 0.25,
 min.group = 3,
 onesided = TRUE,
 uptest = "Amplification", 
 downtest = "Deletion", 
 fc = NA, 
 numcores = NA
){

	if(onesided & !(cndir %in% colnames(fData(cn)))){
		stop ("If onesided is specified, cndir must be in colnames(fData(cn))")
		if(any(!(fData(cn)[, cndir] %in% c(uptest, downtest)))){
			stop("fData(cn)[, cndir] must contain elements from uptest or downtest")
		}
	}
	all_genes<-as.character(fData(gep)[, gepid])
	cat(paste("Running iEDGE differential expression for ", nrow(cn), " alterations: \n\n", sep = ""))

	get_genes<-function(i){
		if (interaction_type == "cis")
			genes.keep<-intersect(all_genes, unique(cisgenes[[i]]))
		else {
			if(is.na(transgenes))
				genes.keep<-setdiff(all_genes, unique(cisgenes[[i]]))
			else 
				genes.keep<-setdiff(transgenes[[i]], unique(cisgenes[[i]]))
		}
		return(genes.keep)

	}

	de.inner<-function(i) { 	#for each alteration
		genes.keep<-genes.keep.list[[i]]
		if (length(genes.keep)<1) res<-NA
		else {

			gep.keep.ind<-match(genes.keep, all_genes)
			gep.keep.ind<-gep.keep.ind[!is.na(gep.keep.ind)]
			gep.keep<-gep[gep.keep.ind,]
			fdat<-fData(gep.keep)

			fdat.add<-fData(cn)[i, ,drop = FALSE]
			fdat.combined<-cbind(fdat, fdat.add)

			if (interaction_type == "cis")	fdat.combined$interaction_type <-"cis"
			else 
				fdat.combined$interaction_type <-"trans"

			fData(gep.keep)<-fdat.combined
			treatment.raw<-suppressWarnings(as.numeric(exprs(cn)[i,]))
			#remove na indices
			rmna<-which(!is.na(treatment.raw))
			gep.keep<-gep.keep[,rmna]
			treatment.raw<-treatment.raw[rmna]

			treatment<-factor(treatment.raw)
			levels(treatment)<-c("control", "case")
			design.treatment<-data.frame(treatment = treatment)

			design.formula<-paste("~0",  paste(colnames(design.treatment), 
				collapse = " + "), 
				sep = " + ")
			design<-model.matrix(data = design.treatment, object = as.formula(design.formula))
			colnames(design)[1:2]<-c("control", "case") #0 always first column if levels are c(0,1)
		

			if (any(apply(design[, c(1,2)], 2, sum) < min.group)){
				cat(paste("less than ", min.group," samples in case or control for ", 
					i, ", cannot run limma", "\n\n", sep = "")
				)
				res<-NA
			}
			else 
				res<-run_limma(eset = gep.keep, design = design, 
					cndir = cndir, onesided = onesided, 
					uptest = uptest, downtest = downtest)
				
		}

		percentProgress<-do.call(sum, lapply(genes.keep.list[1:i], length))/genes.keep.total
		setTxtProgressBar(pb, percentProgress)
		
		return(res)
	}

	n<-nrow(cn)	
	genes.keep.list<-lapply(1:n, function(i) get_genes(i))
	genes.keep.total<-do.call(sum, lapply(genes.keep.list, length))

	if(is.na(numcores)){
		pb<-txtProgressBar(style = 3)
		setTxtProgressBar(pb, 0)
		res<-lapply(1:n, function(i) de.inner(i))
		setTxtProgressBar(pb, 1)
		close(pb)
	} else {
		pb <- txtProgressBar(style = 3)
		setTxtProgressBar(pb, 0)
		res<-mclapply(1:n, function(i) de.inner(i), mc.cores = numcores)
		setTxtProgressBar(pb, 1)
		close(pb)
	}

	res.df<-do.call(rbind, res[!is.na(res)])	
	res.df[, "adj.P.Val.all"]<-p.adjust(res.df$P.Value, method = "fdr")
	
	#reorder columns
	col.first<-c(gepid, "fold.change", "adj.P.Val.all", "t", 
		"high.class", "mean1.unlog", "mean0.unlog", "sd0.unlog", "sd1.unlog")
	col.ord<-c(col.first, setdiff(colnames(res.df), col.first))
	res.df<-res.df[, col.ord]
	
	res.list<-list()
	res.list$full<-res.df
	if(is.na(fc)){	
		res.list$sig<-res.df[res.df[, "adj.P.Val.all"]<fdr.cutoff,]
	} else {
		fc.low<-2^(-log2(fc))
		res.list$sig <- res.df[res.df[, "adj.P.Val.all"]<fdr.cutoff & 
		(res.df[, "fold.change"] > fc | res.df[, "fold.change"] < fc.low),]
	}
	return(res.list)
}

get_sig_split<-function(tab, cnid, gepid){
	res.names<-unique(tab[, cnid])
	res<-lapply(res.names, 
			function(x){
				y<-tab[tab[, cnid] == x,]
				return(unique(as.character(y[, gepid])))
				}#, USE.NAMES = TRUE
				)
	names(res)<-res.names
	return(res)
}

# get_sig_split<-function(tab, cnid, gepid){
# 	res<-sapply(unique(tab[, cnid]), 
# 			function(x){
# 				y<-tab[tab[, cnid] == x,]
# 				return(unique(as.character(y[, gepid])))
# 				}, USE.NAMES = TRUE)
# 	return(res)
# }

get_sig_single<-function(tab, gepid, nm){
	res<-list()
	res[[nm]]<-unique(as.character(tab[, gepid]))
	return(res)
}

#' iEDGE_DE performs differential expression analysis and pathway enrichment and saves to text tables
iEDGE_DE<-function(cn, gep, cisgenes, transgenes,
	header,
	gepid, cnid,  	
	f.dir.out, 
##	gs, #genset for hyper
	fdr.cis.cutoff = 0.25, fdr.trans.cutoff = 0.05, 
	min.group = 3,
	min.drawsize = 3,  
	cndir = "alteration_direction",
	cis.onesided = TRUE, 
	trans.onesided = FALSE,
	fc.cis = NA,
	fc.trans = NA,
#	enrich.heatmap = FALSE,
	... #other parameters in iEDGE_DE_inner
	){

	suppressWarnings(dir.create(f.dir.out, recursive =TRUE))
	cat("Performancing cis analysis...\n")
	res.cis<-iEDGE_DE_inner(gep = gep, 
 		cn = cn,
 		cisgenes = cisgenes, 
 		transgenes = transgenes,
 		gepid = gepid, 
 		cnid = cnid,
 		cndir = cndir,
 		interaction_type = "cis",
 		fdr.cutoff = fdr.cis.cutoff, 
 		onesided = cis.onesided,
 		min.group = min.group,
 		fc = fc.cis,
 		...)


	if(is.na(fc.cis)) fc_cis_header <- ""
	else fc_cis_header<-paste("_fc_", fc.cis, sep = "")
	
	f.out<-paste(f.dir.out, "/", header, "_cis_sig_fdr_", 
		fdr.cis.cutoff, fc_cis_header, ".txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(res.cis$sig, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)

	f.out<-paste(f.dir.out, "/", header, "_cis_full.txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(res.cis$full, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)

	cat("Performancing trans analysis...\n")
	res.trans<-iEDGE_DE_inner(gep = gep, 
 		cn = cn,
 		cisgenes = cisgenes, 
 		transgenes = transgenes,
 		gepid = gepid, 
 		cnid = cnid,
 		cndir = cndir,
 		interaction_type = "trans",
 		fdr.cutoff = fdr.trans.cutoff, 
 		onesided = trans.onesided,
 		min.group = min.group,
 		fc = fc.trans,
 		...)

	res.cis.sig<-res.cis$sig
	res.trans.sig<-res.trans$sig
	res.trans.sig.up<-res.trans.sig[res.trans.sig[, "high.class"] %in% "case",]
	res.trans.sig.dn<-res.trans.sig[res.trans.sig[, "high.class"] %in% "control",]
	
	res.cistrans.sig<-rbind(res.cis.sig, res.trans.sig)
	res.cistrans.sig.up<-rbind(res.cis.sig, res.trans.sig.up)
	res.cistrans.sig.dn<-rbind(res.cis.sig, res.trans.sig.dn)

	if(is.na(fc.trans)) fc_trans_header <- ""
	else fc_trans_header<-paste("_fc_", fc.trans, sep = "")

	f.out<-paste(f.dir.out, "/", header, "_trans_sig_fdr_", 
		fdr.trans.cutoff, fc_trans_header, ".txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(res.trans$sig, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)
	
	cat("Done!\n")


	return(list(cis=res.cis, trans =res.trans#, 
		#hyper = hyper, hyperhm = hyperhm
		))
}

#do enrichment for differential expression
do_de_enrichment<-function(res, header, gep, gepid, cnid, gs, f.dir.out){
	cat("\nPerforming hyperEnrichment analysis...\n")
 	ngenes<-nrow(gep)
 	#make lists of drawns genes
 	res.cis<-res$cis
 	res.trans<-res$trans
	res.cis.sig<-res.cis$sig
	res.trans.sig<-res.trans$sig
	res.trans.sig.up<-res.trans.sig[res.trans.sig[, "high.class"] %in% "case",]
	res.trans.sig.dn<-res.trans.sig[res.trans.sig[, "high.class"] %in% "control",]
	
	res.cistrans.sig<-rbind(res.cis.sig, res.trans.sig)
	res.cistrans.sig.up<-rbind(res.cis.sig, res.trans.sig.up)
	res.cistrans.sig.dn<-rbind(res.cis.sig, res.trans.sig.dn)

	#if(is.na(fc.trans)) fc_trans_header <- ""
	#else fc_trans_header<-paste("_fc_", fc.trans, sep = "")

	drawns<-list()
	drawns[["cis"]]<-get_sig_single(res.cis.sig, gepid, "cis")
	drawns[["trans"]]<-get_sig_single(res.cis.sig, gepid, "trans")
	drawns[["trans.up"]]<-get_sig_single(res.trans.sig.up, gepid, "trans.up")
	drawns[["trans.dn"]]<-get_sig_single(res.trans.sig.dn, gepid, "trans.dn")

	drawns[["cistrans"]]<-get_sig_single(res.cistrans.sig, gepid, "cistrans")
	drawns[["cistrans.up"]]<-get_sig_single(res.cistrans.sig, gepid, "cistrans.up")
	drawns[["cistrans.dn"]]<-get_sig_single(res.cistrans.sig, gepid, "cistrans.dn")

	drawns[["cis_split"]]<-get_sig_split(res.cis.sig, cnid, gepid)
	drawns[["trans_split"]]<-get_sig_split(res.trans.sig, cnid, gepid)
	drawns[["trans.up_split"]]<-get_sig_split(res.trans.sig.up, cnid, gepid)
	drawns[["trans.dn_split"]]<-get_sig_split(res.trans.sig.dn, cnid, gepid)

	drawns[["cistrans_split"]]<-get_sig_split(res.cistrans.sig, cnid, gepid)
	drawns[["cistrans.up_split"]]<-get_sig_split(res.cistrans.sig.up, cnid, gepid)
	drawns[["cistrans.dn_split"]]<-get_sig_split(res.cistrans.sig.dn, cnid, gepid)

 	hyper<-lapply(names(gs), function(i){
 		cat(paste0("Reading geneset:", i, "\n"))
	 	gs.i<-gs[[i]]
	 	gs.file.name<-i
	 	hyper.res<-lapply(names(drawns), function(j){
	 		run_hyperEnrichment_unpruned(ngenes = ngenes, 
	 		gs=gs.i, gs.file.name =gs.file.name, 
	 		drawnList = drawns[[j]], f.dir.out =f.dir.out, 
	 		header = header, header2= j, min.drawsize = min.drawsize, 
	 		verbose = FALSE)
	 		})
	 	names(hyper.res)<-names(drawns)
		return(hyper.res)
	})
	names(hyper)<-names(gs)

	return(list(de.enrich = hyper))
}

calc_sobel_y_z1<-function(x,y,z){

	m1<-lm(t(z) ~ x)
	m2<-lm(t(y) ~ x)
	m3<-lapply(1:nrow(y), function(i){
		lm(t(z) ~ x+y[i, ])
		})

	m2.summary<-summary(m2)
	ny<-nrow(y)

	tau0<-m1$coefficients[2]

	tau1<-lapply(1:ny, function(i) m3[[i]]$coefficients[2])
	beta<-lapply(1:ny, function(i) m3[[i]]$coefficients[3])
	alpha<-m2$coefficients[2,]
	taudiff<-lapply(1:ny, function(i) sign(tau0)*(tau0 - tau1[[i]]))
	sa<-unlist(lapply(1:ny, function(i) m2.summary[[i]]$coefficients[2, "Std. Error"]))

	sb<-lapply(1:ny, function(i) {
			m3_i<-summary(m3[[i]])
			return( m3_i$coefficients[3, "Std. Error"])
		})

	S<-lapply(1:ny, function(i) {
			sqrt(beta[[i]]^2 * sa[i]^2 + alpha[i]^2*sb[[i]]^2)
		})

	Z<-lapply(1:ny, function(i) {
			taudiff[[i]]/S[[i]]
		})

	Pvalue<-lapply(1:ny, function(i) {
			pnorm(-Z[[i]])
		})

	res<-data.frame(yind = vector(), zind = vector(), value = vector())
	for(i in 1:ny){
		res.add<-c(yind = i, zind = 1, 
			taudiff = taudiff[[i]],
			tau0 = tau0,
			tau1 = tau1[[i]],
			tauratio = taudiff[[i]]/(sign(tau0)*tau0),
			Z = Z[[i]],
			S = S[[i]],
			pvalue = Pvalue[[i]])
		res<-rbind(res, res.add)
	}

	colnames(res)<-c("yind", "zind", "taudiff", "tau0", "tau1", "tauratio", "Z", "S", "pvalue")
	res[, "fdr"]<-p.adjust(res[, "pvalue"], method = "fdr")

	return(res)
}

calc_sobel_y1_z<-function(x,y,z){

	m1<-lm(t(z) ~ x)
	m2<-lm(t(y) ~ x)
	m3<-lm(t(z) ~ x+t(y))
	
	m2.summary<-summary(m2)
	nz<-nrow(z)
	tau0 <- m1$coefficients[2,]

	tau1<-m3$coefficients[2,]
	beta<-m3$coefficients[3,]
	alpha<-m2$coefficients[2]

	taudiff<-sign(tau0)*(tau0 - tau1)
	sa<-m2.summary$coefficients[2, "Std. Error"]

	m3.summary<-summary(m3)
	sb<-unlist(lapply(1:nz, function(j) m3.summary[[j]]$coefficients[3, "Std. Error"]))

	S<-sqrt(beta^2 * sa^2 + alpha^2*sb^2)
	Z<-taudiff/S
	
	Pvalue<-pnorm(-Z)
	res<-data.frame(yind = vector(), zind = vector(), value = vector())

	for(j in 1:nz){
		res.add<-c(yind = 1, zind = j, 
			taudiff = taudiff[j],
			tau0 = tau0[j],
			tau1 = tau1[j],
			tauratio = taudiff[j]/(sign(tau0[j])*tau0[j]),
			Z = Z[j],
			S = S[j],
			pvalue = Pvalue[j])
		res<-rbind(res, res.add)
	}
	colnames(res)<-c("yind", "zind", "taudiff", "tau0","tau1", "tauratio", "Z", "S", "pvalue")
	res[, "fdr"]<-p.adjust(res[, "pvalue"], method = "fdr")

	return(res)
}

calc_sobel_y1_z1<-function(x,y,z){
	y<-suppressWarnings(as.numeric(y))
	z<-suppressWarnings(as.numeric(z))
	m1<-lm(z ~ x)
	tau0<-m1$coefficients[2]
	m2<-lm(y ~ x)
	alpha<-m2$coefficients[2]
	m3<-lm(z ~ x + y)
	tau1<-m3$coefficients[2]
	beta<-m3$coefficients[3]
	taudiff<-sign(tau0)*(tau0-tau1)
	
	sa<-summary(m2)$coefficients[2, "Std. Error"]
	sb<-summary(m3)$coefficients[3, "Std. Error"]
		
	S<-sqrt(beta^2 * sa^2 + alpha^2*sb^2)
	Z<-taudiff/S
	pvalue<-pnorm(-Z)

	res<-data.frame(yind = 1, zind = 1, 
			taudiff = taudiff,
			tau0 = tau0,
			tau1 = tau1,
			tauratio = taudiff/(sign(tau0)*tau0),
			Z = Z,
			S = S,
			pvalue = pvalue)

	colnames(res)<-c("yind", "zind", "taudiff", "tau0", "tau1", "tauratio", "Z", "S", "pvalue")

	res[, "fdr"]<-p.adjust(res[, "pvalue"], method = "fdr")

	return(res)
}

calc_sobel_mat<-function(x,y,z){

	m1<-lm(t(z) ~ x)
	m2<-lm(t(y) ~ x)
	m3<-lapply(1:nrow(y), function(i){
		lm(t(z) ~ x+y[i, ])
		})

	m2.summary<-summary(m2)
	ny<-nrow(y)
	nz<-nrow(z)

	#size of z
	tau0<-m1$coefficients[2,]
	tau1<-lapply(1:ny, function(i) m3[[i]]$coefficients[2,])
	beta<-lapply(1:ny, function(i) m3[[i]]$coefficients[3,])
	alpha<-m2$coefficients[2,]
	taudiff<-lapply(1:ny, function(i) sign(tau0) * (tau0 - tau1[[i]]))
	sa<-unlist(lapply(1:ny, function(i) m2.summary[[i]]$coefficients[2, "Std. Error"]))

	sb<-lapply(1:ny, function(i) {
			m3_i<-summary(m3[[i]])
			return(unlist(lapply(1:nz, function(j) m3_i[[j]]$coefficients[3, "Std. Error"])))
		})

	S<-lapply(1:ny, function(i) {
			sqrt(beta[[i]]^2 * sa[i]^2 + alpha[i]^2*sb[[i]]^2)
		})

	Z<-lapply(1:ny, function(i) {
			taudiff[[i]]/S[[i]]
		})

	Pvalue<-lapply(1:ny, function(i) {
		pnorm(-Z[[i]])
		})

	res<-data.frame(yind = vector(), zind = vector(), value = vector())
	for(i in 1:ny){
		for(j in 1:nz){
			res.add<-c(yind = i, zind = j, 
				taudiff = taudiff[[i]][j],
				tau0 = tau0[j],
				tau1 = tau1[[i]][j],
				tauratio = taudiff[[i]][j]/(sign(tau0[j])*tau0[j]),
				Z = Z[[i]][j],
				S = S[[i]][j],
				pvalue = Pvalue[[i]][j])
			res<-rbind(res, res.add)
		}
	}

	colnames(res)<-c("yind", "zind", "taudiff", "tau0", "tau1", "tauratio", "Z", "S", "pvalue")
	res[, "fdr"]<-p.adjust(res[, "pvalue"], method = "fdr")

	return(res)
}

#' @import plyr
subset_group_min<-function(res, by = "trans", metric = "pvalue"){
	res[, "metric"]<-res[, metric]
	res<-ddply(res, by, subset, metric == min(metric))
	return(res[, setdiff(colnames(res), "metric")])
}


#' @import Biobase
prune<-function(f_cis_tab, f_trans_tab, 
	cn, gep,
	alteration_id = "Unique.Name",
	gene_id = "accession", 
	pruning_dir,
	prunecol = "pvalue", prunethres = 0.25,
	gs, ... #other args.file in run_hyperEnrichment_pruned
	) {
	

	res.actual<-list()
	res.sig<-list()

	if(nrow(f_cis_tab) == 0){
		warning("No significant cis genes")
		return(list(all = res.actual, sig = res.sig))
	}
	if(nrow(f_trans_tab) == 0){
		warning("No significant trans genes")
		return(list(all = res.actual, sig = res.sig))
	}

	alt_id<-unique(as.character(f_cis_tab[, alteration_id]))
	cn.fdat<-fData(cn)
	cn.exprs<-exprs(cn)
	ge.fdat<-fData(gep)
	ge.fdat.genes<-as.character(ge.fdat[, gene_id])
	ngenes<-length(ge.fdat.genes)
	ge.exprs<-exprs(gep)

	#set rownames of exprs to gene symbols
	rownames(ge.exprs)<-ge.fdat.genes

	pruning_dir_tables<-paste(pruning_dir, "/tables", sep = "")
	pruning_dir_js<-paste(pruning_dir, "/js", sep = "")

	suppressWarnings(dir.create(pruning_dir_tables, recursive = TRUE))
	suppressWarnings(dir.create(pruning_dir_js, recursive = TRUE))


	hyper<-list()

	#create pruning/hyperEnrichment subdirectories
	suppressWarnings(if(!is.na(gs)){
		pruning_dir_hyper<-paste(pruning_dir, "/hyperEnrichment", sep = "")
		suppressWarnings(dir.create(pruning_dir_hyper, recursive = TRUE))
		
		for(j in names(gs)){
			pruning_dir_hyper_gs<-paste(pruning_dir_hyper, "/", j, sep = "")
			suppressWarnings(dir.create(pruning_dir_hyper_gs))
		}
	})

	cis_summary<-data.frame(c())
	
	for(i in alt_id){

		cat(paste("Running pruning for ", i, "\n"))

		hyper[[i]]<-list()

		i_ind <-which(cn.fdat[, alteration_id] == i)
		print(i_ind)
		x <- suppressWarnings(as.numeric(cn.exprs[i_ind,]))
		cis_genes<-as.character(f_cis_tab[which(f_cis_tab[,alteration_id] == i), gene_id])
		trans_genes<-as.character(f_trans_tab[which(f_trans_tab[,alteration_id] == i), gene_id])
		cis_genes<-intersect(cis_genes, ge.fdat.genes)
		trans_genes<-intersect(trans_genes, ge.fdat.genes)
		cis_genes_n<-length(cis_genes)
		trans_genes_n<-length(trans_genes)

		print("number cis")
		print(cis_genes_n)
		print("number trans")
		print(trans_genes_n)

		if(cis_genes_n >0 & trans_genes_n > 0) {
	
			cis_vec<-ge.exprs[match(cis_genes, ge.fdat.genes),, drop = FALSE]
			trans_vec<-ge.exprs[match(trans_genes, ge.fdat.genes),, drop = FALSE]
			
			res.actual[[i]]<-calc_sobel(x =x, y = cis_vec,z = trans_vec, 
				y.names = rownames(cis_vec), z.names = rownames(trans_vec))

			tab<-res.actual[[i]]
		
			tab<-subset_group_min(tab, by = "trans", metric = "fdr")
			res.sig[[i]]<-tab[tab[, prunecol] < prunethres,]

			suppressWarnings(if(!is.na(gs)){
				
				for(j in names(gs)){ #iterate through geneset compendiums
					gs.curr<-gs[[j]]
					pruning_dir_hyper_gs<-paste(pruning_dir_hyper, "/", j, sep = "")
					cat("Running hyperenrichment...\n")
					hyper[[i]][[j]]<-run_hyperEnrichment_pruned(tab = res.sig[[i]], 
						tab.name = i, gs= gs.curr, ngenes = ngenes,
						f.dir.out = pruning_dir_hyper_gs, ...)
				}

				write_bipartite_JSON(tab = res.sig[[i]], 
						hyper = hyper[[i]], f.dir.out = pruning_dir_js, header = i)

			} 

			else 
				write_bipartite_JSON(tab = res.sig[[i]], 
						f.dir.out = pruning_dir_js, header = i)
			)
			
			#summary of cis genes by prioritization
			for(cis in unique(res.actual[[i]][, "cis"])){
				tot_trans<-sum(res.actual[[i]][, "cis"] %in% cis)
				if(nrow(res.sig[[i]]) == 0){
					mediated_trans<-0
					mediated_trans_weighted<-0
				}
				else {
					mediated_trans<-sum(res.sig[[i]][, "cis"] %in% cis)
					if(mediated_trans>0){
						weights<-res.sig[[i]][res.sig[[i]][, "cis"] %in% cis, "tauratio"]
						weights[weights > 1]<-1
						weights[weights < 0]<-0
						mediated_trans_weighted<-sum(weights)
					} 
					else 
						mediated_trans_weighted<-0
				}

				num_path<-NA
				paths<-NA
				num_path_header<-"num_path"
				paths_header<-"paths"
				names(num_path)<-num_path_header
				names(paths_header)<-paths_header
				
				suppressWarnings(if(!is.na(gs)){
					cis.path<-lapply(names(gs), function(j){
						num_pathway_mediated_trans<-0
						pathway_mediated_trans<-""

						if(cis %in% names(hyper[[i]][[j]][["hyperbyalt"]])){
							hyper_cis<-hyper[[i]][[j]][["hyperbyalt"]][[cis]]
							if(!is.null(hyper_cis))
								num_pathway_mediated_trans<-nrow(hyper_cis)
								if(num_pathway_mediated_trans > 0)
									pathway_mediated_trans<-paste(as.character(hyper_cis[, "category"]),
									 collapse = ",")
						}
						names(num_pathway_mediated_trans)<-paste("numpath_", j, sep = "")
						names(pathway_mediated_trans)<-paste("paths_", j, sep = "")
						return(list(num_pathway_mediated_trans, pathway_mediated_trans))

					})

					num_path<-lapply(cis.path, function(i) i[[1]])
					names(num_path)<-paste(num_path_header, "_", names(gs), sep  = "")
					paths<-lapply(cis.path, function(i) i[[2]])
					names(paths)<-paste(paths_header, "_", names(gs), sep = "")

				})

				cis_summary.add<-data.frame(alteration = i, cis = cis,  
					mediated_trans = mediated_trans, 
					mediated_trans_weighted = mediated_trans_weighted,
					total_trans = tot_trans,
					frac_mediated_trans = mediated_trans/tot_trans, 
					frac_mediated_trans_weighted = mediated_trans_weighted/tot_trans,
					num_path, 
					paths)

				cis_summary<-rbind(cis_summary, cis_summary.add)
			}
			cat("\n")
		}

	}

	
	res.actual.alt<-lapply(names(res.actual), function(i){
		cbind(data.frame(alt = i, res.actual[[i]]))
		})
	names(res.actual.alt)<-names(res.actual)

	res.sig.alt<-lapply(names(res.sig), function(i){
		cbind(data.frame(alt = rep(i, nrow(res.sig[[i]])), res.sig[[i]]))
		})
	names(res.sig.alt)<-names(res.sig)

	res.combined.all<-Reduce(rbind, res.actual.alt)

	write.table(res.combined.all, 
				file = paste(pruning_dir_tables, "/pruned.all.txt", sep = ""),
				col.names = TRUE, row.names = FALSE, sep = "\t")

	res.combined.sig<-Reduce(rbind, res.sig.alt)

	write.table(res.combined.sig, 
				file = paste(pruning_dir_tables, "/pruned.sig.txt", sep = ""),
				col.names = TRUE, row.names = FALSE, sep = "\t")

	if(nrow(cis_summary) > 0){
	cis_summary<-cis_summary[order(cis_summary[,"alteration"],
			-cis_summary[,"frac_mediated_trans_weighted"],decreasing=FALSE),]

	a<-cis_summary[,"alteration"]
	tots<-sapply(unique(a), function(x) sum(a == x))
	ranks<-unlist(sapply(tots, function (x) 1:x))

	cis_summary_full<-cbind(cis_summary, rank = ranks)
	rownames(cis_summary_full)<-NULL
	write.table(cis_summary_full, 
		file = paste(pruning_dir_tables, "/cis_summary.txt", sep = ""),
		col.names = TRUE, row.names = FALSE, sep = "\t")
	
	suppressWarnings(if(!is.na(gs))
		return(list(cis_summary = cis_summary_full, all = res.actual, sig = res.sig, hyper = hyper))
	else 
		return(list(cis_summary = cis_summary_full, all = res.actual, sig = res.sig))
		)
	} else {
		print("no pruned results!")
		return(list(all = res.actual, sig = res.sig))
	}

}

#' @import Biobase
to.eSet<-function(mat, pdat, fdat){
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
