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
	fit2.table<-topTable(fit2, coef=1, adjust.method="BH", number =length(fit2) ,
		sort.by = "none")

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
	fit2.table<-join(x=fit2.table, y=fit2.accessary, by = colnames(fData(eset)), type = "left", match = "first")
	fit2.table$fold.change<-2^fit2.table$logFC
	fit2.table<-fit2.table[order(fit2.table$adj.P.Val, decreasing = FALSE),]
	return(fit2.table)
}

#' iEDGE_combine combines alteration data (e.g. GISTIC2) and gene expression data, using a mapping data frame for matching samples
#' @export
iEDGE_combine<-function(alt, gep, mapping, mapping.cn = "CN", mapping.gep = "GEP", uppercase = TRUE){

	cn <- alt$cn
	cis <- alt$cis

	samples.gep<-colnames(gep)
	samples.cn<-colnames(cn)

	if(uppercase){
		cn_id<-match(toupper(mapping[, mapping.cn]), toupper(samples.cn))
		gep_id<-match(toupper(mapping[, mapping.gep]), toupper(samples.gep))
	} else {
		cn_id<-match(mapping[, mapping.cn], samples.cn)
		gep_id<-match(mapping[, mapping.gep], samples.gep)
	}
	mapping.sub<-data.frame(cn_id, gep_id)
	mapping.sub<-mapping.sub[!apply(mapping.sub, 1, function(i){
		any(is.na(i))
		}),]
	dat<-list()
	dat$cn<-cn[, mapping.sub$cn_id]
	dat$gep<-gep[, mapping.sub$gep_id]
	dat$cisgenes<-cis
	return(dat)
}

#' @import Biobase
#' @export
make_iEDGE<-function(gep, #eset containing log2 gene expression
 cn, #eset containing copy number by sample, must be sorted with respect to columns in gep
 cisgenes, #list of cisgenes in in gep with respect to cn ordered by rows in cn
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
 fc = NA
){

	if(onesided & !(cndir %in% colnames(fData(cn)))){
		stop ("If onesided is specified, cndir must be in colnames(fData(cn))")
		if(any(!(fData(cn)[, cndir] %in% c(uptest, downtest)))){
			stop("fData(cn)[, cndir] must contain elements from uptest or downtest")
		}
	}
	all_genes<-as.character(fData(gep)[, gepid])
	cat(paste("Running iEDGE differential expression for ", nrow(cn), " alterations: \n\n", sep = ""))

	res<-list()
	pb <- txtProgressBar(style = 3)
	setTxtProgressBar(pb, 0)

	n<-nrow(cn)
	for (i in 1:n){ 	#for each alteration
		setTxtProgressBar(pb, i/n)
		if (interaction_type == "cis")
			genes.keep<-intersect(all_genes, unique(cisgenes[[i]]))
		else 
			genes.keep<-setdiff(all_genes, unique(cisgenes[[i]]))

		if (length(genes.keep)<1) res[[i]]<-NA
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
			treatment<-factor(as.numeric(exprs(cn)[i,]))
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
				res[[i]]<-NA
			}
			else 
				res[[i]]<-run_limma(eset = gep.keep, design = design, 
					cndir = cndir, onesided = onesided, 
					uptest = uptest, downtest = downtest)
		}
	}

	setTxtProgressBar(pb, 1)
	close(pb)
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
		res.list$sig <- res.df[res.df[, "adj.P.Val.all"]<fdr.cutoff & (res.df[, "fold.change"] > fc | res.df[, "fold.change"] < fc.low),]
	}
#	res.list$sig<-subset(res.df, adj.P.Val.all < fdr.cutoff)
	return(res.list)
}

#' iEDGE_DE performs differential expression analysis and pathway enrichment and saves to text tables
#' @export
iEDGE_DE<-function(cn, gep, cisgenes,
	header,
	gepid, cnid,  	
	f.dir.out, 
	gs.file, #genset for hyper
	gs.file.name,
	fdr.cis.cutoff = 0.25, fdr.trans.cutoff = 0.05, 
	min.group = 3,
	min.drawsize = 3,  
	cndir = "alteration_direction",
	cis.onesided = TRUE, 
	trans.onesided = FALSE,
	fc.cis = NA,
	fc.trans = NA,
	... #other parameters in make_iEDGE
	){
	dir.create(f.dir.out, recursive =TRUE)
	cat("Performancing cis analysis...\n")
	res.cis<-make_iEDGE(gep = gep, 
 		cn = cn,
 		cisgenes = cisgenes, 
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
	f.out<-paste(f.dir.out, "/", header, "_cis_sig_fdr_", fdr.cis.cutoff, fc_header, ".txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(res.cis$sig, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)

	f.out<-paste(f.dir.out, "/", header, "_cis_full.txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(res.cis$full, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)

	cat("Performancing trans analysis...\n")
	res.trans<-make_iEDGE(gep = gep, 
 		cn = cn,
 		cisgenes = cisgenes, 
 		gepid = gepid, 
 		cnid = cnid,
 		cndir = cndir,
 		interaction_type = "trans",
 		fdr.cutoff = fdr.trans.cutoff, 
 		onesided = trans.onesided,
 		min.group = min.group,
 		fc = fc.trans,
 		...)

	if(is.na(fc.trans)) fc_trans_header <- ""
	else fc_trans_header<-paste("_fc_", fc.trans, sep = "")

	f.out<-paste(f.dir.out, "/", header, "_trans_sig_fdr_", fdr.trans.cutoff, fc_trans_header, ".txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(res.trans$sig, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)
	
	cat("Done!\n")

	cat("\nPerforming hyperEnrichment analysis...\n")

 	ngenes<-nrow(gep)

 	#run hyperenrichment with all sig cis genes in one geneset
 	drawnList<-list(cis.sig = unique(as.character(res.cis$sig[, gepid])))
 	hyper.cis<-run_hyperEnrichment(drawn=drawnList,
	    categories=gs.file,
	    ntotal=ngenes,
	    min.drawsize = min.drawsize, mht = TRUE, verbose = TRUE, order = TRUE)
 	f.out<-paste(f.dir.out, "/", header, "_hyperEnrichment_",gs.file.name,"_cis.txt", sep = "")	
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(hyper.cis, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)
	cat("Done!\n")

	#run hyperenrichment with all sig trans genes in one geneset
	drawnList<-list(trans.sig = unique(as.character(res.trans$sig[, gepid])))
 	hyper.trans<-run_hyperEnrichment(drawn=drawnList,
	    categories=gs.file,
	    ntotal=ngenes,
	    min.drawsize = min.drawsize, mht = TRUE, verbose = TRUE, order = TRUE)
 	f.out<-paste(f.dir.out, "/", header, "_hyperEnrichment_",gs.file.name,"_trans.txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(hyper.trans, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)	
	cat("Done!\n")

	#run hyperenrichment with lists of sig cis genes separated by the alteration it is contained in
	drawnList<-sapply(unique(res.cis$sig[, cnid]), 
		function(x){
			y<-res.cis$sig[res.cis$sig[, cnid] == x,]
			#y<-subset(res.cis$sig, Unique.Name == x)
			return(unique(as.character(y[, gepid])))
			})
	
	names(drawnList)<-unique(res.cis$sig[, cnid])
 	hyper.cis.split<-run_hyperEnrichment(drawn=drawnList,
	    categories=gs.file,
	    ntotal=ngenes,
	    min.drawsize = min.drawsize, mht = TRUE, verbose = TRUE, order = TRUE)

 	f.out<-paste(f.dir.out, "/", header, "_hyperEnrichment_",gs.file.name,"_cis_splitbyalteration.txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(hyper.cis.split, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)		
	cat("Done!\n")

	#run hyperenrichment with lists of sig trans genes separated by the alteration it is contained in
	drawnList<-sapply(unique(res.trans$sig[, cnid]), 
		function(x){
			#y<-subset(res.trans$sig, Unique.Name == x)
			y<-res.trans$sig[res.trans$sig[, cnid] == x,]
			return(unique(as.character(y[, gepid])))
			})	
	names(drawnList)<-unique(res.trans$sig[, cnid])

 	hyper.trans.split<-run_hyperEnrichment(drawn=drawnList,
	    categories=gs.file,
	    ntotal=ngenes,
	    min.drawsize = min.drawsize, mht = TRUE, verbose = TRUE, order = TRUE)

 	f.out<-paste(f.dir.out, "/", header, "_hyperEnrichment_",gs.file.name,"_trans_splitbyalteration.txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(hyper.trans.split, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)

	return(list(cis=res.cis, trans =res.trans, 
		hyper.cis=hyper.cis, hyper.cis.split=hyper.cis.split, 
		hyper.trans=hyper.trans, hyper.trans.split=hyper.trans.split))
}

#' @import infotheo
#' @export
calc_mutinfo<-function(x, #alt
	y,#cis 
	z,#trans
	nbins = 5){
	y<-discretize(y, nbins=nbins)
	z<-discretize(z, nbins=nbins)
	I <- condinformation(x,discretize(z),discretize(y),method="emp")
	return(I)
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
	taudiff<-lapply(1:ny, function(i) tau0 - tau1[[i]])
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

	P.twosided<-lapply(1:ny, function(i) {
			2*pnorm(-abs(Z[[i]]))
		})

	res<-data.frame(yind = vector(), zind = vector(), value = vector())
	for(i in 1:ny){
		res.add<-c(yind = i, zind = j, 
			taudiff = taudiff[[i]][j],
			Z = Z[[i]][j],
			S = S[[i]][j],
			pvalue = P.twosided[[i]][j])
		res<-rbind(res, res.add)
	}

	colnames(res)<-c("yind", "zind", "taudiff", "Z", "S", "pvalue")
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

	taudiff<-tau0 - tau1
	sa<-m2.summary$coefficients[2, "Std. Error"]

	m3.summary<-summary(m3)
	sb<-unlist(lapply(1:nz, function(j) m3.summary[[j]]$coefficients[3, "Std. Error"]))

	S<-sqrt(beta^2 * sa^2 + alpha^2*sb^2)
	Z<-taudiff/S
	
	P.twosided<-2*pnorm(-abs(Z))
	res<-data.frame(yind = vector(), zind = vector(), value = vector())

	for(j in 1:nz){
		res.add<-c(yind = 1, zind = j, 
			taudiff = taudiff[j],
			Z = Z[j],
			S = S[j],
			pvalue = P.twosided[j])
		res<-rbind(res, res.add)
	}
	colnames(res)<-c("yind", "zind", "taudiff", "Z", "S", "pvalue")
	return(res)
}

calc_sobel_y1_z1<-function(x,y,z){
	y<-as.numeric(y)
	z<-as.numeric(z)
	m1<-lm(z ~ x)
	tau0<-m1$coefficients[2]
	m2<-lm(y ~ x)
	alpha<-m2$coefficients[2]
	m3<-lm(z ~ x + y)
	tau1<-m3$coefficients[2]
	beta<-m3$coefficients[3]
	taudiff<-tau0-tau1
	
	sa<-summary(m2)$coefficients[2, "Std. Error"]
	sb<-summary(m3)$coefficients[3, "Std. Error"]
		
	S<-sqrt(beta^2 * sa^2 + alpha^2*sb^2)
	Z<-taudiff/S
	pvalue<-2*pnorm(-abs(Z))

	res<-data.frame(yind = 1, zind = 1, 
			taudiff = taudiff,
			Z = Z,
			S = S,
			pvalue = pvalue)

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

	tau0<-m1$coefficients[2,]

	tau1<-lapply(1:ny, function(i) m3[[i]]$coefficients[2,])
	beta<-lapply(1:ny, function(i) m3[[i]]$coefficients[3,])
	alpha<-m2$coefficients[2,]
	taudiff<-lapply(1:ny, function(i) tau0 - tau1[[i]])
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

	P.twosided<-lapply(1:ny, function(i) {
			2*pnorm(-abs(Z[[i]]))
		})

	res<-data.frame(yind = vector(), zind = vector(), value = vector())
	for(i in 1:ny){
		for(j in 1:nz){
			res.add<-c(yind = i, zind = j, 
				taudiff = taudiff[[i]][j],
				Z = Z[[i]][j],
				S = S[[i]][j],
				pvalue = P.twosided[[i]][j])
			res<-rbind(res, res.add)
		}
	}

	colnames(res)<-c("yind", "zind", "taudiff", "Z", "S", "pvalue")
	return(res)
}

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

run_cmi_hyperenrichment<-function(tab, tab.name, gs.file, ngenes, 
	min.drawsize = 3, 
	hypercol = "fdr", 
	hyperthres = 0.25, 
	f.dir.out){

	drawnList<-list(drawn = unique(tab$trans))
	hyperunion<-run_hyperEnrichment(drawn=drawnList,
		categories=gs.file,
		ntotal=ngenes,
		min.drawsize = min.drawsize, 
		mht = TRUE, 
		verbose = TRUE, 
		order = TRUE, 
		hypercol = hypercol, hyperthres = hyperthres)
	
	if(nrow(hyperunion)>0)
	write.table(hyperunion, file = paste(f.dir.out, "/", tab.name, ".txt", sep = ""), 
		sep = "\t", row.names = FALSE)
		
	cis<-unique(tab$cis)
	hyperbyalt<-lapply(cis, function(j){
			i.sub<-tab[tab$cis %in% j,]
			drawnList <- list(drawn = unique(i.sub$trans))
			hyper<-run_hyperEnrichment(drawn=drawnList,
	    		categories=gs.file,
	    		ntotal=ngenes,
	    		min.drawsize = min.drawsize, 
	    		mht = TRUE, 
	    		verbose = TRUE, 
	    		order = TRUE,
	    		hypercol = hypercol, hyperthres = hyperthres)
			
			if(nrow(hyper)>0)
			write.table(hyper, file = paste(f.dir.out, "/", tab.name, "_", j, ".txt", sep = ""), 
						sep = "\t", row.names = FALSE)
			return(hyper)
			})
	names(hyperbyalt)<-cis

	return(list(hyperunion = hyperunion, hyperbyalt = hyperbyalt))
}


#' @import Biobase
#' @export
prune<-function(f_cis_tab, 
	f_trans_tab, 
	cn, 
	gep,
	alteration_id = "Unique.Name",
	gene_id = "accession",
	seed =7,
	nsamples = 1000, 
	nbins = 5,  
	cmi_dir,
	prunecol = "pvalue", prunethres = 0.25,
	method = "cmi", #cmi or sobel
	... #other args.file in run_cmi_hyperenrichment
	) {
	
	if(!(method %in% c("cmi", "sobel"))){
		stop("method must be one of \"cmi\" or \"sobel\"")
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

	res.actual<-list()
	res.null<-list()
	res.sig<-list()

	cmi_dir_tables<-paste(cmi_dir, "/tables", sep = "")
	cmi_dir_plots<-paste(cmi_dir, "/plots", sep = "")
	cmi_dir_js<-paste(cmi_dir, "/js", sep = "")

	dir.create(cmi_dir_tables, recursive = TRUE)
	dir.create(cmi_dir_plots, recursive = TRUE)

	if(hasArg("gs.file")){
		cmi_dir_hyper<-paste(cmi_dir, "/hyperEnrichment", sep = "")
		dir.create(cmi_dir_hyper, recursive = TRUE)
		hyper<-list()
	}

	p<-list()

	for(i in alt_id){

		cat(paste("alteration: ", i, "\n", sep = " "))	
		i_ind <-which(cn.fdat[, alteration_id] == i)
		x <- as.numeric(cn.exprs[i_ind,])

		cis_genes<-as.character(f_cis_tab[which(f_cis_tab[,alteration_id] == i),gene_id])
		trans_genes<-as.character(f_trans_tab[which(f_trans_tab[,alteration_id] == i),gene_id])
		cis_genes<-intersect(cis_genes, ge.fdat.genes)
		trans_genes<-intersect(trans_genes, ge.fdat.genes)

		cis_genes_n<-length(cis_genes)
		trans_genes_n<-length(trans_genes)


		cat(paste("cis genes: ", cis_genes_n, ", trans genes: ", trans_genes_n, "\n", sep = ""))	
		if (cis_genes_n ==0 | trans_genes_n == 0){

		}
		else {


		cis_vec<-ge.exprs[match(cis_genes, ge.fdat.genes),, drop = FALSE]
		trans_vec<-ge.exprs[match(trans_genes, ge.fdat.genes),, drop = FALSE]

		if(method == "cmi"){

			set.seed(seed)
			sample_cg<-ge.fdat.genes[sample(1:ngenes, nsamples, replace = T)]
			set.seed(seed+10)
			sample_tg<-sample(trans_genes, min(nsamples, trans_genes_n), replace = T)
		
			cis_null_vec<-ge.exprs[match(sample_cg, ge.fdat.genes),, drop = FALSE]
			trans_null_vec<-ge.exprs[match(sample_tg, ge.fdat.genes),, drop = FALSE]

			cat("Running actual model...\n")
			res.actual[[i]]<-lapply(1:nrow(cis_vec),
					function(k){
						y<- cis_vec[k,]
						y.name<-rownames(cis_vec)[k]
						res<-sapply(1:nrow(trans_vec), function(j)
								calc_mutinfo( x = x,y = y, z = trans_vec[j,], 
									nbins = nbins))

						res.df<-data.frame(cis = y.name, trans = rownames(trans_vec), cmi = as.numeric(res))
						return(res.df)	
					})

			res.actual[[i]]<-do.call(rbind, res.actual[[i]])

			cat("Running null model...\n")
			res.null[[i]]<-lapply(1:nrow(cis_null_vec),
					function(k){
						y<- cis_null_vec[k,]
						y.name<-rownames(cis_null_vec)[k]
						res<-sapply(1:nrow(trans_null_vec), function(j)
								calc_mutinfo( x = x,y = y, z = trans_null_vec[j,], 
									nbins = nbins))
						res.df<-data.frame(cis = y.name, trans = rownames(trans_null_vec), cmi = as.numeric(res))
						return(res.df)	
					})

			res.null[[i]]<-do.call(rbind, res.null[[i]])	

			write.table(res.null[[i]], file = paste(cmi_dir_tables, "/null_", i, ".txt", sep = ""),
				col.names = TRUE, row.names = FALSE, sep = "\t")
			
			nullecdf<-ecdf(res.null[[i]]$cmi)
			res.actual[[i]]$pvalue<-nullecdf(res.actual[[i]]$cmi)

			write.table(res.actual[[i]], file = paste(cmi_dir_tables, "/actual_", i, ".txt", sep = ""),
				col.names = TRUE, row.names = FALSE, sep = "\t")

			tab<-res.actual[[i]]
			res.sig[[i]]<-tab[tab[, prunecol] < prunethres,]
			write.table(res.sig[[i]], file = paste(cmi_dir_tables, "/sig_", i, ".txt", sep = ""),
				col.names = TRUE, row.names = FALSE, sep = "\t")

		} else { #method = sobel

			cat("Running actual model...\n")

			res.actual[[i]]<-calc_sobel(x =x,y = cis_vec,z = trans_vec, 
				y.names = rownames(cis_vec), z.names = rownames(trans_vec))

			write.table(res.actual[[i]], file = paste(cmi_dir_tables, "/actual_", i, ".txt", sep = ""),
				col.names = TRUE, row.names = FALSE, sep = "\t")

			tab<-res.actual[[i]]
			res.sig[[i]]<-tab[tab[, prunecol] < prunethres,]
			write.table(res.sig[[i]], file = paste(cmi_dir_tables, "/sig_", i, ".txt", sep = ""),
				col.names = TRUE, row.names = FALSE, sep = "\t")

		}

		if(hasArg("gs.file")){
			cat("Running hyperenrichment...\n")
			hyper[[i]]<-run_cmi_hyperenrichment(tab = res.sig[[i]], tab.name = i, ngenes = ngenes,
			f.dir.out = cmi_dir_hyper, ...)

			cat("Writing cmi js file..\n")
			write_bipartite_JSON(tab = res.sig[[i]], 
				hyper = hyper[[i]], f.dir.out = cmi_dir_js, header = i)

		} else {
			write_bipartite_JSON(tab = res.sig[[i]], 
				f.dir.out = cmi_dir_js, header = i)
		}
		cat("\n")
		}

	}

	if(hasArg("gs.file"))
	return(list(actual = res.actual, null = res.null, sig = res.sig, p = p, hyper = hyper))
	else 
	return(list(actual = res.actual, null = res.null, sig = res.sig, p = p))

}

#' run_iEDGE is the main wrapper for iEDGE, assumes processed data
#' @param iEDGE datalist consisting of cn (alteration), gep (gene expression), and cisgenes (list of genes in each alt)
#' @param outdir output directory
#' @export
run_iEDGE<-function(dat, header, outdir, gs.file = NA, gepid = "SYMBOL", cnid = "Unique.Name", cndir = "alteration_direction",
	fdr.cis.cutoff = 0.25, fdr.trans.cutoff = 0.05, fc.cis = NA, fc.trans = NA, min.drawsize = 3, onesided.cis = TRUE, 
	onesided.trans = FALSE, uptest = "Amplification", downtest = "Deletion", gs.file.name = "c2.cp.v5.0",
	min.group = 2, mutinfo.seed = 7, mutinfo.nsamples = 500, mutinfo.bins = 5,  
	prune.method = "cmi", prune.col = "pvalue", prune.thres = 0.05, hyperthres = 0.25, 
	cis.boxplot = TRUE, trans.boxplot = TRUE, bipartite = TRUE){

	if(is.na(gs.file)){ #use default gmt file c2 cp
		HOMEDIR<-file.path(path.package("iEDGE"), "genesets")
		gs.file<-paste(HOMEDIR, "/c2.cp.v5.0.symbols.gmt", sep =  "")
	}

	cat("Reading genesets..\n")
	c2<-read_gmt(gs.file)$genesets
	c2<-sapply(c2, function(x){
		a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
		return(unique(do.call(c,a)))
	})


	cn<-dat$cn
	gep<-dat$gep
	cisgenes<-dat$cisgenes
	cat(paste("Running iEDGE for data set: ", header, "\n",sep = ""))


	dir.create(outdir, recursive = TRUE)
	base_dir<-paste(outdir,"/", header, sep = "")

	dir.create(base_dir, recursive = TRUE)
	de_dir<-paste(base_dir, "/tables", sep = "")

	cat("Making Differential Expression tables...\n")
	res<-iEDGE_DE(cn, gep, cisgenes,
		header,
		gepid, cnid,  	
		f.dir.out = de_dir, 
		gs.file = c2, #geneset for hyper
		gs.file.name = gs.file.name,
		fdr.cis.cutoff = fdr.cis.cutoff, fdr.trans.cutoff = fdr.trans.cutoff, 
		min.group = min.group,
		min.drawsize = min.drawsize,  
		cis.onesided = onesided.cis, 
		trans.onesided = onesided.trans,
		fc.cis = fc.cis,
		fc.trans = fc.trans,
		uptest = "Amplification",
		downtest = "Deletion"
		)


	cmi_dir<-paste(base_dir, "/cmi", sep = "")
	res.cis.sig<-res[["cis"]][["sig"]]
	res.cis.full<-res[["cis"]][["full"]]
	res.trans.sig<-res[["trans"]][["sig"]]

	if(bipartite){
		cmi<-prune(f_cis_tab =  res.cis.sig, 
			f_trans_tab= res.trans.sig, 
			cn = cn, 
			gep = gep,
			alteration_id = cnid,
			gene_id = gepid,
			seed =mutinfo.seed,
			nsamples = mutinfo.nsamples, 
			cmi_dir = cmi_dir, 
			nbins = mutinfo.bins,
			gs.file = c2,
			prunecol = prune.col, prunethres = prune.thres, 
			method = prune.method,
			min.drawsize = min.drawsize, 
			hypercol = "fdr", 
			hyperthres = hyperthres)
	} else {
		cmi <-NA
	}

	html_dir<-paste(base_dir, "/html", sep = "")
	jsdir<-file.path(path.package("iEDGE"), "javascript")

	iEDGE_UI(cistab = res.cis.sig, cisfulltab = res.cis.full,
		transtab = res.trans.sig, cn = cn, gep = gep, cisgenes = cisgenes,
		outdir = html_dir, jsdir = jsdir, cmi = cmi, cmijsdir = paste(cmi_dir, "/js", sep = ""),
		altid = cnid, geneid = gepid, 
		cis.boxplot = cis.boxplot, trans.boxplot = trans.boxplot, bipartite = bipartite)
}
