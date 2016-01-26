require(plyr)
require(Biobase)

#run limma for specified expression and design matrix
run_limma<-function(eset, design, cndir = "alteration_direction"){
	require(limma)
	fit <- lmFit(eset, design)
	command_str<-paste("makeContrasts(",
	"(", colnames(design)[2] , "-", colnames(design)[1], ")", 
	 ",levels = design)", sep = "")
	contrast.matrix<-eval(parse(text =command_str)) 
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	fit2.table<-topTable(fit2, coef=1, adjust="BH", number =length(fit2) ,
		sort.by = "none")
	if (fit2.table$interaction_type[1] == "cis"){
		p.val.onesided <- sapply(1:nrow(fit2.table), 
			function(i){
				if (fit2.table[, cndir][i] == "Amplification")
					lower.tail<-fit2$t[i] < 0
				else #deletion
					lower.tail <-fit2$t[i] >= 0
				return(pt(q = abs(fit2$t[i]), df = fit2$df.total[i], lower.tail = lower.tail))
			})
		p.val.onesided.adjust<-p.adjust(p.val.onesided, method = "fdr")

		fit2.table$P.Value <- p.val.onesided 
		fit2.table$adj.P.Val <-p.val.onesided.adjust
	}
	fit2.accessary<-run_limma_accessory(eset, design)
	fit2.table<-join(x=fit2.table, y=fit2.accessary, by = colnames(fData(eset)), type = "left",
	match = "first")
	fit2.table$fold.change<-2^fit2.table$logFC
	fit2.table<-fit2.table[order(fit2.table$adj.P.Val, decreasing = FALSE),]
	return(fit2.table)
}

#append additional columns to limma table
run_limma_accessory<-function(eset, design){

	eset.sub.cond1<-eset[,design[,1] == 1]
	eset.sub.cond2<-eset[,design[,2] == 1]

	exprs.cond1<-exprs(eset.sub.cond1)
	exprs.cond2<-exprs(eset.sub.cond2)

	mean.cond1<-apply(exprs.cond1,1, mean)
	sd.cond1<-apply(exprs.cond1,1, sd)
	n.cond1<-rep(ncol(exprs.cond1), nrow(exprs.cond1))
	
	mean.cond2<-apply(exprs.cond2,1, mean)
	sd.cond2<-apply(exprs.cond2,1, sd)
	n.cond2<-rep(ncol(exprs.cond2), nrow(exprs.cond2))

	eset.unlog<-eset
	exprs(eset.unlog)<-2^exprs(eset.unlog) -1

	eset.sub.cond1<-eset.unlog[,design[,1] == 1]
	eset.sub.cond2<-eset.unlog[,design[,2] == 1]

	exprs.cond1.unlog<-exprs(eset.sub.cond1)
	exprs.cond2.unlog<-exprs(eset.sub.cond2)

	mean.cond1.unlog<-apply(exprs.cond1.unlog,1, mean)
	sd.cond1.unlog<-apply(exprs.cond1.unlog,1, sd)
	
	mean.cond2.unlog<-apply(exprs.cond2.unlog,1, mean)
	sd.cond2.unlog<-apply(exprs.cond2.unlog,1, sd)

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

#run iEDGE differential expression for specified gene expression, 
#copy number, cis genes, interaction type, and fdr cutoff
make_iEDGE<-function(gep, #eset containing log2 gene expression
 cn, #eset containing copy number by sample, must be sorted with respect to columns in gep
 cisgenes, #list of cisgenes in in gep with respect to cn ordered by rows in cn
 gepid = "accession", #column name in fData of gep identifying the genes in cisgenes
 cnid = "Unique.Name", #column name in fData of cn identifying the alterations in names(cn)
 cndir, #column name in fData of cn identifying direction of test 
 #(Amplification = high expression in cn = 1 low in 0, Deletion = low in 1 high in 0)
 interaction_type = "cis", #cis or trans
 fdr.cutoff = 0.25
){

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
			
			fdat.add<-fData(cn)[i, ]
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
			if (any(apply(design[, c(1,2)], 2, sum) < 3)){
				cat(paste("less than three samples in case or control for ", 
					alteration_ids[i], ", cannot run limma", "\n\n", sep = "")
				)
				res[[i]]<-NA
			}
			else 
				res[[i]]<-run_limma(eset = gep.keep, design = design, cndir = cndir)
		}
	}

	setTxtProgressBar(pb, 1)
	close(pb)

	res.df<-do.call(rbind, res[!is.na(res)])
	res.df$adj.P.Val.all<-p.adjust(res.df$P.Value, method = "fdr")
	res.list<-list()
	res.list$full<-res.df
	res.list$sig<-subset(res.df, adj.P.Val.all < fdr.cutoff)
	return(res.list)
}


##running iEDGE on cancer cell 2012
run_lymph<-function(){


 	lymph<-readRDS("../data/lymph2012/cancercell2012_all.RDS")

	f.dir.out <- "/Users/amyli/Desktop/git_projects/iEDGE/test2/tables"
	header <- "cancercell2012"
	
	#expression set object for gene expression dataset
	gep<-lymph$gep.log.corrected
	#display only two columns from fData(gep) in limma tables
	fData(gep)<-fData(gep)[, c("accession", "description")]

	#expression set object for alteration dataset
	cn<-lymph$cn.focal.binary

	#list of list representing cis genes in each alteration
	cisgenes<-lymph$gep.cisgenes$genes
	
	gepid<-"accession"
	cnid<-"Unique.Name"
	#make a vector of directions of alteration and add to fData(cn)
	cndirvector<-as.character(sapply(fData(cn)[, cnid], function(x){
		if (grepl("Amp", x)) return("Amplification")
		else return("Deletion")
		}))
	fData(cn)<-data.frame(fData(cn), alteration_direction = cndirvector)
	#display only three columns from fData(cn) in limma tables
	fData(cn)<-fData(cn)[, c("Unique.Name", "Descriptor", "alteration_direction")]
	cndir<-"alteration_direction"
	
	fdr.cis.cutoff<-0.25
	fdr.trans.cutoff<-0.05
	cat("Performancing cis analysis...\n")
	res.cis<-make_iEDGE(gep = gep, 
 		cn = cn,
 		cisgenes = cisgenes, 
 		gepid = gepid, 
 		cnid = cnid,
 		cndir = cndir,
 		interaction_type = "cis",
 		fdr.cutoff = fdr.cis.cutoff)

	f.out<-paste(f.dir.out, "/", header, "_cis_sig_fdr_", fdr.cis.cutoff, ".txt", sep = "")
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
 		fdr.cutoff = fdr.trans.cutoff)
	
	f.out<-paste(f.dir.out, "/", header, "_trans_sig_fdr_", fdr.cis.cutoff, ".txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(res.trans$sig, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)

	f.out<-paste(f.dir.out, "/", header, "_trans_full.txt", sep = "")
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(res.trans$full, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)
	
	cat("Done!\n")



	cat("\nPerforming hyperEnrichment analysis...\n")

 	source("run_hyperEnrichment.R")
 	source("read_gmt.R")
 	c2.file<-"../data/msigdb/c2.cp.v2.5.symbols.gmt"
  	c2<-read_gmt(c2.file)$genesets
 	c2<-sapply(c2, function(x){
    	a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
    	return(unique(do.call(c,a)))
    })


 	#run hyperenrichment with all sig cis genes in one geneset
 	drawnList<-list(cis.sig = unique(as.character(res.cis$sig$accession)))

 	hyper.cis.c2<-run_hyperEnrichment(drawn=drawnList,
    categories=c2,
    ntotal=nrow(lymph$gep.log.corrected),
    min.drawsize = 4, mht = TRUE, verbose = TRUE, order = TRUE)

 	f.out<-paste(f.dir.out, "/", header, "_hyperEnrichment_c2_cis.txt", sep = "")
	
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(hyper.cis.c2, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)
	
	cat("Done!\n")

	#run hyperenrichment with all sig trans genes in one geneset
	drawnList<-list(trans.sig = unique(as.character(res.trans$sig$accession)))

 	hyper.trans.c2<-run_hyperEnrichment(drawn=drawnList,
    categories=c2,
    ntotal=nrow(lymph$gep.log.corrected),
    min.drawsize = 4, mht = TRUE, verbose = TRUE, order = TRUE)

 	f.out<-paste(f.dir.out, "/", header, "_hyperEnrichment_c2_trans.txt", sep = "")
	
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(hyper.trans.c2, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)
	
	cat("Done!\n")

	#run hyperenrichment with lists of sig trans genes separated by the alteration it is contained in

	drawnList<-sapply(unique(res.trans$sig$Unique.Name), 
		function(x){
			y<-subset(res.trans$sig, Unique.Name == x)
			return(unique(as.character(y$accession)))
			})
	
 	hyper.trans.c2.split<-run_hyperEnrichment(drawn=drawnList,
    categories=c2,
    ntotal=nrow(lymph$gep.log.corrected),
    min.drawsize = 4, mht = TRUE, verbose = TRUE, order = TRUE)

 	f.out<-paste(f.dir.out, "/", header, "_hyperEnrichment_c2_trans_splitbyalteration.txt", sep = "")
	
	cat(paste("Writing table to ", f.out, "\n", sep = ""))
	write.table(hyper.trans.c2.split, sep = "\t", col.names = TRUE, row.names = FALSE,
		file = f.out)
	
	cat("Done!\n")

}

	