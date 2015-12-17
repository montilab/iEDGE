require(plyr)
require(Biobase)

run_limma<-function(eset, design){
	require(limma)
	#cat("running limma\n")
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
				if (fit2.table$alteration_direction[i] == "Amplification")
					lower.tail<-fit2$t[i] < 0
				else 
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

make_gistic2ge<-function(x, interaction_type = "cis", 
	region_type = c("focal", "focal_or_arm", "arm", "arm_in_focal")){
	require(plyr)

	gep<-x$gep.log.corrected

	cisgenes<-x$gep.cisgenes$genes

	if (region_type == "focal")
		cn<-x$cn.focal.binary
	else if (region_type == "focal_or_arm")
		cn <-x$cn.focal.or.arm.binary
	else if (region_type == "arm"){
		cn <-x$cn.arm.binary
		cisgenes<-x$gep.cisgenes.arm$genes
	}
	else if (region_type == "arm_in_focal"){
		cn<-x$cn.arm.infocal.binary
		cisgenes<-x$gep.cisgenes.arm.infocal$genes
	}

	all_genes<-as.character(fData(gep)$accession)
	cat(paste("total of ", nrow(cn), " alterations: \n\n", sep = ""))

	alteration_ids<-fData(cn)$Unique.Name
	alteration_desc<-fData(cn)$Descriptor
	res<-list()
	#i: for each alteration
	for (i in 1:nrow(cn)){ 
		cat(paste(i, "\n", sep = ""))
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
			fdat<-fdat[, c("accession", "description")]
			colnames(fdat)<-c("gene_symbol", "gene_description")
			fdat.add<-fData(cn)[i, c("Unique.Name", "Descriptor")]
			fdat.combined<-cbind(fdat, fdat.add)
			if (interaction_type == "cis"){
				fdat.combined$interaction_type <-"cis"
			} else {
				fdat.combined$interaction_type <-"trans"
			}
			if (grepl("Amp", fdat.combined$Unique.Name[1])){
				fdat.combined$alteration_direction <-"Amplification"
			} else{
				fdat.combined$alteration_direction <-"Deletion"
			}
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
					alteration_ids[i], "\n\n", sep = ""))
				res[[i]]<-NA
			}
			else 
				res[[i]]<-run_limma(eset = gep.keep, design = design)
		}
	}
	res.df<-do.call(rbind, res[!is.na(res)])
	res.df$adj.P.Val.all<-p.adjust(res.df$P.Value, method = "fdr")

	res.list<-list()
	res.list$full<-res.df
	if(interaction_type == "cis")
		fdr.cutoff <- 0.25
	else 
		fdr.cutoff <- 0.05
	res.list$sig<-subset(res.df, adj.P.Val.all < fdr.cutoff)

	return(res.list)
}

run_dataset<-function(x, f.dir.out, header){
	cat(paste("Running for ", header, "\n", sep = ""))
	#lymph<-readRDS("/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/cancercell2012/cancercell2012_all.RDS")
	#x<-lymph
	x.de.cis.focal<-make_gistic2ge(x = x, interaction_type = "cis", region_type = c("focal"))
	x.de.trans.focal<-make_gistic2ge(x = x, interaction_type = "trans", region_type = c("focal"))
	x.de.cis.focalorarm<-make_gistic2ge(x = x, interaction_type = "cis", region_type = c("focal_or_arm"))
	x.de.trans.focalorarm<-make_gistic2ge(x = x, interaction_type = "trans", region_type = c("focal_or_arm"))
	
	x.de.cis.arm<-make_gistic2ge(x = x, interaction_type = "cis", region_type = c("arm"))
	x.de.trans.arm<-make_gistic2ge(x = x, interaction_type = "trans", region_type = c("arm"))
	x.de.cis.arm.infocal<-make_gistic2ge(x = x, interaction_type = "cis", region_type = c("arm_in_focal"))
	x.de.trans.arm.infocal<-make_gistic2ge(x = x, interaction_type = "trans", region_type = c("arm_in_focal"))
					
	#require(CBMRtools)
	#f.dir.out<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/cancercell2012"
#full results (Cis only)
	full<-list(cisfocal = x.de.cis.focal$full, 
		cisfocalorarm = x.de.cis.focalorarm$full)
	###save.xlsx(x = full, f.dir = f.dir.out, f.name = paste(header, "_limma_DE_full", sep = ""))

	#full.arm<-list(cis.arm = x.de.cis.arm$full, 
	#	cis.arm.infocal = x.de.cis.arm.infocal$full)	
	#save.xlsx(x = full.arm, f.dir = f.dir.out, f.name = paste(header, "_limma_DE_arm_full", sep = ""))
#sig results
	sig<-list(cisfocal = x.de.cis.focal$sig, 
		cisfocalorarm = x.de.cis.focalorarm$sig, 
		transfocal = x.de.trans.focal$sig,
		transfocalorarm = x.de.trans.focalorarm$sig)
	###save.xlsx(x = sig, f.dir = f.dir.out, f.name = paste(header, "_limma_DE_sig", sep = ""))

	#sig.arm<-list(cis.arm = x.de.cis.arm$sig, 
	#	cis.arm.infocal = x.de.cis.arm.infocal$sig, 
	#	trans.arm = x.de.trans.arm$sig,
	#	trans.arm.infocal = x.de.trans.arm.infocal$sig)
	#save.xlsx(x = sig.arm, f.dir = f.dir.out, f.name = paste(header, "_limma_DE_arm_sig", sep = ""))

	#return(list(full = full, sig = sig, sig.arm = sig.arm))
	return(list(full = full, sig = sig))
}

if (FALSE){
	cat("Running for test data: lymphoma2012")
	lymph <-readRDS("../test/data/cancercell2012_all.RDS")
	f.dir.out <- "../test/tables"
	header <- "cancercell2012"
	res1<-run_dataset(x = lymph, f.dir.out = f.dir.out, header = header)

	unique_name<-"Unique.Name"
	res1.summary<-sapply(names(res1), function(x){
		sapply(names(res1[[x]]), function(y){
			sapply(unique(res1[[x]][[y]][, unique_name]), function(z){
				table_full<-res1[[x]][[y]]
				table_sub<-table_full[table_full[, unique_name] == z,]
				out.file<-paste(f.dir.out, "/", "byalteration", "_", x, "_", y, "_", z, ".txt", sep="")
				write.table(table_sub, file = out.file, quote = FALSE, sep = "\t", 
					row.names = FALSE, col.names = TRUE)
				#df<-data.frame(sig.cutoff = x, analysis.type= y,
				# alteration.name = z, number.genes = nrow(table_sub))
				#return(df)
#				return(list(sig.cutoff = x, analysis.type= y,
#				 alteration.name = z, number.genes = nrow(table_sub)))
			})		
		})
	})

	alterations<-unique(res1[[1]][[1]][, unique_name])

	res1.summary<-sapply(names(res1), function(x){	
		##data frame: row = alteration, col = test
		sapply(names(res1[[x]]), function(y){
			table_full<-res1[[x]][[y]]
			table_sub<-lapply(alterations, function(z) 
				res1[[x]][[y]][res1[[x]][[y]][, unique_name] == z, ])
			n<-sapply(table_sub, function(x) nrow(x))
			return(n)
			})
		})

	res1.summary.cbind<-do.call(cbind, sapply(names(res1.summary), function(x){
		res1.summary.full<-res1.summary[[x]]
		colnames(res1.summary.full)<-paste(x, 
			colnames(res1.summary.full), sep = "_")
		return(res1.summary.full)
		}))
	res1.summary.cbind<-data.frame(alteration_id = alterations, res1.summary.cbind)
	out.file<-paste(f.dir.out, "/", "summarytable.txt", sep = "")
	write.table(res1.summary.cbind, file = out.file, quote = FALSE, sep = "\t", 
					row.names = FALSE, col.names = TRUE)

}

	