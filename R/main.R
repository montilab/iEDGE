sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm,":")           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
}


run_process_data<-function(
	amp.peak.file, amp.region.file, 
	del.peak.file, del.region.file,
	lesions.file, ge.file,
	fmap.file, ge.confounder.file, 
	match.by, fmap.gistic.col, fmap.ge.col, row.ind, ge.format
){

	gistic_ge<-process_data(
		amp.peak.file = amp.peak.file,
		amp.region.file = amp.region.file,
		del.peak.file = del.peak.file,
		del.region.file = del.region.file,
		lesions.file = lesions.file,
		ge.file = ge.file,
		fmap.file = fmap.file,
		ge.confounder.file = ge.confounder.file,
		match.by = match.by,
		fmap.gistic.col = fmap.gistic.col,
		fmap.ge.col = fmap.ge.col, 
		row.ind = row.ind, 
		ge.format = ge.format)
	return(gistic_ge)
}

run_lymphoma2010<-function(){
	source("main.R")
	sourceDir(getwd())

	res.dir="/Users/amyli/Desktop/git_projects/gistic2ge_results/lymphoma2010"
	f.header="lymphoma2010_"
	data.dir="/Users/amyli/Desktop/git_projects/datasets/lymphoma2010"
	gistic.dir=paste(data.dir, "GISTIC_20110620.4amy", sep = "/")
	lesions="all_lesions_file.txt"
	amp.peak="Amp_genes.txt"
	amp.region="Amp_genes_in_region.txt"
	del.peak="Del_genes.txt"
	del.region="Del_genes_in_region.txt"

	amp.peak.file=paste(gistic.dir, amp.peak, sep = "/")
	amp.region.file=paste(gistic.dir, amp.region, sep = "/")
	del.peak.file=paste(gistic.dir, del.peak, sep = "/")
	del.region.file=paste(gistic.dir, del.region, sep = "/")
	lesions.file=paste(gistic.dir, lesions, sep = "/")
	ge.file=paste(data.dir, "lymphoma2010_wold.entrez.mas.gnorm.GS.res", sep = "/")
	fmap.file=paste(data.dir, "snp2express.mapping.lymphoma2010_wold.txt", sep = "/")
	ge.confounder.file=paste(data.dir, "lymphoma2010_wold.03vs10.cls", sep = "/")
	match.by = c("Descriptor", "cytoband")
	fmap.gistic.col = "snp.array.ID"
	fmap.ge.col = "mrna.array.ID"
	row.ind = FALSE
	logged = FALSE

	gistic_ge<-run_process_data(amp.peak.file, amp.region.file, 
		del.peak.file, del.region.file,
		lesions.file, ge.file,
		fmap.file, ge.confounder.file, 
		match.by, fmap.gistic.col, fmap.ge.col, row.ind)

	gistic_ge_amp<-gistic_ge$amp
	gistic_ge_del<-gistic_ge$del

	util.save.rds(gistic_ge_amp, f.dir = res.dir, f.header = f.header, f.name = "gistic_ge_amp")
	util.save.rds(gistic_ge_del, f.dir = res.dir, f.header = f.header, f.name = "gistic_ge_del")


	gistic2ge.full<-do_gistic2ge_all(gistic_ge_amp, gistic_ge_del, logged = logged)
	
	util.save.rds(gistic2ge.full, f.dir = res.dir, f.header = f.header, f.name = "gistic2ge_full")

	params1<-list(cis.fdr = 0.25, trans.fdr = 0.05, sig.FC = NA)
	params2<-list(cis.fdr = 0.25, trans.fdr = 0.05, sig.FC = 1.1)
	params3<-list(cis.fdr = 0.25, trans.fdr = 0.25, sig.FC = 1.1)

	param.list<-list(params1, params2, params3)

	names(param.list)<-sapply(param.list, function(x) 
		paste(names(x), x, sep = "_", collapse = "_"))

	gistic2ge.sig<-list()
	for(i in names(param.list)){
		print(i)
		gistic2ge.sig[[i]]<-do.call(subset_gistic2ge, c(list(x = gistic2ge.full), param.list[[i]]))
		#util.save.xlsx(gistic2ge.sig[[i]], f.dir = res.dir, f.header = f.header,
		#	f.name = paste("gistic2ge_sig", i, sep= "_"))
		util.save.rds(gistic2ge.sig[[i]], f.dir = res.dir, f.header = f.header,
			f.name = paste("gistic2ge_sig", i, sep= "_"))
	}


	c2.file="~/Desktop/monti_lab/year2/DataIntegration/annot/msigdb_v2.5/c2.cp.v2.5.symbols.gmt"
	c3.file="~/Desktop/monti_lab/year2/DataIntegration/annot/msigdb_v2.5/c3.tft.v2.5.symbols.gmt"
 	c2<-read_gmt(c2.file)$genesets
  	c2<-sapply(c2, function(x){
    	a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
    	return(unique(do.call(c,a)))
    })
  	c3<-read_gmt(c3.file)$genesets
  	c3<-sapply(c3, function(x){
    	a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
    	return(unique(do.call(c,a)))
    })



  	cat("running hyperenrichment\n")

    gs.pathway.list<-list(c2 = c2, c3 = c3)
    ntot = nrow(gistic_ge_amp$ge$mat)
   
    hyper.res<-list()
    for(i in c("c2", "c3")){
    	hyper.res[[i]]<-list()
	    for (j in names(gistic2ge.sig)){
			hyper.res[[i]][[j]]<-run_hyperEnrichment_wrapper(
				gs.tab = gistic2ge.sig[[j]], 
				ntot = ntot, 
		  		gs.pathway = gs.pathway.list[[i]])
		}
	}

	hyper.res.sig<-list()
	hyper.fdr = 0.25

	for(i in names(hyper.res)){ 
 		hyper.res.sig[[i]]<-list()
		for (j in names(hyper.res[[i]])){
			hyper.res.sig[[i]][[j]]<-list()
			for (k in names(hyper.res[[i]][[j]])){
				hyper.res.sig[[i]][[j]][[k]]<-subset_hyperEnrichment(hyper.res[[i]][[j]][[k]], 
					hyper.fdr = hyper.fdr)
				util.save.xlsx(hyper.res.sig[[i]][[j]][[k]], f.dir = res.dir, f.header = f.header,
				f.name = paste("hyperEnrichment_sig", i, k, j, "hyper.fdr", hyper.fdr, sep= "_"))   
			}
		}
	}

}

order_gistic<-function(gistic.dir){

	#gistic.dir = "/Users/amyli/Desktop/git_projects/datasets/ricover2015/job.56547077"
	#gistic.dir=paste(data.dir, "job.56547077", sep = "/")
	for (cn.dir in c("amp", "del")){
		tab.file.name <- paste(cn.dir, "_genes.conf_99.txt", sep = "")
		tab.file <- paste(gistic.dir, tab.file.name, sep = "/")

		df<-read.table(tab.file, sep = "\t", fill = TRUE, header = FALSE)
		lesions<-"all_lesions.conf_99.txt"
		lesions.file<-paste(gistic.dir, lesions, sep = "/")
		lesions.obj<-read.table(lesions.file, sep = "\t", header = TRUE)

		lesions.obj<-lesions.obj[!grepl("CN", lesions.obj$Unique.Name),]
		if (cn.dir == "del"){
			lesions.obj<-lesions.obj[grepl("Del", lesions.obj$Unique.Name),]
		} else {
			lesions.obj<-lesions.obj[grepl("Amp", lesions.obj$Unique.Name),]
		}

		lesions.ord<-as.character(lesions.obj$Descriptor)
		lesions.ord<-gsub(" ", "", lesions.ord)

		df.col1<-df$V1
		df.rest<-df[, -1]
		df.rest<-df.rest[, apply(df.rest, 2, function(x) !(all(is.na(x))))]

		df.rest.order<-match(lesions.ord, as.character(t(df.rest[1,])))
		df.rest<-df.rest[, df.rest.order]

		df.combined<-cbind(df.col1, df.rest)
		f.out<-paste(cn.dir, "_genes.conf_99_ordered.txt", sep = "")
		write.table(df.combined, file =  paste(gistic.dir, f.out, sep = "/"),
			sep = "\t", col.names = F, row.names = F)
	}
}




make_gistic_new<-function(gistic.dir){

	#data.dir="/Users/amyli/Desktop/git_projects/datasets/ricover2015"
	#gistic.dir=paste(data.dir, "job.56547077", sep = "/")
	#do region first
	#cn.dir <- "del"


	match_ref<-function(ind){
	 	match_ind<-NA
	 	for (i in 1:length(df.region.genes)){
	 		if (length(intersect(df.region.genes[[i]],ref.df[[ind]])) > (0.6*length(ref.df[[ind]])))
	 			match_ind <- i
	 	}
	 	return(match_ind)
	}

	for(cn.dir in c("amp","del")){
	ref.file<-paste(cn.dir, "_genes.conf_99_ordered.txt", sep = "")
	ref.df<-read.table(paste(gistic.dir, ref.file, sep = "/"), 
		sep = "\t", fill = TRUE, header = FALSE)


	ref.df<-read_cis_genes(file.in.peak = paste(gistic.dir, ref.file, sep = "/"), 
	file.in.region = NA, 
	cn.dir)
	ref.df<-ref.df$genesets.peak

	ref.df<-lapply(ref.df, function(x){
		as.character(sapply(x, function(i) gsub("\\[|\\]", "", sub("\\|.*$", "",i))))
		})


	tab.file.name <- paste("table_", cn.dir, ".conf_99.txt", sep = "")
	tab.file <- paste(gistic.dir, tab.file.name, sep = "/")

	df<-read.table(tab.file, sep = "\t", fill = TRUE, header = TRUE)

	df.region.genes<-sapply(df[, "genes_in_region"], 
				function(x) strsplit(as.character(x), split = ","))

	df.region.genes<-lapply(df.region.genes, function(x){
		as.character(sapply(x, function(i) gsub("\\[|\\]", "", sub("\\|.*$", "",i))))
		})

	df.peak.genes<-sapply(df[, "genes_in_peak"], 
			function(x) strsplit(as.character(x), split = ","))
	
	df.peak.genes<-lapply(df.peak.genes, function(x){
		as.character(sapply(x, function(i) gsub("\\[|\\]", "", sub("\\|.*$", "",i))))
		})
	
	new.ord<-sapply(1:length(ref.df), function(x) match_ref(x))

	df.region.genes.ord<-df.region.genes[new.ord]
	df.peak.genes.ord<-df.peak.genes[new.ord]

	if(cn.dir == "amp") {
		meta.header<-"AmplificationPeak"
	}
	else {
		meta.header<-"DeletionPeak"
	}

	df.meta<-data.frame(index = df[, "index"], 
		id = paste(meta.header, df[,"index"], sep=""))

	df.geneset.peak<-matrix("", 
		nrow =max(sapply(df.peak.genes.ord,length)), 
		ncol =length(df.peak.genes.ord))

	for(i in 1:length(df.peak.genes.ord))
		for(j in 1:length(df.peak.genes.ord[[i]]))
			df.geneset.peak[j,i]<-df.peak.genes.ord[[i]][[j]]
	
	df.all.peak<-rbind(t(df.meta), df.geneset.peak)
	rownames(df.all.peak)[3]<-"genes in wide peak"

	df.geneset.region<-matrix("", 
		nrow =max(sapply(df.region.genes.ord,length)), 
		ncol =length(df.region.genes.ord))

	for(i in 1:length(df.region.genes.ord))
		for(j in 1:length(df.region.genes.ord[[i]]))
			df.geneset.region[j,i]<-df.region.genes.ord[[i]][[j]]

	df.all.region<-rbind(t(df.meta), df.geneset.region)
	rownames(df.all.region)[3]<-"genes in wide peak"

	df.out.peak.name<-paste(cn.dir, "_genes_conf.99_sorted_in_peak_ordered.txt", sep = "")
	write.table(df.all.peak, file = paste(gistic.dir, df.out.peak.name, sep = "/"), 
		sep = "\t", row.names = T, col.names = F)

	df.out.region.name<-paste(cn.dir, "_genes_conf.99_sorted_in_region_ordered.txt", sep = "")
	write.table(df.all.region, file = paste(gistic.dir, df.out.region.name, sep = "/"), 
		sep = "\t", row.names = T, col.names = F)
	}

}


add_gistic_files<-function(gistic.dir){
	order_gistic(gistic.dir)
	make_gistic_new(gistic.dir)
}

make_ricover2015_gistic<-function(){
	data.dir="/Users/amyli/Desktop/git_projects/datasets/ricover2015"
	gistic.dir=paste(data.dir, "job.56547077", sep = "/")

	for(cn.dir in c("amp", "del")){	
		for(boundary in c("peak", "region")){

			tab.file.name <- paste("table_", cn.dir, ".conf_99.txt", sep = "")
			tab.file <- paste(gistic.dir, tab.file.name, sep = "/")

			df<-read.table(tab.file, sep = "\t", fill = TRUE, header = TRUE)
			df<-df[with(df, order(df[, "chromosome"], df[,"region_start"])), ]
			df.col<-paste("genes_in_", boundary, sep = "")
			res.geneset.peak<-sapply(df[, df.col], 
				function(x) strsplit(as.character(x), split = ","))

			if(cn.dir == "amp") meta.header<-"AmplificationPeak"
			else meta.header<-"DeletionPeak"
			df.meta<-data.frame(index = df[, "index"], 
				id = paste(meta.header, df[,"index"], sep=""))

			df.geneset.peak<-matrix("", 
				nrow =max(sapply(res.geneset.peak,length)), 
				ncol =length(res.geneset.peak))
			for(i in 1:length(res.geneset.peak)){
				for(j in 1:length(res.geneset.peak[[i]])){
					df.geneset.peak[j,i]<-res.geneset.peak[[i]][[j]]
				}
			}
			df.all<-rbind(t(df.meta), df.geneset.peak)
			rownames(df.all)[3]<-"genes in wide peak"
			df.out.name<-paste(cn.dir, "_genes_conf.99_sorted_in_", boundary, ".txt", sep = "")
			write.table(df.all, file = paste(gistic.dir, df.out.name, sep = "/"), 
				sep = "\t", row.names = T, col.names = F)
		}
	}
}



run_ricover2015<-function(){
	source("main.R")
	sourceDir(getwd())

	res.dir="/Users/amyli/Desktop/git_projects/gistic2ge_results/ricover2015_new_gistic"
	f.header="ricover2015_"
	data.dir="/Users/amyli/Desktop/git_projects/datasets/ricover2015"
	gistic.dir=paste(data.dir, "job.56547077", sep = "/")	
	lesions="all_lesions.conf_99.txt"
	amp.peak="amp_genes_conf.99_sorted_in_peak_ordered.txt"
	amp.region="amp_genes_conf.99_sorted_in_region_ordered.txt"
	del.peak="del_genes_conf.99_sorted_in_peak_ordered.txt"
	del.region="del_genes_conf.99_sorted_in_region_ordered.txt"


	amp.peak.file=paste(gistic.dir, amp.peak, sep = "/")
	amp.region.file=paste(gistic.dir, amp.region, sep = "/")
	del.peak.file=paste(gistic.dir, del.peak, sep = "/")
	del.region.file=paste(gistic.dir, del.region, sep = "/")
	lesions.file=paste(gistic.dir, lesions, sep = "/")
	ge.file=paste(data.dir, "ricover2014.rma.hgu133plus2hsensgcdf.res", sep = "/")
	

	#fmap<-read.table(fmap.file, sep = "\t", header = T)
	#fmap$snp.array.ID<-gsub("-", ".", fmap$snp.array.ID)
	#fmap.dir.new<-"/Users/amyli/Desktop/git_projects/datasets/ricover2015/ricover2014_mapping_new.txt"
	#write.table(fmap, file = fmap.dir.new, sep = "\t", col.names = T, quote = F)
	#fmap.file<-fmap.dir.new
	fmap.file=paste(data.dir,"ricover2014_mapping_new2.txt", sep = "/")
	ge.confounder.file=paste(data.dir, "ricover2014.batches.cls", sep = "/")
	match.by = c("Unique.Name", "id")
	fmap.gistic.col = "snp.array.ID"
	fmap.ge.col = "mrna.array.ID"
	row.ind = TRUE
	logged = FALSE

	gistic_ge<-run_process_data(amp.peak.file, amp.region.file, 
		del.peak.file, del.region.file,
		lesions.file, ge.file,
		fmap.file, ge.confounder.file, 
		match.by, fmap.gistic.col, fmap.ge.col, row.ind)

	gistic_ge_amp<-gistic_ge$amp
	gistic_ge_del<-gistic_ge$del

	util.save.rds(gistic_ge_amp, f.dir = res.dir, f.header = f.header, f.name = "gistic_ge_amp")
	util.save.rds(gistic_ge_del, f.dir = res.dir, f.header = f.header, f.name = "gistic_ge_del")


	gistic2ge.full<-do_gistic2ge_all(gistic_ge_amp, gistic_ge_del, logged = logged)
	

	util.save.rds(gistic2ge.full, f.dir = res.dir, f.header = f.header, f.name = "gistic2ge_full")

	gistic2ge.full.cis<-gistic2ge.full[grep("cis", names(gistic2ge.full))]
	util.save.xlsx(gistic2ge.full.cis, f.dir = res.dir, f.header = f.header, f.name = "gistic2ge_full_cis")


	params1<-list(cis.fdr = 0.25, trans.fdr = 0.05, sig.FC = NA)
	params2<-list(cis.fdr = 0.25, trans.fdr = 0.05, sig.FC = 1.1)
	params3<-list(cis.fdr = 0.25, trans.fdr = 0.25, sig.FC = 1.1)

	param.list<-list(params1, params2, params3)

	names(param.list)<-sapply(param.list, function(x) 
		paste(names(x), x, sep = "_", collapse = "_"))

	gistic2ge.sig<-list()
	for(i in names(param.list)){
		print(i)
		gistic2ge.sig[[i]]<-do.call(subset_gistic2ge, c(list(x = gistic2ge.full), param.list[[i]]))
		util.save.xlsx(gistic2ge.sig[[i]], f.dir = res.dir, f.header = f.header,
			f.name = paste("gistic2ge_sig", i, sep= "_"))
	}


	c2.file="~/Desktop/monti_lab/year2/DataIntegration/annot/msigdb_v2.5/c2.cp.v2.5.symbols.gmt"
	c3.file="~/Desktop/monti_lab/year2/DataIntegration/annot/msigdb_v2.5/c3.tft.v2.5.symbols.gmt"
 	c2<-read_gmt(c2.file)$genesets
  	c2<-sapply(c2, function(x){
    	a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
    	return(unique(do.call(c,a)))
    })
  	c3<-read_gmt(c3.file)$genesets
  	c3<-sapply(c3, function(x){
    	a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
    	return(unique(do.call(c,a)))
    })


    cat("running hyperenrichment\n")

    gs.pathway.list<-list(c2 = c2, c3 = c3)
    ntot = nrow(gistic_ge_amp$ge$mat)
    hyper.fdr = 0.25
    hyper.res<-list()
    for(i in c("c2", "c3")){
    	hyper.res[[i]]<-list()
	    for (j in names(gistic2ge.sig)){
			hyper.res[[i]][[j]]<-run_hyperEnrichment_wrapper(
				gs.tab = gistic2ge.sig[[j]], 
				ntot = ntot, 
		  		gs.pathway = gs.pathway.list[[i]])
		}
	}

	hyper.res.sig<-list()

	for(i in names(hyper.res)){ 
 		hyper.res.sig[[i]]<-list()
		for (j in names(hyper.res[[i]])){
			hyper.res.sig[[i]][[j]]<-list()
			for (k in names(hyper.res[[i]][[j]])){
				hyper.res.sig[[i]][[j]][[k]]<-subset_hyperEnrichment(hyper.res[[i]][[j]][[k]], 
					hyper.fdr = hyper.fdr)
				util.save.xlsx(hyper.res.sig[[i]][[j]][[k]], f.dir = res.dir, f.header = f.header,
				f.name = paste("hyperEnrichment_sig", i, k, j, "hyper.fdr", hyper.fdr, sep= "_"))   
			}
		}
	}
}


run_TCGA_brca_tn<-function(){
	source("main.R")
	sourceDir(getwd())

	res.dir="/Users/amyli/Desktop/git_projects/gistic2ge_results/TCGA_brca_tn"
	f.header="TCGA_brca_tn_"
	data.dir="/Users/amyli/Desktop/git_projects/datasets/TCGA_brca"
	gistic.dir=paste(data.dir, "gistic2", sep = "/")

	add_gistic_files(gistic.dir)


	lesions="all_lesions.conf_99.txt"
	amp.peak="amp_genes_conf.99_sorted_in_peak_ordered.txt"
	amp.region="amp_genes_conf.99_sorted_in_region_ordered.txt"
	del.peak="del_genes_conf.99_sorted_in_peak_ordered.txt"
	del.region="del_genes_conf.99_sorted_in_region_ordered.txt"

	amp.peak.file=paste(gistic.dir, amp.peak, sep = "/")
	amp.region.file=paste(gistic.dir, amp.region, sep = "/")
	del.peak.file=paste(gistic.dir, del.peak, sep = "/")
	del.region.file=paste(gistic.dir, del.region, sep = "/")
	lesions.file=paste(gistic.dir, lesions, sep = "/")
	ge.file=paste(data.dir, "TCGA_brca_eset_tn_tumor.RDS", sep = "/")
	

	#fmap<-read.table(fmap.file, sep = "\t", header = T)
	#fmap$snp.array.ID<-gsub("-", ".", fmap$snp.array.ID)
	#fmap.dir.new<-"/Users/amyli/Desktop/git_projects/datasets/ricover2015/ricover2014_mapping_new.txt"
	#write.table(fmap, file = fmap.dir.new, sep = "\t", col.names = T, quote = F)
	#fmap.file<-fmap.dir.new
	fmap.file=paste(data.dir,"file_map.txt", sep = "/")
	ge.confounder.file=NA
	match.by = c("Unique.Name", "id")
	fmap.gistic.col = "snp.array.ID"
	fmap.ge.col = "mrna.array.ID"
	row.ind = TRUE
	logged = FALSE
	ge.format = "eset"

	gistic_ge<-run_process_data(amp.peak.file, amp.region.file, 
		del.peak.file, del.region.file,
		lesions.file, ge.file,
		fmap.file, ge.confounder.file, 
		match.by, fmap.gistic.col, fmap.ge.col, row.ind, ge.format)

	gistic_ge_amp<-gistic_ge$amp
	gistic_ge_del<-gistic_ge$del

	util.save.rds(gistic_ge_amp, f.dir = res.dir, f.header = f.header, f.name = "gistic_ge_amp")
	util.save.rds(gistic_ge_del, f.dir = res.dir, f.header = f.header, f.name = "gistic_ge_del")


	gistic2ge.full<-do_gistic2ge_all(gistic_ge_amp, gistic_ge_del, logged = logged)
	
	util.save.rds(gistic2ge.full, f.dir = res.dir, f.header = f.header, f.name = "gistic2ge_full")

	params1<-list(cis.fdr = 0.25, trans.fdr = 0.05, sig.FC = NA)
	params2<-list(cis.fdr = 0.25, trans.fdr = 0.05, sig.FC = 1.1)
	params3<-list(cis.fdr = 0.25, trans.fdr = 0.25, sig.FC = 1.1)

	param.list<-list(params1, params2, params3)

	names(param.list)<-sapply(param.list, function(x) 
		paste(names(x), x, sep = "_", collapse = "_"))

	gistic2ge.sig<-list()
	for(i in names(param.list)){
		print(i)
		gistic2ge.sig[[i]]<-do.call(subset_gistic2ge, c(list(x = gistic2ge.full), param.list[[i]]))
		util.save.xlsx(gistic2ge.sig[[i]], f.dir = res.dir, f.header = f.header,
			f.name = paste("gistic2ge_sig", i, sep= "_"))
	}


	c2.file="~/Desktop/monti_lab/year2/DataIntegration/annot/msigdb_v2.5/c2.cp.v2.5.symbols.gmt"
	c3.file="~/Desktop/monti_lab/year2/DataIntegration/annot/msigdb_v2.5/c3.tft.v2.5.symbols.gmt"
 	c2<-read_gmt(c2.file)$genesets
  	c2<-sapply(c2, function(x){
    	a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
    	return(unique(do.call(c,a)))
    })
  	c3<-read_gmt(c3.file)$genesets
  	c3<-sapply(c3, function(x){
    	a<-lapply(x,function(i) strsplit(i, split="[ |_]///[ |_]")[[1]])
    	return(unique(do.call(c,a)))
    })


  

  	cat("running hyperenrichment\n")

    gs.pathway.list<-list(c2 = c2, c3 = c3)
    ntot = nrow(gistic_ge_amp$ge$mat)
   
    hyper.res<-list()
    for(i in c("c2", "c3")){
    	hyper.res[[i]]<-list()
	    for (j in names(gistic2ge.sig)){
			hyper.res[[i]][[j]]<-run_hyperEnrichment_wrapper(
				gs.tab = gistic2ge.sig[[j]], 
				ntot = ntot, 
		  		gs.pathway = gs.pathway.list[[i]])
		}
	}

	hyper.res.sig<-list()
	hyper.fdr = 0.25

	for(i in names(hyper.res)){ 
 		hyper.res.sig[[i]]<-list()
		for (j in names(hyper.res[[i]])){
			hyper.res.sig[[i]][[j]]<-list()
			for (k in names(hyper.res[[i]][[j]])){
				hyper.res.sig[[i]][[j]][[k]]<-subset_hyperEnrichment(hyper.res[[i]][[j]][[k]], 
					hyper.fdr = hyper.fdr)
				util.save.xlsx(hyper.res.sig[[i]][[j]][[k]], f.dir = res.dir, f.header = f.header,
				f.name = paste("hyperEnrichment_sig", i, k, j, "hyper.fdr", hyper.fdr, sep= "_"))   
			}
		}
	}


}



