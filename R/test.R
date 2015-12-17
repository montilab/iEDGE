##new
some_test<-function(){
	res.dir="/Users/amyli/Desktop/git_projects/gistic2ge_results/ricover2015_new_gistic"
	f.header="ricover2015_"
	data.dir="/Users/amyli/Desktop/git_projects/datasets/ricover2015"
	gistic.dir=paste(data.dir, "job.56547077", sep = "/")	
	lesions="all_lesions.conf_99.txt"
	amp.peak="amp_genes_conf.99_sorted_in_peak.txt"
	amp.region="amp_genes_conf.99_sorted_in_region.txt"
	del.peak="del_genes_conf.99_sorted_in_peak.txt"
	del.region="del_genes_conf.99_sorted_in_region.txt"

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
	fmap.file=paste(data.dir,"ricover2014_mapping_new.txt", sep = "/")
	ge.confounder.file=paste(data.dir, "ricover2014.batches.cls", sep = "/")
	match.by = c("Unique.Name", "id")
	fmap.gistic.col = "snp.array.ID"
	fmap.ge.col = "mrna.array.ID"
	row.ind = TRUE
	logged = TRUE

	cat("reading gistic: amplifications\n")
	gistic.amp.new<-read_gistic(lesions.dir = lesions.file, 
		genes.peak.dir = amp.peak.file, 
		genes.region.dir = amp.region.file,
		match.by =  match.by, 
		cn.dir = "amp")
	cat("reading gistic: deletions\n")
	gistic.del.new<-read_gistic(lesions.dir = lesions.file, 
		genes.peak.dir = del.peak.file, 
		genes.region.dir = del.region.file, 
		match.by =  match.by, 
		cn.dir = "del")
	cat("reading gene expression\n")
	ge<-read_res(f = ge.file, row.ind = row.ind) #this data is unlogged!

	cat("matching ge expression and gistic samples:amplification\n")
	gistic_ge_amp<-match_samples(gistic = gistic.amp.new, 
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


###old
	res.dir="/Users/amyli/Desktop/git_projects/gistic2ge_results/ricover2015_new_gistic"
	f.header="ricover2015_"
	data.dir="/Users/amyli/Desktop/git_projects/datasets/ricover2015"
	gistic.dir=paste(data.dir, "job.50618380", sep = "/")	
	lesions="all_lesions.conf_99.txt"
	amp.peak="amp_genes.conf_99.sorted.txt"
	amp.region="amp_genes.conf_99.sorted_in_region.txt"
	del.peak="del_genes.conf_99.sorted.txt"
	del.region="del_genes.conf_99.sorted_in_region.txt"
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



	#fmap.file=paste(data.dir,"ricover2014_mapping_new.txt", sep = "/")
	#fmap<-read.table(fmap.file, sep = "\t", header = T)
	#fmap.sub<-subset(fmap, snp.array.ID %in% colnames(gistic.amp.new$mat) & 
	#	mrna.array.ID %in% colnames(ge$mat))
	#fmap.dir.new<-"/Users/amyli/Desktop/git_projects/datasets/ricover2015/ricover2014_mapping_new2.txt"
	#write.table(fmap.sub, file = fmap.dir.new, sep = "\t", col.names = T, quote = F)
	#fmap.file<-fmap.dir.new

	ge.confounder.file=paste(data.dir, "ricover2014.batches.cls", sep = "/")
	match.by = c("Descriptor", "cytoband")
	fmap.gistic.col = "snp.array.ID"
	fmap.ge.col = "mrna.array.ID"
	row.ind = TRUE
	logged = TRUE


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
	ge<-read_res(f = ge.file, row.ind = row.ind) #this data is unlogged!

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

}