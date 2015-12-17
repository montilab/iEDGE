require(CBMRtools)
require(Biobase)

read_GISTIC2_focal<-function(all_lesions){
	x<-read.table(all_lesions,sep = "\t", header =T)
	x<-x[, !apply(x, 2, function(i){all(is.na(i))})]
	x.meta.ind<-1:which(colnames(x) == "Amplitude.Threshold")
	x<-x[x$Amplitude.Threshold != "Actual Copy Change Given",]
	#make eset
	fdat<-x[, x.meta.ind]
	for(i in colnames(fdat)){
		fdat[, i]<-gsub(" ", "", fdat[,i])
	}
	dat<-x[, -x.meta.ind]
	pdat<-data.frame(sample_id = colnames(dat))
	eset<-to.eSet(mat = dat, pdat = pdat, fdat = fdat)
	return(eset)
}

read_GISTIC2_broad<-function(broad_lesions){
	x<-read.table(broad_lesions,sep = "\t", header =T)
	x<-x[, !apply(x, 2, function(i){all(is.na(i))})]
	x.meta.ind<-which(colnames(x) == "Chromosome.Arm")
#	x<-x[x$Amplitude.Threshold != "Actual Copy Change Given",]
	#make eset
	fdat<-data.frame(Descriptor = x[, x.meta.ind])
	for(i in colnames(fdat)){
		fdat[, i]<-gsub(" ", "", fdat[,i])
	}
	dat<-x[, -x.meta.ind]
	pdat<-data.frame(sample_id = colnames(dat))
	eset<-to.eSet(mat = dat, pdat = pdat, fdat = fdat)
	return(eset)
#return(x)
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
	for(ind.focal in 1:nrow(eset.focal)){
		name.focal<-fData(eset.focal)$Unique.Name[ind.focal]
		desc.focal<-fData(eset.focal)$Descriptor[ind.focal]
		desc.arm<-get_chr(desc.focal)
		if(grepl("Amplification", fData(eset.focal)$Unique.Name[ind.focal]))
			ind.arm<-which(fData(eset.arm)$Descriptor == desc.arm & grepl("Amplification", fData(eset.arm)$Unique.Name))	
		else
			ind.arm<-which(fData(eset.arm)$Descriptor == desc.arm & grepl("Deletion", fData(eset.arm)$Unique.Name))
		print(ind.arm)
		if(length(ind.arm)!=0){
			##has matching arm level copy number alteration
			x.both<-rbind(exprs(eset.focal[ind.focal,]), exprs(eset.arm[ind.arm,]))
			x.or<-apply(x.both, 2, max) #take max of copy number status
			exprs(eset.or)[ind.focal,]<-x.or
			fData(eset.or)$Unique.Name[ind.focal]<-paste(fData(eset.or)$Unique.Name[ind.focal],
				"_or_arm", sep = "")
		}
	}
	return(eset.or)
}

get_cis_genes<-function(f.genes, eset.focal, direction){

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
	res.genes<-lapply(res.genes, function(x) gsub("\\[|\\]", "", x))
	return(list(meta = res.meta, genes = res.genes))
}
####remove, old function
get_region_genes<-function(cis.genes, f.tab.amp, f.tab.del){
	tab.amp<-read.table(f.tab.amp, sep = "\t", header = TRUE, fill = TRUE)
	tab.del<-read.table(f.tab.del, sep = "\t", header = TRUE, fill = TRUE)

	tab<-rbind(tab.amp, tab.del)

	tab.genes.peak<-strsplit(as.character(tab$genes_in_peak), split = ",")
	tab.genes.region<-strsplit(as.character(tab$genes_in_region), split = ",")

	tab.genes.ind<-sapply(1:length(cis.genes$genes), function(i){
		for(j in 1:length(tab.genes.peak)){
			if(length(cis.genes$genes[[j]]) == length(tab.genes.peak[[j]])){
				if (all(cis.genes$genes[[j]]  == tab.genes.peak[[j]])) return(i)
			}
		} 
		return(NA)
		})
	tab.genes.ord<-tab.genes.region[tab.genes.ind]
	return(tab.genes.ord)
}


if(TRUE){
	f.dir.in<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/wes_7_16_15/job.60875343"
	f.dir.out<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/wes_7_16_15"

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
	saveRDS(eset.focal, file = paste(f.dir.out, "/", "ricover_wes_gistic_all_focal.RDS", sep = ""))
	saveRDS(eset.arm, file = paste(f.dir.out, "/", "ricover_wes_gistic_all_arm.RDS", sep = ""))
	saveRDS(eset.or, file = paste(f.dir.out, "/", "ricover_wes_gistic_all_or.RDS", sep = ""))

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
	saveRDS(eset.focal, file = paste(f.dir.out, "/", "ricover_wes_gistic_all_focal_binary.RDS", sep = ""))
	saveRDS(eset.arm, file = paste(f.dir.out, "/", "ricover_wes_gistic_all_arm_binary.RDS", sep = ""))
	saveRDS(eset.or, file = paste(f.dir.out, "/", "ricover_wes_gistic_all_or_binary.RDS", sep = ""))

	##get cis genes

	eset.focal<-readRDS(paste(f.dir.out, "/", "ricover_wes_gistic_all_focal_binary.RDS", sep = ""))

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

	saveRDS(genes.both, file = paste(f.dir.out, "/", "ricover_wes_gistic_cisgenes.RDS", sep = ""))

	#cis.genes<-readRDS(paste(f.dir.out, "/", "ricover_wes_gistic_cisgenes.RDS", sep = ""))
	#f.tab.amp<-paste(f.dir.in, "/", "table_amp.conf_99.txt", sep = "")
	#f.tab.del<-paste(f.dir.in, "/", "table_del.conf_99.txt", sep = "")
	#region.genes<-get_region_genes(cis.genes, f.tab.amp = f.tab.amp, f.tab.del = f.tab.del)

	#cis.genes.region<-list(meta =cis.genes$meta, genes = region.genes)
	#saveRDS(cis.genes.region, file = paste(f.dir.out, "/", "ricover_wes_gistic_cisgenes_region.RDS", sep = ""))
}
