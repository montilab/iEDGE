#make_boxplot.R
library(CBMRtools)
library(ggplot2)

find_gene_loc<-function(x, gene, is.arm = TRUE){
	if(is.arm == TRUE){
		ind<-which(unlist(lapply(x$gep.cisgenes.arm$genes, 
		function(x) gene %in% x)))
		desc<-fData(x$cn.arm)[ind, c("Descriptor")]
		dir<-fData(x$cn.arm)[ind, c("Unique.Name")]

		#if (grepl("Amp", dir)) dir <- "Amplification"
		#else dir <- "Deletion"
	} else {
		ind<-which(unlist(lapply(x$gep.cisgenes$genes, 
		function(x) gene %in% x)))
		desc<-fData(x$cn.focal)[ind, c("Descriptor")]
		dir<-fData(x$cn.focal)[ind, c("Unique.Name")]	
	}

	dir_short<-as.character(sapply(dir, 
			function(x){
				if (grepl("Amp", x)) return("Amplification")
				else return("Deletion")
				}))

	return(list(gene = gene,
		ind = ind, 
		desc = desc, 
		dir = dir_short))
}

make_boxplot_new<-function(x, loc, opt = c("focal", "arm", "focal_or_arm")){

	gene<-loc$gene
	desc <- loc$desc
	dir <-loc$dir
	ind <- loc$ind
	gep<-x$gep.log.corrected
	print(gene)
	print(desc)
	print(dir)
	print(ind)
	if (opt == "arm") cn<-x$cn.arm
	else if (opt == "focal") cn <- x$cn.focal
	else if (opt == "focal_or_arm") cn <-x$cn.focal.or.arm
	else stop("invalid opts, must be one of focal, arm, focal_or_arm")

	if(grepl("Amp", fData(cn)[ind,"Unique.Name"])){
		cn.dir<-"Amplification"
		cn.dir.short<-"gain"
	}
	else if(grepl("Del", fData(cn)[ind,"Unique.Name"])){
		cn.dir<-"Deletion"
		cn.dir.short <-"loss"
	}
	else {
		stop("direction not found")
	}
	cn.values<-as.factor(exprs(cn[ind,]))
	gep.values<-as.numeric(exprs(gep)[fData(gep)$accession == gene,])
	df<-data.frame(alteration = cn.values, gene_expression = gep.values)
	allmin<-min(df$gene_expression)
	give.n <- function(x){
   		return(c(y = allmin, label = length(x)))
	}
	p1<-ggplot(df, aes(x = alteration, y = gene_expression)) + 
	geom_boxplot() + 
	#facet_wrap(  ~ copy_number, scale = "free_y") + 
	stat_summary(fun.data = give.n, geom = "text") + 
	xlab(paste("Copy number alteration status\n (0 = not altered, 1 = single copy ", 
		cn.dir.short, ", 2 = double copy ", cn.dir.short, ")", sep  = "")) + 
	ylab("Gene Expression") + 
	ggtitle(paste("Copy number alteration", "(", desc, ",", dir, ")", 
		"\n vs gene expression", "(", gene,")", "\n",
		opt,
		sep = ""))
}
###do for cancercell2012

f.dir.in<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/cancercell2012"
f.dir.out<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/cancercell2012/boxplots"

x<-readRDS( file = paste(f.dir.in,"/", "cancercell2012_all.RDS", sep = ""))

loc_TP53_arm<-find_gene_loc(x, gene = "TP53", is.arm = TRUE)
loc_TP53_arm<-lapply(loc_TP53_arm, function(x) {
	n<-length(x)
	x[n]
	})

loc_TP53_focal<-find_gene_loc(x, gene = "TP53", is.arm = FALSE)

loc_BCL2_arm<-find_gene_loc(x, gene = "BCL2", is.arm = TRUE)
loc_BCL2_arm<-lapply(loc_BCL2_arm, function(x) {
	n<-length(x)
	x[1]
	})
#loc_BCL2_focal<-find_gene_loc(x, gene = "BCL2", is.arm = FALSE)
loc_BCL2_focal<-list()
loc_BCL2_focal$gene<-"BCL2"
loc_BCL2_focal$ind<-30
loc_BCL2_focal$desc<- "18q21.32"
loc_BCL2_focal$dir<-"Amplification"

loc_BCL2<-list()

loc_BCL2[["focal"]]<-loc_BCL2_focal
loc_BCL2[["arm"]]<-loc_BCL2_arm
loc_BCL2[["focal_or_arm"]]<-loc_BCL2_focal
p_BCL2<-list()

for(i in c("focal", "arm", "focal_or_arm")){
	p_BCL2[[i]]<-make_boxplot_new(x = x, loc = loc_BCL2[[i]], opt = i)
}

loc_TP53<-list()

loc_TP53[["focal"]]<-loc_TP53_focal
loc_TP53[["arm"]]<-loc_TP53_arm
loc_TP53[["focal_or_arm"]]<-loc_TP53_focal
p_TP53<-list()

for(i in c("focal", "arm", "focal_or_arm")){
	p_TP53[[i]]<-make_boxplot_new(x = x, loc = loc_TP53[[i]], opt = i)
	#ggsave(p_TP53[[i]], file = paste(f.dir.out, "/", "TP53_",opt,".png", sep = ""))

}

##do for ricover
f.dir.in<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/wes_7_16_15"
f.dir.out<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/wes_7_16_15/boxplots"

x<-readRDS( file = paste(f.dir.in,"/", "ricover_wes_all.RDS", sep = ""))

loc_TP53_arm<-find_gene_loc(x, gene = "TP53", is.arm = TRUE)
loc_TP53_arm<-lapply(loc_TP53_arm, function(x) {
	n<-length(x)
	x[n]
	})

loc_TP53_focal<-find_gene_loc(x, gene = "TP53", is.arm = FALSE)

loc_BCL2_arm<-find_gene_loc(x, gene = "BCL2", is.arm = TRUE)
loc_BCL2_arm<-lapply(loc_BCL2_arm, function(x) {
	n<-length(x)
	x[1]
	})
#loc_BCL2_focal<-find_gene_loc(x, gene = "BCL2", is.arm = FALSE)
loc_BCL2_focal<-list()


loc_BCL2_focal$gene<-"BCL2"
loc_BCL2_focal$ind<-30
loc_BCL2_focal$desc<- "18q21.32"
loc_BCL2_focal$dir<-"Amplification"

loc_BCL2<-list()

loc_BCL2[["focal"]]<-loc_BCL2_focal
loc_BCL2[["arm"]]<-loc_BCL2_arm
loc_BCL2[["focal_or_arm"]]<-loc_BCL2_focal
p_BCL2<-list()

for(i in c("focal", "arm", "focal_or_arm")){
	p_BCL2[[i]]<-make_boxplot_new(x = x, loc = loc_BCL2[[i]], opt = i)
}

loc_TP53<-list()

loc_TP53[["focal"]]<-loc_TP53_focal
loc_TP53[["arm"]]<-loc_TP53_arm
loc_TP53[["focal_or_arm"]]<-loc_TP53_focal
p_TP53<-list()

for(i in c("focal", "arm", "focal_or_arm")){
	p_TP53[[i]]<-make_boxplot_new(x = x, loc = loc_TP53[[i]], opt = i)
	ggsave(p_TP53[[i]], file = paste(f.dir.out, "/", "TP53_",opt,".png", sep = ""))

}

make_boxplot<-function(gene, cn, gep, arm = FALSE, cn.focal = NA, wes){
	ind.cn<-which(sapply(wes$gep.cisgenes$genes, function(x) gene %in% x))
	#print(ind.cn)
	if(length(ind.cn) != 1)
		stop("gene not found, or duplicated alteration mapping")
	
	if (arm == TRUE){
		#cn.focal must be defined
		#cn.focal.desc<-fData(cn.focal)$Descriptor[ind.cn]
		#cn.focal.uniquename<-substr(fData(cn.focal)$Unique.Name[ind.cn], 1, 5)
		cn.focal.desc<-fData(cn.focal)$Descriptor
		cn.focal.uniquename<-substr(fData(cn.focal)$Unique.Name, 1, 5)
		ind1<-sapply(fData(cn)$Descriptor, function(x) grepl(x, cn.focal.desc))
		ind2<-sapply(substr(fData(cn)$Unique.Name, 1, 4), function(x) grepl(x, cn.focal.uniquename))
		ind.cn<-which(ind1 & ind2)
	}

	cn.desc<-fData(cn)[ind.cn,"Descriptor"]
	if(grepl("Amp", fData(cn)[ind.cn,"Unique.Name"])){
		cn.dir<-"Amplification"
		cn.dir.short<-"gain"
	}
	else if(grepl("Del", fData(cn)[ind.cn,"Unique.Name"])){
		cn.dir<-"Deletion"
		cn.dir.short <-"loss"
	}
	else {
		stop("direction not found")
	}
	cn.values<-as.factor(exprs(cn[ind.cn,]))
	gep.values<-as.numeric(exprs(gep)[fData(gep)$accession == gene,])
	df<-data.frame(alteration = cn.values, gene_expression = gep.values)
	allmin<-min(df$gene_expression)
	give.n <- function(x){
   		return(c(y = allmin, label = length(x)))
	}
	p1<-ggplot(df, aes(x = alteration, y = gene_expression)) + 
	geom_boxplot() + 
	#facet_wrap(  ~ copy_number, scale = "free_y") + 
	stat_summary(fun.data = give.n, geom = "text") + 
	xlab(paste("Copy number alteration status\n (0 = not altered, 1 = single copy ", 
		cn.dir.short, ", 2 = double copy ", cn.dir.short, ")", sep  = "")) + 
	ylab("Gene Expression") + 
	ggtitle(paste("Copy number alteration", "(", cn.desc, ",", cn.dir, ")", 
		"vs gene expression", "(", gene,")", sep = ""))
	return(list(p1 = p1, ind.cn = ind.cn))
}



make_boxplot_cytoband<-function(cytoband, gene, cn, gep){
	ind.cns<-which(fData(cn)$Descriptor == cytoband)
	print(length(ind.cns))
	#print(ind.cn)
	#if(length(ind.cn) != 1){
	#	stop("gene not found, or duplicated alteration mapping")
	#}
	p<-list()
	for(ind.cn in ind.cns){
	cn.desc<-fData(cn)[ind.cn,"Descriptor"]
	if(grepl("Amp", fData(cn)[ind.cn,"Unique.Name"])){
		cn.dir<-"Amplification"
		cn.dir.short<-"gain"
	}
	else if(grepl("Del", fData(cn)[ind.cn,"Unique.Name"])){
		cn.dir<-"Deletion"
		cn.dir.short <-"loss"
	}
	else {
		stop("direction not found")
	}
	cn.values<-as.factor(exprs(cn[ind.cn,]))
	gep.values<-as.numeric(exprs(gep)[fData(gep)$accession == gene,])
	df<-data.frame(alteration = cn.values, gene_expression = gep.values)
	allmin<-min(df$gene_expression)
	give.n <- function(x){
   		return(c(y = allmin, label = length(x)))
	}
	p1<-ggplot(df, aes(x = alteration, y = gene_expression)) + 
	geom_boxplot() + 
	#facet_wrap(  ~ copy_number, scale = "free_y") + 
    stat_summary(fun.data = give.n, geom = "text") + 
	xlab(paste("Copy number alteration status\n (0 = not altered, 1 = single copy ", 
		cn.dir.short, ", 2 = double copy ", cn.dir.short, ")", sep  = "")) + 
	ylab("Gene Expression") + 
	ggtitle(paste("Copy number alteration", "(", cn.desc,",", cn.dir, ")", 
		"vs gene expression", "(", gene,")", sep = ""))
	p[[paste("ind",ind.cn, sep ="")]]<-p1
	}
	return(p)
}

make_boxplot_ricover_wes<-function(){
f.dir.in<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/wes_7_16_15"
f.dir.out<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/wes_7_16_15/boxplots"

wes<-readRDS( file = paste(f.dir.in,"/", "ricover_wes_all.RDS", sep = ""))


#focal only
cat("BCL2\n")
cat("focal only\n")
gene<-"BCL2" #BCL2
cn<-wes$cn.focal
gep<-wes$gep.log.corrected
cytoband<-"18q21.32"
p1<-make_boxplot_cytoband(cytoband = cytoband, gene = gene, cn, gep)
cn.focal<-as.factor(as.character(exprs(cn)[which(fData(cn)$Descriptor == "18q21.32")[1],]))
ggsave(p1[[1]], file = paste(f.dir.out, "/", "BCL2_focal.png", sep = ""))

##arm only

cat("arm only\n")
cn<-wes$cn.arm
cytoband<-"18q"
p1<-make_boxplot_cytoband(cytoband = cytoband, gene = gene, cn, gep)
p1[[1]]
ggsave(p1[[1]], file = paste(f.dir.out, "/", "BCL2_arm.png", sep = ""))
cn.arm<-as.factor(as.character(exprs(cn)[which(fData(cn)$Descriptor == "18q")[1],]))

##table
cat("cont. table\n")
tab<-as.data.frame.matrix(table(cn.focal,cn.arm))
tab<-cbind(rownames(tab), tab)
tab<-rbind(c("", colnames(tab)[-1]), as.matrix(tab))
tab<-cbind(rep("", nrow(tab)), tab)
tab[2,1]<-"focal"
colnames(tab)<-rep("", ncol(tab))
colnames(tab)[3]<-"arm"

save.xlsx(x=list("Sheet1" = tab), f.dir = f.dir.out,  f.name = "BCL2_contingency_table")

##TP53
gene<-"TP53" 
cat("TP53\n")
cat("focal only\n")
#focal only
cn<-wes$cn.focal
gep<-wes$gep.log.corrected
p1<-make_boxplot(gene, cn, gep, wes=wes)
#ind.cn<-which(sapply(wes$gep.cisgenes$genes, function(x) gene %in% x))
#cn.focal<-as.factor(as.character(exprs(cn)[ind.cn,]))
ggsave(p1$p1, file = paste(f.dir.out, "/", "TP53_focal.png", sep = ""))

##arm only
cat("arm only\n")
cn<-wes$cn.arm
p1<-make_boxplot(gene, cn, gep, wes=wes)
ggsave(p1$p1, file = paste(f.dir.out, "/", "TP53_arm.png", sep = ""))
#cn.arm<-as.factor(as.character(exprs(cn)[ind.cn,]))

##focal or arm 
cat("focal or arm\n")
cn<-wes$cn.focal.or.arm
p1<-make_boxplot( gene, cn, gep, wes=wes)
ggsave(p1$p1, file = paste(f.dir.out, "/", "TP53_focal_or_arm.png", sep = ""))

cat("cont. table\n")
tab<-as.data.frame.matrix(table(cn.focal,cn.arm))
tab<-cbind(rownames(tab), tab)
tab<-rbind(c("", colnames(tab)[-1]), as.matrix(tab))
tab<-cbind(rep("", nrow(tab)), tab)
tab[2,1]<-"focal"
colnames(tab)<-rep("", ncol(tab))
colnames(tab)[3]<-"arm"
save.xlsx(x=list("Sheet1" = tab), f.dir = f.dir.out,  f.name = "TP53_contingency_table")
}



make_boxplot_lymphoma2012<-function(){
f.dir.in<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/cancercell2012"
f.dir.out<-"/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/cancercell2012/boxplots"

wes<-readRDS( file = paste(f.dir.in,"/", "cancercell2012_all.RDS", sep = ""))

#focal only
gene<-"BCL2" #BCL2
cn<-wes$cn.focal
gep<-wes$gep.log.corrected
p1<-make_boxplot(gene, cn, gep, wes =wes)
ggsave(p1$p1, file = paste(f.dir.out, "/", "BCL2_focal.png", sep = ""))
#ind.cn<-which(sapply(wes$gep.cisgenes$genes, function(x) gene %in% x))
#cn.focal<-as.factor(as.character(exprs(cn)[ind.cn,]))
counts.focal<-as.factor(as.character(exprs(cn)[p1$ind.cn,]))

##arm only
cn<-wes$cn.arm
p1<-make_boxplot(gene, cn, gep, arm = TRUE, cn.focal = wes$cn.focal, wes=wes)
ggsave(p1$p1, file = paste(f.dir.out, "/", "BCL2_arm.png", sep = ""))
#cn.arm<-as.factor(as.character(exprs(cn)[p1$ind.cn,]))
counts.arm<-as.factor(as.character(exprs(cn)[p1$ind.cn,]))

##focal or arm
cn<-wes$cn.focal.or.arm
p1<-make_boxplot( gene, cn, gep, wes = wes)
ggsave(p1$p1, file = paste(f.dir.out, "/", "BCL2_focal_or_arm.png", sep = ""))



##table
tab<-as.data.frame.matrix(table(counts.focal,counts.arm))
tab<-cbind(rownames(tab), tab)
tab<-rbind(c("", colnames(tab)[-1]), as.matrix(tab))
tab<-cbind(rep("", nrow(tab)), tab)
tab[2,1]<-"focal"
colnames(tab)<-rep("", ncol(tab))
colnames(tab)[3]<-"arm"

save.xlsx(x=list("Sheet1" = tab), f.dir = f.dir.out,  f.name = "BCL2_contingency_table")

##TP53
gene<-"TP53" 
#focal only
cn<-wes$cn.focal
gep<-wes$gep.log.corrected
p1<-make_boxplot(gene, cn, gep, wes=wes)
#ind.cn<-which(sapply(wes$gep.cisgenes$genes, function(x) gene %in% x))
#cn.focal<-as.factor(as.character(exprs(cn)[ind.cn,]))
ggsave(p1$p1, file = paste(f.dir.out, "/", "TP53_focal.png", sep = ""))
counts.focal<-as.factor(as.character(exprs(cn)[p1$ind.cn,]))


##arm only
cn<-wes$cn.arm
p1<-make_boxplot(gene, cn, gep, arm = TRUE, cn.focal = wes$cn.focal, wes = wes)
ggsave(p1$p1, file = paste(f.dir.out, "/", "TP53_arm.png", sep = ""))
#cn.arm<-as.factor(as.character(exprs(cn)[p1$ind.cn,]))
counts.arm<-as.factor(as.character(exprs(cn)[p1$ind.cn,]))

##focal or arm 
cn<-wes$cn.focal.or.arm
p1<-make_boxplot( gene, cn, gep, wes = wes)
ggsave(p1$p1, file = paste(f.dir.out, "/", "TP53_focal_or_arm.png", sep = ""))


tab<-as.data.frame.matrix(table(counts.focal,counts.arm))
tab<-cbind(rownames(tab), tab)
tab<-rbind(c("", colnames(tab)[-1]), as.matrix(tab))
tab<-cbind(rep("", nrow(tab)), tab)
tab[2,1]<-"focal"
colnames(tab)<-rep("", ncol(tab))
colnames(tab)[3]<-"arm"
save.xlsx(x=list("Sheet1" = tab), f.dir = f.dir.out,  f.name = "TP53_contingency_table")

}
if (TRUE){
cat("Making boxplot for ricover wes\n")
make_boxplot_ricover_wes()
cat("Making boxplot for lymphoma2012\n")
make_boxplot_lymphoma2012()
}