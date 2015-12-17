lymph<-readRDS("/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/cancercell2012/lymphoma2012_all.RDS")

ricover_wes<-readRDS("/Users/amyli/Desktop/git_projects/gistic2ge_results/G2G_results_7_17_15/wes_7_16_15/ricover_wes_all.RDS")

summarize_dataset<-function(x){
	print("gene expression")
	print(dim( x$gep))
	print("cn focal")
	print(dim(x$cn.focal))

}

summarize_dataset(lymph)
summarize_dataset(ricover_wes)