query<-c("AHR", "CYP1B1")




get_corr<-function(mat.query, mat.full){
	n1<-nrow(mat.query)
	n2<-nrow(mat.full)


	res<-matrix(NA, nrow = n1, ncol = n2)
	for(i in 1:n1){
		for(j in 1:n2){
			x<-as.numeric(mat.query[i,])
			y<-as.numeric(mat.full[j,])
			res[i,j]<-cor(x,y)
		}
	}
	rownames(res)<-rownames(mat.query)
	colnames(res)<-rownames(mat.full)
	res<-data.frame(t(res))
	return(res)
}

select_corr<-function(res, thres, dir = c("pos", "neg", "bidirection")){
	if (dir == "pos"){
		cols<-sapply(1:nrow(res), 
		function(x){
			i<-as.numeric(res[x, ])
			return(max(i)> thres)
			})
	} else if (dir == "neg"){
		cols<-sapply(1:nrow(res), 
		function(x){
			i<-as.numeric(res[x, ])
			return(min(i)< (-thres))
			})
	} else{
		cols<-sapply(1:nrow(res), 
		function(x){
			i<-as.numeric(res[x, ])
			return(max(i)> thres | min(i)< (-thres))
			})		
	}

	return(res[cols,])
}


eset<-readRDS("../../datasets/TCGA_brca/TCGA_brca_eset_tn.RDS")
eset<-eset[,  eset$tumor.status == "Tumor"]
mat.full<-log2(exprs(eset) + 1)
mat.sub<-mat.full[query,]

mat.res<-get_corr(mat.sub, mat.full)
mat.res.sig<-select_corr(res =mat.res, thres = 0.4, dir = "bidirection")









