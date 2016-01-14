write_eset_json<-function(eset, fout){
	library(Biobase)
	library(reshape2)
	library(rjson)
	x<-exprs(eset)
	fdat<-fData(eset)
	pdat<-pData(eset)

	x<-as.data.frame(x)
	x<-melt(data.frame(x, id = rownames(x)), id.vars = "id")
	colnames(x)<-c("fdataid", "pdataid", "value")
	x.write<-toJSON(unname(split(x, 1:nrow(x))))

	fdat<-data.frame(fdat, fdataid = rownames(fdat))
	fdat.write<-toJSON(unname(split(fdat, 1:nrow(fdat))))
	
	pdat<-data.frame(pdat, pdataid = rownames(pdat))
	pdat.write<-toJSON(unname(split(pdat, 1:nrow(pdat))))

	x.write<-paste("var exprdata = ", x.write, ";", sep = "")
	fdat.write<-paste("var fdata = ", fdat.write, ";", sep = "")
	pdat.write<-paste("var pdata = ", pdat.write, ";", sep = "")

	all.write<-paste(x.write, fdat.write, pdat.write, sep = "\n")
	write(all.write, file = fout)
}


if (FALSE){
	#test case
	ds<-readRDS("../test/data/cancercell2012_all.RDS")
	eset<-ds$gep
	eset<-eset[1:40, 1:5]
	write_eset_json(eset, "test.js")
}