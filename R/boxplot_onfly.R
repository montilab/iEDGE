#input the dataset object (e.g lymph), 
#testtype (e.g. cn.focal.binary),
#alterationid (e.g. AmplificationPeak1), and gene (e.g.FCGR2B)
#plot boxplot of gene expression for each alteration group

data2html<-function(df){
	head<-"var data = ["
	df.names<-colnames(df)
	df.names.escape<-paste("\"", df.names, "\"", sep = "")
	
	appended.row<-apply(df, 1, function(x) {
		appended<-paste(df.names.escape, x, sep = ":")
		collapsed<-paste(appended, collapse = ",")
		})
	appended.row<-paste("{", appended.row, "}", sep = "")

	body<-paste(appended.row, collapse = ",")
	tail<-"]"
	return(paste(head, body, tail, sep = ""))
}

make_boxplotdataset<-function(
	ds,
	testtype,
	alterationid,
	gene,
	loggep = TRUE
	){

	library(Biobase)
	library(ggplot2)
	altstatus_rowid<-which(fData(ds[[testtype]])[, "Unique.Name"] == alterationid)
	if (length(altstatus_rowid) == 0){
		return(NULL)
	}
	altstatus<-as.numeric(exprs(ds[[testtype]])[altstatus_rowid,])
	altstatus<-as.factor(altstatus)
	gep_rowid<-which(fData(ds[["gep"]])[, "accession"] == gene)
	if (length(gep_rowid) == 0){
		return(NULL)
	}
	gep<-as.numeric(exprs(ds[["gep"]])[gep_rowid,])

	if (loggep){
		gep<-log2(gep+1)
		logind <- "(log2 transformed)"
	}
	df<-data.frame(alteration = altstatus, 
		gene_expression = round(gep, 2), 
		name =  colnames(ds[["gep"]]))
	df$name<-paste("\"", df$name, "\"", sep = "")
	return(df)
}

ds<-lymph
#testtype<-"cn.focal.binary"
#alterationid<-"AmplificationPeak1"
#gene<-"FCGR2B"

if (FALSE){
	lymph <-readRDS("../test/data/cancercell2012_all.RDS")
	f.dir.out <- "../test/figures"
	header <- "cancercell2012"
	if (!file.exists(f.dir.out)){
		dir.create(f.dir.out)
	}


	tab<-list.files("../test/tables")
	tab<-tab[grepl("^byalteration.*", tab)]


	
	for(i in tab){
		if(length(grep("^byalteration.*focal_.*", i)) == 1){
			testtype<-"cn.focal.binary"
		} else if (length(grep("^byalteration.*focalorarm.*", i)) == 1){
			testtype<-"cn.focal.or.arm.binary"
		}

		j <- read.table(paste("../test/tables/", i, sep = ""), 
			sep = "\t",header = T, quote = "")
		for(j.row in 1:nrow(j)){
			alterationid<-as.character(j$Unique.Name[j.row])
			gene<-as.character(j$gene_symbol[j.row])
			print(alterationid)
			print(gene)
			df<-make_boxplotdataset(ds = ds,
			testtype = testtype,
			alterationid = alterationid,
			gene = gene, 
			loggep = TRUE)
			if (!is.null(df)){
				h1<-paste(readLines("../inst/html/boxheader1.html"), 
					collapse = "\n")
				h2<-paste(readLines("../inst/html/boxheader2.html"), 
					collapse = "\n")
				body<-data2html(df)
				boxhtml<-paste(h1, body, h2, sep = "\n")
				gene<-gsub("-|/", ".", gene)
				out.name<- paste(gsub("\\.txt", "", i), "_", gene, ".html", sep = "")
			
				out.name<-paste(f.dir.out, "/",
				 gsub("\\.txt", "", i), "_", gene, ".html",sep = "")
				write(boxhtml, file = out.name)
			}
		}
	}


}