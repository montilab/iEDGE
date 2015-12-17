##same output format as GSA.read.gmt except without the annoying text outputs
#' \code{read_gmt} read from gmt file
#' @param f name of gmt file
#' @param split.char delimiter
read_gmt<-function(f, split.char = '\t'){

	con <- file(f, open = 'r') 

	results.list <- list()
	current.line <- 1
	while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
 		results.list[[current.line]] <- as.character(unlist(strsplit(line, split=split.char)))
 		current.line <- current.line + 1
	} 
	close(con)

	genesets<-sapply(results.list, function(x) x[-(1:2)])
	geneset.names<-unlist(sapply(results.list, function(x) x[1]))
	geneset.descriptions<-unlist(sapply(results.list, function(x) x[2]))
	names(genesets)<-geneset.names
	res<-list(genesets = genesets, geneset.names = geneset.names, geneset.descriptions = geneset.descriptions)
	return(res)
}