
#' \code{util.save.rds} save to .RDS, and optionally .xls (for list of data frames)
#' @param x R object to save
#' @param f.dir directory to save to, default = NA, current working directory
#' @param f.header file header
#' @param f.name file name (appended to file header)
#' @export
util.save.rds<-function(x, f.dir = NA, f.header = "ricover2010_", f.name){
	if (is.na(f.dir)) saveRDS(x, file = f.name)
	else saveRDS(x, file = paste(f.dir, "/", f.header, f.name, ".RDS", sep = ""))
}


#' \code{util.save.xlsx} save to .xlsx  
#' @import XLConnect
#' @param x R object to save, must be a named list of data.frames, row limit =65535 for each data frame
#' @param f.dir directory to save to, default = NA, current working directory
#' @param f.header file header
#' @param f.name file name (appended to file header)
#' @export
util.save.xlsx<-function(x, f.dir = NA, f.header = "ricover2010_", f.name){
	options(java.parameters = "-Xmx3g" )
	require(XLConnect)

	# Load workbook (create if not existing)
	if(is.na(f.dir)) 
		wb.name<-paste(f.name, ".xlsx", sep = "")
	else
		wb.name<-paste(f.dir, "/", f.header, f.name, ".xlsx", sep = "")
	
	if (file.exists(wb.name)){
		file.remove(wb.name)
	}
	
	wb <- loadWorkbook(wb.name, create = TRUE)
	
	for (i in names(x)){ 
		createSheet(wb, name = i)
		writeWorksheet(wb, x[[i]], sheet = i, startRow = 1, startCol = 1)
	}
	saveWorkbook(wb)
}

