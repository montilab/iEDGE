##to be called after make_iEDGE to display interact graphics results


##turns data frame into html code 
##appends js headers
to_table_html<-function(x, 
	x.name, 
	header = "", ## addition js header
	jsdir = "../../inst/javascript"
	){
	library(knitr)
	##add default js header
	x.name.escape<-gsub("\\.", "\\\\.", x.name)

	#made with https://www.datatables.net/download/
	dtcss<-"https://cdn.datatables.net/s/dt/jszip-2.5.0,pdfmake-0.1.18,dt-1.10.10,b-1.1.0,b-colvis-1.1.0,b-flash-1.1.0,b-html5-1.1.0,b-print-1.1.0,cr-1.3.0,fc-3.2.0,fh-3.1.0,r-2.0.0,sc-1.4.0,se-1.1.0/datatables.min.css"
 	dtjs<-"https://cdn.datatables.net/s/dt/jszip-2.5.0,pdfmake-0.1.18,dt-1.10.10,b-1.1.0,b-colvis-1.1.0,b-flash-1.1.0,b-html5-1.1.0,b-print-1.1.0,cr-1.3.0,fc-3.2.0,fh-3.1.0,r-2.0.0,sc-1.4.0,se-1.1.0/datatables.min.js"

	js_header<-paste("<script type=\"text/javascript\" charset=\"utf8\" src=\"", 
		jsdir, 
		"/jquery/dist/jquery.js\"></script>
		<link rel=\"stylesheet\" type=\"text/css\" href=\"", 
		dtcss,
		"\">
		<script type=\"text/javascript\" charset=\"utf8\" src=\"",
		dtjs,
		"\"></script>
		<script type=\"text/javascript\">
    	$(document).ready(function() {
        	$('#", x.name.escape, "').DataTable(", 
        	"{
    		aLengthMenu: [
        	[25, 50, 100, 200, -1],
        	[25, 50, 100, 200, \"All\"]
    		],
    		dom: 'lBfrtip',
    		fixedHeader: true,
    		buttons: [ 'copy',
    			{ extend: 'excelHtml5',
    			title: '", x.name, "',
    			extension: '.xlsx'},
  				{ extend: 'csvHtml5',
    			title: '", x.name, "'},
    			{ extend: 'pdfHtml5',
    			title: '", x.name, "'
    			}],
    		iDisplayLength: -1 }",
     	");} );	
     	</script>", 
		sep = "")

	js_header2<-header
	attrstr<-paste("id=\"", x.name, "\"", " class=\"display\" cellspacing=\"0\" width=\"60%\"", sep = "")
	res<-kable(x, "html", table.attr = attrstr,
	caption = x.name)
	res<-paste(js_header, js_header2, res, sep = "\n")
	return(res)
}

write_byalt_html<-function(tab,#data frame to write
 tabid = "Unique.Name", #column name to split table by
 outdir, #output directory name
 header = "by_alteration_cis"
 ){

	alts<-unique(as.character(tab[, altid]))
	in_table<-lapply(alts, 
		function(x){
			tabsub<-tab[tab[, tabid] == x,]
			tabsub<-as.data.frame(tabsub)
			row.names(tabsub)<-NULL
			return(tabsub)
		})	
	names(in_table)<-paste(header, "_", alts, sep = "")

	##addboxplots js reference
	in_table_headers1<-"<script type=\"text/javascript\" charset=\"utf8\" src=\"../../inst/javascript/addboxplot.js\"></script>"
	in_table_headers2<-"<script type=\"text/javascript\" charset=\"utf8\" src=\"../../inst/javascript/addgenecard.js\"></script>"
	in_table_headers<-paste(in_table_headers1, in_table_headers2, sep = "\n")
	

	in_table_headers<-rep(in_table_headers, length(in_table))

 	jsdir <- "../../inst/javascript"
 	jsdir <-rep(jsdir, length(in_table))

	##write by-alteration DE tables to html
	res<-mapply(to_table_html, in_table, names(in_table), in_table_headers, jsdir)
	out_dir<-"../test2/html"
	sapply(names(res), 
		function(x){
		write(res[[x]], file = paste(outdir, "/", x, ".html", sep = ""))
		})

	cat("\nDone\n")
	
}


write_byalt_boxplot<-function(tab,#data frame to write
 tabid = "Unique.Name", #column name to split table by
 geneid = "accession", #column name for gene id
 outdir,
 header = "byalteration_cis",
 cn,
 gep
 ){

 	source("boxplot_onfly.R")
 	if (!file.exists(outdir)){
		dir.create(outdir)
	}
	for(j.row in 1:nrow(tab)){

		alterationid<-as.character(tab[, tabid][j.row])
		gene<-as.character(tab[, geneid][j.row])
		df<-make_boxplotdataset(cn = cn,
		gep = gep,
		testtype = testtype,
		alterationid = alterationid,
		gene = gene, 
		tabid = tabid,
		geneid = geneid,
		loggep = FALSE)
		if (!is.null(df)){
			h1<-paste(readLines("../inst/html/boxheader1.html"), 
				collapse = "\n")
			h2<-paste(readLines("../inst/html/boxheader2.html"), 
				collapse = "\n")
			body<-data2html(df)
			boxhtml<-paste(h1, body, h2, sep = "\n")
			gene<-gsub("-|/", ".", gene)
			out.name<- paste(outdir, "/", header, "_", alterationid, "_", gene, ".html", sep = "")
			write(boxhtml, file = out.name)
		}
	}
}



indir<-"../test2/tables"
header<-"cancercell2012"
cisfdr<-0.25
transfdr<-0.05
outdir<-"../test2/html"

cistab<-read.table(paste(indir, "/", header, "_cis_sig_fdr_",cisfdr, ".txt",
 sep = ""), header= TRUE)

transtab<-read.table(paste(indir, "/", header, "_trans_sig_fdr_", transfdr, ".txt",
 sep = ""), header= TRUE)

altid<-"Unique.Name"

##writing by alteration DE tables
write_byalt_html(cistab,#data frame to write
 tabid = altid, #column name to split table by
 outdir = outdir, #output directory name
 header = "byalteration_cis"
 )

write_byalt_html(transtab,#data frame to write
 tabid = altid, #column name to split table by
 outdir = outdir, #output directory name
 header = "byalteration_trans"
 )
##writing by alteration and by gene boxplots

lymph<-readRDS("../data/lymph2012/cancercell2012_all.RDS")
write_byalt_boxplot(cistab,#data frame to write
 tabid = altid, #column name to split table by
 geneid = "accession", #column name for gene id
 outdir = "../test2/figures",
 header = "byalteration_cis",
 cn = lymph[["cn.focal.binary"]],
 gep = lymph[["gep.log.corrected"]]
 )
write_byalt_boxplot(transtab,#data frame to write
 tabid = altid, #column name to split table by
 geneid = "accession", #column name for gene id
 outdir = "../test2/figures",
 header = "byalteration_trans",
 cn = lymph[["cn.focal.binary"]],
 gep = lymph[["gep.log.corrected"]]
 )


alterations<-unique(as.character(cistab[, altid]))

summarytab<-lapply(alterations, function(x){
	numcis <- length(which(cistab[, tabid] == x))
	numtrans <- length(which(transtab[, tabid] == x))

	return(data.frame(alteration_id = x, 
		cis = numcis, 
		trans = numtrans))
	})

summarytab<-do.call(rbind, summarytab)

row.names(summarytab)<-NULL
jsdir<-"../../inst/javascript"
addlinksheader <-paste("<script type=\"text/javascript\" charset=\"utf8\" src=\"",
jsdir,"/addlinks.js\"></script>", sep = "")


#add links
summarytable.html<-to_table_html(x = summarytab, 
	x.name = "summarytable", 
	header = addlinksheader,
	jsdir = jsdir)

write(summarytable.html, file = paste(outdir, "/", "summarytable", ".html", sep = ""))





