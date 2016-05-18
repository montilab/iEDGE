library(XML)

###
# scripts for boxplot
##

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

make_boxplotdataset<-function(cn,
	gep,
	testtype,
	alterationid,
	gene,
	tabid = "Unique.Name",
	geneid = "accession",
	loggep = TRUE
	){

	library(Biobase)
	library(ggplot2)
	altstatus_rowid<-which(fData(cn)[, tabid] == alterationid)
	if (length(altstatus_rowid) == 0){
		return(NULL)
	}
	altstatus<-as.numeric(exprs(cn)[altstatus_rowid,])
	altstatus<-as.factor(altstatus)
	gep_rowid<-which(fData(gep)[, geneid] == gene)
	if (length(gep_rowid) == 0){
		return(NULL)
	}

	emat<-exprs(gep)
	ematrow<-as.numeric(emat[gep_rowid,])

	if (loggep){
		emat<-emat+(-1*min(emat))
		ematrow<-as.numeric(emat[gep_rowid,])
		ematrow<-log2(ematrow+1)
		logind <- "(log2 transformed)"
	}
	df<-data.frame(alteration = altstatus, 
		gene_expression = round(ematrow, 2), 
		name =  colnames(gep))
	df$name<-paste("\"", df$name, "\"", sep = "")
	return(df)
}



##turns data frame into html code 
##appends js headers
to_table_html<-function(x, 
	x.name, 
	header = "", ## addition js header
	jsdir = "../../inst/javascript",
	hidecol = FALSE
	){
	library(knitr)
	##add default js header
	x.name.escape<-gsub("\\.", "\\\\.", x.name)

	#made with https://www.datatables.net/download/
	dtcss<-"https://cdn.datatables.net/s/dt/jszip-2.5.0,pdfmake-0.1.18,dt-1.10.10,b-1.1.0,b-colvis-1.1.0,b-flash-1.1.0,b-html5-1.1.0,b-print-1.1.0,cr-1.3.0,fc-3.2.0,fh-3.1.0,r-2.0.0,sc-1.4.0,se-1.1.0/datatables.min.css"
 	dtjs<-"https://cdn.datatables.net/s/dt/jszip-2.5.0,pdfmake-0.1.18,dt-1.10.10,b-1.1.0,b-colvis-1.1.0,b-flash-1.1.0,b-html5-1.1.0,b-print-1.1.0,cr-1.3.0,fc-3.2.0,fh-3.1.0,r-2.0.0,sc-1.4.0,se-1.1.0/datatables.min.js"

 	if(hidecol == TRUE){
	js_header<-paste("<head>
	<script src=\"http://ajax.aspnetcdn.com/ajax/jQuery/jquery-1.12.0.min.js\"></script>
	</head>
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

    		        columnDefs: [
            {
                targets: [ 0 ],
                visible: false,
                searchable: false
            }
        ],	
    		iDisplayLength: -1 }",
     	");} );	
     	</script>", 
		sep = "")
	} else {
		js_header<-paste("<head>
	<script src=\"http://ajax.aspnetcdn.com/ajax/jQuery/jquery-1.12.0.min.js\"></script>
	</head>
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
	}

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
 header = "by_alteration_cis",
 boxplot_link = "./boxplots", 
 hidecol
 ){

	alts<-unique(as.character(tab[, tabid]))
	in_table<-lapply(alts, 
		function(x){
			tabsub<-tab[tab[, tabid] == x,]
			tabsub<-as.data.frame(tabsub)
			row.names(tabsub)<-NULL
			return(tabsub)
		})	

	names(in_table)<-paste(header, "_", alts, sep = "")

	##addboxplots js reference

	boxplot_link<-paste("<script> var dir_out = \"", boxplot_link, "\"; </script>",sep = "")

	in_table_headers1<-"<script type=\"text/javascript\" charset=\"utf8\" src=\"addboxplot.js\"></script>"
	in_table_headers2<-"<script type=\"text/javascript\" charset=\"utf8\" src=\"addgenecard.js\"></script>"
	in_table_headers<-paste(boxplot_link,in_table_headers1, in_table_headers2, sep = "\n")
	
	in_table_headers<-rep(in_table_headers, length(in_table))

 	jsdir <- "../../inst/javascript"
 	jsdir <-rep(jsdir, length(in_table))

	##write by-alteration DE tables to html
	res<-mapply(to_table_html, in_table, names(in_table), in_table_headers, jsdir, hidecol)
	out_dir<-"../test2/html"
	sapply(names(res), 
		function(x){
		write(res[[x]], file = paste(outdir, "/", x, ".html", sep = ""))
		})
	
}

write_byalt_boxplot<-function(tab,#data frame to write
 tabid = "Unique.Name", #column name to split table by
 geneid = "accession", #column name for gene id
 outdir,
 header = "byalteration_cis",
 cn,
 gep, 
 basedir
 ){

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
			h1<-paste(readLines(paste(basedir, "/boxheader1.html", sep = "")), 
				collapse = "\n")
			h2<-paste(readLines(paste(basedir, "/boxheader2.html", sep = "")), 
				collapse = "\n")
			body<-data2html(df)
			boxhtml<-paste(h1, body, h2, sep = "\n")
			gene<-gsub("-|/", ".", gene)
			out.name<- paste(outdir, "/", header, "_", alterationid, "_", gene, ".html", sep = "")
			write(boxhtml, file = out.name)
		}
	}
}


###
# scripts for bipartite html
##

make_bipartite_html<-function(f.dir.in, f.dir.out, header = "", headeradd = "", 
	template){

	dir.create(f.dir.out, recursive = TRUE)
	#template <- "../inst/javascript/template_bipartite.html"
	f.files.in<-list.files(f.dir.in)
	for(i in f.files.in){
		#copy json files to out directory
		i.base<-gsub(paste(header, "_", sep = ""), "", i)
		if(grepl(".js", i.base)){
			i.ext<-".js"
		} else{
			i.ext<-".html"
		}

		i.basename<-gsub(i.ext, "", i.base)
		file.copy(from=paste(f.dir.in, "/",i, sep =""), to= paste(f.dir.out,"/", i.basename, headeradd, i.ext, sep = "") , 
	          overwrite = TRUE, recursive = FALSE, 
	          copy.mode = TRUE)
		doc.html = htmlTreeParse(template,
		           useInternal = TRUE)
		temp.src = "temp.js"
		new.src = paste(i.basename, headeradd, i.ext,sep = "")
		#replace src file name of json data
		nodes<-getNodeSet(doc.html, "//script")

		lapply(nodes, function(n) {	
			nodesrc = xmlAttrs(n)
			if("src" %in% names(nodesrc)){
				if(nodesrc["src"] == temp.src){		
					xmlAttrs(n)["src"]<-new.src
				}
			}
			return(n);
		})

		f.out  = paste(f.dir.out, "/", gsub(".js","",i.base), headeradd, ".html", sep = "")
		print(f.out)
		#save to html for each JSON input file
		
		saveXML(doc.html, file = f.out)
	}
}

make_iEDGE_ui<-function(cistab, transtab, cn, gep, cisgenes,
	outdir, jsdir, cmijsdir, altid = "Unique.Name", geneid = "accession"){

	dir.create(outdir)

	jsfiles<-list.files(jsdir, full.names = TRUE)

	for(js in jsfiles){
		file.copy(from=js, to=outdir, 
		          overwrite = TRUE, recursive = FALSE, 
		          copy.mode = TRUE)
	}

	cat("Writing by alteration cis DE tables...\n")
	##writing by alteration DE tables
	if(nrow(cistab)>0)
	write_byalt_html(cistab,#data frame to write
	 tabid = altid, #column name to split table by
	 outdir = outdir, #output directory name
	 header = "byalteration_cis", 
	 hidecol = FALSE
	 )

	cat("Writing by alteration trans DE tables...\n")
	if(nrow(transtab)>0)
	write_byalt_html(transtab,#data frame to write
	 tabid = altid, #column name to split table by
	 outdir = outdir, #output directory name
	 header = "byalteration_trans",
	 hidecol = FALSE
	 )
	##writing by alteration and by gene boxplots

	cat("Writing cis boxplots...\n")
	outdirfig<-paste(outdir, "/", "boxplots", sep = "")
	dir.create(outdirfig)
	write_byalt_boxplot(cistab,#data frame to write
	 tabid = altid, #column name to split table by
	 geneid = geneid, #column name for gene id
	 outdir = outdirfig,
	 header = "byalteration_cis",
	 cn = cn,
	 gep = gep, 
	 basedir = outdir
	 )

	cat("Writing trans boxplots...\n")
	write_byalt_boxplot(transtab,#data frame to write
	 tabid = altid, #column name to split table by
	 geneid = geneid, #column name for gene id
	 outdir = outdirfig,
	 header = "byalteration_trans",
	 cn = cn,
	 gep = gep, 
	 basedir = outdir
	 )


	cat("Writing summary table...\n")
	#alterations<-unique(as.character(cistab[, altid]))
	alterations<-unique(as.character(fData(cn)[, altid]))
	summarytab<-lapply(alterations, function(x){

		numcis <- length(which(cistab[, altid] == x))
		numtrans <- length(which(transtab[, altid] == x))
		
		ind.cis<-which(fData(cn)[,altid] == x)
		cytoband<-fData(cn)$Descriptor[ind.cis]
		cislist <- cisgenes[[ind.cis]]
		cislist <- paste(cislist, collapse = ",")
		numAlt<-length(which(exprs(cn)[ind.cis,] == 1))
		numNormal<-length(which(exprs(cn)[ind.cis,] == 0))

		numbipartitecis<-0
		numbipartitetrans<-0

		if(x %in% names(cmi$sig)){
			numbipartitecis<-length(unique(cmi$sig[[x]]$cis))
			numbipartitetrans<-length(unique(cmi$sig[[x]]$trans))
		}

		return(data.frame(alteration_id = x, 
			cytoband = cytoband,
			cis = numcis, 
			trans = numtrans, 
			bipartite= paste("(",numbipartitecis, "/", numbipartitetrans,")", sep = ""),
			samples_altered = numAlt,
			samples_normal = numNormal,
			genes_in_alteration = cislist)
		)
		})

	summarytab<-do.call(rbind, summarytab)
	summarytab<-data.frame(index = 1:nrow(summarytab), summarytab)

	row.names(summarytab)<-NULL
	#jsdir<-"../inst/javascript"

	addlinksheader <-paste("<script type=\"text/javascript\" charset=\"utf8\" src=\"addlinks.js\"></script>", sep = "")

	#add links
	cat("Adding links...\n")
	summarytable.html<-to_table_html(x = summarytab, 
		x.name = "summarytable", 
		header = addlinksheader,
		jsdir = jsdir,
		hidecol = TRUE)
	write(summarytable.html, file = paste(outdir, "/", "index.html", sep = ""))

	##writing bipartite
	cat("Writing bipartite graphs...\n")
	make_bipartite_html(f.dir.in = cmijsdir, 
		f.dir.out = paste(outdir, "/bipartiteplots", sep = ""), 
		header = "", headeradd = "", 
		template = paste(outdir, "/", "template_bipartite.html", sep = ""))

}



