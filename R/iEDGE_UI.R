to.hex<-function(x){
	cols<-col2rgb(x)
	red<-cols[1]
	green<-cols[2]
	blue<-cols[3]
	return(rgb(red, green, blue, maxColorValue = 255))
}

#' @import RColorBrewer
#' @export
get_colors<-function(){
	colors1<-brewer.pal(8,"Set1")
	qual_col_pals<-brewer.pal.info[brewer.pal.info$category == 'qual',]
	colors2<-unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
	set.seed(1)
	colors3<-as.character(sapply(unique(sample(colors(),500)), to.hex))
	colorsall<-union(colors1, union(colors2, colors3))
	return(colorsall)
}

cis_order<-function(dat){
	cat("Performing cis ordering...\n")
	cis_in<-unique(as.character(dat[,1]))
	trans_in<-unique(as.character(dat[,2]))
	n_edges<-suppressWarnings(as.numeric(sapply(cis_in, function(i){
		sum(dat[,1] == i)
		})))

	res<-cis_in[order(n_edges, decreasing = TRUE)]
	return(res)
}

trans_order<-function(dat){
	cis_in<-unique(as.character(dat[,1]))
	trans_in<-unique(as.character(dat[,2]))	
	n<-length(trans_in)
	cat(paste("Performing trans ordering on ", n, " sets...\n", sep =""))
	if(n<3) return(trans_in)
	mat<-matrix(0,n,n)
	dat.2<-dat[,2]
	dat.1<-dat[,1]
	in_i<-list()

	#in-edges for each right node
	for(i in 1:n) {
		in_i[[i]]<-as.character(dat.1[dat.2 %in% trans_in[i]])
	}
	for(i in 1:n){
		for(j in 1:i){
			a<-intersect(in_i[[i]], in_i[[j]])
			b<-union(in_i[[i]], in_i[[j]])
			mat[i,j]<-length(a)/length(b)
			mat[j,i]<-mat[i,j]
		}
	}
	mat_dist<-1-mat
	mat_dist_obj<-as.dist(mat_dist)
	trans_ord<-hclust(mat_dist_obj, method = "complete")
	trans_ord<-order.dendrogram(as.dendrogram(trans_ord))
	res<-trans_in[order(trans_ord)]
	return(res)
}

write_bipartite_JSON<-function(tab, hyper, f.dir.out, header){
	suppressWarnings(dir.create(f.dir.out))
	ntab<-nrow(tab)

	if(ntab>=1){	
		colors<-get_colors()
		colors.max<-length(colors)

		ncis<-length(unique(tab$cis))
		ntrans<-length(unique(tab$trans))

		cis_ord<-cis_order(tab)
		trans_ord<-trans_order(tab)	
		cis_df<-data.frame(genes=cis_ord, 
			numEdges=sapply(cis_ord, 
				function(i){
					sum(tab$cis %in% i)
				}), 
			color = rep(colors, ceiling(ncis/colors.max))[1:ncis])
		
		trans_df<-data.frame(genes=trans_ord, 
			numEdges=sapply(trans_ord, 
				function(i){
					sum(tab$trans %in% i)
				}),
			color = rep(colors, ceiling(ntrans/colors.max))[1:ntrans])

		cis.write<-write_JSON_df(df = cis_df, df.name = "var cis")
		trans.write<-write_JSON_df(df = trans_df, df.name = "var trans")
		f.tab.ord<-tab[order(match(tab$trans, trans_ord)),]
		edges.write<-write_JSON_df(df = f.tab.ord, df.name = "var edges")

		if(hasArg(hyper)){

			hu<-lapplyn(hyper, function(i) get_hyper_wrapper(i[["hyperunion"]]))
			hu.gs<-lapplyn(hu, function(i) {
				res<-i[["hypergsets"]]
				if(is.null(res)) return(vector())
				else return(res)
				})
			hu.edges<-lapplyn(hu, function(i) i[["hyperedges"]])

			ha<-lapplyn(hyper, function(i){	
					lapplyn(i[["hyperbyalt"]], function(j)
						get_hyper_wrapper(j))					
				})			

			ha.gs<-lapplyn(ha, 
					function(i) lapplyn(i, function(k){
						res<-k[["hypergsets"]]
						if(is.null(res)) return(vector())
						else return(res)
					}))

			ha.edges<-lapplyn(ha, 
					function(i) lapplyn(i, function(k) k[["hyperedges"]]))

			hyper.union.hypergsets<-write_JSON_list(df = hu.gs, df.name = "hypergsets")
			hyper.union.hyperedges<-write_JSON_list(df = hu.edges, df.name = "hyperedges")
			hyper.union.hypergsets.cis<-write_JSON_list(df = ha.gs, df.name = "hypergsetscis")
			hyper.union.hyperedges.cis<-write_JSON_list(df = ha.edges, df.name = "hyperedgescis")

			all.write<-paste(cis.write, trans.write, edges.write, 
				hyper.union.hypergsets, hyper.union.hyperedges, 
				hyper.union.hypergsets.cis,
				hyper.union.hyperedges.cis,
				sep = "\n")
			write(all.write, file = paste(f.dir.out, "/", header, ".js", sep = ""), append = FALSE)
			
		} else {
			all.write<-paste(cis.write, trans.write, edges.write,
				sep = "\n")
			write(all.write, file = paste(f.dir.out, "/", header, ".js", sep = ""), append = FALSE)
		}
	}
}

#' @import jsonlite
write_JSON_list<-function(df, df.name){
	res.header<-paste(df.name, " = ", sep = "")
	res.body<-toJSON(df)
	res<-paste(res.header, res.body, ";", "\n", sep = "")
	return(res)
}

write_JSON_df<-function(df, df.name){
	if(is.null(df)) return("")
	res.header<-paste(df.name, " = ", sep = "")
	res<-paste(res.header, "[", paste("[", apply(df, 1, function(i){
				paste(paste("\"", i, "\"", sep =""), collapse =",")}), "]", collapse = ","), "]", sep = "")
	res<-paste(res, ";", "\n", sep = "")
	return(res)
}

get_hyper_edges<-function(df){
	edges<-do.call(rbind, lapply(1:nrow(df), function(x){
		transstr<-as.character(df$hits[x])
		trans<-as.character(gsub(" ", "", strsplit(transstr, split = ",")[[1]]))
		gsets<-as.character(df$category[x])
		return(data.frame(trans = trans, enrichsets =rep(gsets, length(trans))))
		}))
	#enrich_ord<-trans_order(edges)
	return(edges)
}

get_hyper_wrapper<-function(f.tab.hyper){
	colors<-get_colors()
	colors.max<-length(colors)
	if(nrow(f.tab.hyper) > 0){
		f.tab.hyper.edges<-get_hyper_edges(f.tab.hyper)
		nhyper<-length(unique(f.tab.hyper.edges$enrichsets))
		hyper_ord<-trans_order(f.tab.hyper.edges)
		hyper_df<-data.frame(genes=hyper_ord, 
		numEdges=sapply(hyper_ord, 
			function(i){
				sum(f.tab.hyper.edges$enrich %in% i)
			}),
		color = rep(colors, ceiling(nhyper/colors.max))[1:nhyper], 
		pval = suppressWarnings(as.numeric(sapply(hyper_ord, 
			function(i){
				f.tab.hyper[which(f.tab.hyper$category %in% i),"pval"]
			}))),
		fdr = suppressWarnings(as.numeric(sapply(hyper_ord, 
			function(i){
				f.tab.hyper[which(f.tab.hyper$category %in% i),"fdr"]
			}))))
		colnames(hyper_df)<-NULL

	colnames(f.tab.hyper.edges)<-NULL
	return(list(hypergsets = hyper_df, hyperedges = f.tab.hyper.edges))

	} else {
		return(list(hypergsets = NULL, hyperedges = NULL))
	}	
}

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

# writes json for pathway enrichment heatmap
#' @import jsonlite RColorBrewer
generate_heatmap_json<-function(header, gs.names, inputdir, outdir) {

	#TODO
	gs.names.base<-gsub(".gmt", "", basename(gs.names))

	finnames <- paste0(header, "_hyperEnrichment_", gs.names.base, "_cistrans_split.txt")

	names(finnames)<-gs.names.base
	datoutlist<-list()

	for(i in names(finnames)){

		finname<-finnames[i]
		finfull<-paste0(inputdir, "/", finname)
		dat<-read.table(finfull, header = T)

		##colorscales
		colors<-rev(brewer.pal(4, "Reds"))
		colorbreaks<-c(0.01, 0.05, 0.25)
		thres<-colorbreaks[length(colorbreaks)-1]

		##heatmapData
		heatmapData<-dat[, c("set", "category",  "fdr")]
		if(nrow(heatmapData)<1) return(NULL)
		rowLabs<-as.character(unique(heatmapData$category))
		colLabs<-as.character(unique(heatmapData$set))

		##remove rows and columns with all in-significant FDR
		rowLabsFilter<-sapply(rowLabs, function(i){
			cond<-any(heatmapData[heatmapData$category %in% i, "fdr"]<thres)
			if(cond) return(TRUE) else return(FALSE)
			})

		colLabsFilter<-sapply(colLabs, function(i){
			cond<-any(heatmapData[heatmapData$set %in% i, "fdr"]<thres)
			if(cond) return(TRUE) else return(FALSE)
			})

		heatmapData<-heatmapData[heatmapData$category %in% rowLabs[rowLabsFilter],]
		heatmapData<-heatmapData[heatmapData$set %in% colLabs[colLabsFilter],]

		rowLabs<-as.character(unique(heatmapData$category))
		colLabs<-as.character(unique(heatmapData$set))

		nR<-length(rowLabs)
		nC<-length(colLabs)

		mat<-sapply(rowLabs, function(i){
				sapply(colLabs, function(j){
					num<-as.numeric(heatmapData[heatmapData$category %in% i & heatmapData$set %in% j, 
						"fdr"])
				})
			})


		if(nR < 3) {
			row.ord<-1:nR
		} else {
			row.d<-dist(t(mat), method = "euclidean")
			row.ord<-order.dendrogram(as.dendrogram(hclust(row.d, method = "ward.D")))
		}
		if(nC < 3){
			col.ord<-1:nC
		} else {
			col.d<-dist(mat, method = "euclidean")
			col.ord<-order.dendrogram(as.dendrogram(hclust(col.d, method = "ward.D")))
		}

		rowLabs<-rowLabs[row.ord]
		colLabs<-colLabs[col.ord]
		
		heatmapData$x<-sapply(heatmapData$set, function(i){
			which(colLabs %in% i)
			})

		heatmapData$y<-sapply(heatmapData$category, function(i){
			which(rowLabs %in% i)
			})

		datout<-list(colors = colors, 
			colorbreaks = colorbreaks, 
			rowLabs = rowLabs,
			colLabs = colLabs,
			heatmapData = heatmapData)

		datoutlist[[i]]<-datout
	}

	dat.json<-toJSON(datoutlist, pretty = TRUE)
	meta.json<-toJSON(list(header = header), pretty = TRUE)
	foutname<-"combined_cistrans_split.json"
	foutfull<-paste0(outdir, "/", foutname)

	write_wrapper<-function(headername = "data", content, metacontent, fout){
		write("dataset = \n", append = FALSE, file = fout)
		write(content, append = TRUE, file = fout)
		write(";\n", append = TRUE, file= fout)
		write("heatmap_metadata = \n", append = TRUE, file = fout)
		write(metacontent, append = TRUE, file = fout)
		write(";", append = TRUE, file= fout)
	}

	write_wrapper(content = dat.json, metacontent=meta.json ,fout = foutfull)
}

##turns data frame into html code 
##appends js headers
#' @import knitr
to_table_html<-function(x, 
	x.name, 
	header = "", ## addition js header
	jsdir = "../../inst/javascript",
	hidecol = FALSE
	){

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
    		fixedHeader: false,
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
	res<-kable(x, "html", table.attr = attrstr)#caption = x.name)
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
	if(!is.na(boxplot_link)){
		in_table_headers1<-"<script type=\"text/javascript\" charset=\"utf8\" src=\"addboxplot.js\"></script>"
	} else {
		boxplot_link<-""
		in_table_headers1<-""
	}
	boxplot_link<-paste("<script> var dir_out = \"", boxplot_link, "\"; </script>", sep = "")
	in_table_headers2<-"<script type=\"text/javascript\" charset=\"utf8\" src=\"addgenecard.js\"></script>"
	in_table_headers<-paste(boxplot_link, in_table_headers1, in_table_headers2, sep = "\n")
	in_table_headers<-rep(in_table_headers, length(in_table))

 	jsdir <- "../../inst/javascript"
 	jsdir <-rep(jsdir, length(in_table))

	##write by-alteration DE tables to html
	res<-mapply(to_table_html, in_table, names(in_table), in_table_headers, jsdir, hidecol)
	sapply(names(res), 
		function(x){
		write(res[[x]], file = paste(outdir, "/", x, ".html", sep = ""))
		})
	
}

#' @import Biobase
make_boxplotdataset<-function(cn,
	gep,
	testtype,
	alterationid,
	gene,
	tabid = "Unique.Name",
	geneid = "accession",
	loggep = TRUE
	){
	altstatus_rowid<-which(fData(cn)[, tabid] == alterationid)
	if (length(altstatus_rowid) == 0){
		return(NULL)
	}
	altstatus<-suppressWarnings(as.numeric(exprs(cn)[altstatus_rowid,]))
	altstatus<-as.factor(altstatus)
	gep_rowid<-which(fData(gep)[, geneid] == gene)
	if (length(gep_rowid) == 0){
		return(NULL)
	}

	emat<-exprs(gep)
	ematrow<-suppressWarnings(as.numeric(emat[gep_rowid,]))

	if (loggep){
		emat<-emat+(-1*min(emat))
		ematrow<-suppressWarnings(as.numeric(emat[gep_rowid,]))
		ematrow<-log2(ematrow+1)
		logind <- "(log2 transformed)"
	}
	df<-data.frame(alteration = altstatus, 
		gene_expression = round(ematrow, 2), 
		name =  colnames(gep))
	df$name<-paste("\"", df$name, "\"", sep = "")
	return(df)
}


write_byalt_boxplot<-function(tab,#data frame to write
 tabid = "Unique.Name", #column name to split table by
 geneid = "accession", #column name for gene id
 outdir,
 header = "byalteration_cis",
 cn,
 gep, 
 basedir,
 testtype
 ){

 	if (!file.exists(outdir)){
		suppressWarnings(dir.create(outdir))
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
#' @import XML
make_bipartite_html<-function(f.dir.in, f.dir.out, 
	header = "", 
	headeradd = "", 
	template){

	suppressWarnings(dir.create(f.dir.out, recursive = TRUE))
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

		doc.html = htmlTreeParse(template, useInternalNodes = TRUE)
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
		#save to html for each JSON input file
		saveXML(doc.html, file = f.out)
	}
}

get_summary<-function(x, cistab, cisfulltab, transtab, altid, altdesc, cn, cisgenes,pruning){
	numcis <- 0
	numtrans <- 0
	if (nrow(cistab)>=1) numcis <- length(which(cistab[, altid] == x))
	if (nrow(transtab)>=1) numtrans <- length(which(transtab[, altid] == x))
	
	ind.cis<-which(fData(cn)[,altid] == x)[1]
	if(is.na(altdesc)) descriptor<-NA
	else descriptor<-fData(cn)[, altdesc][ind.cis]
	
	cislist <- cisgenes[[ind.cis]]
	cislist <- paste(cislist, collapse = ",")

	numAlt<-length(which(exprs(cn)[ind.cis,] == 1))
	numNormal<-length(which(exprs(cn)[ind.cis,] == 0))

	numbipartitecis<-0
	numbipartitetrans<-0

	numcisfull<-0
	if (nrow(cisfulltab)>=1)
	numcisfull<-length(which(cisfulltab[, altid] == x))
	numbipart<-NA
	suppressWarnings(if(!is.na(pruning)){
		numbipartitecis<-0
		numbipartitetrans<-0
		pruning.sig<-pruning[["sig"]]
		pruning.actual<-pruning[["actual"]]
		if(x %in% names(pruning.sig)){
			numbipartitecis<-length(unique(pruning.sig[[x]][,"cis"]))
			numbipartitetrans<-length(unique(pruning.sig[[x]][,"trans"]))
		}
		numbipart<-paste("(",numbipartitecis, "/", numbipartitetrans,")", sep = "")	
	})
	res<-data.frame(alteration_id = x, 
		alteration_description = descriptor,
		cis = paste("(", numcis, "/", numcisfull, ")", sep = ""), 
		trans = numtrans, 
		bipartite= numbipart,
		num_altered = numAlt,
		num_normal = numNormal,
		genes_in_alteration = cislist)
	return(res)
}

#Takes as input our summary table with links
#Outputs our full index file
to_main_index<-function(existingHTML, jsdir) {
	indexupperfname <- paste0(jsdir, '/indexupper.html')
	indexhtmlupper<- readChar(indexupperfname, file.info(indexupperfname)$size)
	indexlowerfname <- paste0(jsdir, '/indexlower.html')
	indexhtmllower <- readChar(indexlowerfname, file.info(indexlowerfname)$size)
	return(paste0(indexhtmlupper, "<br></br>", existingHTML, indexhtmllower));
}

#' iEDGE_UI makes the user interface for iEDGE reports
#' @export
iEDGE_UI<-function(cistab, cisfulltab, transtab, cn, gep, cisgenes,
	outdir, jsdir = file.path(path.package("iEDGE"), "javascript"), pruning, pruningjsdir,
	altid = "Unique.Name", altdesc = "Descriptor", geneid = "accession", 
	cis.boxplot = TRUE, trans.boxplot = TRUE, bipartite = TRUE, heatmap = TRUE, heatmapdatadir = NA, 
	header = NA, gs.names = NA){

	##TODO: if heatmap json was already generated by heatmap parameter is FALSE, remove heatmap file
	if (heatmap) {
		generate_heatmap_json(header, gs.names, heatmapdatadir, outdir) 
	}
	suppressWarnings(dir.create(outdir))

	jsfiles<-list.files(jsdir, full.names = TRUE)

	for(js in jsfiles){
		file.copy(from=js, to=outdir, 
		          overwrite = TRUE, recursive = FALSE, 
		          copy.mode = TRUE)
	}

	cat("Writing by alteration cis full DE tables...\n")
	boxplot_link_cis <- NA
	boxplot_link_trans <- NA

	if(cis.boxplot){
		boxplot_link_cis <- "./boxplots"
		cat("Writing cis boxplots...\n")
		outdirfig<-paste(outdir, "/", "boxplots", sep = "")
		suppressWarnings(dir.create(outdirfig))
		write_byalt_boxplot(cisfulltab,#data frame to write
		 tabid = altid, #column name to split table by
		 geneid = geneid, #column name for gene id
		 outdir = outdirfig,
		 header = "byalteration_cis",
		 cn = cn,
		 gep = gep, 
		 basedir = outdir
		 )
	}

	if(trans.boxplot){
		boxplot_link_trans <- "./boxplots"
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
	} 

	cat("Writing by alteration cis DE tables...\n")
	if(nrow(cisfulltab)>0)
		write_byalt_html(cisfulltab,#data frame to write
		 tabid = altid, #column name to split table by
		 outdir = outdir, #output directory name
		 header = "byalteration_cis", 
		 boxplot_link = boxplot_link_cis,
		 hidecol = FALSE
	 )

	cat("Writing by alteration trans DE tables...\n")
	if(nrow(transtab)>0)
		write_byalt_html(transtab,#data frame to write
		 tabid = altid, #column name to split table by
		 outdir = outdir, #output directory name
		 header = "byalteration_trans",
		 boxplot_link = boxplot_link_trans,
		 hidecol = FALSE
	 )

	##writing by alteration and by gene boxplots
	cat("Writing summary table...\n")
	alterations<-unique(as.character(fData(cn)[, altid]))

	if(bipartite){
		summarytab<-lapply(alterations, 
			function(x) get_summary(x, cistab, cisfulltab, transtab, altid, altdesc, cn, cisgenes, pruning = pruning))
		addlinksheader <-paste("<script type=\"text/javascript\" charset=\"utf8\" src=\"addlinks.js\"></script>", sep = "")

	} else {
		summarytab<-lapply(alterations, 
			function(x) get_summary(x, cistab, cisfulltab, transtab, altid, altdesc, cn, cisgenes, pruning = NA))
		addlinksheader <-paste("<script type=\"text/javascript\" charset=\"utf8\" src=\"addlinksnobipartite.js\"></script>", sep = "")
	}
	
	summarytab<-do.call(rbind, summarytab)
	summarytab<-data.frame(index = 1:nrow(summarytab), summarytab)
	summarytab<-summarytab[, !apply(summarytab, 2, function(i) all(is.na(i)))]
	row.names(summarytab)<-NULL
	
	#add links
	cat("Adding links to main table...\n")
	summarytable.html<-to_table_html(x = summarytab, 
		x.name = "summarytable", 
		header = addlinksheader,
		jsdir = jsdir,
		hidecol = TRUE)
	cat("Adding main table to main index...\n")
	summarytable.html<-to_main_index(summarytable.html, jsdir)
	write(summarytable.html, file = paste(outdir, "/", "index.html", sep = ""))

	##writing bipartite
	if(bipartite){
		cat("Writing bipartite graphs...\n")
		make_bipartite_html(f.dir.in = pruningjsdir, 
			f.dir.out = paste(outdir, "/bipartiteplots", sep = ""), 
			header = "", headeradd = "", 
			template = paste(outdir, "/", "template_bipartite.html", sep = ""))
	}
	cat(paste0("Done constructing html report at directory: ", outdir, "\n"))
	browseURL(paste0(outdir, "/index.html"))
	return(summarytab)
}
