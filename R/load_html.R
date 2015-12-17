
##turns data frame into html code 
##appends js headers
to_table_html<-function(x, x.name, 
	header = ""## addition js header
	){

	library(knitr)

	##add default js header
	x.name.escape<-gsub("\\.", "\\\\.", x.name)
	js_header<-paste("<script type=\"text/javascript\" charset=\"utf8\" src=\"./bower_components/jquery/dist/jquery.js\"></script>

<link rel=\"stylesheet\" type=\"text/css\" href=\"./bower_components/datatables.net-dt/css/jquery.dataTables.css\">

<script type=\"text/javascript\" charset=\"utf8\" src=\"./bower_components/datatables.net/js/jquery.dataTables.js\"></script>

<script type=\"text/javascript\">
    $(document).ready(function() {
        $('#", x.name.escape, "').DataTable();
    } );
</script>", sep = "")

	js_header2<-header

	attrstr<-paste("id=\"", x.name, "\"", sep = "")
	res<-kable(x, "html", table.attr = attrstr,
	caption = x.name)
	res<-paste(js_header, js_header2, res, sep = "\n")
	return(res)
}

##main script

##load by-alteration DE tables
in_dir<-"../test/tables"
in_files<-list.files(path = in_dir, pattern = "^byalteration*", full.names = FALSE)
in_table<-lapply(in_files, function(x){
 	res<-read.table(paste(in_dir, "/", x, sep = ""), sep = "\t", header = TRUE, quote = "")
 	return(res)	
})
names(in_table)<-gsub("\\.txt", "", in_files)

##write by-alteration DE tables to html
res<-mapply(to_table_html, in_table, names(in_table))
out_dir<-"../test/html"
sapply(names(res), 
	function(x){
		write(res[[x]], file = paste(out_dir, "/", x, ".html", sep = ""))
		})

#write summary DE table to html
summarytable<-read.table(paste(in_dir, "/", "summarytable.txt", sep = ""), header = TRUE, quote = "")
	addinksheader <-"<script type=\"text/javascript\" charset=\"utf8\" src=\"../../R/addlinks.js\"></script>
"
#add links
summarytable.html<-to_table_html(x = summarytable, x.name = "summarytable", header = addlinksheader)
write(summarytable.html, file = paste(out_dir, "/", "summarytable", ".html", sep = ""))


