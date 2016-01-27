
// add links from summarytable.html to open pages for individual alteration-based reports
function addGeneCardToDoc(tableName){
    	// Our columns

        var str0 = "\#";
        var str1 = str0.concat(tableName);
        var str2 = str1.concat(" thead tr th:eq(0)");
        column_name  = $(str2);
  
        //selection: element of second column
        //row_name: text of element of first column
        function addGeneCard(selection, row_name, tableName) {
        	//console.log(selection);
        	row_name = row_name.replace(/\ /gi,'');
            
            //text of first column
        	prev_value = selection.text();
       
        	new_text = "<a href=\"" + "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + row_name +
             "&keywords=" + row_name + "\">" + prev_value + "</a>";

        	selection.html(new_text);
        }

        rows = $(str1.concat(" tbody tr"));
        rows.each( function(index, element){
        	fst = $(element).children(0);
        	fst_name = fst.eq(0).text();
            snd = fst.eq(1);
            addGeneCard(snd , fst_name, tableName);
        });
    }

$(document).ready( function() {
    var htmlName = document.location.href.match(/[^\/]+$/)[0];
    var tableName = htmlName.replace("\.html", "");
    addGeneCardToDoc(tableName);
});

$(document).on("click", function() {
    var htmlName = document.location.href.match(/[^\/]+$/)[0];
    var tableName = htmlName.replace("\.html", "");
    addGeneCardToDoc(tableName);
});
$(document).change( function(){
    var htmlName = document.location.href.match(/[^\/]+$/)[0];
    var tableName = htmlName.replace("\.html", "");
    addGeneCardToDoc(tableName);
});


