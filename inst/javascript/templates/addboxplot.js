

// add links from summarytable.html to open pages for individual alteration-based reports
function addBoxplotToDoc(tableName, dir_out){
    	// Our columns

        var str0 = "\#";
        var str1 = str0.concat(tableName);
        var str2 = str1.concat(" thead tr th:eq(0)");
        column_name  = $(str2);
 
        function addBoxplot(selection, row_name, tableName, dir_out) {
        	
        	row_name = row_name.replace(/\ /gi,'');
    
        	prev_value = selection.text();
           
            link_name = tableName + "_" + prev_value + ".html";

            link_name = link_name.replace(/ /g, "");
            link_name = link_name.replace(/\//g, ".");
            link_name = link_name.replace(/-/g, ".");

        	new_text = "<a href=\"" + dir_out + "/" + link_name + "\">" + prev_value + "</a>";

        	selection.html(new_text);
        }

        rows = $(str1.concat(" tbody tr"));
        //rows = $("#summarytable tbody tr");

        rows.each( function(index, element){
        	fst = $(element).children(0);
        	r_name = fst.eq(0);
            addBoxplot(r_name, r_name.text(), tableName, dir_out);
        });
    }

$(document).ready( function() {
    var htmlName = document.location.href.match(/[^\/]+$/)[0];
    var tableName = htmlName.replace("\.html", "");
    addBoxplotToDoc(tableName, dir_out);
});

$(document).on("click", function() {
    var htmlName = document.location.href.match(/[^\/]+$/)[0];
    var tableName = htmlName.replace("\.html", "");
    addBoxplotToDoc(tableName, dir_out);
});
$(document).change( function(){
    var htmlName = document.location.href.match(/[^\/]+$/)[0];
    var tableName = htmlName.replace("\.html", "");
    addBoxplotToDoc(tableName, dir_out);
});


