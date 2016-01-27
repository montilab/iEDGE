
// add links from summarytable.html to open pages for individual alteration-based reports
function addBoxplotToDoc(tableName){
    	// Our columns

        var str0 = "\#";
        var str1 = str0.concat(tableName);
        var str2 = str1.concat(" thead tr th:eq(0)");
        column_name  = $(str2);
       // column_name = $("#summarytable thead tr th:gt(0)");
        function addBoxplot(selection, row_name, tableName) {
        	
        	row_name = row_name.replace(/\ /gi,'');
    
        	prev_value = selection.text();
        
            dir_out = "../../test2/figures/";
           // link_name = "summarytable.html"
           
            link_name = tableName + "_" + prev_value + ".html";

            //replace special symbols not allowed in file names
            link_name = link_name.replace(/ /g, "");
            link_name = link_name.replace(/\//g, ".");
            link_name = link_name.replace(/-/g, ".");

        //    console.log(link_name)
        	new_text = "<a href=\"" + dir_out + "/" + link_name + "\">" + prev_value + "</a>";

        	selection.html(new_text);
        }

        rows = $(str1.concat(" tbody tr"));
        //rows = $("#summarytable tbody tr");

        rows.each( function(index, element){
        	fst = $(element).children(0);
        	r_name = fst.eq(0);
            addBoxplot(r_name, r_name.text(), tableName);
        });
    }

//var htmlName = document.location.href.match(/[^\/]+$/)[0];
//var tableName = htmlName.replace("\.html", "");
//console.log(tableName);
//alert(tableName)




$(document).ready( function() {
    var htmlName = document.location.href.match(/[^\/]+$/)[0];
    var tableName = htmlName.replace("\.html", "");
    addBoxplotToDoc(tableName);
});

$(document).on("click", function() {
    var htmlName = document.location.href.match(/[^\/]+$/)[0];
    var tableName = htmlName.replace("\.html", "");
    addBoxplotToDoc(tableName);
});
$(document).change( function(){
    var htmlName = document.location.href.match(/[^\/]+$/)[0];
    var tableName = htmlName.replace("\.html", "");
    addBoxplotToDoc(tableName);
});


