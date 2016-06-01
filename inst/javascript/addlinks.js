function addLinksToDoc(){
    	// Our columns
        column_name = $("#summarytable thead tr th:gt(0)");
        function addLink(selection, column_name, row_name) {
        	column_name = column_name.replace(/\ /gi,'');
        	row_name = row_name.replace(/\ /gi,'');
        	link_name = "byalteration_" + column_name + "_" + row_name + ".html";
        	prev_value = selection.text();

        	if(prev_value == 0) { // no links for empty tables
        		return;
        	}
        	new_text = "<a href=\"" + link_name + "\">" + prev_value + "</a>";
        	//console.log(new_text);
        	selection.html(new_text);
        }

        function addLinkBipartite(selection, column_name, row_name) {
            column_name = column_name.replace(/\ /gi,'');
            row_name = row_name.replace(/\ /gi,'');
            link_name = "./bipartiteplots/" + row_name + ".html";
            prev_value = selection.text();
            prev_value = prev_value.replace(/ /gi, '');
         
            if(prev_value == 0) { // no links for empty tables
                return;
            }
            if(prev_value == "(0/0)") { // no links for empty tables
     
                return;
            }
            new_text = "<a href=\"" + link_name + "\">" + prev_value + "</a>";
            //console.log(new_text);
            selection.html(new_text);
        }

        rows = $("#summarytable tbody tr");
        rows.each( function(index, element){
        	fst = $(element).children(0);
        	rest = fst.nextAll();
        	//rest contains all the elements we want to add links to
        	//console.log(rest.size());
        	r_name = fst.eq(0).text();
        	for (var i = 0; i < rest.size(); i++) {
        		c_name = column_name.eq(i).text();
                c_name = c_name.replace(/ /gi, '');
        		elt = rest.eq(i);
        		//console.log([r_name, c_name, elt.text()]);
                if(c_name == 'bipartite'){
                    addLinkBipartite(elt, c_name, r_name);
                    
                } else if (c_name == 'cis' || c_name == 'trans' || c_name == 'num_genes_in_alteration'){
        		  addLink(elt, c_name, r_name);
                }
        	}
        });
    }

$(document).ready( function() {
    	addLinksToDoc();
});
$(document).on("click", addLinksToDoc);
$(document).change(addLinksToDoc);
