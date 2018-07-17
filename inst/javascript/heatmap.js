$(document).ready(function() {
	
	var init = function(dataName){

		var data = dataset[dataName];
		var colors = data.colors;
		var colorbreaks = data.colorbreaks;
		var heatmapData = data.heatmapData;
		var rowLabs = data.rowLabs;
		var colLabs = data.colLabs;

		var nRow = rowLabs.length;
		var nCol = colLabs.length;

		var maxWidth = 5000;
		var maxHeight = 5000;

		var heatmapWidth = Math.min(nCol*12, maxWidth);
		var heatmapHeight = Math.min(nRow*12, maxHeight);

		var paddingDim = [320, 150];

		//set color gradient for the heatmap
		var colorRange = d3.scale.threshold().domain(colorbreaks).range(colors);

		//set color legend
		//concatenate white label to end to have same number of elements as labels
		var colorsNew = colors.concat(['#FFFFFF']) 
		var legendTileDim = [60, 30];
		var legendY = 60;
		var legendLabels = [0].concat(colorbreaks).concat([1])

		//create panel for color legends
		var legendPanel = d3.select(".legend-panel")		
			.attr("width", heatmapWidth+paddingDim[0])
			.attr("height", 70)
			.append("g")
			.attr("transform", "translate(" + 0 + "," + 0 + ")")

		var legend = legendPanel.selectAll("rect")
			.data(colorsNew, function(each) {return each;})
			.enter()
			.append("rect")
			.attr("x", function(d, index) { return (legendTileDim[0] * index + (heatmapWidth/2)); })
			.attr("y", legendY)
			.attr("width", legendTileDim[0])
			.attr("height", legendTileDim[1])
			.attr("class", "legendTile")
			.style("fill", function(d, index) { return colorsNew[index]; });

		var legendText = legendPanel.selectAll("text")
			.data(legendLabels, function(each){return each;})
			.enter()
			.append("text")
			.text(function(d, index) { return legendLabels[index]; })
			.attr("x", function(d, index) { return (legendTileDim[0] * index) + (heatmapWidth/2); })
			.attr("y", legendY - 10)
			.attr("class", "legendText")
			.style("text-anchor", "middle")
			.style("alignment-baseline", "bottom");

		var heatMap = d3.select(".heat-map")
			.attr("width", heatmapWidth+paddingDim[0])
			.attr("height", heatmapHeight+paddingDim[1])
			.attr("class", "heat-map")
			.append("g")
			.attr("transform", "translate(" + paddingDim[0] + "," + 0 + ")")

		//add blank label to the beginning so labels appear between ticks instead of on ticks
		var rowLabsNew = [""].concat(rowLabs);
		var colLabsNew = [""].concat(colLabs);

		//set up axis scales
		var xScale = d3.scale.ordinal().domain(colLabsNew).rangePoints([0, heatmapWidth]);
		var yScale = d3.scale.ordinal().domain(rowLabsNew).rangePoints([0, heatmapHeight]);

		var xAxis = d3.svg.axis()
			.scale(xScale)
			.orient("bottom");

		var yAxis = d3.svg.axis()
			.scale(yScale)
			.orient("left");
					
		//slide text from aligning on tick mark to aligning between tick marks
		var tileWidth = heatmapWidth/nCol;
		var tileHeight = heatmapHeight/nRow;

		//space between axis lines and text
		var tickToLabel = 5; 

		//set tooltip for heatmap tiles, default invisible
		var tooltip = d3.select("#wrapper")
			.append("g")
			.attr("class", "tooltip")
			.style("opacity", 0);

		//draw x and y axis and labels
		heatMap.append("g")
			.attr("class", "x axis")
			.call(xAxis)
			.attr("transform", "translate(" + 0 + ", " + heatmapHeight + ")")
			.selectAll("text")
			.attr("transform", "translate(" + tileWidth/2 + ", " + tickToLabel + ") rotate(90)")
			.style("text-anchor", "start")
			.style("alignment-baseline", "bottom");

		heatMap.append("g")
			.attr("class", "y axis")
			.call(yAxis)
			.selectAll("text")
			.attr("transform", "translate(" + tickToLabel + ", " + -1*tileHeight/2 + ")")
			.style("alignment-baseline", "bottom");

		//draw heatmap tiles
		var rectangles = heatMap.selectAll("rect.stats")
			.data(heatmapData)
			.enter()
			.append("rect")
			.attr("class", "stats")
			.attr("x", function(item) {
				return tileWidth * (item.x-1);
			})
			.attr("y", function(item) {
				return tileHeight * (item.y-1);
			})
			.attr("width", tileWidth)
			.attr("height", tileHeight)
			.style("fill", function(each) {
				return colorRange(each.fdr);
			})
			.on("mouseover", function(each) {
				tooltip.transition()
				.style("opacity", 0.7);
				tooltip.html(each.set + "<br/>" + each.category + "<br/>" + "FDR:" + each.fdr)
				.style("left", (d3.event.pageX) - 70 + "px")
				.style("top", (d3.event.pageY) - 90 + "px");
			})
			.on("mouseout", function(each) {
				tooltip.transition()
				.style("opacity", 0);
			});
	}

	var remove = function(){
		d3.selectAll("g")
		.remove();
		}

	//create heatmap
	var datasetItems = Object.keys(dataset);

	var datasetMenu = d3.select("#menu")
	datasetMenu.append("select")
		.selectAll("option")
		.data(datasetItems)
		.enter()
		.append("option")
		.attr("value", function(d, index){
			return datasetItems[index];
		})
		.text(function(d, index){
			return datasetItems[index];
		})

	init(datasetItems[0]);

 	datasetMenu.on('change', function(){
 		
 		//clear existing panels
 		remove()

 		var selected = d3.select(this)
            .select("select")
            .property("value")
    
    	//initialize new panels
        init(selected)
    });

    $("#heatmap-title").html("Pathway enrichment of cis-or-trans signatures for dataset: " + heatmap_metadata.header[0]);
});