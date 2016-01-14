	

 d3.csv("/Users/amyli/Desktop/git_projects/iEDGE/test/html/prices.csv", function(data)
    {
     console.log(data[0]);
    })


  var visualization = d3plus.viz()
    .container("#viz")
    .data(data)
    .type("box")
    //.id("name")
    .x("month")
    .y("price")
    .time("month")
    .ui([{ 
        "label": "Visualization Type",
        "method": "type", 
        "value": ["scatter","box"]
      }])
    .draw()
