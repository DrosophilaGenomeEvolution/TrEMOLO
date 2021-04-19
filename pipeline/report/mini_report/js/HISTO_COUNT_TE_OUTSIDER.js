var max_nb_x = 20;
var indice   = 1;
var size_x   = 0;

// Initialize the plot with the first dataset
//update(data1)
var db_TE = [0, 0, 0, 0];
var db = data;
db_TE[0] = db;

const div = d3.select("body").append("div")
    .attr("class", "tooltip")         
    .style("opacity", 0)
    .style("position", "absolute")
    .style("color", "black")
    .style("background-color", "white")
    .style("font-size", "1.5em")
    .style("padding", "6px")
    .style("border-radius", "10px")
    .style("box-shadow", "5px 5px 5px black")
    
//solidMDFXE755@mail

//set the dimensions and margins of the graph
var margin = {top: 10, right: 30, bottom: 90, left: 40},
    width = 820 - margin.left - margin.right,
    height = 720 - margin.top - margin.bottom;

// append the svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");


var x = d3.scaleBand()
  .range([ 0, width ])
  .padding(0.2);
var xAxis = svg.append("g")
  .attr("transform", "translate(0," + height + ")")

// Initialize the Y axis
var y = d3.scaleLinear()
  .range([ height, 0]);

var yAxis = svg.append("g")
  .attr("class", "myYaxis")


function update(data) {
  data.forEach(d => d.value = +d.value);

  // Update the X axis
  x.domain(data.map(function(d) { return d.name; }))
  xAxis.call(d3.axisBottom(x)).selectAll("text")
     .attr("transform", "translate(-10,0)rotate(-45)")
     .style("text-anchor", "end")
     .style("font-size", "1.4em");
     
  max_y = +d3.max(data, function(d) { return parseInt(d.value) });
  console.log("max", parseInt(max_y), max_y)

  // // Update the Y axis
  y.domain([0, parseInt(max_y) ]);
  yAxis.transition().duration(1000).call(d3.axisLeft(y))

  // Create the u variable
  var u = svg.selectAll("rect")
    .data(data)

  u
    .enter()
    .append("rect") // Add a new rect for each new elements
    
    .merge(u) // get the already existing elements as well
    .on("mouseover", function(d) {
            div.style("display", "inherit")
            div.transition()        
                .duration(200)      
                .style("opacity", .9);
            div.html("<stron>TE</strong> : " + d.name + " <br> value : " + d.value)
                .style("left", (d3.event.pageX + 10) + "px")     
                .style("top", (d3.event.pageY - 50) + "px")
        })
        .on("mousemove", function(d) {
            div.html("TE : " + d.name + " <br> value : " + d.value)
                .style("left", (d3.event.pageX + 10) + "px")     
                .style("top", (d3.event.pageY - 50) + "px")
        })
        .on("mouseout", function(d) {
            div.style("opacity", 0)
                .style("-webkit-user-select", "none") // For Webkit
                .style("-khtml-user-select", "none")
                .style("-moz-user-select", "none") // For Mozilla
                .style("display", "none")

        })
    .transition() // and apply changes to all of them
    .duration(800)
      .attr("x", function(d) { return x(d.name); })
      .attr("y", function(d) { return y(d.value); })
      .attr("width", x.bandwidth())
      .attr("height", function(d) { return height - y(d.value); })
      .attr("fill", "#FF0000")//005EC7//E0583B

  // If less group in the new dataset, I delete the ones not in use anymore
  u
    .exit()
    .remove()
}

function search(name){
    for(i in db_TE[0]){
        
        if( db_TE[0][i]["name"] == name ){
            data = db_TE[0].slice(parseInt(i), parseInt(i) + max_nb_x)
            update(data);
        }
    }
}

// function changes (type){
//     if (type == "LTR") {
//         db = db_TE[1]
//     }
//     else if (type == "DNA") {
//         db = db_TE[2]
//     }
//     else if (type == "NLR") {
//         db = db_TE[3]
//     }
//     else{
//         db = db_TE[0]
//     }

//     size_x   = db.length
//     //document.getElementById("nb").value = "NB-TE=" + size_x.toString()
//     d3.select("#nb").attr("data-to", size_x.toString())
//     d3.select("#nb").html(size_x.toString())

//     document.getElementById("nb").value = "NB-TE=" + size_x.toString()
//     indice   = 0
//     max_nb_x = 20
//     data     = db.slice(indice, indice + max_nb_x)
//     update(data)
//     size = parseInt(size_x) - parseInt(max_nb_x)
//     //console.log("size", size , "size_x", size_x, "max_nb_x", max_nb_x)
//     document.getElementById("slider").max = size.toString();
//     document.getElementById("slider").value = "0"
//     //d3.select("#slider").attr("max", size_x.toString())
// }

function next() {
    console.log("next")

    if( indice + max_nb_x <= size_x - max_nb_x ){
        indice += max_nb_x;
    }
    else{
        indice = size_x - max_nb_x;
    }
    data = db.slice(indice, indice + max_nb_x)
    update(data)
}


function pred() {
    console.log("pred")

    if( indice - max_nb_x >= 0 ){
        indice -= max_nb_x;
    }
    else{
        indice = 0;
    } 

    data = db.slice(indice, indice + max_nb_x)
    update(data)
}

function position(value){

    if( value + max_nb_x >= 0 || value + max_nb_x <= size_x - max_nb_x ){
        indice = value;
    }
    else if( value + max_nb_x < 0 ){
        indice = 0;
    } 
    else if( value + max_nb_x > size_x - max_nb_x){
        indice = size_x - max_nb_x;
    }

    console.log("value", value, "indice", indice)
    data = db.slice(indice, indice + max_nb_x)
    update(data)
}

d3.select("#nBin").on("input", function() {
    max_nb_x = +this.value
    data     = db.slice(indice, indice + max_nb_x)
    size_x   = db.length;
    size     = parseInt(size_x) - parseInt(max_nb_x)

    update(data)
    document.getElementById("slider").max   = size.toString();
    document.getElementById("slider").value = "0"
});


d3.select("#slider").on("input", function() {
    position(+this.value)
});


d3.select("#search").on("input", function() {
    search(this.value)
});

size_x   = data.length;
size     = parseInt(size_x) - parseInt(max_nb_x)

update(data.slice(1, 30))
document.getElementById("slider").max   = size.toString();
document.getElementById("slider").value = "0"
