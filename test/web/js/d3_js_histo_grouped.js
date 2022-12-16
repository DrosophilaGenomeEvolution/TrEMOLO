var max_nb_x = 20;
var indice   = 0;
var db       = data.slice(indice, indice + max_nb_x)
var size_x   = data.length


// set the dimensions and margins of the graph
var margin = {top: 10, right: 30, bottom: 100, left: 50},
    width  = 720 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

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

// ----------------
// Create a tooltip
// ----------------
var tooltip = d3.select("#my_dataviz")
  .append("div")
  .style("opacity", 0)
  .attr("class", "tooltip")
  .style("position", "relative")
  .style("background-color", "white")
  .style("border", "solid")
  .style("border-width", "1px")
  .style("border-radius", "5px")
  .style("padding", "10px")

  // Three function that change the tooltip when user hover / move / leave a cell
  var mouseover = function(d) {
    var subgroupName = d3.select(this.parentNode).datum().key;
    var subgroupValue = d.data[subgroupName];
    tooltip
        .html("subgroup: " + subgroupName + "<br>" + "Value: " + subgroupValue)
        .style("opacity", 1)
  }
  var mousemove = function(d) {
    // tooltip
    //   .style("left", (d3.mouse(this)[0]+90) + "px") // It is important to put the +90: other wise the tooltip is exactly where the point is an it creates a weird effect
    //   .style("top", (d3.mouse(this)[1]) + "px")
  }

  var mouseleave = function(d) {
    tooltip.style("opacity", 0)
  }

  var color = d3.scaleOrdinal()
    .range(['#1A86DA','#4DB1FF','#24E1F3'])

  //stack the data? --> stack per subgroup
  var stackedData = d3.stack();
  var subgroups = ""
  var groups = ""

  svg.append("g").attr("class", "th");

function cal_max(data){
  max=0
  for(var i= 0; i < data.length; i++)
  {
    sum=0;
       //console.log(data[i]);
       for(var e in data[i])
       {
          if(e!="group" ){
             sum+= parseInt(data[i][e]);
          }
      }
      if (max<sum) {
        max=sum
      }
  }
  return max
}


function update(data, columns) {
  // console.log(data)
  // console.log(columns)
  //console.log(typeof data)

  // List of subgroups = header of the csv files = soil condition here
  subgroups = columns
  //console.log(subgroups)
  // List of groups = species here = value of the first column called group -> I show them on the X axis
  groups = d3.map(data, function(d){return(d.group)}).keys()
  //console.log(groups)

  //console.log(groups, "goupe", x)

  x.domain(groups)
  xAxis.call(d3.axisBottom(x)).selectAll("text")
     .attr("transform", "translate(-10,0)rotate(-45)")
     .style("text-anchor", "end")
     .style("font-size", "1.4em");

  //calcul max for update
  max = cal_max(data);
  y.domain([0, max+10]);
  yAxis.transition().duration(1000).call(d3.axisLeft(y))

  // //stack the data? --> stack per subgroup
  stackedData = d3.stack()
    .keys(subgroups)
    (data)

  color.domain(subgroups)

  svg.select("g.th").remove();
  svg.append("g").attr("class", "th");
  var u = svg.select("g.th")
    .selectAll("g")
    .data(stackedData)
  // Show the bars
  
    u.enter().append("g").merge(u)
      .attr("fill", function(d) { return color(d.key); })
      .selectAll("rect")
      // enter a second time = loop subgroup per subgroup to add all rectangles
      .data(function(d) { return d; })
      .enter().append("rect").merge(u)
      .on("mouseover", mouseover)
      .on("mousemove", mousemove)
      .on("mouseleave", mouseleave)
      .transition() // and apply changes to all of them
    .duration(800)
        .attr("x", function(d) { return x(d.data.group); })
        .attr("y", function(d) { return y(d[1]); })
        .attr("height", function(d) { return y(d[0]) - y(d[1]); })
        .attr("width",x.bandwidth())
        //.attr("stroke", "grey")

    u
    .exit()
    .remove()
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

    //console.log("value", value, "indice", indice)
    db = data.slice(indice, indice + max_nb_x)
    update(db, data.columns.slice(1))
}

update(db, data.columns.slice(1))
//console.log(data)

size = parseInt(size_x) - parseInt(max_nb_x)
document.getElementById("slider").max   = size.toString();
document.getElementById("slider").value = "0"
if (size<0) {
    size = 0;
}
$("#slider").closest(".container").closest(".container").find(".box-minmax span:nth-child(2)").text(size.toString());


d3.select("#slider").on("input", function() {
    position(+this.value)
});



//https://codepen.io/mayuMPH/pen/ZjxGEY
var rangeSliders = document.getElementsByClassName ("rs-range");
var rangeBullets = document.getElementsByClassName("rs-label");

Array.prototype.filter.call(rangeSliders, function(rangeSlider){
    rangeSlider.addEventListener("input", showSliderValue, false);
});

Array.prototype.filter.call(rangeSliders, function(rangeSlider){
    rangeBullet = rangeSlider.closest(".range-slider").children[0];

    rangeBullet.innerHTML  = rangeSlider.value;
    var bulletPosition     = (rangeSlider.value /rangeSlider.max);
    rangeBullet.style.left = (bulletPosition * 578) + "px";
});

function showSliderValue(event) {
    rangeSlider = event.target;
    rangeBullet = rangeSlider.closest(".range-slider").children[0];

    rangeBullet.innerHTML  = rangeSlider.value;
    var bulletPosition     = (rangeSlider.value /rangeSlider.max);
    rangeBullet.style.left = (bulletPosition * 578) + "px";
}

