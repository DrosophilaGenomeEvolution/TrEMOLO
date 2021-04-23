var max_nb_x = 20;
var indice   = 1;
var size_x   = 0;

// Initialize the plot with the first dataset
//update(data1)
var db_TE = [0, 0, 0, 0];
var db = 0;

d3.csv("https://mites.infinity-atom.fr/db.csv", function(data_all) {
    db_TE[0] = data_all;
    db       = db_TE[0];
    data     = db.slice(1, 30)
    size_x   = db.length;
    size     = parseInt(size_x) - parseInt(max_nb_x)

    update(data)
    document.getElementById("slider").max   = size.toString();
    document.getElementById("slider").value = "0"

});

d3.csv("https://mites.infinity-atom.fr/db_LTR.csv", function(data_all) {
    db_TE[1] = data_all;
});

d3.csv("https://mites.infinity-atom.fr/db_DNA.csv", function(data_all) {
    db_TE[2] = data_all;
});

d3.csv("https://mites.infinity-atom.fr/db_NLR.csv", function(data_all) {
    db_TE[3] = data_all;
});


const div = d3.select("body").append("div")
    .attr("class", "tooltip")         
    .style("opacity", 0)
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

d3.csv("https://mites.infinity-atom.fr/db_NLR.csv", function(data_all) {
    max_y = d3.max(data_all, function(d) { return parseInt(d.value) });
    all = data_all.map(function(d) { return parseInt(d.value); })

    //console.log(max_y)
});

function update(data) {
  data.forEach(d => d.value = +d.value);

  // Update the X axis
  x.domain(data.map(function(d) { return d.name; }))
  xAxis.call(d3.axisBottom(x)).selectAll("text")
     .attr("transform", "translate(-10,0)rotate(-45)")
     .style("text-anchor", "end")
     .style("font-size", "1.4em");
     
  max_y = +d3.max(data, function(d) { return parseInt(d.value) });
  console.log("max", parseInt(max_y) + 100, max_y)

  // // Update the Y axis
  y.domain([0, parseInt(max_y) + 100 ]);
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
      .attr("fill", "#0E0009")//005EC7//E0583B

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

function changes (type){
    if (type == "LTR") {
        db = db_TE[1]
    }
    else if (type == "DNA") {
        db = db_TE[2]
    }
    else if (type == "NLR") {
        db = db_TE[3]
    }
    else{
        db = db_TE[0]
    }

    size_x   = db.length
    //document.getElementById("nb").value = "NB-TE=" + size_x.toString()
    d3.select("#nb").attr("data-to", size_x.toString())
    d3.select("#nb").html(size_x.toString())

    document.getElementById("nb").value = "NB-TE=" + size_x.toString()
    indice   = 0
    max_nb_x = 20
    data     = db.slice(indice, indice + max_nb_x)
    update(data)
    size = parseInt(size_x) - parseInt(max_nb_x)
    //console.log("size", size , "size_x", size_x, "max_nb_x", max_nb_x)
    document.getElementById("slider").max = size.toString();
    document.getElementById("slider").value = "0"
    //d3.select("#slider").attr("max", size_x.toString())
}

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


(function ($) {
    $.fn.countTo = function (options) {
        options = options || {};
        
        return $(this).each(function () {
            // set options for current element
            var settings = $.extend({}, $.fn.countTo.defaults, {
                from:            $(this).data('from'),
                to:              $(this).data('to'),
                speed:           $(this).data('speed'),
                refreshInterval: $(this).data('refresh-interval'),
                decimals:        $(this).data('decimals')
            }, options);
            
            // how many times to update the value, and how much to increment the value on each update
            var loops = Math.ceil(settings.speed / settings.refreshInterval),
                increment = (settings.to - settings.from) / loops;
            
            // references & variables that will change with each update
            var self = this,
                $self = $(this),
                loopCount = 0,
                value = settings.from,
                data = $self.data('countTo') || {};
            
            $self.data('countTo', data);
            
            // if an existing interval can be found, clear it first
            if (data.interval) {
                clearInterval(data.interval);
            }
            data.interval = setInterval(updateTimer, settings.refreshInterval);
            
            // initialize the element with the starting value
            render(value);
            
            function updateTimer() {
                value += increment;
                loopCount++;
                
                render(value);
                
                if (typeof(settings.onUpdate) == 'function') {
                    settings.onUpdate.call(self, value);
                }
                
                if (loopCount >= loops) {
                    // remove the interval
                    $self.removeData('countTo');
                    clearInterval(data.interval);
                    value = settings.to;
                    
                    if (typeof(settings.onComplete) == 'function') {
                        settings.onComplete.call(self, value);
                    }
                }
            }
            
            function render(value) {
                var formattedValue = settings.formatter.call(self, value, settings);
                $self.html(formattedValue);
            }
        });
    };
    
    $.fn.countTo.defaults = {
        from: 0,               // the number the element should start at
        to: 0,                 // the number the element should end at
        speed: 1000,           // how long it should take to count between the target numbers
        refreshInterval: 100,  // how often the element should be updated
        decimals: 0,           // the number of decimal places to show
        formatter: formatter,  // handler for formatting the value before rendering
        onUpdate: null,        // callback method for every time the element is updated
        onComplete: null       // callback method for when the element finishes updating
    };
    
    function formatter(value, settings) {
        return value.toFixed(settings.decimals);
    }
}(jQuery));





$('.btn_type_TE').on("click", function ($) {
  $('.count-number').data('countToOptions', {
    formatter: function (value, options) {
        console.log(value.toFixed(options.decimals).replace(/\B(?=(?:\d{3})+(?!\d))/g, ','), "bbbb")
      return value.toFixed(options.decimals).replace(/\B(?=(?:\d{3})+(?!\d))/g, ',');
    }
  });
  
  // start all the timers
  $('.timer').each(count);  
  
  function count(options) {
    var $this = $(this);
    options = $.extend({}, options || {}, $this.data('countToOptions') || {});
    console.log(options, "options")
    $this.countTo(options);
  }
})

//
jQuery(function ($) {
  // custom formatting example
  $('.count-number').data('countToOptions', {
    formatter: function (value, options) {
      return value.toFixed(options.decimals).replace(/\B(?=(?:\d{3})+(?!\d))/g, ',');
    }
  });
  
  // start all the timers
  $('.timer').each(count);  
  
  function count(options) {
    var $this = $(this);
    options = $.extend({}, options || {}, $this.data('countToOptions') || {});
    $this.countTo(options);
  }

});



