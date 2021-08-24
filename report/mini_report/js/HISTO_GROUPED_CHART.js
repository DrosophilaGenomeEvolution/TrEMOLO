// Themes begin
am4core.useTheme(am4themes_animated);
//am4core.useTheme(am4themes_dataviz);


// Create chart instance
var chart = am4core.create("chartdiv", am4charts.XYChart3D);

// var data = [{"group": "Max-element_FIVE_PRIME", 'TSD KO': "5", "TSD OK": "1"}, {"group": "Transpac_FIVE_PRIME", 'TSD KO': "1", "TSD OK": "0"}, {"group": "springer_FIVE_PRIME", 'TSD KO': "1", "TSD OK": "0"}, {"group": "gypsy_FIVE_PRIME", 'TSD KO': "1", "TSD OK": "0"}, {"group": "HMS-Beagle_FIVE_PRIME", 'TSD KO': "1", "TSD OK": "0"}, {"group": "roo_FIVE_PRIME", 'TSD KO': "2", "TSD OK": "0"}, {"group": "412_FIVE_PRIME", 'TSD KO': "2", "TSD OK": "0"}, {"group": "ZAM_FIVE_PRIME", 'TSD KO': "0", "TSD OK": "2"}, {"group": "gtwin_FIVE_PRIME", 'TSD KO': "2", "TSD OK": "0"}, {"group": "blood_FIVE_PRIME", 'TSD KO': "2", "TSD OK": "0"}, {"group": "mdg3_FIVE_PRIME", 'TSD KO': "6", "TSD OK": "0"}]

// //data["columns"] = ["group",  "TSD OK", 'TSD KO']


// Add data
function add_series( categorie, valueY) {

    var series = chart.series.push(new am4charts.ColumnSeries3D());
    series.dataFields.valueY = valueY;
    series.dataFields.categoryX = "group";
    series.name = categorie;
    series.clustered = false;
    series.columns.template.tooltipText = " {name} : [bold]{valueY}[/]";
    series.columns.template.fillOpacity = 0.9;

    return series;
}


function init(data) {

    chart.data = data;

    // Create axes
    var categoryAxis = chart.xAxes.push(new am4charts.CategoryAxis());
    categoryAxis.dataFields.category = "group";
    categoryAxis.renderer.grid.template.location = 0;
    categoryAxis.renderer.minGridDistance = 30;

    var valueAxis = chart.yAxes.push(new am4charts.ValueAxis());
    valueAxis.title.text = "COUNT";
    valueAxis.renderer.labels.template.adapter.add("text", function(text) {
      return text + "";
    });

    // Configure axis label
    var label = categoryAxis.renderer.labels.template;
    label.truncate = true;
    label.maxWidth = 200;
    label.tooltipText = "{category}";

    categoryAxis.events.on("sizechanged", function(ev) {
      var axis = ev.target;
      var cellWidth = axis.pixelWidth / (axis.endIndex - axis.startIndex);
      if (cellWidth < axis.renderer.labels.template.maxWidth) {
        axis.renderer.labels.template.rotation = -45;
        axis.renderer.labels.template.horizontalCenter = "right";
        axis.renderer.labels.template.verticalCenter = "middle";
      }
      else {
        axis.renderer.labels.template.rotation = 0;
        axis.renderer.labels.template.horizontalCenter = "middle";
        axis.renderer.labels.template.verticalCenter = "top";
      }
    });

    chart.legend = new am4charts.Legend();
    
    //Create Series
    for (const key in chart.data["columns"]) {
        if(chart.data["columns"][key] != "group"){
            series_list.push(add_series(chart.data["columns"][key], chart.data["columns"][key]));
        }
    }
}

var series_list = [];
init(data);

function update(data) {
    chart.data = data;
    chart.reinit();
}

//Range Gestion
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


//**********
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
    db = data.slice(indice, indice + max_nb_x);
    update(db);
}

var max_nb_x = 15;
var indice   = 0;
var db       = data.slice(indice, indice + max_nb_x)
var size_x   = data.length
update(db);

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
