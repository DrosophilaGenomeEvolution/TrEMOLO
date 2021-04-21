// var BioCircosGenome = [      // Configure your own genome here.
//      ["2L" , 23513712],
//      ["2R" , 25286936],
//      ["3L" , 28110227],
//      ["3R" , 32079331],
//      ["4" , 1348131],
//    ["X" , 23542271]
// ];


for (var key in dico_histo) {
    $("#pet-select").append("<option value=\"" + key + "\">" + key + "</option>")
}

var BioCircos01 = null;

function init(name_TE) {
    $("#biocircos").html("");

    var var_histo = null
    if (name_TE != "ALL") {
      var_histo = dico_histo[name_TE];
    }
    else{
      var_histo = HISTOGRAM;
    }

    BioCircos01 = new BioCircos(BACKGROUND01,var_histo,BioCircosGenome,{       // Initialize BioCircos.js with "BioCircosGenome" and Main configuration
       target : "biocircos",                              // Main configuration "target"
       svgWidth : 900,                                  // Main configuration "svgWidth"
       svgHeight : 600,                                 // Main configuration "svgHeight"
       innerRadius: 246,
       outerRadius: 252,
       zoom : true,
       genomeBorder : {
          display : true,
          borderColor : "#000",
          borderSize : 0.5
       },
       ticks : {
          display : true,
          len : 5,
          color : "#000",
          textSize : 10,
          textColor : "#000",
          scale : 2000000
       },
       genomeLabel : {
          display : true,
          textSize : 15,
          textColor : "#000",
          dx : 0.028,
          dy : "-0.55em"
       },
       genomeFillColor: ["#999999"],
       HISTOGRAMMouseEvent: true,
       HISTOGRAMMouseOutDisplay: true,
       HISTOGRAMMouseOutAnimationTime: 500,
       HISTOGRAMMouseOutColor: "gold",
       HISTOGRAMMouseOutOpacity: 1,
       HISTOGRAMMouseOutStrokeColor: "red",
       HISTOGRAMMouseOutStrokeWidth: 0,
       HISTOGRAMMouseOverDisplay: true,
       HISTOGRAMMouseOverColor: "red",
       HISTOGRAMMouseOverOpacity: 1,
       HISTOGRAMMouseOverStrokeColor: "#F26223",
       HISTOGRAMMouseOverStrokeWidth: 5,
    });

  BioCircos01.draw_genome(BioCircos01.genomeLength);
}


init("ALL");


$(document).on('change', '#pet-select', function(){
  init($("#pet-select").val());
});


