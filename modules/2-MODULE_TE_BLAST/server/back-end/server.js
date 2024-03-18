const http = require("http");
const fs = require('fs');
const querystring = require('node:querystring');
const path = require('path');
const yaml = require('js-yaml')



const host = 'localhost';
let port = 8030;
let port_client = 8031;
const url = require('url');


try {
    const yamlData = yaml.load(fs.readFileSync('../config.yaml', 'utf8'));
    port = yamlData.PORT_SERVER;
    port_client = yamlData.PORT_CLIENT;
} 
catch (error) {
    console.error(error);
}

console.log("port:", port);

// Gobal
var index_TE_INFOS = {};
var index_DEPTH = {};

function TE_in_selected(dico, TE_selected) {
    if ( TE_selected.length == 0 ) {
        return true;
    }

    for(const index in TE_selected){
        if( dico["sseqid"] == TE_selected[index] ){
            return true;
        }
    }
    return false;
}


function read_data_html(data, chrom, TE_selected, MODE_P, generation) {
    var seqs   = "";
    var TE     = "";
    var query  = "";
    var id     = "";
    var index_chrom = 0;

    // other
    var class_a = "";

    for (const index in data["data"]) {
                
        var dico = data["data"][index];
        if ( TE_in_selected(dico, TE_selected) && ( chrom == dico["chrom"] || chrom == undefined ) ) {

            // 2L_RaGOO_RaGOO:<INS>:547321:547321:sniffles.INS.806:3:IMPRECISE:2
            if (id != data["data"][index]["qseqid"].split(":")[4]) {
                start = new Intl.NumberFormat("de-DE").format(data["data"][index]["qseqid"].split(":")[2]);
                end   = new Intl.NumberFormat("de-DE").format(data["data"][index]["qseqid"].split(":")[3]);

                seqs += `<div class="id-separate" data-locus="${data["data"][index]["qseqid"].split(":")[0] + ":" + start + "-" + end}" data-index="${index_chrom}" style="position:relative">
                            <div class="dot" style="position:absolute; left:-10px"></div>
                            <div class="dot" style="position:absolute; left:240px"></div>
                            <div class="dot" style="position:absolute; left:490px"></div>
                        </div>`;
                id    = data["data"][index]["qseqid"].split(":")[4];
                index_chrom += 1;
            }

            if (query != data["data"][index]["qseqid"]) {
                seqs += `<div class="query-separate"></div>`;
                query = data["data"][index]["qseqid"];
            }


            if (TE != data["data"][index]["sseqid"]) {
                seqs += `<div class="te-name"><a>${data["data"][index]["sseqid"]}</a></div>`;
                TE = data["data"][index]["sseqid"];
            }

            if( index_TE_INFOS[generation] != undefined && index_TE_INFOS[generation][id] != undefined ){
                class_a = "TE_infos";
            }
            else{
                class_a = "";
            }
            
            
            //color
            for(pos in dico["positions"]){
                seqs += `<div class="leftGr grpos">`;
                const red = dico["positions"][pos]["strand"] == 1 ? "darkblue" : "lightcoral";

                if ("no_match_send" in dico["positions"][pos]) {
                    seqs += `<div style="width: ${dico["positions"][pos]["no_match_size"][MODE_P] - 4}px; padding-left: ${dico["positions"][pos]["no_match_size"][MODE_P] - 4}px; margin-left: ${dico["positions"][pos]["posq"][MODE_P]}px;" class="h4 white"></div>`;
                    seqs += `<div class="red ${red}"> <a class="${class_a}" href="#" qseqid="${dico["qseqid"]}" sseqid="${dico["sseqid"]}" strand="${dico["positions"][pos]["strand"]}" sstart="${dico["positions"][pos]["sstart"]}" send="${dico["positions"][pos]["send"]}" qstart="${dico["positions"][pos]["qstart"]}" qend="${dico["positions"][pos]["qend"]}" query_seq_size="${dico["positions"][pos]["query_seq_size"]}" subject_seq_size="${dico["positions"][pos]["subject_seq_size"]}" gapopen="${dico["positions"][pos]["gapopen"]}" mismatch="${dico["positions"][pos]["mismatch"]}" evalue="${dico["positions"][pos]["evalue"]}" bitscore="${dico["positions"][pos]["bitscore"]}" pident="${dico["positions"][pos]["pident"]}" index="${index}" pos="${pos}" style="width: ${dico["positions"][pos]["size"][MODE_P]}px; "></a> </div>`;
                }
                else{
                    if("posq" in dico["positions"][pos]){
                        seqs += `<div class="red ${red}" style="margin-left: ${dico["positions"][pos]["posq"][MODE_P]}px;"> <a class="${class_a}" href="#" qseqid="${dico["qseqid"]}" sseqid="${dico["sseqid"]}" strand="${dico["positions"][pos]["strand"]}" sstart="${dico["positions"][pos]["sstart"]}" send="${dico["positions"][pos]["send"]}" qstart="${dico["positions"][pos]["qstart"]}" qend="${dico["positions"][pos]["qend"]}" query_seq_size="${dico["positions"][pos]["query_seq_size"]}" subject_seq_size="${dico["positions"][pos]["subject_seq_size"]}" index="${index}" pos="${pos}" gapopen="${dico["positions"][pos]["gapopen"]}" mismatch="${dico["positions"][pos]["mismatch"]}" evalue="${dico["positions"][pos]["evalue"]}" bitscore="${dico["positions"][pos]["bitscore"]}" pident="${dico["positions"][pos]["pident"]}" style="width: ${dico["positions"][pos]["size"][MODE_P]}px; "></a> </div>`
                    }
                    else{
                        seqs += `<div class="red ${red}"> <a class="${class_a}" href="#" qseqid="${dico["qseqid"]}" sseqid="${dico["sseqid"]}" strand="${dico["positions"][pos]["strand"]}" style="width: ${dico["positions"][pos]["size"][MODE_P]}px;"></a> </div>`
                    }
                }

                if("rest_size" in dico["positions"][pos] ){
                    seqs += `<div style="width: ${dico["positions"][pos]["rest_size"][MODE_P] - 4}px; padding-left: ${dico["positions"][pos]["rest_size"][MODE_P] - 4}px;" class="h4 white"></div>`;
                }
                seqs += `</div>`;
            }
        }   
    }

    return seqs;
}




function read_data_html_with_index_TE(data, chrom, TE_selected, MODE_P, generation) {
    var seqs   = "";
    var TE     = "";
    var query  = "";
    var id     = "";
    var index_chrom = 0;

    console.log("read with index TE")

    // rawdata    = fs.readFileSync(`index/index_TE_${generation}.json`);
    // index_TE   = JSON.parse(rawdata);

    index_TE = dico_generation_data[generation]["index_TE"];

    for(const i_TE in TE_selected){
        TE = TE_selected[i_TE];
        for(const i_indexTE in index_TE[TE]){

            var index = index_TE[TE][i_indexTE];

            var dico = data["data"][index];
            if ( ( chrom == dico["chrom"] || chrom == undefined ) ) {

                // 2L_RaGOO_RaGOO:<INS>:547321:547321:sniffles.INS.806:3:IMPRECISE:2
                if (id != data["data"][index]["qseqid"].split(":")[4]) {
                    start = new Intl.NumberFormat("de-DE").format(data["data"][index]["qseqid"].split(":")[2]);
                    end   = new Intl.NumberFormat("de-DE").format(data["data"][index]["qseqid"].split(":")[3]);

                    seqs += `<div class="id-separate" data-locus="${data["data"][index]["qseqid"].split(":")[0] + ":" + start + "-" + end}" data-index="${index_chrom}" style="position:relative">
                                <div class="dot" style="position:absolute; left:-10px"></div>
                                <div class="dot" style="position:absolute; left:240px"></div>
                                <div class="dot" style="position:absolute; left:490px"></div>
                            </div>`;
                    id    = data["data"][index]["qseqid"].split(":")[4];
                    index_chrom += 1;
                }

                if (query != data["data"][index]["qseqid"]) {
                    seqs += `<div class="query-separate"></div>`;
                    query = data["data"][index]["qseqid"];
                }


                if (TE != data["data"][index]["sseqid"]) {
                    seqs += `<div class="te-name"><a>${data["data"][index]["sseqid"]}</a></div>`;
                    TE = data["data"][index]["sseqid"];
                }
                
                if( index_TE_INFOS[generation] != undefined && index_TE_INFOS[generation][id] != undefined ){
                    class_a = "TE_infos";
                }
                else{
                    
                    class_a = "";
                }

                // console.log(generation, id, index_TE_INFOS[generation][id])

                //color
                for(pos in dico["positions"]){
                    seqs += `<div class="leftGr grpos">`;
                    const red = dico["positions"][pos]["strand"] == 1 ? "darkblue" : "lightcoral";

                    if ("no_match_send" in dico["positions"][pos]) {
                        seqs += `<div style="width: ${dico["positions"][pos]["no_match_size"][MODE_P] - 4}px; padding-left: ${dico["positions"][pos]["no_match_size"][MODE_P] - 4}px; margin-left: ${dico["positions"][pos]["posq"][MODE_P]}px;" class="h4 white"></div>`;
                        seqs += `<div class="red ${red}"> <a class="${class_a}" href="#" qseqid="${dico["qseqid"]}" sseqid="${dico["sseqid"]}" strand="${dico["positions"][pos]["strand"]}" sstart="${dico["positions"][pos]["sstart"]}" send="${dico["positions"][pos]["send"]}" qstart="${dico["positions"][pos]["qstart"]}" qend="${dico["positions"][pos]["qend"]}" query_seq_size="${dico["positions"][pos]["query_seq_size"]}" subject_seq_size="${dico["positions"][pos]["subject_seq_size"]}" gapopen="${dico["positions"][pos]["gapopen"]}" mismatch="${dico["positions"][pos]["mismatch"]}" evalue="${dico["positions"][pos]["evalue"]}" bitscore="${dico["positions"][pos]["bitscore"]}" pident="${dico["positions"][pos]["pident"]}" index="${index}" pos="${pos}" style="width: ${dico["positions"][pos]["size"][MODE_P]}px; "></a> </div>`;
                    }
                    else{
                        if("posq" in dico["positions"][pos]){
                            seqs += `<div class="red ${red}" style="margin-left: ${dico["positions"][pos]["posq"][MODE_P]}px;"> <a class="${class_a}" href="#" qseqid="${dico["qseqid"]}" sseqid="${dico["sseqid"]}" strand="${dico["positions"][pos]["strand"]}" sstart="${dico["positions"][pos]["sstart"]}" send="${dico["positions"][pos]["send"]}" qstart="${dico["positions"][pos]["qstart"]}" qend="${dico["positions"][pos]["qend"]}" query_seq_size="${dico["positions"][pos]["query_seq_size"]}" subject_seq_size="${dico["positions"][pos]["subject_seq_size"]}" index="${index}" pos="${pos}" gapopen="${dico["positions"][pos]["gapopen"]}" mismatch="${dico["positions"][pos]["mismatch"]}" evalue="${dico["positions"][pos]["evalue"]}" bitscore="${dico["positions"][pos]["bitscore"]}" pident="${dico["positions"][pos]["pident"]}" style="width: ${dico["positions"][pos]["size"][MODE_P]}px; "></a> </div>`
                        }
                        else{
                            seqs += `<div class="red ${red}"> <a class="${class_a}" href="#" qseqid="${dico["qseqid"]}" sseqid="${dico["sseqid"]}" strand="${dico["positions"][pos]["strand"]}" style="width: ${dico["positions"][pos]["size"][MODE_P]}px;"></a> </div>`
                        }
                    }

                    if("rest_size" in dico["positions"][pos] ){
                        seqs += `<div style="width: ${dico["positions"][pos]["rest_size"][MODE_P] - 4}px; padding-left: ${dico["positions"][pos]["rest_size"][MODE_P] - 4}px;" class="h4 white"></div>`;
                    }
                    seqs += `</div>`;
                }
            }   
        }
    }
    console.log("send seqs...")
    return seqs;
}


function getDataChrom(data, chrom, TE_selected) {
    var id = ""
    var index_chrom = 0;
    data_chrom = [];

    
    for (const index in data["data"]) {
        
        var dico = data["data"][index];
        if ( TE_in_selected(dico, TE_selected) && ( chrom == dico["chrom"] || chrom == '' ) ) {
            if (id != data["data"][index]["qseqid"].split(":")[4]) {
                start = new Intl.NumberFormat("de-DE").format(data["data"][index]["qseqid"].split(":")[2]);
                end   = new Intl.NumberFormat("de-DE").format(data["data"][index]["qseqid"].split(":")[3]);

                id    = data["data"][index]["qseqid"].split(":")[4];

                if( chrom != ''){
                    data_chrom.push(`${data["data"][index]["qseqid"].split(":")[0] + ":" + start + "-" + end}`);
                    index_chrom += 1;
                }
                else{
                    data_chrom = [];
                }
            }
        }

    }   
    return data_chrom;
}


function onlyUnique(value, index, self) {
    return self.indexOf(value) === index;
}


function create_chrom_index(generation) {
    console.log("BEGIN | create_chrom_index : ", generation);

    rawdata    = fs.readFileSync(`data/data_${generation}_sorted_by_position.json`);
    data       = JSON.parse(rawdata);

    pred_chrom = ""
    dico_index_chrom = {}
    dico_index_chrom_interval = {}
    dico_index_TE = {}
    for (const index in data["data"]) {
        
        var dico  = data["data"][index];
        var chrom = dico["chrom"];
        var TE    = dico["sseqid"];

        //chrom
        if( dico_index_chrom[chrom] === undefined ){
            dico_index_chrom[chrom] = [parseInt(index)];
        }
        else if(chrom in dico_index_chrom){
            dico_index_chrom[chrom].push(parseInt(index));
        }

        //chrom interval
        if( dico_index_chrom_interval[chrom] === undefined ){
            dico_index_chrom_interval[chrom] = [parseInt(index)];
            if( pred_chrom != "" ){
                dico_index_chrom_interval[pred_chrom].push(parseInt(index-1));
            }
            pred_chrom = chrom
        }
        else if(chrom in dico_index_chrom_interval){
            if(pred_chrom != chrom){

                dico_index_chrom_interval[pred_chrom].push(parseInt(index-1));
                dico_index_chrom_interval[chrom].push(parseInt(index));
                pred_chrom = chrom
            }
        }

        //TE
        if( dico_index_TE[TE] === undefined ){
            dico_index_TE[TE] = [parseInt(index)];
        }
        else if(TE in dico_index_TE){
            dico_index_TE[TE].push(parseInt(index));
        }

    }

    dico_index_chrom_interval[pred_chrom].push(parseInt(data["data"].length-1));

    //chrom
    fs.writeFileSync(`index/index_chrom_${generation}.json`, JSON.stringify(dico_index_chrom), err => {
        if (err) {
            console.error(err);
        }
    });

    //chrom interval
    fs.writeFileSync(`index/index_chrom_interval_${generation}.json`, JSON.stringify(dico_index_chrom_interval), err => {
        if (err) {
            console.error(err);
        }
    });

    //TE
    fs.writeFileSync(`index/index_TE_${generation}.json`, JSON.stringify(dico_index_TE), err => {
        if (err) {
            console.error(err);
        }
    });

    console.log("END | create_chrom_index : ", generation)
}


// Comparaison du bitscore pour le trie
function compare_bitscore(s1, s2) {
    max_s1 = Math.max.apply(Math, s1["positions"].map(function(pos) { return pos.bitscore; }));
    max_s2 = Math.max.apply(Math, s2["positions"].map(function(pos) { return pos.bitscore; }));
    
    if (max_s1 > max_s2) {
        // console.log("max_s1", max_s1, "max_s2", max_s2, "-1", s1.qseqid, s2.qseqid)
        return -1;
    }
    if (max_s1 < max_s2) {
        // console.log("max_s1", max_s1, "max_s2", max_s2, "1", s1.qseqid, s2.qseqid)
        return 1;
    }
    // console.log("max_s1", max_s1, "max_s2", max_s2, "0", s1.qseqid, s2.qseqid)
    return 0;
}



// création du data_set trié par bitscore
function create_data_sorted_by_bitscore(generation) {

    console.log("BEGIN | create_data_sorted_by_bitscore : ", generation);
    rawdata    = fs.readFileSync(`data/data_${generation}.json`);
    data       = JSON.parse(rawdata);

    data["data"] = data["data"].sort(compare_bitscore);
    
    //data_sorted_by_bitscore
    fs.writeFileSync(`data/data_${generation}_sorted_by_bitscore.json`, JSON.stringify(data), err => {
      if (err) {
        console.error(err);
      }
    });
    console.log("END | create_data_sorted_by_bitscore : ", generation);
}

function create_index_TE_INFOS(generation){

    console.log("BEGIN | create_index_TE_INFOS", generation)

    data = fs.readFileSync(`TE_INFOS/TE_INFOS_${generation}.json`)

    const jsonArray = JSON.parse(data);
    const index_TE_INFOS = {};

    jsonArray.forEach((obj, index) => {
      const id = obj["TE|ID"].split("|")[1];
      index_TE_INFOS[id] = index;
    });

    // console.log('Index :', index_TE_INFOS);
    console.log("END | create_index_TE_INFOS", generation)
    return index_TE_INFOS;
}


function create_index_DEPTH(generation){

    console.log("BEGIN | create_index_DEPTH", generation)

    data = fs.readFileSync(`DEPTH/DEPTH_${generation}.json`)

    const jsonArray = JSON.parse(data);
    const index_DEPTH = {};

    jsonArray.forEach((obj, index) => {
        const id = obj["info_TE"].split(":")[4];
        index_DEPTH[id] = index;
    });

    // console.log('Index :', index_DEPTH);
    console.log("END | create_index_DEPTH", generation)
    return index_DEPTH;
}

// Comparaison des position pour le trie
function compare_position(s1, s2) {
    chrom_s1 = s1["qseqid"].split(":")[0];
    chrom_s2 = s2["qseqid"].split(":")[0];

    start_s1 = parseInt(s1["qseqid"].split(":")[2]);
    start_s2 = parseInt(s2["qseqid"].split(":")[2]);

    if (chrom_s1 > chrom_s2) {
        return 1;
    }
    else if (chrom_s1 < chrom_s2) {
        return -1;
    }
    else {
        if (start_s1 > start_s2) {
            return 1;
        }
        
        if (start_s1 < start_s2) {
            return -1;
        }
        return 0;
    }
}

// création du data_set trié par bitscore
function create_data_sorted_by_position(generation) {

    console.log("BEGIN | create_data_sorted_by_position : ", generation);
    rawdata    = fs.readFileSync(`data/data_${generation}.json`);
    data       = JSON.parse(rawdata);

    data["data"] = data["data"].sort(compare_position);
    
    //data_sorted_by_bitscore
    fs.writeFileSync(`data/data_${generation}_sorted_by_position.json`, JSON.stringify(data), err => {
      if (err) {
        console.error(err);
      }
    });
    console.log("END | create_data_sorted_by_position : ", generation);
}

async function create_index(generation){

    console.log('\nSTEP-1')
    await fs.exists(`data/data_${generation}_sorted_by_bitscore.json`, function (doesExist) {  
        if (doesExist) {  
            console.log('\nfile sorted by bitscore founded', generation);  
        } else {  
            console.log('\nfile sorted by bitscore Not found', generation); 
            create_data_sorted_by_bitscore(generation);
        }  
    });

    console.log('\nSTEP-2')
    await fs.exists(`data/data_${generation}_sorted_by_position.json`, function (doesExist) {  
        if (doesExist) {  
            console.log('\nfile sorted by position founded', generation);  
        } else {  
            console.log('\nfile sorted by position Not found', generation); 
            create_data_sorted_by_position(generation);
        }  
    });

    console.log('\nSTEP-3')
    await fs.exists(`index/index_TE_${generation}.json`, function (doesExist) {  
        if (doesExist) {  
            console.log('\nfile index chrom founded', generation);  
        } else {  
            console.log('\nfile index  Not found', generation); 
            create_chrom_index(generation); 
        }  
    });

    console.log('\nSTEP-4')
    await fs.exists(`TE_INFOS/TE_INFOS_${generation}.json`, function (doesExist) {  
        if (doesExist) {  
            console.log('\nfile TE INFOS founded', generation);
            index_TE_INFOS[generation] = create_index_TE_INFOS(generation);
            console.log("size index_TE_INFOS", Object.keys(index_TE_INFOS[generation]).length)
        } else {  
            console.log('\nfile TE INFOS Not found', generation); 
        }  
    });

    console.log('\nSTEP-5')
    await fs.exists(`DEPTH/DEPTH_${generation}.json`, function (doesExist) {  
        if (doesExist) {  
            console.log('\nfile DEPTH existe', generation);
            index_DEPTH[generation] = create_index_DEPTH(generation);
            console.log("size index_DEPTH", Object.keys(index_DEPTH[generation]).length)
        } else {  
            console.log('\nfile DEPTH Not found', generation); 
        }  
    });

}


var dico_generation_data = {}
var dico_generation_index = {}
const requestListener = function (req, res) {
    
    // console.log(req.headers.origin)
    res.setHeader("Access-Control-Allow-Origin", "*");
    // headers = req.headers['access-control-request-headers'];
    
    res.setHeader("Access-Control-Allow-Headers", "content-type");

    let date = new Date().toJSON();
    console.log(`\n[${date}]`);

    console.log("req-url:", req.url)

    if( req.url.split('?').length > 1 ){
        
        generation = querystring.parse(req.url.split('?')[1])['gen'];
        if ( dico_generation_data[generation] === undefined ) {
            var rawdata = fs.readFileSync(`data/data_${generation}_sorted_by_bitscore.json`);
            var data    = JSON.parse(rawdata);

            dico_generation_data[generation] = {bitscore: data}
            
            var rawdata2 = fs.readFileSync(`data/data_${generation}_sorted_by_position.json`);
            var data2    = JSON.parse(rawdata2);

            dico_generation_data[generation]["position"] = data2;

            rawdata3   = fs.readFileSync(`index/index_TE_${generation}.json`);
            index_TE   = JSON.parse(rawdata3);

            dico_generation_data[generation]["index_TE"] = index_TE;

            console.log("data loaded. gen =", generation)
        }

        console.log("request:", req.url.split('?')[0], "params:", req.url.split('?')[1])

        switch (req.url.split('?')[0]) {
            case "/chrom":
                console.log("Asking chrom, gen =", generation)
                res.setHeader("Content-Type", "application/json");
                res.writeHead(200);

                rawdata_CHROM    = fs.readFileSync(`index/index_chrom_${generation}.json`);
                data_CHROM       = JSON.parse(rawdata_CHROM);
                keys_CHROM       = Object.keys(data_CHROM)

                console.log("keys_CHROM", keys_CHROM)

                //GET Chrom
                var array_Chrom = [];
                for (const index in keys_CHROM) {
                     array_Chrom.push(keys_CHROM[index]);
                }
                
                res.end(JSON.stringify(array_Chrom));
                break
            case "/TE":
                res.setHeader("Content-Type", "application/json");
                res.writeHead(200);

                rawdata_TE    = fs.readFileSync(`index/index_TE_${generation}.json`);
                data_TE       = JSON.parse(rawdata_TE);

                name = querystring.parse(req.url.split('?')[1])['name']

                keys_TE = Object.keys(data_TE)
                console.log("name TE:", name)

                //GET TE
                var array_TE = [];
                for (const index in keys_TE) {
                     array_TE.push(keys_TE[index]);
                }

                if (name == undefined) {
                    array_TE = array_TE.sort();
                }
                else if (typeof(name) == 'string') {
                    array_TE = array_TE.filter(word => word.toLowerCase().includes(name.toLowerCase())).sort()
                }
                else if (Array.isArray(name)) {                    
                    array_TE = array_TE.filter(word => 
                        {
                            for(const n in name){
                                if(word.toLowerCase().includes(name[n].toLowerCase())){
                                    return true;
                                }
                            }
                            return false;
                        }).sort();
                }

                res.end(JSON.stringify(array_TE.slice(0, 20)));
                break
            case "/data_chrom":
                res.setHeader("Content-Type", "application/json");
                res.writeHead(200);

                sort_b = querystring.parse(req.url.split('?')[1])['sort_bitscore'];
                chrom  = querystring.parse(req.url.split('?')[1])['chrom']
                TE     = querystring.parse(req.url.split('?')[1])['TE[]']
                MODE_P = querystring.parse(req.url.split('?')[1])['MODE_P'] != undefined ? querystring.parse(req.url.split('?')[1])['MODE_P'] : 1;
                
                if(chrom != undefined && chrom != ""){
                    console.log("/data_chrom == load data... gen =", generation)

                    data = sort_b == "true" ? dico_generation_data[generation]["bitscore"] : dico_generation_data[generation]["position"];

                    if(Array.isArray(TE)){
                        data_chrom = getDataChrom(data, chrom, TE);
                    }
                    else if(typeof(TE) == 'string'){
                        data_chrom = getDataChrom(data, chrom, [TE]);
                    }
                    else{
                        data_chrom = getDataChrom(data, chrom, []);
                    }

                    res.end(JSON.stringify(data_chrom));
                }
                else{
                    res.end(JSON.stringify({}));
                }
                break
            case "/data":
                res.setHeader("Content-Type", "text/html");
                res.writeHead(200);

                sort_b = querystring.parse(req.url.split('?')[1])['sort_bitscore'];
                chrom  = querystring.parse(req.url.split('?')[1])['chrom']
                TE     = querystring.parse(req.url.split('?')[1])['TE[]']
                MODE_P = querystring.parse(req.url.split('?')[1])['MODE_P'] != undefined && querystring.parse(req.url.split('?')[1])['MODE_P'] != "" ? querystring.parse(req.url.split('?')[1])['MODE_P'] : 1;
                
                if(chrom != undefined && chrom != ""){
                    console.log("/data == load data... gen =", generation)

                    data = sort_b == "true" ? dico_generation_data[generation]["bitscore"] : dico_generation_data[generation]["position"];

                    console.log("data loaded. gen =", generation)
                
                    if(Array.isArray(TE)){
                        if(TE.length > 0 && sort_b != "true"){
                            html = read_data_html_with_index_TE(data, chrom, TE, MODE_P, generation);
                        }
                        else{
                            html = read_data_html(data, chrom, TE, MODE_P, generation);
                        }
                    }
                    else if(typeof(TE) == 'string'){
                        if(TE != "" && sort_b != "true"){
                            html = read_data_html_with_index_TE(data, chrom, [TE], MODE_P, generation);
                        }
                        else{
                            html = read_data_html(data, chrom, [TE], MODE_P, generation);
                        }
                    }
                    else{
                        html = read_data_html(data, chrom, [], MODE_P, generation);
                    }

                    res.end(html);
                }
                else{
                    res.end("");
                }
                break
            case "/TE_INFOS":
                res.setHeader("Content-Type", "application/json");
                res.writeHead(200);

                ID = querystring.parse(req.url.split('?')[1])['ID'];

                console.log("ID:", ID)
                data = fs.readFileSync(`TE_INFOS/TE_INFOS_${generation}.json`)

                const jsonArrayTE_INFOS = JSON.parse(data);

                dataDEPTH = fs.readFileSync(`DEPTH/DEPTH_${generation}.json`)

                const jsonArrayDEPTH = JSON.parse(dataDEPTH);

                console.log(generation, Object.keys(index_TE_INFOS[generation]).length, jsonArrayTE_INFOS[index_TE_INFOS[generation][ID]])
                console.log("DEPTH:", JSON.stringify(jsonArrayDEPTH[index_DEPTH[generation][ID]]))

                res.end(JSON.stringify({"gen": generation, "ID": ID, "TE_INFOS": jsonArrayTE_INFOS[index_TE_INFOS[generation][ID]], "DEPTH": jsonArrayDEPTH[index_DEPTH[generation][ID]]}));

                break
        }
    }
    else{
        res.writeHead(500);
        res.end(JSON.stringify({'error':1}))
    }
    date = new Date().toJSON();
    console.log(`--[${date}]`);
    console.log("...end")
}

const server = http.createServer(requestListener);

server.listen(port, async () => {
    const directoryPath = path.join(__dirname, 'data');
    console.log("directoryPath:", directoryPath)

    await create_index("GEN");

    fs.readdir(directoryPath, function (err, files) {
        //handling error
        if (err) {
            console.log('Unable to scan directory: ' + err);
            return 
        } 
        //listing all files using forEach
        files.forEach(function (file) {
            // Do whatever you want to do with the file
            console.log(`Server is running on http://${host}:${port_client}?gen=${file.replace("data_", "").replace(".json", "").replace("_sorted_by_bitscore", "")}`); 
        });
    });
    
    date = new Date().toJSON();
    console.log(`[${date}]`);
    console.log("Begin...")
});



