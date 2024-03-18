const http = require("http");
const fs = require('fs');
const querystring = require('node:querystring');
const yaml = require('js-yaml')

const host = 'localhost';
let port = 8031;
let port_server = 8030;

const url = require('url');


try {
    const yamlData = yaml.load(fs.readFileSync('../config.yaml', 'utf8'));
    port = yamlData.PORT_CLIENT;
    port_server = yamlData.PORT_SERVER;
} 
catch (error) {
    console.error(error);
}

console.log("port:", port, port_server);

// const requestListener = function (req, res) {
    
//     res.setHeader("Access-Control-Allow-Origin", "*");
//     res.setHeader("Access-Control-Allow-Headers", "*");

        
//     switch (req.url.split('?')[0]) {
//         case "/":
//             res.setHeader("Content-Type", "text/html");
//             res.writeHead(200);

//             html = fs.readFileSync(`client.html`);

//             res.write(html);  
//             res.end();  
//             break
//     }

// }

requestListener = function (req, res) {
   
    generation = querystring.parse(req.url.split('?')[1])['gen'];

    generation = generation == undefined ? "GEN" : generation;

    console.log("generation:", generation);

    const nameFile = "./client_vuetify3.html";

    fs.readFile(nameFile, 'utf-8', function (err, html) {
        if (err) {
            throw err; 
        }

        res.writeHeader(200, {"Content-Type": "text/html"}); 
        html = html.replaceAll('%%GENERATION%%', generation);
        html = html.replaceAll('%%PORT%%', port_server);
        res.write(html);  
        res.end();
    });
}

const server = http.createServer(requestListener);
server.listen(port, host, () => {
    console.log(`Server is running on http://${host}:${port}`);
});



