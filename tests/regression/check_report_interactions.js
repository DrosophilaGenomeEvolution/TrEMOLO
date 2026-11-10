// Execute the actual report startup/renderers with a DOM double and real report data.
const fs=require('node:fs'), vm=require('node:vm'), assert=require('node:assert/strict');
const path=require('node:path');
const root=path.resolve(__dirname,'../..');
if(!process.argv[2]) throw new Error('Usage: node tests/regression/check_report_interactions.js /path/REPORT');
const reportDir=path.resolve(process.argv[2]);
const nodes={};
class Element{
constructor(){this.children=[];this.dataset={};this.attributes={};this.style={};this.value='';this.listeners={};}
get firstElementChild(){return this.children[0];}
removeChild(c){this.children=this.children.filter(x=>x!==c);}
replaceChildren(...c){this.children=c;}
appendChild(c){this.children.push(c);return c;}
append(...c){this.children.push(...c);}
setAttribute(k,v){this.attributes[k]=v;}
addEventListener(k,f){(this.listeners[k]||=[]).push(f);}
after(c){nodes[c.id]=c;}
}
const qmd=fs.readFileSync(root+'/report/quarto/report.qmd','utf8');
for(const match of qmd.matchAll(/\bid="([^"]+)"/g)) nodes[match[1]]=new Element();
nodes['tremolo-report-data']={textContent:fs.readFileSync(path.join(reportDir,'report-data.json'),'utf8')};
nodes['trm-heatmap-window'].value='1000000';
nodes['trm-depth-guide'].value='10';nodes['trm-depth-frequency-guide'].value='50';
nodes['trm-timing-order'].value='duration';
const document={getElementById:id=>nodes[id]||null,createElement:()=>new Element(),createElementNS:()=>new Element(),title:''};
vm.runInNewContext(fs.readFileSync(root+'/report/quarto/report.js','utf8'),{document,Intl,Map,Set,window:{innerWidth:1280,innerHeight:900}}, {timeout:30000});
assert.ok(nodes['trm-detector-note'].textContent.includes('reported TE calls'));
assert.ok(nodes['trm-genome-heatmap'].children.length>0);
const catalogue=JSON.parse(nodes['tremolo-report-data'].textContent).resident_te;
if(catalogue && catalogue.available){
assert.ok(nodes['trm-structure-preview-note'].textContent.includes('relation pairs retained'));
for(const id of ['trm-structure-identity','trm-structure-coverage']) nodes[id].value='80';
nodes['trm-structure-identity'].listeners.input[0]();
assert.ok(nodes['trm-structure-preview-note'].textContent.includes('relation pairs retained'));
nodes['trm-structure-reset'].listeners.click[0]();
assert.ok(nodes['trm-structure-preview-note'].textContent.includes('relation pairs retained'));
}
nodes['trm-filter-family'].value='__no_matching_family__';
nodes['trm-filter-family'].listeners.change[0]();
assert.match(nodes['trm-heatmap-impact'].textContent,/filter affected/);
console.log('Complete report initialization and interaction checks passed');
