const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
class Element {
  constructor() { this.children=[]; this.attributes={}; this.style={}; this.listeners={}; this.value=''; this.offsetWidth=300; this.offsetHeight=140; }
  replaceChildren(...children) { this.children=children; }
  appendChild(child) { this.children.push(child); return child; }
  append(...children) { this.children.push(...children); }
  setAttribute(k,v) { this.attributes[k]=v; }
  addEventListener(k,f) { (this.listeners[k] ||= []).push(f); }
  getBoundingClientRect() { return {left:50,bottom:80}; }
}
const nodes={};
for(const id of ['trm-detector-table-body','trm-detector-note','trm-raw-sv-summary','trm-genome-heatmap','trm-heatmap-window','trm-heatmap-tooltip','trm-heatmap-impact','trm-heatmap-maximum','trm-structure-identity','trm-structure-coverage','trm-structure-identity-value','trm-structure-coverage-value','trm-structure-preview','trm-structure-preview-note']) nodes[id]=new Element();
const calls=[
  {event_id:'sniffles.INS.1',event_type:'INS',source:'OUTSIDER',chrom:'chr1',anchor:10},
  {event_id:'TrEMOLO.INS.1',event_type:'INS_DEL',source:'OUTSIDER',chrom:'chr1',anchor:20},
  {event_id:'HARD.2',event_type:'HARD',source:'OUTSIDER',chrom:'chr1',anchor:75},
  {event_id:'Assemblytics_w_1',event_type:'DELETION',source:'INSIDER',chrom:'chr2',anchor:30},
  {event_id:'custom',event_type:'DEL',source:'UNKNOWN',chrom:'chr2',anchor:100},
];
const copies=[{target:'GENOME',copy_id:'a',identity:90,coverage:80},{target:'GENOME',copy_id:'b',identity:70,coverage:40},{target:'REFERENCE',copy_id:'a',identity:20,coverage:10}];
const relations=[{target:'GENOME',source_copy_id:'a',target_copy_id:'b',relation:'nested_candidate'},{target:'REFERENCE',source_copy_id:'a',target_copy_id:'a',relation:'same_family_overlap'}];
nodes['tremolo-report-data']={textContent:JSON.stringify({calls,chromosome_summary:[{chrom:'chr1',length:100},{chrom:'chr2',length:100}],sv_counts_available:true,sv_counts:{INS:40,DEL:10},resident_te:{copies,relations,thresholds:{},summary:{relation_types:{nested_candidate:1,same_family_overlap:1}}}})};
const document={getElementById:id=>nodes[id],createElement:()=>new Element(),createElementNS:()=>new Element()};
let source=fs.readFileSync(path.join(__dirname,'../../report/quarto/report.js'),'utf8');
source=source.slice(0,source.lastIndexOf('  document.title ='))+'globalThis.api={detectorCounts,detectorEvent,renderDetectorSummary,buildHeatmap,renderGenomeHeatmap,filteredStructureRelations,renderStructurePreview,state};})();';
const context={document,Intl,Map,Set,window:{innerWidth:1000,innerHeight:800}};
vm.runInNewContext(source,context);
const api=context.api;
assert.equal(api.detectorEvent(calls[1]),'Insertions');
const counts=api.detectorCounts(calls);
assert.equal(counts.get('TrEMOLO').Insertions,1);
assert.equal(counts.get('TrEMOLO').HARD,1);
assert.equal(counts.get('Other / unknown').Deletions,1);
api.renderDetectorSummary();
assert.equal(nodes['trm-detector-table-body'].children.length,4);
assert.equal(nodes['trm-raw-sv-summary'].children[0].children[1].textContent,'50');
let model=api.buildHeatmap(50);
assert.equal(model.chromosomes.length,2);
assert.equal(model.chromosomes[0].bins[0].all,2);
assert.equal(model.chromosomes[1].bins[1].all,1); // terminal insertion goes in final window
const maximum=model.maximum;
api.state.filtered=[api.state.filtered[0]];
model=api.buildHeatmap(50);
assert.equal(model.maximum,maximum);
assert.equal(model.chromosomes[0].bins[0].kept,1);
assert.equal(model.chromosomes[1].bins[1].kept,0);
nodes['trm-heatmap-window'].value='50';
api.renderGenomeHeatmap();
const cells=nodes['trm-genome-heatmap'].children[0].children.filter(n=>String(n.attributes.class).includes('trm-heatmap-cell'));
assert.equal(cells.length,4);
assert.ok(cells[0].attributes.class.includes('affected'));
assert.ok(cells[3].attributes.class.includes('affected'));
cells[0].listeners.pointerenter[0]({clientX:990,clientY:790});
assert.equal(nodes['trm-heatmap-tooltip'].hidden,false);
assert.match(nodes['trm-heatmap-tooltip'].children[3].textContent,/OUTSIDER: 1 selected \/ 2/);
assert.ok(parseInt(nodes['trm-heatmap-tooltip'].style.left)<700);
cells[0].listeners.keydown[0]({key:'Escape'});
assert.equal(nodes['trm-heatmap-tooltip'].hidden,true);
api.state.filtered=[];api.renderGenomeHeatmap();
assert.match(nodes['trm-heatmap-impact'].textContent,/4 \/ 4 windows filter affected/);
assert.equal(api.filteredStructureRelations(0,0).length,2);
assert.equal(api.filteredStructureRelations(65,20).length,1);
assert.equal(api.filteredStructureRelations(80,20).length,0); // both endpoints must pass
nodes['trm-structure-identity'].value='65';nodes['trm-structure-coverage'].value='20';
api.renderStructurePreview();
assert.match(nodes['trm-structure-preview-note'].textContent,/1 \/ 2 relation pairs retained/);
assert.equal(nodes['trm-structure-preview'].children[1].children[1].textContent,'1 / 1');
console.log('Detector, heatmap, tooltip and structure-preview checks passed');
