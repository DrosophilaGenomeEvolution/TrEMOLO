// Exercise the real SVG renderer, guide counts and missing-depth states.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
class Element {
  constructor() { this.children = []; this.attributes = {}; this.value = ''; }
  replaceChildren(...items) { this.children = items; }
  appendChild(child) { this.children.push(child); return child; }
  append(...items) { this.children.push(...items); }
  setAttribute(name, value) { this.attributes[name] = value; }
}
const nodes = {};
for (const id of ['trm-frequency-depth-chart', 'trm-depth-scale', 'trm-depth-guide', 'trm-depth-frequency-guide', 'trm-depth-summary']) nodes[id]=new Element();
nodes['trm-depth-scale'].value='linear';
nodes['trm-depth-guide'].value='10';
nodes['trm-depth-frequency-guide'].value='50';
const calls = [[2,100],[100,100],[2,10],[100,10]].map(([depth,freq],i) => ({
  family:'TE',event_id:`e${i}`,source:i%2?'INSIDER':'OUTSIDER',chrom:'chr1',start:i,event_type:'INS',
  display_frequency:freq,frequency_depth:depth,frequency_support_reads:depth*freq/100,
  frequency_depth_status:'matched',frequency_depth_method:'test',frequency_depth_table:'raw'
}));
calls.push({...calls[0],frequency_depth:0,frequency_depth_status:'zero_depth'});
calls.push({...calls[0],frequency_depth:null,frequency_depth_status:'unavailable'});
calls.push({...calls[0],display_frequency:null,frequency_depth_status:'frequency_unavailable'});
nodes['tremolo-report-data']={textContent:JSON.stringify({calls})};
const document={getElementById:id=>nodes[id],createElement:()=>new Element(),createElementNS:()=>new Element()};
let source=fs.readFileSync(path.join(__dirname,'../../report/quarto/report.js'),'utf8');
source=source.slice(0,source.lastIndexOf('  document.title ='))+'globalThis.api={renderFrequencyDepthChart,state};})();';
const context={document,Intl,Map,Set};
vm.runInNewContext(source,context);
const api=context.api;
const points=()=>nodes['trm-frequency-depth-chart'].children[0].children.filter(c=>String(c.attributes.class).includes('trm-frequency-depth-point'));
api.renderFrequencyDepthChart();
assert.equal(points().length,4);
assert.match(nodes['trm-depth-summary'].textContent,/2 frequency estimates omitted/);
assert.match(nodes['trm-depth-summary'].textContent,/1 without frequency/);
assert.match(nodes['trm-depth-summary'].textContent,/high frequency \/ low depth: 1/i);
assert.ok(Number(points()[0].attributes.cx)<Number(points()[1].attributes.cx));
assert.match(points()[0].children[0].textContent,/depth 2 reads; TE-present support 2/);
const linearX=points()[0].attributes.cx;
nodes['trm-depth-scale'].value='log';api.renderFrequencyDepthChart();
assert.notEqual(points()[0].attributes.cx,linearX);
assert.ok(points().every(p=>Number.isFinite(Number(p.attributes.cx))));
nodes['trm-depth-guide'].value='200';api.renderFrequencyDepthChart();
assert.match(nodes['trm-depth-summary'].textContent,/high frequency \/ low depth: 2/i);
api.state.filtered=[calls[1]];api.renderFrequencyDepthChart();
assert.equal(points().length,1);
api.state.filtered=[calls[4]];api.renderFrequencyDepthChart();
assert.match(nodes['trm-frequency-depth-chart'].children[0].textContent,/No positive, matching depth/);
console.log('Frequency/depth plot checks passed');
