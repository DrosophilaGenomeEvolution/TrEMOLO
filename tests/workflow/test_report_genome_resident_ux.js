const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
class Element {
  constructor() { this.children=[]; this.attributes={}; this.style={}; }
  replaceChildren(...children) { this.children=children; }
  appendChild(child) { this.children.push(child); return child; }
  append(...children) { this.children.push(...children); }
  setAttribute(k,v) { this.attributes[k]=v; }
}
const nodes={};
for (const id of ['trm-genome-chart','trm-resident-full-value','trm-resident-identity-value','trm-resident-partial-value','trm-resident-calibration-summary','trm-resident-calibration-context','trm-resident-discovery-floor','trm-resident-calibration-chart','bars']) nodes[id]=new Element();
const calls=[['long',50],['short',25],['unknown',30]].map(([chrom,anchor])=>({chrom,anchor,family:'TE',event_id:chrom,source:'OUTSIDER'}));
const copies=[[90,90],[40,90],[10,90],[90,60]].map(([coverage,identity])=>({coverage,identity}));
nodes['tremolo-report-data']={textContent:JSON.stringify({calls,chromosome_summary:[{chrom:'long',length:100},{chrom:'short',length:50},{chrom:'unknown',length:null}],resident_te:{copies,thresholds:{full_length_coverage:80,high_confidence_pident:80,partial_coverage:20,min_pident:65,min_aligned_bp:80}}})};
const document={getElementById:id=>nodes[id],createElement:()=>new Element(),createElementNS:()=>new Element()};
let source=fs.readFileSync(path.join(__dirname,'../../report/quarto/report.js'),'utf8');
source=source.slice(0,source.lastIndexOf('  document.title ='))+'globalThis.api={renderGenomeChart,renderResidentCalibration,renderSimpleBars,state,residentState};})();';
const context={document,Intl,Map,Set};vm.runInNewContext(source,context);
const api=context.api;
api.renderGenomeChart();
const tracks=()=>nodes['trm-genome-chart'].children[0].children.filter(n=>String(n.attributes.class).startsWith('trm-genome-line'));
let lines=tracks();
assert.equal((lines[0].attributes.x2-lines[0].attributes.x1)/(lines[1].attributes.x2-lines[1].attributes.x1),2);
assert.match(lines[0].children[0].textContent,/100 bp/);
assert.match(lines[2].attributes.class,/unknown/);
const shortEnd=lines[1].attributes.x2;
api.state.filtered=[calls[1]];api.renderGenomeChart();
assert.equal(tracks()[0].attributes.x2,shortEnd); // scale stays stable after filtering
api.renderResidentCalibration();
const cards=nodes['trm-resident-calibration-summary'].children;
assert.equal(cards.length,3);
assert.equal(cards[0].children[1].textContent,'1');
assert.equal(cards[1].children[1].textContent,'2'); // high coverage/low identity stays partial
assert.equal(cards[2].children[1].textContent,'1');
assert.match(cards[1].children[2].textContent,/50% of reported copies/);
assert.match(nodes['trm-resident-discovery-floor'].textContent,/65%/);
api.residentState.partialCoverage=50;api.renderResidentCalibration();
assert.equal(nodes['trm-resident-calibration-summary'].children[1].children[1].textContent,'1');
api.renderSimpleBars('bars',[{label:'same_family_overlap_candidate',value:2}]);
assert.equal(nodes.bars.children[0].children[0].title,'same_family_overlap_candidate');
console.log('Genome scale, resident summary and tooltip checks passed');
