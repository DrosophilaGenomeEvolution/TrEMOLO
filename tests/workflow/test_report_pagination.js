// Exercise the actual report renderers without requiring browser dependencies.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");
const path = require("node:path");

class Element {
  constructor() { this.children = []; this.dataset = {}; this.listeners = {}; this.attributes = {}; }
  replaceChildren() { this.children = []; }
  appendChild(child) { this.children.push(child); return child; }
  setAttribute(key, value) { this.attributes[key] = value; }
  addEventListener(key, fn) { (this.listeners[key] ||= []).push(fn); }
  after(child) { nodes[child.id] = child; }
}
const nodes = {};
const document = {
  getElementById: (id) => nodes[id] || null,
  createElement: () => new Element(),
};
const candidates = Array.from({ length: 251 }, (_, i) => ({
  event_id: `event${i}`, source: "OUTSIDER", chrom: "chr1", start: i, end: i + 1,
  candidate_te: "A", reported_te: "A", assignment: "reported_primary", evidence_channels: ["direct"],
  evidence_fraction: 0.5, candidate_group_id: `group${i}`,
}));
nodes["tremolo-report-data"] = { textContent: JSON.stringify({
  calls: candidates.map((c) => ({ ...c, family: "A" })),
  ambiguous_calls: { available: true, candidates },
}) };
for (const id of ["trm-call-table-body", "trm-filter-count", "trm-page-status", "trm-page-previous", "trm-page-next",
  "trm-ambiguous-call-note", "trm-ambiguous-call-table-body", "trm-download-ambiguous-calls",
  "trm-download-candidate-evidence", "trm-ambiguous-page-status", "trm-ambiguous-page-previous", "trm-ambiguous-page-next"]) nodes[id] = new Element();
let source = fs.readFileSync(path.join(__dirname, "../../report/quarto/report.js"), "utf8");
// Skip unrelated chart initialization; retain the real renderer implementations.
source = source.slice(0, source.lastIndexOf("  document.title =")) +
  "globalThis.testReport = { renderTable, renderAmbiguousCalls, state, ambiguousState };})();";
const context = { document, Intl, Map, Set };
vm.runInNewContext(source, context);
const api = context.testReport;
api.renderTable();
api.renderAmbiguousCalls();
assert.equal(nodes["trm-call-table-body"].children.length, 50);
assert.equal(nodes["trm-ambiguous-call-table-body"].children.length, 50);

for (const prefix of ["trm-page", "trm-ambiguous-page"]) {
  const button = nodes[`${prefix}-numbers`].children.find((el) => el.textContent === 5);
  button.listeners.click[0]();
  assert.equal(nodes[`${prefix}-status`].textContent, "Page 5 / 6");
}
assert.equal(nodes["trm-call-table-body"].children[0].children[4].textContent, "event200");
assert.equal(nodes["trm-ambiguous-call-table-body"].children[0].children[1].textContent, "event200");
api.ambiguousState.page = 5;
api.renderAmbiguousCalls();
assert.equal(nodes["trm-ambiguous-call-table-body"].children.length, 1);
assert.equal(nodes["trm-ambiguous-page-next"].disabled, true);
assert.equal(nodes["trm-download-ambiguous-calls"].listeners.click.length, 1);
assert.equal(nodes["trm-download-candidate-evidence"].listeners.click.length, 1);
api.state.filtered = [];
api.renderTable();
assert.equal(nodes["trm-page-status"].textContent, "No result");
assert.equal(nodes["trm-page-next"].disabled, true);
console.log("Report pagination checks passed");
