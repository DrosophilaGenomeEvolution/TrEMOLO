// Run the actual module script and generated local chunks without browser packages.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");

const output = process.argv[2];
const qmd = fs.readFileSync(path.join(output, "report.qmd"), "utf8");
const bootstrap = qmd.split('<script id="structure-report-data" type="application/json">')[1].split("</script>")[0];
const source = fs.readFileSync(path.join(__dirname, "../../modules/2-MODULE_TE_BLAST/report.js"), "utf8");

function browser(fail = () => false) {
  const requests = [];
  class Element {
    constructor(tag = "div") {
      this.tagName = tag;
      this.children = [];
      this.listeners = {};
      this.attributes = {};
      this.dataset = {};
      this.style = { setProperty() {} };
      this.value = "";
      this.textContent = "";
    }
    appendChild(child) {
      if (child.tagName === "fragment") this.children.push(...child.children);
      else this.children.push(child);
      return child;
    }
    append(...children) { children.forEach((child) => this.appendChild(child)); }
    replaceChildren(...children) { this.children = []; this.append(...children); }
    setAttribute(key, value) { this.attributes[key] = value; }
    addEventListener(key, fn) { (this.listeners[key] ||= []).push(fn); }
    remove() {}
    scrollIntoView() {}
    fire(key) { (this.listeners[key] || []).forEach((fn) => fn({ key: "Enter", preventDefault() {} })); }
  }
  const nodes = {};
  for (const match of qmd.matchAll(/id="([^"]+)"/g)) nodes[match[1]] = new Element();
  nodes["structure-report-data"].textContent = bootstrap;
  nodes["trm-structure-sort"].value = "signal";
  nodes["trm-structure-min-components"].value = "0";
  let context;
  const document = {
    getElementById: (id) => nodes[id],
    createElement: (tag) => new Element(tag),
    createDocumentFragment: () => new Element("fragment"),
    head: { appendChild(script) {
      requests.push(script.src);
      queueMicrotask(() => {
        if (fail(script.src)) { script.onerror(); return; }
        vm.runInContext(fs.readFileSync(path.join(output, script.src), "utf8"), context);
        script.onload();
      });
    } },
  };
  context = vm.createContext({ document, window: {}, Intl, setTimeout, clearTimeout,
    Option: class extends Element { constructor(text, value) { super("option"); this.textContent = text; this.value = value; } },
  });
  vm.runInContext(source, context);
  return { nodes, requests };
}

async function until(condition) {
  const deadline = Date.now() + 30000;
  while (Date.now() < deadline) {
    if (condition()) return;
    await new Promise((resolve) => setTimeout(resolve, 5));
  }
  throw new Error("Report did not reach the expected state");
}
function button(node, text) {
  const found = node.children.find((child) => child.tagName === "button" && child.textContent === text);
  assert.ok(found, `Missing button: ${text}`);
  return found;
}
function jump(node, page) {
  node.children.find((child) => child.tagName === "input").value = String(page);
  button(node, "Go to page").fire("click");
}

(async () => {
  const { nodes, requests } = browser();
  const gallery = nodes["trm-structure-gallery"];
  const pages = nodes["trm-structure-pages"];
  await until(() => gallery.children.length === 50);
  if (process.argv.includes("--scale")) {
    const count = JSON.parse(bootstrap).summary.queries;
    assert.match(nodes["trm-structure-filter-note"].textContent, new RegExp(`^${new Intl.NumberFormat().format(count)} /`));
    assert.equal(gallery.children.length, 50);
    nodes["trm-structure-search"].value = `event${count}`;
    nodes["trm-structure-search"].fire("input");
    await until(() => /^1 \/ /.test(nodes["trm-structure-filter-note"].textContent) && gallery.children.length === 1);
    assert.match(gallery.children[0].children[0].children[0].children[1].textContent, /No TE signal/);
    console.log(`Scale browser harness passed: ${count} sequences, 50 cards per page, global search including no-match queries; heap ${Math.round(process.memoryUsage().heapUsed / 1024 / 1024)} MiB`);
    return;
  }
  assert.equal(nodes["trm-event-body"].children.length, 50);
  assert.equal(requests.filter((p) => p.includes("evidence-")).length, 1);
  assert.match(nodes["trm-structure-filter-note"].textContent, /251 \/ 251/);
  assert.match(gallery.children[0].children[0].children[0].children.at(-1).textContent,
               /Source genomic alignment: secondary · SAM flag: 256/);
  const initialDetails = gallery.children[0].children.at(-1);
  initialDetails.open = true;
  initialDetails.fire("toggle");
  const rawEvidence = initialDetails.children[1].children.at(-1);
  const rawTable = rawEvidence.children[1];
  const rawBody = rawTable.children[0].children[1];
  const rawPages = rawTable.children[1];
  assert.equal(rawBody.children.length, 50);
  button(rawPages, "Next").fire("click");
  assert.equal(rawBody.children.length, 10);
  assert.equal(button(rawPages, "Next").disabled, true);

  jump(pages, 5);
  await until(() => gallery.children.length === 50);
  assert.match(gallery.children[0].children[0].children[0].children[0].textContent, /event201$/);
  assert.equal(button(pages, "5").attributes["aria-current"], "page");
  button(pages, "Next").fire("click");
  await until(() => gallery.children.length === 1);
  assert.equal(button(pages, "Next").disabled, true);

  // Filters operate on the entire index, including queries outside the current page.
  nodes["trm-structure-search"].value = "event251";
  nodes["trm-structure-search"].fire("input");
  await until(() => /^1 \/ 251/.test(nodes["trm-structure-filter-note"].textContent) && gallery.children.length === 1);
  assert.match(nodes["trm-structure-filter-note"].textContent, /^1 \/ 251/);
  const cardDetails = gallery.children[0].children.at(-1);
  cardDetails.open = true;
  cardDetails.fire("toggle");
  assert.equal(cardDetails.children[1].dataset.rendered, "yes");
  assert.equal(cardDetails.children[1].children[2].children[0].children[1].children.length, 3);

  jump(nodes["trm-event-pages"], 5);
  assert.equal(nodes["trm-event-body"].children[0].children[1].textContent, "event201");
  jump(nodes["trm-event-pages"], 6);
  assert.equal(nodes["trm-event-body"].children.length, 1);
  nodes["trm-event-body"].children[0].fire("click");
  await until(() => gallery.children.length === 1);
  assert.equal(nodes["trm-structure-search"].value, "event251");

  nodes["trm-structure-search"].value = "does-not-exist";
  nodes["trm-structure-search"].fire("input");
  await until(() => gallery.children.length === 0 && nodes["trm-structure-status"].textContent === "");
  assert.equal(gallery.children.length, 0);
  assert.equal(nodes["trm-structure-empty"].hidden, false);
  nodes["trm-structure-reset"].fire("click");
  await until(() => gallery.children.length === 50);

  const missingEvidence = browser((file) => file.includes("evidence-"));
  await until(() => missingEvidence.nodes["trm-structure-status"].attributes.role === "alert");
  assert.match(missingEvidence.nodes["trm-structure-status"].textContent, /Keep report-assets/);
  const missingIndex = browser((file) => file.includes("index-"));
  await until(() => missingIndex.nodes["trm-structure-status"].attributes.role === "alert");
  assert.equal(missingIndex.nodes["trm-structure-filters"].disabled, true);
  console.log("Offline structure report loading, pagination, filters and error handling passed");
})().catch((error) => { console.error(error); process.exitCode = 1; });
