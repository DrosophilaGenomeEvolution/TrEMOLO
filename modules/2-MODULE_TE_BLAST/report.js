(() => {
  "use strict";
  const data = JSON.parse(document.getElementById("structure-report-data").textContent);
  const structures = data.structures || [];
  const events = data.events || [];
  const components = data.components || [];
  const matches = data.matches || [];
  const $ = (id) => document.getElementById(id);
  const ns = "http://www.w3.org/2000/svg";
  const integer = new Intl.NumberFormat();
  const decimal = new Intl.NumberFormat(undefined, { maximumFractionDigits: 2 });
  const colors = ["#087f78", "#d66a26", "#5868b2", "#b04f72", "#6b8e23", "#7651a1", "#2f7ebc", "#b58a21"];
  const colorByTe = new Map();
  [...new Set(matches.map((match) => match.subject_te))].sort().forEach((te, index) => colorByTe.set(te, colors[index % colors.length]));
  const matchesByQuery = new Map();
  const componentsByStructure = new Map();
  matches.forEach((match) => {
    if (!matchesByQuery.has(match.query_id)) matchesByQuery.set(match.query_id, []);
    matchesByQuery.get(match.query_id).push(match);
  });
  components.forEach((component) => {
    if (!componentsByStructure.has(component.structure_id)) componentsByStructure.set(component.structure_id, []);
    componentsByStructure.get(component.structure_id).push(component);
  });

  function svgElement(name, attributes = {}, text = "") {
    const node = document.createElementNS(ns, name);
    Object.entries(attributes).forEach(([key, value]) => node.setAttribute(key, value));
    if (text) node.textContent = text;
    return node;
  }

  function renderSummary() {
    const summary = data.summary || {};
    const values = [
      [summary.queries || 0, "insertion sequences"],
      [summary.hsps || 0, "raw HSPs retained in output"],
      [summary.retained_matches || 0, "normalized retained matches"],
      [summary.components || 0, "provisional components"],
      [summary.final_events || 0, "reported final events"],
    ];
    values.forEach(([value, label]) => {
      const card = document.createElement("div");
      card.className = "trm-structure-card";
      const strong = document.createElement("strong");
      strong.textContent = integer.format(value);
      const span = document.createElement("span");
      span.textContent = label;
      card.append(strong, span);
      $("trm-structure-summary").appendChild(card);
    });
  }

  function populateFilters() {
    [...new Set(structures.map((row) => row.chrom))].sort().forEach((value) => $("trm-structure-chrom").appendChild(new Option(value, value)));
    [...new Set(structures.map((row) => row.classification))].sort().forEach((value) => $("trm-structure-class").appendChild(new Option(value.replaceAll("_", " "), value)));
    ["trm-structure-chrom", "trm-structure-class", "trm-structure-final", "trm-structure-min-components"].forEach((id) => $(id).addEventListener("input", renderCatalogue));
  }

  function renderEvents() {
    const body = $("trm-event-body");
    const candidates = events.filter((event) => Number(event.multi_component_queries) > 0);
    candidates.forEach((event) => {
      const tr = document.createElement("tr");
      [`${event.chrom}:${integer.format(event.locus_start)}`, event.event_id, event.reported_te === "." ? "—" : event.reported_te, event.query_count, event.queries_with_components, event.multi_component_queries, event.multi_family_queries, event.maximum_component_count, event.classification.replaceAll("_", " ")].forEach((value) => {
        const td = document.createElement("td");
        td.textContent = value;
        tr.appendChild(td);
      });
      tr.addEventListener("click", () => {
        const selected = structures.filter((row) => row.event_id === event.event_id).sort((a, b) => Number(b.component_count) - Number(a.component_count))[0];
        if (selected) renderStructure(selected);
      });
      body.appendChild(tr);
    });
    if (!candidates.length) {
      const tr = document.createElement("tr");
      const td = document.createElement("td");
      td.colSpan = 9;
      td.textContent = "No event has more than one provisional component at the current thresholds.";
      tr.appendChild(td);
      body.appendChild(tr);
    }
  }

  function filteredStructures() {
    const chrom = $("trm-structure-chrom").value;
    const classification = $("trm-structure-class").value;
    const finalCall = $("trm-structure-final").value;
    const minimum = Math.max(0, Number($("trm-structure-min-components").value) || 0);
    return structures.filter((row) => (!chrom || row.chrom === chrom) && (!classification || row.classification === classification) && (!finalCall || row.final_call === finalCall) && Number(row.component_count) >= minimum);
  }

  function renderCatalogue() {
    const rows = filteredStructures();
    $("trm-structure-filter-note").textContent = `${integer.format(rows.length)} insertion structures match the filters.`;
    const body = $("trm-structure-body");
    body.replaceChildren();
    rows.forEach((row) => {
      const tr = document.createElement("tr");
      [`${row.chrom}:${integer.format(row.locus_start)}`, row.event_id, integer.format(row.query_length), row.reported_te === "." ? "—" : row.reported_te, row.retained_matches, row.component_count, row.alternative_matches, row.classification.replaceAll("_", " ")].forEach((value) => {
        const td = document.createElement("td");
        td.textContent = value;
        tr.appendChild(td);
      });
      tr.addEventListener("click", () => renderStructure(row));
      body.appendChild(tr);
    });
    if (rows.length) renderStructure(rows[0]);
  }

  function parseSegments(value) {
    if (!value || value === ".") return [];
    return value.split(";").map((segment) => segment.split("-").map(Number));
  }

  function renderStructure(structure) {
    const queryMatches = (matchesByQuery.get(structure.query_id) || []).slice().sort((a, b) => Number(a.query_start) - Number(b.query_start) || Number(b.bitscore_sum) - Number(a.bitscore_sum));
    const queryComponents = (componentsByStructure.get(structure.structure_id) || []).slice().sort((a, b) => Number(a.rank) - Number(b.rank));
    $("trm-structure-title").textContent = `${structure.event_id} · ${structure.chrom}:${integer.format(structure.locus_start)}`;
    $("trm-structure-note").textContent = `${structure.classification.replaceAll("_", " ")} · ${structure.component_count} proposed components and ${structure.alternative_matches} retained overlapping alternatives. Reported TE: ${structure.reported_te === "." ? "not a final call" : structure.reported_te}.`;
    renderDiagram(structure, queryMatches, queryComponents);
    renderMatches(queryMatches);
  }

  function renderDiagram(structure, queryMatches, queryComponents) {
    const svg = $("trm-structure-svg");
    svg.replaceChildren();
    const left = 170, right = 1160, queryY = 55;
    const length = Math.max(1, Number(structure.query_length));
    const scale = (value) => left + Number(value) / length * (right - left);
    svg.appendChild(svgElement("rect", { x: left, y: queryY, width: right - left, height: 22, rx: 5, class: "trm-st-query" }));
    svg.appendChild(svgElement("text", { x: 16, y: queryY + 17, class: "trm-st-label" }, `Insertion (${integer.format(length)} bp)`));
    [0, 0.25, 0.5, 0.75, 1].forEach((fraction) => {
      const x = scale(fraction * length);
      svg.appendChild(svgElement("line", { x1: x, x2: x, y1: queryY + 24, y2: 500, class: "trm-st-axis", opacity: 0.2 }));
      svg.appendChild(svgElement("text", { x, y: 96, "text-anchor": "middle", class: "trm-st-label" }, integer.format(fraction * length)));
    });
    queryComponents.forEach((component, index) => {
      const y = 125 + index * 44;
      svg.appendChild(svgElement("text", { x: 16, y: y + 16, class: "trm-st-label" }, `C${component.rank} ${component.te_name} ${component.strand}`));
      parseSegments(component.query_segments).forEach(([start, end]) => {
        const rectangle = svgElement("rect", { x: scale(start), y, width: Math.max(2, scale(end) - scale(start)), height: 22, rx: 3, fill: colorByTe.get(component.te_name), class: "trm-st-component" });
        rectangle.appendChild(svgElement("title", {}, `${component.te_name}: ${start}–${end} bp`));
        svg.appendChild(rectangle);
      });
    });
    const matchStart = 145 + Math.max(queryComponents.length, 1) * 44;
    queryMatches.slice(0, 12).forEach((match, index) => {
      const y = matchStart + index * 20;
      svg.appendChild(svgElement("text", { x: 16, y: y + 4, class: "trm-st-label" }, match.subject_te));
      parseSegments(match.query_segments).forEach(([start, end]) => svg.appendChild(svgElement("line", { x1: scale(start), x2: scale(end), y1: y, y2: y, class: match.status === "retained" ? "trm-st-match-retained" : "trm-st-match-rejected" })));
    });
    if (queryMatches.length > 12) svg.appendChild(svgElement("text", { x: 16, y: matchStart + 250, class: "trm-st-label" }, `+ ${queryMatches.length - 12} matches in table`));
  }

  function renderMatches(rows) {
    const body = $("trm-match-body");
    body.replaceChildren();
    rows.forEach((row) => {
      const tr = document.createElement("tr");
      [row.subject_te, `${integer.format(row.query_start)}–${integer.format(row.query_end)}`, row.subject_strand, `${decimal.format(row.weighted_identity)}%`, `${decimal.format(row.consensus_coverage)}%`, decimal.format(row.bitscore_sum), row.status, row.assignment.replaceAll("_", " "), row.filter_reasons === "." ? "—" : row.filter_reasons].forEach((value, index) => {
        const td = document.createElement("td");
        td.textContent = value;
        if (index === 6 && row.status === "rejected") td.className = "trm-st-rejected";
        if (index === 7 && row.assignment === "reported_primary") td.className = "trm-st-primary";
        tr.appendChild(td);
      });
      body.appendChild(tr);
    });
  }

  renderSummary();
  renderEvents();
  populateFilters();
  renderCatalogue();
})();
