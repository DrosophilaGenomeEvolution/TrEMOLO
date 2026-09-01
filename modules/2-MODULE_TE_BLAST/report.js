(() => {
  "use strict";

  const data = JSON.parse(document.getElementById("structure-report-data").textContent);
  const structures = data.structures || [];
  const events = data.events || [];
  const components = data.components || [];
  const matches = data.matches || [];
  const hsps = data.hsps || [];
  const $ = (id) => document.getElementById(id);
  const integer = new Intl.NumberFormat();
  const decimal = new Intl.NumberFormat(undefined, { maximumFractionDigits: 2 });
  const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: "base" });
  const matchesByQuery = groupBy(matches, "query_id");
  const hspsByQuery = groupBy(hsps, "query_id");
  const componentsByStructure = groupBy(components, "structure_id");
  const familiesByQuery = new Map();

  structures.forEach((structure) => {
    const families = new Set();
    if (structure.reported_te && structure.reported_te !== ".") families.add(structure.reported_te);
    (matchesByQuery.get(structure.query_id) || []).forEach((match) => families.add(match.subject_te));
    (componentsByStructure.get(structure.structure_id) || []).forEach((component) => families.add(component.te_name));
    familiesByQuery.set(structure.query_id, [...families].sort(collator.compare));
  });

  function groupBy(rows, field) {
    const grouped = new Map();
    rows.forEach((row) => {
      if (!grouped.has(row[field])) grouped.set(row[field], []);
      grouped.get(row[field]).push(row);
    });
    return grouped;
  }

  function element(name, className, text) {
    const node = document.createElement(name);
    if (className) node.className = className;
    if (text !== undefined && text !== null) node.textContent = text;
    return node;
  }

  function appendCell(row, value, className = "") {
    const cell = element("td", className, value);
    row.appendChild(cell);
    return cell;
  }

  function label(value) {
    return String(value || "").replaceAll("_", " ");
  }

  function percent(value) {
    return `${decimal.format(Number(value) || 0)}%`;
  }

  function scientific(value) {
    const number = Number(value);
    if (!Number.isFinite(number)) return "—";
    if (number === 0) return "0";
    return number.toExponential(2);
  }

  function teColor(te) {
    let hash = 0;
    for (const character of te) hash = ((hash << 5) - hash + character.codePointAt(0)) | 0;
    const hue = Math.abs(hash) % 360;
    return `hsl(${hue} 52% 43%)`;
  }

  function parseSegments(value) {
    if (!value || value === ".") return [];
    return value.split(";").map((segment) => segment.split("-").map(Number));
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
    values.forEach(([value, description]) => {
      const card = element("div", "trm-structure-card");
      card.append(element("strong", "", integer.format(value)), element("span", "", description));
      $("trm-structure-summary").appendChild(card);
    });
  }

  function appendOptions(id, values) {
    const select = $(id);
    values.forEach((value) => select.appendChild(new Option(label(value), value)));
  }

  function populateFilters() {
    appendOptions("trm-structure-chrom", [...new Set(structures.map((row) => row.chrom))].sort(collator.compare));
    appendOptions(
      "trm-structure-te",
      [...new Set(structures.flatMap((row) => familiesByQuery.get(row.query_id) || []))].sort(collator.compare),
    );
    appendOptions("trm-structure-class", [...new Set(structures.map((row) => row.classification))].sort(collator.compare));
    [
      "trm-structure-chrom",
      "trm-structure-te",
      "trm-structure-class",
      "trm-structure-final",
      "trm-structure-min-components",
      "trm-structure-search",
      "trm-structure-sort",
    ].forEach((id) => $(id).addEventListener("input", renderLandscape));
    $("trm-structure-reset").addEventListener("click", resetFilters);
  }

  function resetFilters() {
    ["trm-structure-chrom", "trm-structure-te", "trm-structure-class", "trm-structure-final", "trm-structure-search"]
      .forEach((id) => { $(id).value = ""; });
    $("trm-structure-min-components").value = "0";
    $("trm-structure-sort").value = "signal";
    renderLandscape();
  }

  function genomicOrder(first, second) {
    return collator.compare(first.chrom, second.chrom)
      || Number(first.locus_start) - Number(second.locus_start)
      || Number(first.locus_end) - Number(second.locus_end)
      || collator.compare(first.event_id, second.event_id);
  }

  function signalOrder(first, second) {
    const firstMatches = (matchesByQuery.get(first.query_id) || []).length;
    const secondMatches = (matchesByQuery.get(second.query_id) || []).length;
    return Number(second.component_count) - Number(first.component_count)
      || Number(second.alternative_matches) - Number(first.alternative_matches)
      || Number(second.retained_matches) - Number(first.retained_matches)
      || secondMatches - firstMatches
      || Number(second.final_call === "yes") - Number(first.final_call === "yes")
      || genomicOrder(first, second);
  }

  function filteredStructures() {
    const chrom = $("trm-structure-chrom").value;
    const te = $("trm-structure-te").value;
    const classification = $("trm-structure-class").value;
    const finalCall = $("trm-structure-final").value;
    const minimum = Math.max(0, Number($("trm-structure-min-components").value) || 0);
    const search = $("trm-structure-search").value.trim().toLocaleLowerCase();
    const order = $("trm-structure-sort").value;
    const rows = structures.filter((row) => {
      if (chrom && row.chrom !== chrom) return false;
      if (te && !(familiesByQuery.get(row.query_id) || []).includes(te)) return false;
      if (classification && row.classification !== classification) return false;
      if (finalCall && row.final_call !== finalCall) return false;
      if (Number(row.component_count) < minimum) return false;
      if (search) {
        const haystack = [
          row.chrom, row.locus_start, row.locus_end, row.event_id, row.query_id,
          row.structure_id, row.reported_te, row.classification,
          ...(familiesByQuery.get(row.query_id) || []),
        ].join(" ").toLocaleLowerCase();
        if (!haystack.includes(search)) return false;
      }
      return true;
    });
    if (order === "position") rows.sort(genomicOrder);
    else if (order === "components") rows.sort((a, b) => Number(b.component_count) - Number(a.component_count) || signalOrder(a, b));
    else if (order === "length") rows.sort((a, b) => Number(b.query_length) - Number(a.query_length) || genomicOrder(a, b));
    else rows.sort(signalOrder);
    return rows;
  }

  function badge(text, modifier = "") {
    return element("span", `trm-structure-badge ${modifier}`.trim(), text);
  }

  function structureTitle(structure, queryComponents, queryMatches) {
    const componentFamilies = [...new Set(queryComponents.map((component) => component.te_name))];
    if (structure.reported_te && structure.reported_te !== ".") {
      const distinctEvidence = componentFamilies.filter((family) => family !== structure.reported_te);
      return distinctEvidence.length
        ? `${structure.reported_te} · evidence: ${componentFamilies.join(" + ")}`
        : structure.reported_te;
    }
    if (componentFamilies.length) return componentFamilies.join(" + ");
    if (queryMatches.length) {
      const strongest = queryMatches.slice().sort((a, b) => Number(b.bitscore_sum) - Number(a.bitscore_sum))[0];
      return `Candidate: ${strongest.subject_te}`;
    }
    return "No TE signal";
  }

  function createStructureCard(structure) {
    const queryMatches = (matchesByQuery.get(structure.query_id) || []).slice()
      .sort((a, b) => Number(a.query_start) - Number(b.query_start) || Number(b.bitscore_sum) - Number(a.bitscore_sum));
    const queryComponents = (componentsByStructure.get(structure.structure_id) || []).slice()
      .sort((a, b) => Number(a.rank) - Number(b.rank));
    const queryHsps = hspsByQuery.get(structure.query_id) || [];
    const card = element("article", `trm-structure-item trm-class-${structure.classification}`);
    card.id = `structure-${structure.structure_id}`;
    card.setAttribute("role", "listitem");

    const header = element("header", "trm-structure-item-header");
    const identity = element("div", "trm-structure-identity");
    identity.append(
      element("p", "trm-structure-kicker", `${structure.chrom}:${integer.format(structure.locus_start)}–${integer.format(structure.locus_end)} · ${structure.event_id}`),
      element("h3", "", structureTitle(structure, queryComponents, queryMatches)),
      element("p", "trm-structure-query", `${integer.format(structure.query_length)} bp insertion sequence · ${structure.query_id}`),
    );
    const badges = element("div", "trm-structure-badges");
    badges.appendChild(badge(label(structure.classification), `trm-badge-${structure.classification}`));
    if (structure.final_call === "yes") badges.appendChild(badge("reported call", "trm-badge-final"));
    badges.appendChild(badge(`${structure.component_count} component${Number(structure.component_count) === 1 ? "" : "s"}`));
    if (Number(structure.alternative_matches)) badges.appendChild(badge(`${structure.alternative_matches} alternative${Number(structure.alternative_matches) === 1 ? "" : "s"}`, "trm-badge-alternative"));
    header.append(identity, badges);
    card.append(header, renderStructureMap(structure, queryMatches, queryComponents));

    const details = element("details", "trm-structure-details");
    const summary = element(
      "summary",
      "",
      `Evidence details · ${integer.format(queryMatches.length)} normalized match${queryMatches.length === 1 ? "" : "es"} · ${integer.format(queryHsps.length)} raw HSP${queryHsps.length === 1 ? "" : "s"}`,
    );
    const detailBody = element("div", "trm-structure-detail-body");
    details.append(summary, detailBody);
    details.addEventListener("toggle", () => {
      if (details.open && !detailBody.dataset.rendered) {
        renderDetails(detailBody, structure, queryMatches, queryHsps);
        detailBody.dataset.rendered = "yes";
      }
    });
    card.appendChild(details);
    return card;
  }

  function evidenceRow(structure, row, kind, index) {
    const isComponent = kind === "component";
    const te = isComponent ? row.te_name : row.subject_te;
    const strand = isComponent ? row.strand : row.subject_strand;
    const coverage = row.consensus_coverage;
    const identity = row.weighted_identity;
    const wrapper = element("div", `trm-structure-map-row trm-map-${kind}`);
    const rowLabel = element("div", "trm-structure-map-label");
    const prefix = isComponent ? `C${row.rank}` : kind === "alternative" ? "Alt" : "Candidate";
    rowLabel.append(
      element("strong", "", `${prefix} · ${te} ${strand}`),
      element("span", "", `${percent(coverage)} consensus · ${percent(identity)} identity`),
    );
    const track = element("div", "trm-structure-track");
    track.setAttribute("aria-label", `${te} coordinates on insertion sequence`);
    const length = Math.max(1, Number(structure.query_length));
    parseSegments(row.query_segments).forEach(([start, end]) => {
      const segment = element("span", "trm-structure-segment");
      const left = Math.max(0, Math.min(100, 100 * start / length));
      const width = Math.max(0.35, Math.min(100 - left, 100 * (end - start) / length));
      segment.style.left = `${left}%`;
      segment.style.width = `${width}%`;
      segment.style.setProperty("--trm-segment-color", teColor(te));
      if (width >= 12) segment.textContent = te;
      const status = isComponent ? "provisional component" : label(row.status);
      segment.title = `${te} · ${integer.format(start)}–${integer.format(end)} bp · ${strand} · ${status} · ${percent(coverage)} consensus coverage · ${percent(identity)} identity`;
      track.appendChild(segment);
    });
    wrapper.append(rowLabel, track);
    wrapper.style.setProperty("--trm-row-index", index);
    return wrapper;
  }

  function renderStructureMap(structure, queryMatches, queryComponents) {
    const map = element("section", "trm-structure-map");
    const axis = element("div", "trm-structure-axis");
    axis.append(element("span", "", "0"), element("strong", "", "Insertion sequence"), element("span", "", `${integer.format(structure.query_length)} bp`));
    map.appendChild(axis);

    if (queryComponents.length) {
      queryComponents.forEach((component, index) => map.appendChild(evidenceRow(structure, component, "component", index)));
      const componentMatchIds = new Set(queryComponents.map((component) => component.match_id));
      queryMatches
        .filter((match) => match.status === "retained" && !componentMatchIds.has(match.match_id))
        .slice(0, 3)
        .forEach((match, index) => map.appendChild(evidenceRow(structure, match, "alternative", queryComponents.length + index)));
    } else if (queryMatches.length) {
      queryMatches
        .slice()
        .sort((a, b) => Number(b.bitscore_sum) - Number(a.bitscore_sum))
        .slice(0, 3)
        .forEach((match, index) => map.appendChild(evidenceRow(structure, match, "rejected", index)));
      if (queryMatches.length > 3) map.appendChild(element("p", "trm-structure-more", `+ ${queryMatches.length - 3} additional candidates in evidence details`));
    } else {
      const empty = element("div", "trm-structure-map-empty", "No alignment to the TE database");
      empty.appendChild(element("span", "", "The insertion sequence remains visible and can be retained or excluded with filters."));
      map.appendChild(empty);
    }
    return map;
  }

  function detailFacts(structure) {
    const facts = [
      ["Structure ID", structure.structure_id],
      ["Query ID", structure.query_id],
      ["Reported TE", structure.reported_te === "." ? "Not a final call" : structure.reported_te],
      ["Interpretation", label(structure.interpretation)],
      ["Retained matches", structure.retained_matches],
      ["Provisional components", structure.component_count],
      ["Overlapping alternatives", structure.alternative_matches],
    ];
    const grid = element("dl", "trm-structure-facts");
    facts.forEach(([term, value]) => {
      const fact = document.createElement("div");
      fact.append(element("dt", "", term), element("dd", "", value));
      grid.appendChild(fact);
    });
    return grid;
  }

  function tableWithHeader(columns) {
    const wrap = element("div", "trm-structure-table-wrap trm-detail-table-wrap");
    const table = element("table", "trm-structure-table");
    const head = document.createElement("thead");
    const row = document.createElement("tr");
    columns.forEach((column) => row.appendChild(element("th", "", column)));
    head.appendChild(row);
    const body = document.createElement("tbody");
    table.append(head, body);
    wrap.appendChild(table);
    return { wrap, body };
  }

  function renderDetails(target, structure, queryMatches, queryHsps) {
    target.appendChild(detailFacts(structure));
    target.appendChild(element("h4", "", "Normalized TE matches"));
    if (queryMatches.length) {
      const matchTable = tableWithHeader([
        "TE", "Query segments", "Strand", "Query coverage", "Consensus coverage",
        "Identity", "HSPs", "Bitscore", "E-value", "Status", "Assignment", "Filter reason",
      ]);
      queryMatches.forEach((row) => {
        const tr = document.createElement("tr");
        [
          row.subject_te,
          row.query_segments,
          row.subject_strand,
          percent(row.query_coverage),
          percent(row.consensus_coverage),
          percent(row.weighted_identity),
          row.hsp_count,
          decimal.format(row.bitscore_sum),
          scientific(row.best_evalue),
          label(row.status),
          label(row.assignment),
          row.filter_reasons === "." ? "—" : label(row.filter_reasons),
        ].forEach((value) => appendCell(tr, value));
        if (row.status === "rejected") tr.className = "trm-st-rejected-row";
        matchTable.body.appendChild(tr);
      });
      target.appendChild(matchTable.wrap);
    } else {
      target.appendChild(element("p", "trm-structure-muted", "No normalized match was produced for this sequence."));
    }

    const raw = element("details", "trm-raw-evidence");
    raw.appendChild(element("summary", "", `Raw BLAST HSPs · ${queryHsps.length}`));
    if (queryHsps.length) {
      const hspTable = tableWithHeader([
        "TE", "Query interval", "Consensus interval", "Strand", "Aligned bp",
        "Identity", "Bitscore", "E-value", "Retained", "Filter reason",
      ]);
      queryHsps.forEach((row) => {
        const tr = document.createElement("tr");
        [
          row.subject_te,
          `${integer.format(row.query_start)}–${integer.format(row.query_end)}`,
          `${integer.format(row.subject_start)}–${integer.format(row.subject_end)}`,
          row.subject_strand,
          integer.format(row.alignment_length),
          percent(row.pident),
          decimal.format(row.bitscore),
          scientific(row.evalue),
          row.retained,
          row.filter_reasons === "." ? "—" : label(row.filter_reasons),
        ].forEach((value) => appendCell(tr, value));
        hspTable.body.appendChild(tr);
      });
      raw.appendChild(hspTable.wrap);
    } else {
      raw.appendChild(element("p", "trm-structure-muted", "No raw HSP for this insertion sequence."));
    }
    target.appendChild(raw);
  }

  function renderLandscape() {
    const rows = filteredStructures();
    $("trm-structure-filter-note").textContent = `${integer.format(rows.length)} / ${integer.format(structures.length)} insertion structures visible`;
    $("trm-structure-empty").hidden = rows.length !== 0;
    const fragment = document.createDocumentFragment();
    rows.forEach((row) => fragment.appendChild(createStructureCard(row)));
    $("trm-structure-gallery").replaceChildren(fragment);
  }

  function focusEvent(eventId) {
    ["trm-structure-chrom", "trm-structure-te", "trm-structure-class", "trm-structure-final"]
      .forEach((id) => { $(id).value = ""; });
    $("trm-structure-min-components").value = "0";
    $("trm-structure-search").value = eventId;
    $("trm-structure-sort").value = "components";
    renderLandscape();
    $("trm-structure-gallery").scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function renderEvents() {
    const body = $("trm-event-body");
    const candidates = events
      .filter((event) => Number(event.multi_component_queries) > 0)
      .sort((a, b) => Number(b.multi_family_queries) - Number(a.multi_family_queries) || Number(b.multi_component_queries) - Number(a.multi_component_queries) || genomicOrder(a, b));
    candidates.forEach((event) => {
      const tr = document.createElement("tr");
      [
        `${event.chrom}:${integer.format(event.locus_start)}`,
        event.event_id,
        event.reported_te === "." ? "—" : event.reported_te,
        event.query_count,
        event.queries_with_components,
        event.multi_component_queries,
        event.multi_family_queries,
        event.maximum_component_count,
        label(event.classification),
      ].forEach((value) => appendCell(tr, value));
      tr.tabIndex = 0;
      tr.setAttribute("role", "button");
      tr.title = `Show all structures for ${event.event_id}`;
      tr.addEventListener("click", () => focusEvent(event.event_id));
      tr.addEventListener("keydown", (keyboardEvent) => {
        if (keyboardEvent.key === "Enter" || keyboardEvent.key === " ") {
          keyboardEvent.preventDefault();
          focusEvent(event.event_id);
        }
      });
      body.appendChild(tr);
    });
    if (!candidates.length) {
      const tr = document.createElement("tr");
      const td = appendCell(tr, "No event has more than one provisional component at the current thresholds.");
      td.colSpan = 9;
      body.appendChild(tr);
    }
  }

  renderSummary();
  populateFilters();
  renderLandscape();
  renderEvents();
})();
