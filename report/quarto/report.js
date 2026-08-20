(() => {
  "use strict";

  const node = document.getElementById("tremolo-report-data");
  if (!node) return;
  const data = JSON.parse(node.textContent);
  const calls = data.calls || [];
  const PAGE_SIZE = 50;
  const EVENT_TYPE_COLORS = [
    "#087f78", "#6f5bd3", "#ed8a3b", "#c44f70", "#3d7db7",
    "#a06b20", "#6b8e23", "#8b5e83", "#4f737b", "#b24b3e",
  ];
  const state = {
    filtered: calls.slice(),
    page: 0,
    frequencyPositionChrom: null,
    hiddenFrequencyPositionTypes: new Set(),
  };
  const $ = (id) => document.getElementById(id);
  const formatInteger = new Intl.NumberFormat("en-US", { maximumFractionDigits: 0 });
  const formatDecimal = new Intl.NumberFormat("en-US", { maximumFractionDigits: 4 });

  function valueOrDash(value, suffix = "") {
    return value === null || value === undefined ? "—" : `${formatDecimal.format(value)}${suffix}`;
  }

  function appendTextCell(row, value) {
    const cell = document.createElement("td");
    cell.textContent = value === null || value === undefined || value === "" ? "—" : String(value);
    row.appendChild(cell);
    return cell;
  }

  function selectOptions(element, values, allLabel) {
    element.replaceChildren();
    const all = document.createElement("option");
    all.value = "";
    all.textContent = allLabel;
    element.appendChild(all);
    values.forEach((value) => {
      const option = document.createElement("option");
      option.value = value;
      option.textContent = value;
      element.appendChild(option);
    });
  }

  function unique(field) {
    return [...new Set(calls.map((call) => call[field]))].sort((a, b) => String(a).localeCompare(String(b)));
  }

  function initializeFilters() {
    selectOptions($("trm-filter-source"), unique("source"), "All sources");
    selectOptions($("trm-filter-type"), unique("event_type"), "All event types");
    selectOptions($("trm-filter-chrom"), unique("chrom"), "All chromosomes");
    selectOptions($("trm-filter-family"), unique("family"), "All TE families");
    ["source", "type", "chrom", "family", "tsd", "frequency"].forEach((key) => {
      $(`trm-filter-${key}`).addEventListener("change", applyFilters);
    });
    $("trm-filter-search").addEventListener("input", applyFilters);
    $("trm-reset-filters").addEventListener("click", resetFilters);
    $("trm-download-calls").addEventListener("click", downloadCalls);
    $("trm-page-previous").addEventListener("click", () => changePage(-1));
    $("trm-page-next").addEventListener("click", () => changePage(1));
    $("trm-frequency-position-chrom").addEventListener("change", (event) => {
      state.frequencyPositionChrom = event.target.value;
      state.hiddenFrequencyPositionTypes.clear();
      renderFrequencyPositionChart();
    });
  }

  function applyFilters() {
    const source = $("trm-filter-source").value;
    const type = $("trm-filter-type").value;
    const chrom = $("trm-filter-chrom").value;
    const family = $("trm-filter-family").value;
    const tsd = $("trm-filter-tsd").value;
    const frequency = $("trm-filter-frequency").value;
    const search = $("trm-filter-search").value.trim().toLocaleLowerCase();
    state.filtered = calls.filter((call) => {
      if (source && call.source !== source) return false;
      if (type && call.event_type !== type) return false;
      if (chrom && call.chrom !== chrom) return false;
      if (family && call.family !== family) return false;
      if (tsd === "present" && call.tsd === null) return false;
      if (tsd === "missing" && call.tsd !== null) return false;
      if (frequency === "present" && call.display_frequency === null) return false;
      if (frequency === "missing" && call.display_frequency !== null) return false;
      if (search) {
        const haystack = [call.family, call.event_id, call.tremolo_id, call.chrom, call.event_type]
          .join(" ").toLocaleLowerCase();
        if (!haystack.includes(search)) return false;
      }
      return true;
    });
    state.page = 0;
    state.hiddenFrequencyPositionTypes.clear();
    renderSelection();
  }

  function resetFilters() {
    ["source", "type", "chrom", "family", "tsd", "frequency"].forEach((key) => {
      $(`trm-filter-${key}`).value = "";
    });
    $("trm-filter-search").value = "";
    applyFilters();
  }

  function summaryCard(label, value, detail) {
    const card = document.createElement("article");
    card.className = "trm-card";
    const labelNode = document.createElement("span");
    labelNode.className = "trm-card-label";
    labelNode.textContent = label;
    const valueNode = document.createElement("strong");
    valueNode.className = "trm-card-value";
    valueNode.textContent = value;
    const detailNode = document.createElement("span");
    detailNode.className = "trm-card-detail";
    detailNode.textContent = detail;
    card.append(labelNode, valueNode, detailNode);
    return card;
  }

  function renderCards() {
    const target = $("trm-summary-cards");
    target.replaceChildren();
    const selected = state.filtered;
    const families = new Set(selected.map((call) => call.family)).size;
    const chromosomes = new Set(selected.map((call) => call.chrom)).size;
    const tsd = selected.filter((call) => call.tsd !== null).length;
    const frequency = selected.filter((call) => call.display_frequency !== null);
    const insider = selected.filter((call) => call.source === "INSIDER").length;
    const outsider = selected.filter((call) => call.source === "OUTSIDER").length;
    target.append(
      summaryCard("Detected calls", formatInteger.format(selected.length), `${insider} INSIDER · ${outsider} OUTSIDER`),
      summaryCard("TE families", formatInteger.format(families), "distinct family labels"),
      summaryCard("Chromosomes", formatInteger.format(chromosomes), "with at least one selected call"),
      summaryCard("TSD confirmed", formatInteger.format(tsd), selected.length ? `${(100 * tsd / selected.length).toFixed(1)}% of calls` : "no selected calls"),
      summaryCard("Frequency available", formatInteger.format(frequency.length), frequency.length ? `mean ${valueOrDash(frequency.reduce((sum, call) => sum + call.display_frequency, 0) / frequency.length, "%")}` : "no estimate")
    );
  }

  function emptyChart(target, message) {
    const empty = document.createElement("div");
    empty.className = "trm-empty";
    empty.textContent = message;
    target.replaceChildren(empty);
  }

  function renderFamilyChart() {
    const target = $("trm-family-chart");
    const grouped = new Map();
    state.filtered.forEach((call) => {
      if (!grouped.has(call.family)) grouped.set(call.family, { family: call.family, INSIDER: 0, OUTSIDER: 0, UNKNOWN: 0 });
      grouped.get(call.family)[call.source] += 1;
    });
    const rows = [...grouped.values()]
      .map((row) => ({ ...row, total: row.INSIDER + row.OUTSIDER + row.UNKNOWN }))
      .sort((a, b) => b.total - a.total || a.family.localeCompare(b.family))
      .slice(0, 25);
    if (!rows.length) return emptyChart(target, "No family in the current selection");
    const maximum = Math.max(...rows.map((row) => row.total));
    target.replaceChildren();
    rows.forEach((item) => {
      const row = document.createElement("div");
      row.className = "trm-bar-row";
      const label = document.createElement("span");
      label.className = "trm-bar-label";
      label.textContent = item.family;
      label.title = item.family;
      const track = document.createElement("span");
      track.className = "trm-bar-track";
      track.style.width = `${Math.max(4, 100 * item.total / maximum)}%`;
      [["INSIDER", "trm-bar-insider"], ["OUTSIDER", "trm-bar-outsider"], ["UNKNOWN", "trm-bar-missing"]].forEach(([source, className]) => {
        if (!item[source]) return;
        const bar = document.createElement("span");
        bar.className = `trm-bar ${className}`;
        bar.style.width = `${100 * item[source] / item.total}%`;
        bar.title = `${source}: ${item[source]}`;
        track.appendChild(bar);
      });
      const value = document.createElement("span");
      value.className = "trm-bar-value";
      value.textContent = item.total;
      row.append(label, track, value);
      target.appendChild(row);
    });
  }

  function renderFrequencyChart() {
    const target = $("trm-frequency-chart");
    const values = state.filtered.map((call) => call.display_frequency).filter((value) => value !== null);
    if (!values.length) return emptyChart(target, "No frequency estimate in the current selection");
    const bins = Array(10).fill(0);
    values.forEach((value) => {
      const index = Math.min(9, Math.max(0, Math.floor(Number(value) / 10)));
      bins[index] += 1;
    });
    const maximum = Math.max(...bins, 1);
    const histogram = document.createElement("div");
    histogram.className = "trm-histogram";
    bins.forEach((count, index) => {
      const bin = document.createElement("div");
      bin.className = "trm-bin";
      const countNode = document.createElement("span");
      countNode.className = "trm-bin-count";
      countNode.textContent = count;
      const bar = document.createElement("span");
      bar.className = "trm-bin-bar";
      bar.style.height = `${Math.max(2, 100 * count / maximum)}%`;
      const label = document.createElement("span");
      label.textContent = `${index * 10}–${(index + 1) * 10}`;
      bin.append(countNode, bar, label);
      histogram.appendChild(bin);
    });
    target.replaceChildren(histogram);
  }

  function svgElement(name, attributes = {}) {
    const element = document.createElementNS("http://www.w3.org/2000/svg", name);
    Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, value));
    return element;
  }

  function eventTypeColorMap() {
    const types = [...new Set(calls.map((call) => call.event_type))]
      .sort((a, b) => String(a).localeCompare(String(b)));
    return new Map(types.map((eventType, index) => [
      eventType,
      EVENT_TYPE_COLORS[index % EVENT_TYPE_COLORS.length],
    ]));
  }

  function compactPosition(value) {
    if (Math.abs(value) >= 1e9) return `${formatDecimal.format(value / 1e9)}G`;
    if (Math.abs(value) >= 1e6) return `${formatDecimal.format(value / 1e6)}M`;
    if (Math.abs(value) >= 1e3) return `${formatDecimal.format(value / 1e3)}k`;
    return formatInteger.format(value);
  }

  function renderFrequencyPositionLegend(types, colors) {
    const target = $("trm-frequency-position-legend");
    target.replaceChildren();
    types.forEach((eventType) => {
      const button = document.createElement("button");
      const hidden = state.hiddenFrequencyPositionTypes.has(eventType);
      button.type = "button";
      button.className = `trm-legend-button${hidden ? " trm-legend-button-hidden" : ""}`;
      button.setAttribute("aria-pressed", hidden ? "false" : "true");
      button.title = `${hidden ? "Show" : "Hide"} ${eventType}`;
      const marker = document.createElement("i");
      marker.style.background = colors.get(eventType);
      const label = document.createElement("span");
      label.textContent = eventType;
      button.append(marker, label);
      button.addEventListener("click", () => {
        if (state.hiddenFrequencyPositionTypes.has(eventType)) {
          state.hiddenFrequencyPositionTypes.delete(eventType);
        } else {
          state.hiddenFrequencyPositionTypes.add(eventType);
        }
        renderFrequencyPositionChart();
      });
      target.appendChild(button);
    });
  }

  function renderFrequencyPositionChart() {
    const target = $("trm-frequency-position-chart");
    const selector = $("trm-frequency-position-chrom");
    const frequencyCalls = state.filtered.filter((call) => call.display_frequency !== null);
    const chromosomes = [...new Set(frequencyCalls.map((call) => call.chrom))]
      .sort((a, b) => String(a).localeCompare(String(b)));

    if (!chromosomes.length) {
      selector.replaceChildren();
      selector.disabled = true;
      $("trm-frequency-position-legend").replaceChildren();
      state.frequencyPositionChrom = null;
      return emptyChart(target, "No frequency estimate in the current selection");
    }

    selector.disabled = false;
    if (!chromosomes.includes(state.frequencyPositionChrom)) {
      const globallySelectedChromosome = $("trm-filter-chrom").value;
      state.frequencyPositionChrom = chromosomes.includes(globallySelectedChromosome)
        ? globallySelectedChromosome
        : chromosomes[0];
      state.hiddenFrequencyPositionTypes.clear();
    }
    selectOptions(selector, chromosomes, "");
    selector.removeChild(selector.firstElementChild);
    selector.value = state.frequencyPositionChrom;

    const chromosomeCalls = frequencyCalls.filter(
      (call) => call.chrom === state.frequencyPositionChrom
    );
    const types = [...new Set(chromosomeCalls.map((call) => call.event_type))]
      .sort((a, b) => String(a).localeCompare(String(b)));
    const colors = eventTypeColorMap();
    renderFrequencyPositionLegend(types, colors);
    const visibleCalls = chromosomeCalls.filter(
      (call) => !state.hiddenFrequencyPositionTypes.has(call.event_type)
    );
    if (!visibleCalls.length) {
      return emptyChart(target, "All event types are hidden for this chromosome");
    }

    const width = 1100;
    const height = 420;
    const margin = { top: 22, right: 25, bottom: 54, left: 72 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;
    const chromosomeInfo = (data.chromosome_summary || []).find(
      (item) => item.chrom === state.frequencyPositionChrom
    );
    const observedMaximum = Math.max(...chromosomeCalls.map((call) => call.start), 1);
    const xMaximum = Math.max(chromosomeInfo && chromosomeInfo.length || 0, observedMaximum, 1);
    const x = (value) => margin.left + plotWidth * Math.min(1, Math.max(0, value / xMaximum));
    const y = (value) => margin.top + plotHeight * (1 - Math.min(100, Math.max(0, value)) / 100);
    const svg = svgElement("svg", {
      viewBox: `0 0 ${width} ${height}`,
      role: "img",
      "aria-label": `TE frequency by position on ${state.frequencyPositionChrom}`,
    });

    [0, 25, 50, 75, 100].forEach((value) => {
      const gridY = y(value);
      svg.appendChild(svgElement("line", {
        x1: margin.left, y1: gridY, x2: width - margin.right, y2: gridY,
        class: "trm-scatter-grid",
      }));
      const label = svgElement("text", {
        x: margin.left - 12, y: gridY + 4,
        class: "trm-scatter-axis-label", "text-anchor": "end",
      });
      label.textContent = `${value}%`;
      svg.appendChild(label);
    });
    for (let index = 0; index <= 5; index += 1) {
      const value = xMaximum * index / 5;
      const tickX = x(value);
      svg.appendChild(svgElement("line", {
        x1: tickX, y1: margin.top, x2: tickX, y2: height - margin.bottom,
        class: "trm-scatter-grid trm-scatter-grid-vertical",
      }));
      const label = svgElement("text", {
        x: tickX, y: height - margin.bottom + 22,
        class: "trm-scatter-axis-label", "text-anchor": "middle",
      });
      label.textContent = compactPosition(value);
      svg.appendChild(label);
    }

    const xTitle = svgElement("text", {
      x: margin.left + plotWidth / 2, y: height - 8,
      class: "trm-scatter-axis-title", "text-anchor": "middle",
    });
    xTitle.textContent = `Position on ${state.frequencyPositionChrom} (bp)`;
    const yTitle = svgElement("text", {
      x: 17, y: margin.top + plotHeight / 2,
      class: "trm-scatter-axis-title", "text-anchor": "middle",
      transform: `rotate(-90 17 ${margin.top + plotHeight / 2})`,
    });
    yTitle.textContent = "Frequency (%)";
    svg.append(xTitle, yTitle);

    visibleCalls.forEach((call) => {
      const point = svgElement("circle", {
        cx: x(call.start),
        cy: y(call.display_frequency),
        r: 5,
        fill: colors.get(call.event_type),
        class: "trm-frequency-position-point",
        tabindex: "0",
      });
      const title = svgElement("title");
      title.textContent = [
        `${call.family} · ${call.event_id}`,
        `${call.event_type} · ${call.source} · strand ${call.strand}`,
        `${call.chrom}:${formatInteger.format(call.start)}`,
        `frequency ${valueOrDash(call.display_frequency, "%")}`,
      ].join("\n");
      point.appendChild(title);
      svg.appendChild(point);
    });
    target.replaceChildren(svg);
  }

  function renderGenomeChart() {
    const target = $("trm-genome-chart");
    const selectedChroms = new Set(state.filtered.map((call) => call.chrom));
    const chromosomes = data.chromosome_summary.filter((item) => selectedChroms.has(item.chrom));
    if (!chromosomes.length) return emptyChart(target, "No genomic position in the current selection");
    const width = 1100;
    const left = 150;
    const right = 25;
    const rowHeight = 34;
    const height = 35 + chromosomes.length * rowHeight;
    const svg = svgElement("svg", { viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": "TE calls along chromosomes" });
    chromosomes.forEach((item, index) => {
      const y = 25 + index * rowHeight;
      const label = svgElement("text", { x: 5, y: y + 4, class: "trm-genome-label" });
      label.textContent = item.chrom;
      svg.appendChild(label);
      svg.appendChild(svgElement("line", { x1: left, y1: y, x2: width - right, y2: y, class: "trm-genome-line" }));
      const chromosomeCalls = state.filtered.filter((call) => call.chrom === item.chrom);
      const inferredLength = Math.max(...chromosomeCalls.map((call) => call.anchor), 1);
      const length = item.length || inferredLength;
      chromosomeCalls.forEach((call) => {
        const x = left + (width - left - right) * Math.min(1, Math.max(0, call.anchor / length));
        const point = svgElement("circle", {
          cx: x,
          cy: y,
          r: 4.2,
          fill: call.source === "INSIDER" ? "#6f5bd3" : call.source === "OUTSIDER" ? "#ed8a3b" : "#a8b4b9",
          class: "trm-call-point",
        });
        const title = svgElement("title");
        title.textContent = `${call.family} · ${call.event_id} · ${call.chrom}:${formatInteger.format(call.anchor)}`;
        point.appendChild(title);
        svg.appendChild(point);
      });
    });
    target.replaceChildren(svg);
  }

  function renderSimpleBars(targetId, rows) {
    const target = $(targetId);
    const total = rows.reduce((sum, row) => sum + row.value, 0);
    if (!total) return emptyChart(target, "No call in the current selection");
    target.replaceChildren();
    rows.forEach((item) => {
      const row = document.createElement("div");
      row.className = "trm-bar-row";
      const label = document.createElement("span");
      label.className = "trm-bar-label";
      label.textContent = item.label;
      const track = document.createElement("span");
      track.className = "trm-bar-track";
      const bar = document.createElement("span");
      bar.className = `trm-bar ${item.className || "trm-bar-neutral"}`;
      bar.style.width = `${100 * item.value / total}%`;
      track.appendChild(bar);
      const value = document.createElement("span");
      value.className = "trm-bar-value";
      value.textContent = item.value;
      row.append(label, track, value);
      target.appendChild(row);
    });
  }

  function renderEvidence() {
    const confirmed = state.filtered.filter((call) => call.tsd !== null).length;
    renderSimpleBars("trm-tsd-summary", [
      { label: "Confirmed", value: confirmed, className: "trm-bar-neutral" },
      { label: "Not found", value: state.filtered.length - confirmed, className: "trm-bar-missing" },
    ]);
    const strands = { "+": 0, "-": 0, ".": 0 };
    state.filtered.forEach((call) => { strands[call.strand] = (strands[call.strand] || 0) + 1; });
    renderSimpleBars("trm-strand-summary", [
      { label: "Sense (+)", value: strands["+"] || 0 },
      { label: "Antisense (−)", value: strands["-"] || 0, className: "trm-bar-outsider" },
      { label: "Unknown", value: strands["."] || 0, className: "trm-bar-missing" },
    ]);
  }

  function renderTable() {
    const body = $("trm-call-table-body");
    body.replaceChildren();
    const pages = Math.max(1, Math.ceil(state.filtered.length / PAGE_SIZE));
    state.page = Math.min(state.page, pages - 1);
    const start = state.page * PAGE_SIZE;
    state.filtered.slice(start, start + PAGE_SIZE).forEach((call) => {
      const row = document.createElement("tr");
      const sourceCell = document.createElement("td");
      const badge = document.createElement("span");
      badge.className = `trm-badge trm-badge-${call.source}`;
      badge.textContent = call.source;
      sourceCell.appendChild(badge);
      row.appendChild(sourceCell);
      appendTextCell(row, call.chrom);
      appendTextCell(row, formatInteger.format(call.anchor));
      appendTextCell(row, call.family);
      appendTextCell(row, call.event_id);
      appendTextCell(row, call.event_type);
      appendTextCell(row, call.strand);
      appendTextCell(row, call.tsd);
      appendTextCell(row, valueOrDash(call.identity, "%"));
      appendTextCell(row, valueOrDash(call.coverage, "%"));
      appendTextCell(row, valueOrDash(call.te_size));
      appendTextCell(row, valueOrDash(call.display_frequency, "%"));
      appendTextCell(row, valueOrDash(call.sv_size));
      body.appendChild(row);
    });
    $("trm-page-status").textContent = state.filtered.length ? `Page ${state.page + 1} / ${pages}` : "No result";
    $("trm-page-previous").disabled = state.page === 0;
    $("trm-page-next").disabled = state.page >= pages - 1;
    $("trm-filter-count").textContent = `${formatInteger.format(state.filtered.length)} / ${formatInteger.format(calls.length)} calls`;
  }

  function changePage(delta) {
    state.page += delta;
    renderTable();
    document.querySelector(".trm-table-wrap").scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function downloadCalls() {
    const columns = ["source", "chrom", "start", "end", "anchor", "family", "event_id", "event_type", "strand", "tsd", "identity", "coverage", "te_size", "display_frequency", "sv_size", "tremolo_id"];
    const lines = [columns.join("\t")];
    state.filtered.forEach((call) => {
      lines.push(columns.map((column) => call[column] === null || call[column] === undefined ? "." : String(call[column]).replace(/[\t\r\n]/g, " ")).join("\t"));
    });
    const blob = new Blob([`${lines.join("\n")}\n`], { type: "text/tab-separated-values;charset=utf-8" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "TrEMOLO.filtered-calls.tsv";
    link.click();
    setTimeout(() => URL.revokeObjectURL(link.href), 1000);
  }

  function renderProximity() {
    const groups = data.proximity_groups || [];
    const note = $("trm-proximity-note");
    note.textContent = `${groups.length} candidate group${groups.length === 1 ? "" : "s"} contain at least two calls within a complete ${data.report.locus_window_bp} bp anchor span. These groups are not merged in TE_INFOS.bed.`;
    const body = $("trm-proximity-table-body");
    body.replaceChildren();
    groups.slice(0, 100).forEach((group) => {
      const row = document.createElement("tr");
      appendTextCell(row, group.candidate_group);
      appendTextCell(row, group.chrom);
      appendTextCell(row, `${formatInteger.format(group.start)}–${formatInteger.format(group.end)} (${group.span} bp)`);
      appendTextCell(row, group.calls);
      appendTextCell(row, group.families.join(", "));
      appendTextCell(row, group.sources.join(", "));
      appendTextCell(row, group.event_ids.join(", "));
      body.appendChild(row);
    });
    if (!groups.length) {
      const row = document.createElement("tr");
      const cell = document.createElement("td");
      cell.colSpan = 7;
      cell.textContent = "No nearby-call candidate in this run.";
      row.appendChild(cell);
      body.appendChild(row);
    }
  }

  function keyValues(values, preferredKeys = []) {
    const list = document.createElement("dl");
    list.className = "trm-key-values";
    const keys = preferredKeys.length ? preferredKeys.filter((key) => key in values) : Object.keys(values);
    keys.forEach((key) => {
      const term = document.createElement("dt");
      term.textContent = key;
      const description = document.createElement("dd");
      description.textContent = values[key];
      list.append(term, description);
    });
    return list;
  }

  function renderContext() {
    const mapping = $("trm-mapping-stats");
    const mappingKeys = ["raw total sequences", "reads mapped", "reads unmapped", "non-primary alignments", "bases mapped (cigar)", "average length", "error rate"];
    if (Object.keys(data.mapping_stats).length) mapping.replaceChildren(keyValues(data.mapping_stats, mappingKeys));
    else emptyChart(mapping, "Mapping statistics are not available for this configuration");

    const sv = $("trm-sv-summary");
    if (Object.keys(data.sv_counts).length) sv.replaceChildren(keyValues(data.sv_counts));
    else emptyChart(sv, "OUTSIDER structural-variant calls are not available");

    const provenance = $("trm-provenance");
    provenance.replaceChildren();
    const metadata = document.createElement("div");
    metadata.className = "trm-input-card";
    const heading = document.createElement("strong");
    heading.textContent = "Report contract";
    metadata.appendChild(heading);
    const contract = {
      "schema version": data.schema_version,
      "TE_INFOS SHA-256": data.report.te_infos_sha256,
      "work directory": data.report.work_directory || "—",
      "enabled sources": (data.report.enabled_sources || []).join(", ") || "—",
      "candidate proximity window": `${data.report.locus_window_bp} bp`,
    };
    if (data.report.author) contract.author = data.report.author;
    metadata.appendChild(keyValues(contract));
    provenance.appendChild(metadata);
    const inputs = (data.input_manifest && data.input_manifest.inputs) || {};
    Object.entries(inputs).forEach(([name, item]) => {
      const card = document.createElement("div");
      card.className = "trm-input-card";
      const title = document.createElement("strong");
      title.textContent = name;
      const path = document.createElement("code");
      path.textContent = item.path || "—";
      const checksum = document.createElement("div");
      checksum.className = "trm-caption";
      checksum.style.textAlign = "left";
      checksum.textContent = item.sha256 ? `SHA-256 ${item.sha256}` : "checksum unavailable";
      card.append(title, path, checksum);
      provenance.appendChild(card);
    });
  }

  function renderSelection() {
    renderCards();
    renderFamilyChart();
    renderFrequencyChart();
    renderGenomeChart();
    renderFrequencyPositionChart();
    renderEvidence();
    renderTable();
  }

  document.title = data.report.title || document.title;
  $("trm-report-title").textContent = data.report.title || "Transposable-element landscape";
  const sourceText = (data.report.enabled_sources || []).join(" + ") || "No enabled evidence source";
  $("trm-report-subtitle").textContent = `${sourceText} · ${data.summary.calls} legacy-compatible calls · data schema ${data.schema_version}`;
  initializeFilters();
  renderSelection();
  renderProximity();
  renderContext();
})();
