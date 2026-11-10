(() => {
  "use strict";

  const node = document.getElementById("tremolo-report-data");
  if (!node) return;
  const data = JSON.parse(node.textContent);
  const calls = data.calls || [];
  const resident = data.resident_te || { available: false, summary: {}, copies: [], thresholds: {} };
  const residentCopies = resident.copies || [];
  const ambiguous = data.ambiguous_calls || { available: false, summary: {}, candidates: [] };
  const ambiguousCandidates = ambiguous.candidates || [];
  const candidateEvidence = ambiguous.evidence || [];
  const candidateEvidenceByGroup = new Map();
  candidateEvidence.forEach((item) => {
    if (!candidateEvidenceByGroup.has(item.candidate_group_id)) candidateEvidenceByGroup.set(item.candidate_group_id, []);
    candidateEvidenceByGroup.get(item.candidate_group_id).push(item);
  });
  const PAGE_SIZE = 50;
  const ambiguousState = { page: 0 };
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
  const residentState = {
    filtered: [],
    page: 0,
    fullCoverage: Number(resident.thresholds.full_length_coverage || 80),
    fullIdentity: Number(resident.thresholds.high_confidence_pident || 80),
    partialCoverage: Number(resident.thresholds.partial_coverage || 20),
  };
  const $ = (id) => document.getElementById(id);
  const formatInteger = new Intl.NumberFormat("en-US", { maximumFractionDigits: 0 });
  const formatDecimal = new Intl.NumberFormat("en-US", { maximumFractionDigits: 4 });
  const formatPercent = new Intl.NumberFormat("en-US", { maximumFractionDigits: 1 });

  function valueOrDash(value, suffix = "") {
    return value === null || value === undefined ? "—" : `${formatDecimal.format(value)}${suffix}`;
  }

  function appendTextCell(row, value) {
    const cell = document.createElement("td");
    cell.textContent = value === null || value === undefined || value === "" ? "—" : String(value);
    cell.title = cell.textContent;
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

  function renderPageNumbers(prefix, current, pages, navigate) {
    const status = $(`${prefix}-status`);
    let numbers = $(`${prefix}-numbers`);
    if (!numbers) {
      numbers = document.createElement("span");
      numbers.id = `${prefix}-numbers`;
      numbers.className = "trm-page-numbers";
      numbers.setAttribute("aria-label", "Page numbers");
      status.after(numbers);
    }
    numbers.replaceChildren();
    const start = Math.max(0, Math.min(current - 3, pages - 7));
    const visible = new Set([0, pages - 1]);
    for (let page = start; page < Math.min(pages, start + 7); page += 1) visible.add(page);
    let previous = -1;
    [...visible].sort((a, b) => a - b).forEach((page) => {
      if (page > previous + 1) {
        const gap = document.createElement("span");
        gap.textContent = "…";
        numbers.appendChild(gap);
      }
      const button = document.createElement("button");
      button.type = "button";
      button.className = "trm-button trm-button-secondary";
      button.textContent = page + 1;
      button.setAttribute("aria-label", `Go to page ${page + 1}`);
      if (page === current) button.setAttribute("aria-current", "page");
      button.addEventListener("click", () => navigate(page));
      numbers.appendChild(button);
      previous = page;
    });
    if (pages > 9) {
      const label = document.createElement("label");
      label.textContent = "Go to page ";
      const input = document.createElement("input");
      input.type = "number";
      input.min = 1;
      input.max = pages;
      input.value = current + 1;
      input.setAttribute("aria-label", "Page number");
      input.addEventListener("change", () => {
        const page = Number(input.value);
        if (Number.isInteger(page) && page >= 1 && page <= pages) navigate(page - 1);
        else input.value = current + 1;
      });
      label.appendChild(input);
      numbers.appendChild(label);
    }
  }

  function initializeFilters() {
    selectOptions($("trm-filter-source"), unique("source"), "All sources");
    selectOptions($("trm-filter-type"), unique("event_type"), "All event types");
    selectOptions($("trm-filter-chrom"), unique("chrom"), "All chromosomes");
    selectOptions($("trm-filter-family"), unique("family"), "All TE families");
    ["source", "type", "chrom", "family", "tsd", "frequency"].forEach((key) => {
      $(`trm-filter-${key}`).addEventListener("change", applyFilters);
    });
    ["trm-depth-scale", "trm-depth-guide", "trm-depth-frequency-guide"].forEach((id) => {
      $(id).addEventListener("input", renderFrequencyDepthChart);
    });
    $("trm-heatmap-window").addEventListener("change", renderGenomeHeatmap);
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
        depthDetails(call),
      ].join("\n");
      point.appendChild(title);
      svg.appendChild(point);
    });
    target.replaceChildren(svg);
  }

  function depthDetails(call) {
    if (call.frequency_depth === null || call.frequency_depth === undefined) {
      return `depth unavailable (${call.frequency_depth_status || "unavailable"})`;
    }
    return `depth ${valueOrDash(call.frequency_depth)} reads; TE-present support ${valueOrDash(call.frequency_support_reads)}; ${call.frequency_depth_method}; ${call.frequency_depth_table}; ${call.frequency_depth_status}`;
  }

  function renderFrequencyDepthChart() {
    const target = $("trm-frequency-depth-chart");
    const eligible = state.filtered.filter((call) => call.frequency_depth_status === "matched"
      && Number.isFinite(call.frequency_depth) && call.frequency_depth > 0
      && Number.isFinite(call.display_frequency) && call.display_frequency >= 0 && call.display_frequency <= 100);
    const withFrequency = state.filtered.filter((call) => Number.isFinite(call.display_frequency));
    const omitted = withFrequency.length - eligible.length;
    const numberOrDefault = (id, fallback) => {
      const raw = $(id).value;
      return raw !== "" && Number.isFinite(Number(raw)) ? Number(raw) : fallback;
    };
    const depthGuide = Math.max(1, numberOrDefault("trm-depth-guide", 10));
    const frequencyGuide = Math.min(100, Math.max(0, numberOrDefault("trm-depth-frequency-guide", 50)));
    const quadrants = [0, 0, 0, 0];
    eligible.forEach((call) => {
      const highFrequency = call.display_frequency >= frequencyGuide;
      const highDepth = call.frequency_depth >= depthGuide;
      quadrants[(highFrequency ? 0 : 2) + (highDepth ? 0 : 1)] += 1;
    });
    $("trm-depth-summary").textContent = `${eligible.length} plotted · ${omitted} frequency estimates omitted (missing/zero/invalid depth or inconsistent evidence) · ${state.filtered.length - withFrequency.length} without frequency. Guides: frequency ≥ ${frequencyGuide}%, depth ≥ ${depthGuide} reads. High frequency / high depth: ${quadrants[0]}; high frequency / low depth: ${quadrants[1]}; low frequency / high depth: ${quadrants[2]}; low frequency / low depth: ${quadrants[3]}.`;
    if (!eligible.length) return emptyChart(target, "No positive, matching depth for the frequency estimates in this selection");
    const width = 1100, height = 420;
    const margin = { top: 22, right: 25, bottom: 54, left: 72 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;
    const logarithmic = $("trm-depth-scale").value === "log";
    let observedMaximum = depthGuide;
    eligible.forEach((call) => { observedMaximum = Math.max(observedMaximum, call.frequency_depth); });
    const minimum = logarithmic ? 1 : 0;
    const maximum = Math.max(2, observedMaximum * 1.05);
    const transform = logarithmic ? Math.log10 : (value) => value;
    const x = (value) => margin.left + plotWidth * (transform(value) - transform(minimum)) / (transform(maximum) - transform(minimum));
    const y = (value) => margin.top + plotHeight * (1 - value / 100);
    const svg = svgElement("svg", { viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": "TE frequency versus local depth" });
    [0, 25, 50, 75, 100].forEach((value) => {
      svg.appendChild(svgElement("line", { x1: margin.left, x2: width - margin.right, y1: y(value), y2: y(value), class: "trm-scatter-grid" }));
      const label = svgElement("text", { x: margin.left - 12, y: y(value) + 4, class: "trm-scatter-axis-label", "text-anchor": "end" });
      label.textContent = `${value}%`;
      svg.appendChild(label);
    });
    for (let i = 0; i <= 5; i += 1) {
      const value = logarithmic ? 10 ** (Math.log10(maximum) * i / 5) : maximum * i / 5;
      const tick = x(value);
      svg.appendChild(svgElement("line", { x1: tick, x2: tick, y1: margin.top, y2: height - margin.bottom, class: "trm-scatter-grid trm-scatter-grid-vertical" }));
      const label = svgElement("text", { x: tick, y: height - margin.bottom + 22, class: "trm-scatter-axis-label", "text-anchor": "middle" });
      label.textContent = formatDecimal.format(value);
      svg.appendChild(label);
    }
    svg.appendChild(svgElement("line", { x1: x(depthGuide), x2: x(depthGuide), y1: margin.top, y2: height - margin.bottom, class: "trm-depth-guide" }));
    svg.appendChild(svgElement("line", { x1: margin.left, x2: width - margin.right, y1: y(frequencyGuide), y2: y(frequencyGuide), class: "trm-depth-guide" }));
    const xTitle = svgElement("text", { x: margin.left + plotWidth / 2, y: height - 8, class: "trm-scatter-axis-title", "text-anchor": "middle" });
    xTitle.textContent = `Frequency depth (reads; ${logarithmic ? "logarithmic" : "linear"} scale)`;
    const yTitle = svgElement("text", { x: 17, y: margin.top + plotHeight / 2, class: "trm-scatter-axis-title", "text-anchor": "middle", transform: `rotate(-90 17 ${margin.top + plotHeight / 2})` });
    yTitle.textContent = "TE frequency (%)";
    svg.append(xTitle, yTitle);
    eligible.forEach((call) => {
      const point = svgElement("circle", { cx: x(call.frequency_depth), cy: y(call.display_frequency), r: 5,
        fill: call.source === "INSIDER" ? "var(--trm-insider)" : "var(--trm-outsider)", class: "trm-frequency-depth-point trm-frequency-position-point", tabindex: "0" });
      const title = svgElement("title");
      title.textContent = `${call.family} · ${call.event_id}\n${call.source} · ${call.event_type} · ${call.chrom}:${call.start}\nfrequency ${valueOrDash(call.display_frequency, "%")}\n${depthDetails(call)}`;
      point.appendChild(title);
      svg.appendChild(point);
    });
    target.replaceChildren(svg);
  }

  function detectorLabel(call) {
    if (call.source === "OUTSIDER" && /^sniffles(?:[. _-]|\d)/i.test(call.event_id || "")) return "Sniffles";
    if (call.source === "OUTSIDER" && /^(?:TrEMOLO|HARD|SOFT)[. _-]/i.test(call.event_id || "")) return "TrEMOLO";
    if (call.source === "INSIDER" && /^Assemblytics/i.test(call.event_id || "")) return "Assemblytics";
    return "Other / unknown";
  }

  function detectorEvent(call) {
    const type = String(call.event_type || "").toUpperCase();
    if (type.startsWith("SOFT")) return "SOFT";
    if (type.startsWith("HARD")) return "HARD";
    if (type.startsWith("INS") || type.includes("EXPANSION")) return "Insertions";
    if (type.startsWith("DEL") || type.includes("CONTRACTION")) return "Deletions";
    return "Other";
  }

  function detectorCounts(rows) {
    const counts = new Map();
    rows.forEach((call) => {
      const label = detectorLabel(call);
      if (!counts.has(label)) counts.set(label, { Insertions: 0, Deletions: 0, SOFT: 0, HARD: 0, Other: 0, Total: 0 });
      counts.get(label)[detectorEvent(call)] += 1;
      counts.get(label).Total += 1;
    });
    return counts;
  }

  function renderDetectorSummary() {
    const baseline = detectorCounts(calls), selected = detectorCounts(state.filtered);
    const body = $("trm-detector-table-body");
    body.replaceChildren();
    ["Sniffles", "TrEMOLO", "Assemblytics", "Other / unknown"].forEach((detector) => {
      const all = baseline.get(detector);
      if (!all) return;
      const row = document.createElement("tr");
      appendTextCell(row, detector);
      ["Insertions", "Deletions", "SOFT", "HARD", "Other", "Total"].forEach((kind) => {
        appendTextCell(row, `${formatInteger.format((selected.get(detector) || {})[kind] || 0)} / ${formatInteger.format(all[kind])}`);
      });
      body.appendChild(row);
    });
    $("trm-detector-note").textContent = `${formatInteger.format(state.filtered.length)} selected / ${formatInteger.format(calls.length)} reported TE calls. Each cell shows selected / all; final call rows, not unique loci. SOFT/HARD channels are separate from insertion/deletion columns.`;
    const raw = $("trm-raw-sv-summary");
    raw.replaceChildren();
    if (!data.sv_counts_available) {
      raw.textContent = "Raw Sniffles VCF counts unavailable in this report.";
    } else {
      const counts = data.sv_counts || {};
      const total = Object.values(counts).reduce((sum, count) => sum + count, 0);
      raw.append(summaryCard("Raw Sniffles VCF", formatInteger.format(total), "all SV records before TE annotation; independent of call filters"),
        summaryCard("Raw INS / DEL", `${formatInteger.format(counts.INS || 0)} / ${formatInteger.format(counts.DEL || 0)}`, "Do not add raw SV counts to reported TE calls"));
    }
  }

  const MAGMA = ["#000004", "#140e36", "#3b0f70", "#641a80", "#8c2981", "#b73779", "#de4968", "#f7705c", "#fe9f6d", "#fecf92", "#fcfdbf"];
  function magmaColor(value, maximum) {
    if (!value) return "var(--trm-heatmap-empty)";
    const position = Math.min(1, value / Math.max(1, maximum)) * (MAGMA.length - 1);
    const low = Math.floor(position), high = Math.min(MAGMA.length - 1, low + 1), fraction = position - low;
    const rgb = [1, 3, 5].map((offset) => Math.round(parseInt(MAGMA[low].slice(offset, offset + 2), 16) * (1 - fraction) + parseInt(MAGMA[high].slice(offset, offset + 2), 16) * fraction));
    return `rgb(${rgb.join(",")})`;
  }

  function buildHeatmap(requestedSize) {
    const byChrom = new Map();
    calls.forEach((call) => {
      if (!byChrom.has(call.chrom)) byChrom.set(call.chrom, []);
      byChrom.get(call.chrom).push(call);
    });
    const lengths = new Map((data.chromosome_summary || []).map((item) => [item.chrom, item.length]));
    const spans = [];
    byChrom.forEach((rows, chrom) => {
      let observed = 1;
      rows.forEach((call) => { observed = Math.max(observed, call.anchor + 1); });
      spans.push({ chrom, length: lengths.get(chrom) || observed, known: Boolean(lengths.get(chrom)), rows });
    });
    spans.sort((a, b) => a.chrom.localeCompare(b.chrom, undefined, { numeric: true }));
    let size = Math.max(1, requestedSize);
    // Bound DOM size on large assemblies; disclose the actual window size in UI.
    while (spans.reduce((sum, item) => sum + Math.ceil(item.length / size), 0) > 20000 && spans.length < 20000) size *= 2;
    const selected = new Set(state.filtered);
    let maximum = 1;
    const chromosomes = spans.map((item) => {
      const bins = Array.from({ length: Math.ceil(item.length / size) }, (_, i) => ({
        chrom: item.chrom, start: i * size, end: Math.min(item.length, (i + 1) * size),
        baseline: { INSIDER: 0, OUTSIDER: 0, UNKNOWN: 0 }, selected: { INSIDER: 0, OUTSIDER: 0, UNKNOWN: 0 }, all: 0, kept: 0,
      }));
      item.rows.forEach((call) => {
        // A terminal insertion at chromosome length belongs to the final window.
        if (!Number.isFinite(call.anchor) || call.anchor < 0 || call.anchor > item.length) return;
        const bin = bins[Math.min(bins.length - 1, Math.floor(call.anchor / size))];
        const source = call.source === "INSIDER" || call.source === "OUTSIDER" ? call.source : "UNKNOWN";
        bin.baseline[source] += 1; bin.all += 1;
        if (selected.has(call)) { bin.selected[source] += 1; bin.kept += 1; }
        maximum = Math.max(maximum, bin.all);
      });
      return { ...item, bins };
    });
    return { chromosomes, size, maximum };
  }

  function showWindowTooltip(event, bin, point) {
    const tooltip = $("trm-heatmap-tooltip");
    tooltip.replaceChildren();
    const heading = document.createElement("strong");
    heading.textContent = `${bin.chrom}:${formatInteger.format(bin.start)}–${formatInteger.format(bin.end)} bp`;
    tooltip.appendChild(heading);
    const total = document.createElement("p");
    total.textContent = `${bin.kept} selected / ${bin.all} reported calls · ${bin.all - bin.kept} excluded by filters`;
    tooltip.appendChild(total);
    ["INSIDER", "OUTSIDER", "UNKNOWN"].forEach((source) => {
      if (source === "UNKNOWN" && !bin.baseline.UNKNOWN) return;
      const line = document.createElement("div");
      line.textContent = `${source}: ${bin.selected[source]} selected / ${bin.baseline[source]} before filtering`;
      tooltip.appendChild(line);
    });
    tooltip.hidden = false;
    const bounds = point.getBoundingClientRect();
    const x = Number.isFinite(event.clientX) ? event.clientX : bounds.left;
    const y = Number.isFinite(event.clientY) ? event.clientY : bounds.bottom;
    tooltip.style.left = `${Math.max(8, Math.min(window.innerWidth - tooltip.offsetWidth - 12, x + 14))}px`;
    tooltip.style.top = `${Math.max(8, Math.min(window.innerHeight - tooltip.offsetHeight - 12, y + 14))}px`;
  }

  function renderGenomeHeatmap() {
    const target = $("trm-genome-heatmap");
    $("trm-heatmap-tooltip").hidden = true;
    const requested = Number($("trm-heatmap-window").value) || 1000000;
    const model = buildHeatmap(requested);
    if (!model.chromosomes.length) {
      $("trm-heatmap-impact").textContent = "No reported calls in this analysis.";
      return emptyChart(target, "No call windows available");
    }
    let affected = 0, windows = 0, maximumLength = 1;
    model.chromosomes.forEach((item) => {
      maximumLength = Math.max(maximumLength, item.length);
      item.bins.forEach((bin) => { windows += 1; if (bin.kept !== bin.all) affected += 1; });
    });
    $("trm-heatmap-impact").textContent = `${affected} / ${windows} windows filter affected · ${windows - affected} unchanged · actual window: ${formatInteger.format(model.size)} bp${model.size !== requested ? " (increased to bound display size)" : ""}. Chromosomes with reported calls remain visible, including filtered-out regions.`;
    $("trm-heatmap-impact").setAttribute("data-filtered", affected > 0 ? "true" : "false");
    $("trm-heatmap-maximum").textContent = `${model.maximum} calls / window`;
    const width = 1100, left = 190, right = 90, rowHeight = 38, top = 32;
    const height = top + model.chromosomes.length * rowHeight + 10;
    const x = (value) => left + (width - left - right) * value / maximumLength;
    const svg = svgElement("svg", { viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": "TE call counts per chromosome window with filter impact" });
    for (let i = 0; i <= 4; i += 1) {
      const tick = svgElement("text", { x: x(maximumLength * i / 4), y: 16, class: "trm-scatter-axis-label", "text-anchor": "middle" });
      tick.textContent = `${compactPosition(maximumLength * i / 4)} bp`; svg.appendChild(tick);
    }
    model.chromosomes.forEach((item, index) => {
      const y = top + rowHeight * index;
      const label = svgElement("text", { x: 5, y: y + 16, class: "trm-genome-label" });
      label.textContent = item.chrom.length > 25 ? `${item.chrom.slice(0, 22)}…` : item.chrom;
      const labelTitle = svgElement("title"); labelTitle.textContent = `${item.chrom}${item.known ? "" : " — observed span; chromosome length unavailable"}`;
      label.appendChild(labelTitle); svg.appendChild(label);
      item.bins.forEach((bin) => {
        const changed = bin.kept !== bin.all;
        const point = svgElement("rect", { x: x(bin.start), y, width: Math.max(0.25, x(bin.end) - x(bin.start)), height: 23, rx: 2,
          fill: magmaColor(bin.kept, model.maximum), class: `trm-heatmap-cell${changed ? " trm-heatmap-affected" : ""}`, tabindex: "0",
          "aria-label": `${bin.chrom}:${bin.start}-${bin.end}: ${bin.kept} selected of ${bin.all}; ${changed ? "filter affected" : "unchanged"}`,
          "aria-describedby": "trm-heatmap-tooltip" });
        point.addEventListener("pointerenter", (event) => showWindowTooltip(event, bin, point));
        point.addEventListener("pointermove", (event) => showWindowTooltip(event, bin, point));
        point.addEventListener("focus", (event) => showWindowTooltip(event, bin, point));
        point.addEventListener("pointerleave", () => { $("trm-heatmap-tooltip").hidden = true; });
        point.addEventListener("blur", () => { $("trm-heatmap-tooltip").hidden = true; });
        point.addEventListener("keydown", (event) => { if (event.key === "Escape") $("trm-heatmap-tooltip").hidden = true; });
        svg.appendChild(point);
      });
    });
    target.replaceChildren(svg);
  }

  function filteredStructureRelations(identity, coverage) {
    const eligible = new Set(residentCopies.filter((copy) => copy.identity >= identity && copy.coverage >= coverage)
      .map((copy) => JSON.stringify([copy.target, copy.copy_id])));
    return (resident.relations || []).filter((relation) => eligible.has(JSON.stringify([relation.target, relation.source_copy_id]))
      && eligible.has(JSON.stringify([relation.target, relation.target_copy_id])));
  }

  function renderStructurePreview() {
    const identity = Number($("trm-structure-identity").value) || 0;
    const coverage = Number($("trm-structure-coverage").value) || 0;
    $("trm-structure-identity-value").textContent = `${identity}%`;
    $("trm-structure-coverage-value").textContent = `${coverage}%`;
    const target = $("trm-structure-preview");
    target.replaceChildren();
    if (!Array.isArray(resident.relations)) return emptyChart(target, "Rebuild the report to enable relation-level threshold preview.");
    const rows = filteredStructureRelations(identity, coverage);
    const counts = new Map(); rows.forEach((row) => counts.set(row.relation, (counts.get(row.relation) || 0) + 1));
    $("trm-structure-preview-note").textContent = `${rows.length} / ${resident.relations.length} relation pairs retained. Both copies must pass identity and consensus coverage; no regrouping or re-annotation. Independent of table filters.`;
    ["same_family_overlap", "nested_candidate", "overlap_candidate"].forEach((type) => {
      const total = (resident.summary.relation_types || {})[type] || 0;
      target.appendChild(summaryCard(type, `${formatInteger.format(counts.get(type) || 0)} / ${formatInteger.format(total)}`, "preview / all reported relations"));
    });
  }

  function renderGenomeChart() {
    const target = $("trm-genome-chart");
    const selectedChroms = new Set(state.filtered.map((call) => call.chrom));
    const chromosomes = data.chromosome_summary.filter((item) => selectedChroms.has(item.chrom));
    if (!chromosomes.length) return emptyChart(target, "No genomic position in the current selection");
    const width = 1100, left = 190, right = 150, rowHeight = 42;
    const height = 52 + chromosomes.length * rowHeight;
    // Keep one absolute scale across chromosomes and filter changes. Unknown
    // lengths use the observed span, explicitly labelled rather than inferred as fact.
    const spans = new Map((data.chromosome_summary || []).map((item) => [item.chrom, item.length || 1]));
    calls.forEach((call) => spans.set(call.chrom, Math.max(spans.get(call.chrom) || 1, call.anchor || 0)));
    let maximum = 1;
    spans.forEach((span) => { maximum = Math.max(maximum, span); });
    const x = (bp) => left + (width - left - right) * bp / maximum;
    const svg = svgElement("svg", { viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": "TE calls along chromosomes on a shared base-pair scale" });
    for (let i = 0; i <= 4; i += 1) {
      const value = maximum * i / 4;
      svg.appendChild(svgElement("line", { x1: x(value), x2: x(value), y1: 26, y2: height - 12, class: "trm-scatter-grid trm-scatter-grid-vertical" }));
      const tick = svgElement("text", { x: x(value), y: 14, "text-anchor": "middle", class: "trm-scatter-axis-label" });
      tick.textContent = `${compactPosition(value)} bp`;
      svg.appendChild(tick);
    }
    chromosomes.forEach((item, index) => {
      const y = 45 + index * rowHeight;
      const label = svgElement("text", { x: 5, y: y + 4, class: "trm-genome-label" });
      label.textContent = item.chrom.length > 25 ? `${item.chrom.slice(0, 22)}…` : item.chrom;
      const labelTitle = svgElement("title");
      labelTitle.textContent = item.chrom;
      label.appendChild(labelTitle);
      svg.appendChild(label);
      const span = item.length || spans.get(item.chrom);
      const track = svgElement("line", { x1: left, y1: y, x2: x(span), y2: y, class: `trm-genome-line${item.length ? "" : " trm-genome-length-unknown"}`, tabindex: "0" });
      const lengthText = item.length ? `${formatInteger.format(item.length)} bp` : "Length unavailable; showing observed span";
      const trackTitle = svgElement("title");
      trackTitle.textContent = `${item.chrom}: ${lengthText}`;
      track.appendChild(trackTitle);
      svg.appendChild(track);
      const size = svgElement("text", { x: width - right + 14, y: y + 4, class: "trm-genome-label" });
      size.textContent = item.length ? `${formatDecimal.format(item.length / 1e6)} Mb` : "Length unknown";
      const sizeTitle = svgElement("title");
      sizeTitle.textContent = lengthText;
      size.appendChild(sizeTitle);
      svg.appendChild(size);
      state.filtered.filter((call) => call.chrom === item.chrom).forEach((call) => {
        const point = svgElement("circle", { cx: x(Math.max(0, call.anchor)), cy: y, r: 4.2,
          fill: call.source === "INSIDER" ? "#6f5bd3" : call.source === "OUTSIDER" ? "#ed8a3b" : "#a8b4b9", class: "trm-call-point", tabindex: "0" });
        const title = svgElement("title");
        title.textContent = `${call.family} · ${call.event_id} · ${call.chrom}:${formatInteger.format(call.anchor)}\nChromosome: ${lengthText}`;
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
      label.title = item.label;
      label.tabIndex = 0;
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

  function genomicAlignmentSummary(call) {
    if (call.source === "INSIDER") return "Not applicable (assembly call)";
    const counts = Object.entries(call.genomic_alignment_counts || {});
    return counts.length
      ? counts.map(([status, count]) => `${status.replaceAll("_", " ")}: ${count}`).join("; ")
      : "Unknown (no source metadata)";
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
      appendTextCell(row, genomicAlignmentSummary(call));
      appendTextCell(row, valueOrDash(call.frequency_depth)).title = depthDetails(call);
      body.appendChild(row);
    });
    $("trm-page-status").textContent = state.filtered.length ? `Page ${state.page + 1} / ${pages}` : "No result";
    $("trm-page-previous").disabled = state.page === 0;
    $("trm-page-next").disabled = state.page >= pages - 1;
    renderPageNumbers("trm-page", state.page, pages, (page) => {
      state.page = page;
      renderTable();
    });
    $("trm-filter-count").textContent = `${formatInteger.format(state.filtered.length)} / ${formatInteger.format(calls.length)} calls`;
  }

  function changePage(delta) {
    state.page += delta;
    renderTable();
    document.querySelector(".trm-table-wrap").scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function downloadCalls() {
    const columns = ["source", "chrom", "start", "end", "anchor", "family", "event_id", "event_type", "strand", "tsd", "identity", "coverage", "te_size", "display_frequency", "sv_size", "tremolo_id", "source_genomic_alignments", "frequency_depth", "frequency_support_reads", "frequency_depth_method", "frequency_depth_table", "frequency_depth_status"];
    const lines = [columns.join("\t")];
    state.filtered.forEach((call) => {
      lines.push(columns.map((column) => column === "source_genomic_alignments" ? genomicAlignmentSummary(call) : call[column] === null || call[column] === undefined ? "." : String(call[column]).replace(/[\t\r\n]/g, " ")).join("\t"));
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

  function renderAmbiguousCalls() {
    const note = $("trm-ambiguous-call-note");
    const body = $("trm-ambiguous-call-table-body");
    const download = $("trm-download-ambiguous-calls");
    const evidenceDownload = $("trm-download-candidate-evidence");
    body.replaceChildren();
    if (!ambiguous.available) {
      note.textContent = "The normalized ambiguous-call table is not available in this report.";
      download.disabled = true;
      evidenceDownload.disabled = true;
      $("trm-ambiguous-page-previous").disabled = true;
      $("trm-ambiguous-page-next").disabled = true;
      $("trm-ambiguous-page-status").textContent = "No result";
      return;
    }
    const summary = ambiguous.summary || {};
    note.textContent = ambiguousCandidates.length
      ? `${formatInteger.format(summary.ambiguous_calls || 0)} reported variable calls retain ${formatInteger.format(summary.candidate_rows || 0)} TE candidates. Unique calls are intentionally omitted.`
      : "No reported variable TE call has more than one distinct candidate label. Unique calls are intentionally omitted.";
    const pages = Math.max(1, Math.ceil(ambiguousCandidates.length / PAGE_SIZE));
    ambiguousState.page = Math.max(0, Math.min(ambiguousState.page, pages - 1));
    const start = ambiguousState.page * PAGE_SIZE;
    ambiguousCandidates.slice(start, start + PAGE_SIZE).forEach((candidate) => {
      const row = document.createElement("tr");
      appendTextCell(row, candidate.source);
      appendTextCell(row, candidate.event_id);
      appendTextCell(row, `${candidate.chrom}:${formatInteger.format(candidate.start)}–${formatInteger.format(candidate.end)}`);
      appendTextCell(row, candidate.reported_te);
      appendTextCell(row, candidate.candidate_te);
      const assignment = document.createElement("td");
      const badge = document.createElement("span");
      badge.className = `trm-badge trm-candidate-${candidate.assignment}`;
      badge.textContent = candidate.assignment === "reported_primary" ? "reported primary" : "alternative";
      assignment.appendChild(badge);
      row.appendChild(assignment);
      appendTextCell(row, `${candidate.candidate_rank} / ${candidate.candidate_count}`);
      appendTextCell(row, candidate.evidence_count);
      appendTextCell(row, valueOrDash(100 * candidate.evidence_fraction, "%"));
      appendTextCell(row, candidate.evidence_channels.join(", "));
      const relationshipLabels = {
        same_sequence_annotation_ambiguity: "Competing annotations on one sequence",
        same_sequence_multiple_regions: "Different regions on one sequence",
        different_sequence_evidence: "Different insertion sequences",
        mixed_sequence_evidence: "Mixed sequence evidence",
        unresolved_family_candidates: "Unresolved sequence relationship",
      };
      appendTextCell(row, relationshipLabels[candidate.ambiguity_type] || candidate.ambiguity_type);
      appendTextCell(row, [candidate.query_count, candidate.sequence_count, candidate.read_count, candidate.vcf_query_count].map((value) => value === null || value === undefined ? "?" : value).join(" / "));
      appendTextCell(row, candidate.provenance_status || "unavailable");
      const evidenceCell = document.createElement("td");
      const details = document.createElement("details");
      const summaryNode = document.createElement("summary");
      summaryNode.textContent = "Inspect raw HSPs";
      details.appendChild(summaryNode);
      details.addEventListener("toggle", () => {
        if (!details.open || details.dataset.rendered) return;
        details.dataset.rendered = "yes";
        const records = (candidateEvidenceByGroup.get(candidate.candidate_group_id) || []).filter((item) => item.candidate_te === candidate.candidate_te);
        if (!records.length) {
          const note = document.createElement("p");
          note.textContent = "Sequence provenance unavailable for this candidate.";
          details.appendChild(note);
        }
        records.forEach((item) => {
          const line = document.createElement("p");
          line.textContent = `${item.query_id} · read ${item.read_id} · ${item.sequence_origin} · source genomic alignment ${item.genomic_alignment_status || (item.sequence_origin === "vcf_reported" ? "not_applicable" : "unknown")} (SAM flag ${item.genomic_alignment_flag ?? "."}) · query [${item.query_start}, ${item.query_end}) · TE [${item.te_start}, ${item.te_end}) ${item.strand} · identity ${item.identity}% · insertion ${item.insertion_length} bp · SHA256 ${item.insertion_sha256}`;
          details.appendChild(line);
        });
      });
      evidenceCell.appendChild(details);
      row.appendChild(evidenceCell);
      body.appendChild(row);
    });
    download.disabled = !ambiguousCandidates.length;
    if (!download.dataset.bound) {
      download.addEventListener("click", downloadAmbiguousCalls);
      download.dataset.bound = "yes";
    }
    evidenceDownload.disabled = !candidateEvidence.length;
    if (!evidenceDownload.dataset.bound) {
      evidenceDownload.dataset.bound = "yes";
      evidenceDownload.addEventListener("click", () => {
      const columns = ["candidate_group_id", "event_id", "evidence_channel", "query_id", "read_id", "sequence_origin", "insertion_sha256", "insertion_length", "candidate_te", "query_start", "query_end", "te_start", "te_end", "strand", "identity", "aligned_bp", "evalue", "bitscore", "relation_types", "genomic_alignment_flag", "genomic_alignment_status"];
      const content = [columns.join("\t"), ...candidateEvidence.map((item) => columns.map((column) => candidateTsvCell(item[column])).join("\t"))].join("\n") + "\n";
      const link = document.createElement("a");
      link.href = URL.createObjectURL(new Blob([content], { type: "text/tab-separated-values;charset=utf-8" }));
      link.download = "TE_CALL_CANDIDATE_EVIDENCE.tsv";
      link.click();
      URL.revokeObjectURL(link.href);
      });
    }
    $("trm-ambiguous-page-status").textContent = ambiguousCandidates.length
      ? `Page ${ambiguousState.page + 1} / ${pages}` : "No result";
    $("trm-ambiguous-page-previous").disabled = ambiguousState.page === 0;
    $("trm-ambiguous-page-next").disabled = ambiguousState.page >= pages - 1;
    renderPageNumbers("trm-ambiguous-page", ambiguousState.page, pages, (page) => {
      ambiguousState.page = page;
      renderAmbiguousCalls();
    });
  }

  function candidateTsvCell(value) {
    const text = value === null || value === undefined ? "." : String(value).replace(/[\t\r\n]/g, " ");
    return text.includes('"') ? '"' + text.replaceAll('"', '""') + '"' : text;
  }

  function downloadAmbiguousCalls() {
    const columns = [
      "candidate_group_id", "source", "chrom", "start", "end", "event_id",
      "tremolo_id", "event_type", "reported_te", "candidate_te", "assignment",
      "candidate_rank", "candidate_count", "evidence_count", "evidence_fraction",
      "evidence_channels", "ambiguity_type",
      "query_count", "sequence_count", "read_count", "vcf_query_count", "provenance_status", "relation_types",
    ];
    const lines = [columns.join("\t")];
    ambiguousCandidates.forEach((candidate) => {
      lines.push(columns.map((column) => {
        const value = column === "evidence_channels"
          ? candidate.evidence_channels.join(";")
          : candidate[column];
        return candidateTsvCell(value);
      }).join("\t"));
    });
    const blob = new Blob([`${lines.join("\n")}\n`], { type: "text/tab-separated-values;charset=utf-8" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "TrEMOLO.TE_CALL_CANDIDATES.tsv";
    link.click();
    setTimeout(() => URL.revokeObjectURL(link.href), 1000);
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

  function residentTier(copy) {
    if (copy.coverage >= residentState.fullCoverage && copy.identity >= residentState.fullIdentity) return "full_length";
    if (copy.coverage >= residentState.partialCoverage) return "partial";
    return "degraded_relic";
  }

  function residentCandidates(copy) {
    if (!copy.candidates || !copy.candidates.length) return copy.te_candidates.join(", ");
    return copy.candidates.map((candidate) =>
      `${candidate.te} [${candidate.assignment}; cov ${formatDecimal.format(candidate.coverage)}%; id ${formatDecimal.format(candidate.identity)}%]`
    ).join(" · ");
  }

  function initializeResidentFilters() {
    const uniqueResident = (field) => [...new Set(residentCopies.map((copy) => copy[field]))]
      .sort((a, b) => String(a).localeCompare(String(b)));
    selectOptions($("trm-resident-filter-target"), uniqueResident("target"), "All targets");
    selectOptions($("trm-resident-filter-chrom"), uniqueResident("chrom"), "All chromosomes");
    selectOptions($("trm-resident-filter-status"), uniqueResident("status"), "All configured tiers");
    ["target", "chrom", "status", "ambiguity"].forEach((key) => {
      $(`trm-resident-filter-${key}`).addEventListener("change", applyResidentFilters);
    });
    $("trm-resident-filter-search").addEventListener("input", applyResidentFilters);
    $("trm-resident-reset-filters").addEventListener("click", () => {
      ["target", "chrom", "status"].forEach((key) => { $(`trm-resident-filter-${key}`).value = ""; });
      $("trm-resident-filter-ambiguity").value = "multi_te";
      $("trm-resident-filter-search").value = "";
      applyResidentFilters();
    });
    $("trm-resident-download").addEventListener("click", downloadResidentCopies);
    $("trm-resident-page-previous").addEventListener("click", () => changeResidentPage(-1));
    $("trm-resident-page-next").addEventListener("click", () => changeResidentPage(1));
  }

  function applyResidentFilters() {
    const target = $("trm-resident-filter-target").value;
    const chrom = $("trm-resident-filter-chrom").value;
    const status = $("trm-resident-filter-status").value;
    const ambiguity = $("trm-resident-filter-ambiguity").value;
    const search = $("trm-resident-filter-search").value.trim().toLocaleLowerCase();
    residentState.filtered = residentCopies.filter((copy) => {
      if (target && copy.target !== target) return false;
      if (chrom && copy.chrom !== chrom) return false;
      if (status && copy.status !== status) return false;
      if (ambiguity === "multi_te" && copy.te_candidates.length <= 1) return false;
      if (ambiguity === "ambiguous" && copy.candidate_count <= 1) return false;
      if (ambiguity === "unique" && copy.candidate_count > 1) return false;
      if (search) {
        const haystack = [copy.copy_id, copy.target, copy.chrom, copy.primary_te, ...copy.te_candidates]
          .join(" ").toLocaleLowerCase();
        if (!haystack.includes(search)) return false;
      }
      return true;
    });
    residentState.page = 0;
    renderResidentTable();
  }

  function renderResidentCards() {
    const summary = resident.summary || {};
    const tiers = summary.tiers || {};
    const target = $("trm-resident-cards");
    target.replaceChildren(
      summaryCard("Resident copies", formatInteger.format(summary.copies || 0), Object.entries(summary.targets || {}).map(([name, count]) => `${name} ${count}`).join(" · ")),
      summaryCard("Primary families", formatInteger.format(summary.families || 0), "distinct primary assignments"),
      summaryCard("Full length", formatInteger.format(tiers.full_length || 0), "configured classification"),
      summaryCard("Partial / relic", `${formatInteger.format(tiers.partial || 0)} / ${formatInteger.format(tiers.degraded_relic || 0)}`, "configured classification"),
      summaryCard("Copies with multiple candidate matches", formatInteger.format(summary.ambiguous_copies || 0), "Several candidate matches, possibly assigned to the same TE family"),
      summaryCard("Copies with multiple candidate families", formatInteger.format(summary.multi_te_copies || 0), "Subset of copies with multiple matches; alternative annotations, not confirmed multi-TE structures")
    );
  }

  function initializeResidentThresholds() {
    const controls = [
      ["full", "fullCoverage"], ["identity", "fullIdentity"], ["partial", "partialCoverage"],
    ];
    controls.forEach(([id, field]) => {
      const input = $(`trm-resident-${id}-threshold`);
      input.value = residentState[field];
      input.addEventListener("input", () => {
        residentState[field] = Number(input.value);
        renderResidentCalibration();
      });
    });
    $("trm-resident-reset-thresholds").addEventListener("click", () => {
      residentState.fullCoverage = Number(resident.thresholds.full_length_coverage || 80);
      residentState.fullIdentity = Number(resident.thresholds.high_confidence_pident || 80);
      residentState.partialCoverage = Number(resident.thresholds.partial_coverage || 20);
      controls.forEach(([id, field]) => { $(`trm-resident-${id}-threshold`).value = residentState[field]; });
      renderResidentCalibration();
    });
  }

  function renderResidentCalibration() {
    $("trm-resident-full-value").textContent = `${residentState.fullCoverage}%`;
    $("trm-resident-identity-value").textContent = `${residentState.fullIdentity}%`;
    $("trm-resident-partial-value").textContent = `${residentState.partialCoverage}%`;
    const counts = { full_length: 0, partial: 0, degraded_relic: 0 };
    residentCopies.forEach((copy) => { counts[residentTier(copy)] += 1; });
    const summary = $("trm-resident-calibration-summary");
    summary.replaceChildren();
    const tiers = [
      ["full_length", "Full-length", `Coverage ≥ ${residentState.fullCoverage}% and identity ≥ ${residentState.fullIdentity}%`],
      ["partial", "Partial", `Coverage ≥ ${residentState.partialCoverage}%; full-length criteria not met`],
      ["degraded_relic", "Degraded / relic", `Coverage < ${residentState.partialCoverage}%; full-length criteria not met`],
    ];
    tiers.forEach(([key, label, rule]) => {
      const percentage = residentCopies.length ? 100 * counts[key] / residentCopies.length : 0;
      const card = summaryCard(label, formatInteger.format(counts[key]), `${formatPercent.format(percentage)}% of reported copies · ${rule}`);
      card.className += ` trm-tier-card trm-tier-${key}`;
      summary.appendChild(card);
    });
    $("trm-resident-calibration-context").textContent = `Classification preview across ${formatInteger.format(residentCopies.length)} reported copies. Counts use the current sliders and all reported copies, independently of the resident table filters.`;
    $("trm-resident-discovery-floor").textContent = `Discovery filters: identity ≥ ${resident.thresholds.min_pident}% · aligned length ≥ ${resident.thresholds.min_aligned_bp} bp. These are fixed pipeline settings; moving the sliders cannot recover evidence excluded during discovery.`;

    const target = $("trm-resident-calibration-chart");
    if (!residentCopies.length) return emptyChart(target, "No resident copy to calibrate");
    const width = 1000;
    const height = 500;
    const margin = { top: 25, right: 25, bottom: 55, left: 65 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;
    const x = (value) => margin.left + plotWidth * Math.min(100, Math.max(0, value)) / 100;
    const y = (value) => margin.top + plotHeight * (1 - Math.min(100, Math.max(0, value)) / 100);
    const svg = svgElement("svg", { viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": "Resident TE coverage and identity calibration" });
    [0, 20, 40, 60, 80, 100].forEach((value) => {
      svg.appendChild(svgElement("line", { x1: x(value), y1: margin.top, x2: x(value), y2: height - margin.bottom, class: "trm-scatter-grid trm-scatter-grid-vertical" }));
      svg.appendChild(svgElement("line", { x1: margin.left, y1: y(value), x2: width - margin.right, y2: y(value), class: "trm-scatter-grid" }));
      const xLabel = svgElement("text", { x: x(value), y: height - margin.bottom + 22, class: "trm-scatter-axis-label", "text-anchor": "middle" });
      xLabel.textContent = `${value}%`;
      svg.appendChild(xLabel);
      const yLabel = svgElement("text", { x: margin.left - 10, y: y(value) + 4, class: "trm-scatter-axis-label", "text-anchor": "end" });
      yLabel.textContent = `${value}%`;
      svg.appendChild(yLabel);
    });
    const stride = Math.max(1, Math.ceil(residentCopies.length / 2500));
    residentCopies.filter((copy, index) => index % stride === 0).forEach((copy) => {
      const colors = { full_length: "#087f78", partial: "#ed8a3b", degraded_relic: "#a8b4b9" };
      const point = svgElement("circle", { cx: x(copy.coverage), cy: y(copy.identity), r: 3.2, fill: colors[residentTier(copy)], class: "trm-resident-point" });
      const title = svgElement("title");
      title.textContent = `${copy.primary_te} · ${copy.copy_id}\n${copy.target}:${copy.chrom}:${copy.start}-${copy.end}\ncoverage ${copy.coverage}% · identity ${copy.identity}%`;
      point.appendChild(title);
      svg.appendChild(point);
    });
    svg.appendChild(svgElement("line", { x1: x(residentState.fullCoverage), y1: margin.top, x2: x(residentState.fullCoverage), y2: height - margin.bottom, class: "trm-threshold-line" }));
    svg.appendChild(svgElement("line", { x1: margin.left, y1: y(residentState.fullIdentity), x2: width - margin.right, y2: y(residentState.fullIdentity), class: "trm-threshold-line" }));
    svg.appendChild(svgElement("line", { x1: x(residentState.partialCoverage), y1: margin.top, x2: x(residentState.partialCoverage), y2: height - margin.bottom, class: "trm-threshold-line trm-threshold-line-secondary" }));
    const xTitle = svgElement("text", { x: margin.left + plotWidth / 2, y: height - 8, class: "trm-scatter-axis-title", "text-anchor": "middle" });
    xTitle.textContent = "Consensus coverage";
    const yTitle = svgElement("text", { x: 16, y: margin.top + plotHeight / 2, class: "trm-scatter-axis-title", transform: `rotate(-90 16 ${margin.top + plotHeight / 2})`, "text-anchor": "middle" });
    yTitle.textContent = "Identity";
    svg.append(xTitle, yTitle);
    target.replaceChildren(svg);
  }

  function renderResidentEvidence() {
    const fragmentRows = Object.entries((resident.summary || {}).fragment_status || {}).map(([label, value]) => ({
      label: label.replace("rejected:", "Rejected · "), value,
      className: label === "accepted" ? "trm-bar-neutral" : "trm-bar-missing",
    }));
    renderSimpleBars("trm-resident-fragment-summary", fragmentRows);
    const relationRows = Object.entries((resident.summary || {}).relation_types || {}).map(([label, value]) => ({ label, value, className: "trm-bar-outsider" }));
    renderSimpleBars("trm-resident-relation-summary", relationRows);
  }

  function renderResidentTable() {
    const body = $("trm-resident-table-body");
    body.replaceChildren();
    const pages = Math.max(1, Math.ceil(residentState.filtered.length / PAGE_SIZE));
    residentState.page = Math.min(residentState.page, pages - 1);
    const start = residentState.page * PAGE_SIZE;
    residentState.filtered.slice(start, start + PAGE_SIZE).forEach((copy) => {
      const row = document.createElement("tr");
      appendTextCell(row, copy.target);
      appendTextCell(row, `${copy.chrom}:${formatInteger.format(copy.start)}–${formatInteger.format(copy.end)} (${copy.strand})`);
      appendTextCell(row, copy.primary_te);
      const candidates = appendTextCell(row, residentCandidates(copy));
      candidates.className = "trm-candidate-cell";
      appendTextCell(row, valueOrDash(copy.coverage, "%"));
      appendTextCell(row, valueOrDash(copy.identity, "%"));
      appendTextCell(row, copy.aligned_bp);
      appendTextCell(row, copy.fragments);
      appendTextCell(row, copy.status);
      appendTextCell(row, copy.structure);
      body.appendChild(row);
    });
    $("trm-resident-page-status").textContent = residentState.filtered.length ? `Page ${residentState.page + 1} / ${pages}` : "No result";
    $("trm-resident-page-previous").disabled = residentState.page === 0;
    $("trm-resident-page-next").disabled = residentState.page >= pages - 1;
    $("trm-resident-filter-count").textContent = `${formatInteger.format(residentState.filtered.length)} / ${formatInteger.format(residentCopies.length)} copies`;
  }

  function changeResidentPage(delta) {
    residentState.page += delta;
    renderResidentTable();
  }

  function downloadResidentCopies() {
    const columns = ["target", "chrom", "start", "end", "strand", "copy_id", "primary_te", "te_candidates", "coverage", "identity", "aligned_bp", "fragments", "candidate_count", "status", "structure"];
    const lines = [columns.join("\t")];
    residentState.filtered.forEach((copy) => {
      lines.push(columns.map((column) => {
        const value = column === "te_candidates" ? copy.te_candidates.join(";") : copy[column];
        return value === null || value === undefined ? "." : String(value).replace(/[\t\r\n]/g, " ");
      }).join("\t"));
    });
    const blob = new Blob([`${lines.join("\n")}\n`], { type: "text/tab-separated-values;charset=utf-8" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "TrEMOLO.resident-TE-copies.tsv";
    link.click();
    setTimeout(() => URL.revokeObjectURL(link.href), 1000);
  }

  function renderResidentSection() {
    if (!resident.available) {
      $("trm-resident-content").hidden = true;
      const note = $("trm-resident-unavailable");
      note.hidden = false;
      note.textContent = "Resident TE annotation is not available in this run. Enable CHOICE.PIPELINE.TE_GENOME to populate this independent assembly catalog.";
      return;
    }
    renderResidentCards();
    initializeResidentThresholds();
    initializeResidentFilters();
    renderResidentCalibration();
    renderResidentEvidence();
    ["trm-structure-identity", "trm-structure-coverage"].forEach((id) => $(id).addEventListener("input", renderStructurePreview));
    $("trm-structure-reset").addEventListener("click", () => {
      $("trm-structure-identity").value = 0;
      $("trm-structure-coverage").value = 0;
      renderStructurePreview();
    });
    renderStructurePreview();
    applyResidentFilters();
  }

  function initializeTimings() {
    const timings = (data.timings || {}).steps || [];
    const branch = $("trm-timing-branch");
    const order = $("trm-timing-order");
    if (!branch || !order) return;
    selectOptions(branch, [...new Set(timings.map((item) => item.branch))].sort(), "All branches");
    function selected() {
      return timings.filter((item) => !branch.value || item.branch === branch.value).sort((a, b) =>
        order.value === "step" ? a.step.localeCompare(b.step) : b.seconds - a.seconds || a.step.localeCompare(b.step));
    }
    function duration(seconds) {
      if (seconds < 60) return `${formatDecimal.format(seconds)} s`;
      if (seconds < 3600) return `${formatDecimal.format(seconds / 60)} min`;
      return `${formatDecimal.format(seconds / 3600)} h`;
    }
    function render() {
      const rows = selected();
      const target = $("trm-timing-chart");
      target.replaceChildren();
      $("trm-timing-summary").textContent = `${rows.length} measured jobs · sum of job durations: ${duration(rows.reduce((sum, item) => sum + item.seconds, 0))}`;
      if (!rows.length) { emptyChart(target, "No execution-time measurements available."); return; }
      const maximum = Math.max(...rows.map((item) => item.seconds), 0.001);
      rows.forEach((item) => {
        const row = document.createElement("div");
        row.className = "trm-bar-row";
        const label = document.createElement("span");
        label.className = "trm-bar-label";
        label.textContent = item.step;
        label.title = `${item.branch} · ${item.step} · measurement ${item.measurement}`;
        const track = document.createElement("div");
        track.className = "trm-bar-track";
        const bar = document.createElement("span");
        bar.className = `trm-bar ${item.branch === "INSIDER" ? "trm-bar-insider" : item.branch === "OUTSIDER" ? "trm-bar-outsider" : "trm-bar-neutral"}`;
        bar.style.width = `${100 * item.seconds / maximum}%`;
        bar.title = `${item.seconds} seconds`;
        track.appendChild(bar);
        const value = document.createElement("span");
        value.className = "trm-bar-value";
        value.textContent = duration(item.seconds);
        row.append(label, track, value);
        target.appendChild(row);
      });
    }
    branch.addEventListener("change", render);
    order.addEventListener("change", render);
    $("trm-download-timings").addEventListener("click", () => {
      const columns = ["step", "branch", "measurement", "seconds", "benchmark"];
      const lines = [columns.join("\t"), ...selected().map((item) => columns.map((key) => String(item[key]).replace(/[\t\r\n]/g, " ")).join("\t"))];
      const url = URL.createObjectURL(new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values;charset=utf-8" }));
      const link = document.createElement("a");
      link.href = url; link.download = "TrEMOLO.execution-times.tsv"; link.click();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    });
    render();
  }

  function renderIntegration() {
    const integration = data.integration || {};
    const section = $("trm-integration-section");
    if (!section || !integration.available) return;
    section.hidden = false;
    $("trm-integration-summary").textContent = `${integration.complete ? "Complete" : "PARTIAL"} reconstruction: ${integration.integrated} of ${integration.expected_insertions} eligible insertions integrated; ${integration.rejected} rejected; ${integration.excluded} non-insertion events excluded.`;
    const rejected = $("trm-integration-rejections");
    rejected.replaceChildren();
    (integration.rejections || []).forEach((item) => {
      const line = document.createElement("li");
      line.textContent = `${item.event_id}: ${item.reason}`;
      rejected.appendChild(line);
    });
    rejected.hidden = !(integration.rejections || []).length;
  }

  function renderSelection() {
    renderCards();
    renderDetectorSummary();
    renderGenomeHeatmap();
    renderFamilyChart();
    renderFrequencyChart();
    renderGenomeChart();
    renderFrequencyPositionChart();
    renderFrequencyDepthChart();
    renderEvidence();
    renderTable();
  }

  document.title = data.report.title || document.title;
  $("trm-report-title").textContent = data.report.title || "Transposable-element landscape";
  const sourceText = (data.report.enabled_sources || []).join(" + ") || "No enabled evidence source";
  $("trm-report-subtitle").textContent = `${sourceText} · ${data.summary.calls} reported calls · data schema ${data.schema_version}`;
  initializeFilters();
  renderSelection();
  renderProximity();
  $("trm-ambiguous-page-previous").addEventListener("click", () => {
    ambiguousState.page -= 1;
    renderAmbiguousCalls();
  });
  $("trm-ambiguous-page-next").addEventListener("click", () => {
    ambiguousState.page += 1;
    renderAmbiguousCalls();
  });
  renderAmbiguousCalls();
  renderResidentSection();
  renderContext();
  initializeTimings();
  renderIntegration();
})();
