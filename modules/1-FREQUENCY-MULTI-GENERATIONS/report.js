(() => {
  "use strict";
  const data = JSON.parse(document.getElementById("population-report-data").textContent);
  const observations = data.observations || [];
  const trajectories = data.trajectories || [];
  const samples = data.samples || [];
  const byAllele = new Map();
  observations.forEach((row) => {
    if (!byAllele.has(row.allele_id)) byAllele.set(row.allele_id, []);
    byAllele.get(row.allele_id).push(row);
  });
  const $ = (id) => document.getElementById(id);
  const ns = "http://www.w3.org/2000/svg";
  const number = new Intl.NumberFormat(undefined, { maximumFractionDigits: 4 });
  const integer = new Intl.NumberFormat();
  const colors = ["#087f78", "#d66a26", "#5868b2", "#b04f72", "#6b8e23", "#7651a1", "#2f7ebc", "#b58a21"];
  const colorByTimepoint = new Map(samples.map((sample, index) => [sample.timepoint_label, colors[index % colors.length]]));

  function svgElement(name, attributes = {}, text = "") {
    const node = document.createElementNS(ns, name);
    Object.entries(attributes).forEach(([key, value]) => node.setAttribute(key, value));
    if (text) node.textContent = text;
    return node;
  }

  function renderSummary() {
    const summary = data.summary || {};
    const values = [
      [summary.samples || 0, "samples"],
      [summary.loci || 0, "provisional loci"],
      [summary.alleles || 0, "alleles"],
      [summary.observation_status?.observed || 0, "quantified observations"],
      [summary.observation_status?.missing || 0, "missing combinations"],
    ];
    const root = $("trm-pop-summary");
    values.forEach(([value, label]) => {
      const card = document.createElement("div");
      card.className = "trm-pop-card";
      const strong = document.createElement("strong");
      strong.textContent = integer.format(value);
      const span = document.createElement("span");
      span.textContent = label;
      card.append(strong, span);
      root.appendChild(card);
    });
  }

  function populateFilters() {
    const chroms = [...new Set(trajectories.map((row) => row.chrom))].sort();
    chroms.forEach((chrom) => $("trm-pop-chrom").appendChild(new Option(chrom, chrom)));
    [...new Set(trajectories.map((row) => row.te_family))].sort().forEach((family) => $("trm-pop-family").appendChild(new Option(family, family)));
    [...new Set(trajectories.map((row) => row.trend))].sort().forEach((trend) => $("trm-pop-trend").appendChild(new Option(trend.replaceAll("_", " "), trend)));
    ["trm-pop-chrom", "trm-pop-family", "trm-pop-trend", "trm-pop-min-observed"].forEach((id) => $(id).addEventListener("input", renderSelection));
  }

  function selectedTrajectories() {
    const chrom = $("trm-pop-chrom").value;
    const family = $("trm-pop-family").value;
    const trend = $("trm-pop-trend").value;
    const minimum = Math.max(1, Number($("trm-pop-min-observed").value) || 1);
    return trajectories.filter((row) => row.chrom === chrom && (!family || row.te_family === family) && (!trend || row.trend === trend) && Number(row.observed_samples) >= minimum);
  }

  function drawAxes(svg, xTicks, xScale, xLabel) {
    const left = 78, right = 1170, top = 24, bottom = 370;
    [0, 0.25, 0.5, 0.75, 1].forEach((value) => {
      const y = bottom - value * (bottom - top);
      svg.append(svgElement("line", { x1: left, x2: right, y1: y, y2: y, class: "trm-pop-grid" }));
      svg.append(svgElement("text", { x: left - 12, y: y + 5, "text-anchor": "end", class: "trm-pop-label" }, number.format(value)));
    });
    xTicks.forEach(([value, label]) => {
      const x = xScale(value);
      svg.append(svgElement("line", { x1: x, x2: x, y1: top, y2: bottom, class: "trm-pop-grid" }));
      svg.append(svgElement("text", { x, y: bottom + 24, "text-anchor": "middle", class: "trm-pop-label" }, label));
    });
    svg.append(svgElement("line", { x1: left, x2: right, y1: bottom, y2: bottom, class: "trm-pop-axis" }));
    svg.append(svgElement("line", { x1: left, x2: left, y1: top, y2: bottom, class: "trm-pop-axis" }));
    svg.append(svgElement("text", { x: (left + right) / 2, y: 420, "text-anchor": "middle", class: "trm-pop-label" }, xLabel));
    svg.append(svgElement("text", { x: 18, y: 200, transform: "rotate(-90 18 200)", "text-anchor": "middle", class: "trm-pop-label" }, "Allele frequency"));
  }

  function renderScatter(rows) {
    const svg = $("trm-pop-scatter");
    svg.replaceChildren();
    if (!rows.length) return;
    const left = 78, right = 1170, top = 24, bottom = 370;
    const maximum = Number(data.chromosome_lengths?.[rows[0].chrom]) || Math.max(...rows.map((row) => Number(row.end)), 1);
    const xScale = (value) => left + Number(value) / maximum * (right - left);
    drawAxes(svg, [0, 0.25, 0.5, 0.75, 1].map((fraction) => [fraction * maximum, integer.format(fraction * maximum)]), xScale, "Genomic position (bp)");
    rows.forEach((trajectory) => {
      (byAllele.get(trajectory.allele_id) || []).filter((row) => row.observation_status === "observed").forEach((row) => {
        const point = svgElement("circle", {
          cx: xScale(trajectory.start),
          cy: bottom - Number(row.frequency) * (bottom - top),
          r: 5.5,
          fill: colorByTimepoint.get(row.timepoint) || "#667881",
          class: "trm-pop-point",
        });
        point.appendChild(svgElement("title", {}, `${trajectory.te_family} · ${row.sample_id} · ${number.format(row.frequency)}`));
        point.addEventListener("click", () => renderTrajectory(trajectory));
        svg.appendChild(point);
      });
    });
  }

  function renderTable(rows) {
    const body = $("trm-pop-trajectory-body");
    body.replaceChildren();
    rows.slice(0, 1000).forEach((row) => {
      const tr = document.createElement("tr");
      [row.chrom, integer.format(row.start), row.te_family, row.strand, row.observed_samples, row.missing_samples, row.observed_timepoints, row.slope_per_timepoint, row.r_squared, row.trend.replaceAll("_", " ")].forEach((value) => {
        const td = document.createElement("td");
        td.textContent = value === "." ? "—" : value;
        tr.appendChild(td);
      });
      tr.addEventListener("click", () => renderTrajectory(row));
      body.appendChild(tr);
    });
  }

  function renderTrajectory(trajectory) {
    const rows = (byAllele.get(trajectory.allele_id) || []).slice().sort((a, b) => samples.find((item) => item.sample_id === a.sample_id).timepoint - samples.find((item) => item.sample_id === b.sample_id).timepoint);
    $("trm-pop-selected-title").textContent = `${trajectory.te_family} · ${trajectory.chrom}:${integer.format(trajectory.start)}`;
    $("trm-pop-selected-note").textContent = `${trajectory.trend.replaceAll("_", " ")} · ${trajectory.observed_samples} observed, ${trajectory.missing_samples} missing. Missing observations are shown in the table but are not drawn as zero.`;
    const body = $("trm-pop-observation-body");
    body.replaceChildren();
    rows.forEach((row) => {
      const tr = document.createElement("tr");
      [row.sample_id, row.timepoint, row.replicate, row.frequency === "." ? "—" : number.format(row.frequency), row.observation_status.replaceAll("_", " "), row.filter].forEach((value, index) => {
        const td = document.createElement("td");
        td.textContent = value;
        if (index === 4) td.className = `trm-pop-status-${row.observation_status}`;
        tr.appendChild(td);
      });
      body.appendChild(tr);
    });
    const svg = $("trm-pop-trajectory");
    svg.replaceChildren();
    const observed = rows.filter((row) => row.observation_status === "observed").map((row) => ({ ...row, numericTimepoint: samples.find((item) => item.sample_id === row.sample_id).timepoint }));
    if (!observed.length) return;
    const left = 78, right = 1170, top = 24, bottom = 370;
    const minTime = Math.min(...samples.map((sample) => sample.timepoint));
    const maxTime = Math.max(...samples.map((sample) => sample.timepoint));
    const span = Math.max(maxTime - minTime, 1);
    const xScale = (value) => left + (Number(value) - minTime) / span * (right - left);
    drawAxes(svg, samples.map((sample) => [sample.timepoint, sample.timepoint_label]), xScale, "Timepoint");
    const points = observed.map((row) => `${xScale(row.numericTimepoint)},${bottom - Number(row.frequency) * (bottom - top)}`).join(" ");
    if (observed.length > 1) svg.appendChild(svgElement("polyline", { points, class: "trm-pop-line" }));
    observed.forEach((row) => svg.appendChild(svgElement("circle", { cx: xScale(row.numericTimepoint), cy: bottom - Number(row.frequency) * (bottom - top), r: 7, fill: colorByTimepoint.get(row.timepoint) || "#087f78", class: "trm-pop-point" })));
  }

  function renderSelection() {
    const rows = selectedTrajectories();
    $("trm-pop-filter-note").textContent = `${integer.format(rows.length)} allele trajectories match the filters. Up to 1,000 rows are displayed in the catalogue.`;
    renderScatter(rows);
    renderTable(rows);
    if (rows.length) renderTrajectory(rows[0]);
  }

  renderSummary();
  populateFilters();
  renderSelection();
})();
