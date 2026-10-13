const $ = (id) => document.getElementById(id);
const COLORS = { A: "#39f0cf", T: "#7299ff", C: "#ffb86b", G: "#d28cff" };
let latestAnalysis = null;
let toastTimer;

function showToast(message) {
  const toast = $("toast");
  toast.textContent = message;
  toast.classList.add("show");
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => toast.classList.remove("show"), 2600);
}
function sequenceInput() { return $("sequence").value; }
function updateCount() {
  const seq = sequenceInput().replace(/\s/g, "");
  $("baseCount").textContent = `${seq.length.toLocaleString()} ${seq.length === 1 ? "base" : "bases"}`;
}
const CODONS = {
  UUU:"F",UUC:"F",UUA:"L",UUG:"L",UCU:"S",UCC:"S",UCA:"S",UCG:"S",UAU:"Y",UAC:"Y",UAA:"*",UAG:"*",UGU:"C",UGC:"C",UGA:"*",UGG:"W",
  CUU:"L",CUC:"L",CUA:"L",CUG:"L",CCU:"P",CCC:"P",CCA:"P",CCG:"P",CAU:"H",CAC:"H",CAA:"Q",CAG:"Q",CGU:"R",CGC:"R",CGA:"R",CGG:"R",
  AUU:"I",AUC:"I",AUA:"I",AUG:"M",ACU:"T",ACC:"T",ACA:"T",ACG:"T",AAU:"N",AAC:"N",AAA:"K",AAG:"K",AGU:"S",AGC:"S",AGA:"R",AGG:"R",
  GUU:"V",GUC:"V",GUA:"V",GUG:"V",GCU:"A",GCC:"A",GCA:"A",GCG:"A",GAU:"D",GAC:"D",GAA:"E",GAG:"E",GGU:"G",GGC:"G",GGA:"G",GGG:"G"
};
function cleanDNA(value) {
  if (typeof value !== "string") throw new Error("Enter a DNA sequence as text.");
  const lines = value.split(/\r?\n/), headers = lines.filter(line => line.trimStart().startsWith(">"));
  if (headers.length > 1) throw new Error("Use one FASTA record at a time.");
  const sequence = lines.filter(line => !line.trimStart().startsWith(">" )).join("").toUpperCase().replace(/\s/g, "");
  if (!sequence) throw new Error("Enter a DNA sequence.");
  const invalid = [...new Set(sequence.replace(/[ATCG]/g, ""))].sort();
  if (invalid.length) throw new Error(`DNA sequence contains invalid character(s): ${invalid.join(", ")}`);
  if (sequence.length > 200000) throw new Error("For this demo, sequences are limited to 200,000 bases.");
  return sequence;
}
function localAnalyze(value) {
  const sequence = cleanDNA(value), counts = { A: 0, T: 0, C: 0, G: 0 };
  for (const base of sequence) counts[base]++;
  const percentages = Object.fromEntries("ATCG".split("").map(base => [base, counts[base] / sequence.length * 100]));
  const complement = sequence.replace(/[ATCG]/g, base => ({ A:"T", T:"A", C:"G", G:"C" })[base]);
  const rna = sequence.replace(/T/g, "U"), protein = [];
  for (let i = 0; i + 2 < rna.length; i += 3) {
    const amino = CODONS[rna.slice(i, i + 3)];
    if (!amino || amino === "*") break;
    protein.push(amino);
  }
  return { sequence, length: sequence.length, counts, percentages,
    gc_content: (counts.G + counts.C) / sequence.length * 100,
    complement, reverse_complement: [...complement].reverse().join(""), rna, protein: protein.join("") };
}
function localOrfs(value) {
  const sequence = cleanDNA(value), complement = sequence.replace(/[ATCG]/g, base => ({ A:"T", T:"A", C:"G", G:"C" })[base]);
  const orfs = [];
  for (const [strand, dna] of [["+", sequence], ["-", [...complement].reverse().join("")]]) {
    const rna = dna.replace(/T/g, "U");
    for (let frame = 0; frame < 3; frame++) {
      let starts = [];
      for (let i = frame; i + 2 < rna.length; i += 3) {
        const codon = rna.slice(i, i + 3);
        if (codon === "AUG") starts.push(i);
        else if (["UAA", "UAG", "UGA"].includes(codon)) {
          for (const start of starts) {
            const region = dna.slice(start, i + 3);
            orfs.push({ strand, frame: frame + 1, start: start + 1, end: i + 3,
              length: i + 3 - start, sequence: region, protein: localAnalyze(region).protein });
          }
          starts = [];
        }
      }
    }
  }
  orfs.sort((a, b) => a.strand.localeCompare(b.strand) || a.start - b.start || a.frame - b.frame);
  return { orfs };
}
function localApi(path, payload) {
  if (path === "/api/analyze") return localAnalyze(payload.sequence);
  if (path === "/api/motif") {
    const sequence = cleanDNA(payload.sequence), motif = cleanDNA(payload.motif), positions = [];
    let start = 0, index;
    while ((index = sequence.indexOf(motif, start)) >= 0) { positions.push(index + 1); start = index + 1; }
    return { motif, positions, count: positions.length };
  }
  if (path === "/api/gc-window") {
    const sequence = cleanDNA(payload.sequence), window = Number(payload.window), step = Number(payload.step);
    if (!Number.isInteger(window) || window < 1 || window > sequence.length) throw new Error("Window size must be a whole number between 1 and the sequence length.");
    if (!Number.isInteger(step) || step < 1) throw new Error("Step size must be a whole number of at least 1.");
    const gcPrefix = [0];
    for (const base of sequence) gcPrefix.push(gcPrefix[gcPrefix.length - 1] + (base === "G" || base === "C" ? 1 : 0));
    const points = [];
    for (let start = 0; start <= sequence.length - window; start += step) {
      const gc = gcPrefix[start + window] - gcPrefix[start];
      points.push({ position: start + 1, gc: gc / window * 100 });
    }
    return { window, step, points };
  }
  if (path === "/api/compare") {
    const sequence1 = cleanDNA(payload.sequence1), sequence2 = cleanDNA(payload.sequence2);
    if (sequence1.length !== sequence2.length) throw new Error("For this simple comparison, both sequences must have the same length.");
    const mismatch_positions = [];
    for (let i = 0; i < sequence1.length; i++) if (sequence1[i] !== sequence2[i]) mismatch_positions.push(i + 1);
    const mismatches = mismatch_positions.length;
    return { sequence1, sequence2, length: sequence1.length, mismatches, mismatch_positions,
      matches: sequence1.length - mismatches, identity: (sequence1.length - mismatches) / sequence1.length * 100 };
  }
  if (path === "/api/orfs") return localOrfs(payload.sequence);
  throw new Error("This analysis is not available in offline mode.");
}
async function api(path, payload) {
  if (location.protocol === "file:") return localApi(path, payload);
  const response = await fetch(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload)
  });
  const data = await response.json().catch(() => ({}));
  if (!response.ok) {
    if (response.status === 404 || response.status === 405) return localApi(path, payload);
    throw new Error(data.error || "The request could not be completed.");
  }
  return data;
}
function showError(message) {
  $("sequenceError").textContent = message;
  $("sequenceError").hidden = false;
}
function clearError() {
  $("sequenceError").hidden = true;
  $("sequenceError").textContent = "";
}
function setBusy(button, busy, label) {
  button.disabled = busy;
  if (busy) {
    button.dataset.oldText = button.innerHTML;
    button.textContent = label || "Working…";
  } else if (button.dataset.oldText) {
    button.innerHTML = button.dataset.oldText;
    delete button.dataset.oldText;
  }
}
function svgNode(tag, attrs = {}, text = "") {
  const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
  Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
  if (text !== "") node.textContent = text;
  return node;
}
function drawBars(counts, percentages) {
  const host = $("compositionChart");
  host.replaceChildren();
  const svg = svgNode("svg", { viewBox: "0 0 480 220", role: "img", "aria-label": "Bar chart of nucleotide percentages" });
  const left = 42, top = 12, chartW = 415, chartH = 160;
  [0, 25, 50, 75, 100].forEach((tick) => {
    const y = top + chartH - tick / 100 * chartH;
    svg.append(svgNode("line", { x1: left, y1: y, x2: 466, y2: y, stroke: "#243746", "stroke-dasharray": "3 5" }));
    svg.append(svgNode("text", { x: left - 10, y: y + 4, fill: "#728b9b", "font-size": 10, "text-anchor": "end", "font-family": "DM Mono, monospace" }, String(tick)));
  });
  ["A", "T", "C", "G"].forEach((base, i) => {
    const x = left + 38 + i * 99;
    const h = percentages[base] / 100 * chartH;
    svg.append(svgNode("rect", { x, y: top + chartH - h, width: 48, height: Math.max(0, h), rx: 5, fill: COLORS[base], opacity: .88 }));
    svg.append(svgNode("text", { x: x + 24, y: 192, fill: "#c1d0d9", "font-size": 12, "text-anchor": "middle", "font-family": "DM Mono, monospace" }, base));
    svg.append(svgNode("text", { x: x + 24, y: Math.max(10, top + chartH - h - 9), fill: COLORS[base], "font-size": 10, "text-anchor": "middle", "font-family": "DM Mono, monospace" }, `${percentages[base].toFixed(1)}%`));
  });
  host.append(svg);
  const legend = $("compositionLegend");
  legend.replaceChildren();
  ["A", "T", "C", "G"].forEach(base => {
    const item = document.createElement("span");
    const swatch = document.createElement("i");
    swatch.style.background = COLORS[base];
    item.append(swatch, document.createTextNode(`${base} · ${counts[base]}`));
    legend.append(item);
  });
}
function drawBaseStats(result) {
  const host = $("baseStats");
  host.replaceChildren();
  ["A", "T", "C", "G"].forEach(base => {
    const row = document.createElement("div");
    row.className = "base-stat";
    row.style.setProperty("--base-color", COLORS[base]);
    const letter = document.createElement("span");
    letter.className = "base-letter";
    letter.textContent = base;
    const bar = document.createElement("div");
    bar.className = "base-bar";
    const fill = document.createElement("span");
    fill.style.setProperty("--width", `${result.percentages[base]}%`);
    bar.append(fill);
    const value = document.createElement("span");
    value.className = "base-value";
    value.textContent = `${result.counts[base]} · ${result.percentages[base].toFixed(1)}%`;
    row.append(letter, bar, value);
    host.append(row);
  });
}
function drawLineChart(host, points) {
  host.replaceChildren();
  if (!points.length) { host.textContent = "No windows to display."; return; }
  const width = 480, height = 170, left = 36, right = 12, top = 12, bottom = 28;
  const chartW = width - left - right, chartH = height - top - bottom;
  const svg = svgNode("svg", { viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": "GC percentage by window start position" });
  [0, 50, 100].forEach(tick => {
    const y = top + chartH - tick / 100 * chartH;
    svg.append(svgNode("line", { x1: left, y1: y, x2: width - right, y2: y, stroke: "#243746", "stroke-dasharray": "3 5" }));
    svg.append(svgNode("text", { x: left - 7, y: y + 3, fill: "#728b9b", "font-size": 9, "text-anchor": "end" }, `${tick}`));
  });
  // Keep SVG rendering responsive when a small step produces many windows.
  const stride = Math.max(1, Math.ceil(points.length / 1000));
  const plotted = points.filter((_, i) => i === 0 || i === points.length - 1 || i % stride === 0);
  const maxPos = Math.max(1, plotted.length - 1);
  const coords = plotted.map((p, i) => ({
    x: left + i / maxPos * chartW,
    y: top + chartH - p.gc / 100 * chartH
  }));
  const path = coords.map((p, i) => `${i ? "L" : "M"} ${p.x} ${p.y}`).join(" ");
  svg.append(svgNode("path", { d: path, fill: "none", stroke: "#39f0cf", "stroke-width": 2.5, "stroke-linejoin": "round", "stroke-linecap": "round" }));
  if (coords.length <= 120) coords.forEach(p => svg.append(svgNode("circle", { cx: p.x, cy: p.y, r: 2.2, fill: "#a0fff0" })));
  const first = points[0].position, last = points[points.length - 1].position;
  svg.append(svgNode("text", { x: left, y: height - 7, fill: "#728b9b", "font-size": 9 }, `Pos ${first}`));
  svg.append(svgNode("text", { x: width - right, y: height - 7, fill: "#728b9b", "font-size": 9, "text-anchor": "end" }, `Pos ${last}`));
  host.append(svg);
}
function addText(parent, tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  node.textContent = text;
  parent.append(node);
  return node;
}
async function runAnalysis() {
  const button = $("analyzeBtn");
  clearError();
  setBusy(button, true, "Analyzing…");
  try {
    const result = await api("/api/analyze", { sequence: sequenceInput() });
    const currentSequence = sequenceInput().toUpperCase().replace(/\s/g, "");
    if (currentSequence !== result.sequence) return;
    latestAnalysis = result;
    $("results").hidden = false;
    $("lengthValue").textContent = result.length.toLocaleString();
    $("gcValue").textContent = `${result.gc_content.toFixed(2)}%`;
    $("atValue").textContent = `${(100 - result.gc_content).toFixed(2)}%`;
    drawBars(result.counts, result.percentages);
    drawBaseStats(result);
    $("complementValue").textContent = result.complement || "—";
    $("reverseValue").textContent = result.reverse_complement || "—";
    $("rnaValue").textContent = result.rna || "—";
    $("proteinValue").textContent = result.protein || "(no complete amino acids before a stop codon)";
    $("results").scrollIntoView({ behavior: "smooth", block: "start" });
    showToast("Sequence analysis complete.");
  } catch (error) {
    showError(error.message);
  } finally {
    setBusy(button, false);
  }
}
function requireAnalysis() {
  if (!latestAnalysis) {
    showToast("Analyze a valid DNA sequence first.");
    $("sequence").focus();
    return false;
  }
  return true;
}
async function runMotif() {
  if (!requireAnalysis()) return;
  const out = $("motifOutput");
  out.textContent = "Searching…";
  try {
    const result = await api("/api/motif", { sequence: sequenceInput(), motif: $("motifInput").value });
    out.replaceChildren();
    addText(out, "div", "", `${result.count} occurrence${result.count === 1 ? "" : "s"} found`);
    addText(out, "div", "muted", result.positions.length ? `1-based positions: ${result.positions.join(", ")}` : "No matches in this sequence.");
  } catch (error) { out.textContent = error.message; }
}
async function runGC() {
  if (!requireAnalysis()) return;
  const button = $("gcBtn");
  setBusy(button, true, "Plotting…");
  try {
    const result = await api("/api/gc-window", {
      sequence: sequenceInput(),
      window: Number($("windowInput").value),
      step: Number($("stepInput").value)
    });
    drawLineChart($("gcChart"), result.points);
    $("gcOutput").textContent = `${result.points.length} windows · ${result.window}-base window · ${result.step}-base step`;
  } catch (error) { $("gcOutput").textContent = error.message; }
  finally { setBusy(button, false); }
}
async function runCompare() {
  if (!requireAnalysis()) return;
  const out = $("compareOutput");
  out.textContent = "Comparing…";
  try {
    const result = await api("/api/compare", { sequence1: sequenceInput(), sequence2: $("compareInput").value });
    out.replaceChildren();
    addText(out, "div", "", `${result.identity.toFixed(2)}% identity`);
    addText(out, "div", "muted", `${result.matches} matching · ${result.mismatches} mismatching positions`);
    if (result.mismatch_positions.length) addText(out, "div", "muted", `Mismatch positions: ${result.mismatch_positions.join(", ")}`);
  } catch (error) { out.textContent = error.message; }
}
async function runORFs() {
  if (!requireAnalysis()) return;
  const button = $("orfBtn"), out = $("orfOutput");
  setBusy(button, true, "Scanning…");
  out.textContent = "Scanning both strands…";
  try {
    const result = await api("/api/orfs", { sequence: sequenceInput() });
    out.replaceChildren();
    if (!result.orfs.length) {
      addText(out, "div", "", "No start-to-stop ORFs found.");
    } else {
      addText(out, "div", "", `${result.orfs.length} candidate ORF(s) found`);
      result.orfs.slice(0, 20).forEach((orf, i) => {
        addText(out, "div", "muted", `${i + 1}. Strand ${orf.strand} · frame ${orf.frame} · positions ${orf.start}–${orf.end} · ${orf.length} bases`);
        addText(out, "div", "", orf.protein ? `Protein: ${orf.protein}` : "Protein: —");
      });
      if (result.orfs.length > 20) addText(out, "div", "muted", `Showing first 20 of ${result.orfs.length} candidates.`);
    }
    addText(out, "div", "muted", "Coordinates refer to the strand being scanned; these are candidates, not proof of gene function.");
  } catch (error) { out.textContent = error.message; }
  finally { setBusy(button, false); }
}
async function downloadReport() {
  if (!requireAnalysis()) return;
  const button = $("reportBtn");
  setBusy(button, true, "Preparing report…");
  try {
    let blob;
    if (location.protocol === "file:") {
      const result = localAnalyze(sequenceInput());
      const rows = [["Metric", "Value"], ["Sequence", result.sequence], ["Length (bases)", result.length], ["GC content (%)", result.gc_content]];
      for (const base of "ATCG") rows.push([`${base} count`, result.counts[base]], [`${base} percentage`, result.percentages[base]]);
      rows.push(["Complement", result.complement], ["Reverse complement", result.reverse_complement], ["RNA transcript", result.rna], ["Protein (translation from first base)", result.protein]);
      const escapeCSV = value => `"${String(value).replace(/"/g, '""')}"`;
      blob = new Blob([rows.map(row => row.map(escapeCSV).join(",")).join("\r\n")], { type: "text/csv;charset=utf-8" });
    } else {
      const response = await fetch("/api/report", {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sequence: sequenceInput() })
      });
      if (!response.ok) {
        if (response.status !== 404 && response.status !== 405) {
          const data = await response.json().catch(() => ({}));
          throw new Error(data.error || "Could not generate report.");
        }
        const result = localAnalyze(sequenceInput());
        const rows = [["Metric", "Value"], ["Sequence", result.sequence], ["Length (bases)", result.length], ["GC content (%)", result.gc_content]];
        for (const base of "ATCG") rows.push([`${base} count`, result.counts[base]], [`${base} percentage`, result.percentages[base]]);
        rows.push(["Complement", result.complement], ["Reverse complement", result.reverse_complement], ["RNA transcript", result.rna], ["Protein (translation from first base)", result.protein]);
        const escapeCSV = value => `"${String(value).replace(/"/g, '""')}"`;
        blob = new Blob([rows.map(row => row.map(escapeCSV).join(",")).join("\r\n")], { type: "text/csv;charset=utf-8" });
      } else {
        blob = await response.blob();
      }
    }
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url; link.download = "genescope_report.csv";
    document.body.append(link); link.click(); link.remove();
    // Let the browser begin consuming the blob before releasing its object URL.
    setTimeout(() => URL.revokeObjectURL(url), 1000);
    showToast("CSV report downloaded.");
  } catch (error) { showToast(error.message); }
  finally { setBusy(button, false); }
}
$("sequence").addEventListener("input", () => {
  updateCount();
  latestAnalysis = null;
  $("results").hidden = true;
});
$("fastaFile").addEventListener("change", async event => {
  const file = event.target.files && event.target.files[0];
  if (!file) return;
  try {
    $("sequence").value = await file.text();
    $("sequence").dispatchEvent(new Event("input", { bubbles: true }));
    clearError();
    showToast(`Loaded ${file.name}. Click Analyze sequence to continue.`);
  } catch {
    showError("Could not read that file. Try a plain-text FASTA file.");
  } finally {
    event.target.value = "";
  }
});
$("sampleBtn").addEventListener("click", () => {
  $("sequence").value = "ATGAAATAGGCCATGCGTAA";
  updateCount(); clearError(); showToast("Sample sequence loaded.");
});
$("randomBtn").addEventListener("click", () => {
  const bases = "ATCG", length = 120, values = new Uint8Array(length);
  if (window.crypto && window.crypto.getRandomValues) window.crypto.getRandomValues(values);
  else for (let i = 0; i < length; i++) values[i] = Math.floor(Math.random() * 256);
  $("sequence").value = Array.from(values, value => bases[value & 3]).join("");
  latestAnalysis = null;
  $("results").hidden = true;
  updateCount(); clearError();
  showToast("Generated 120 bases of synthetic random DNA.");
});
$("clearBtn").addEventListener("click", () => {
  $("sequence").value = ""; latestAnalysis = null; updateCount(); clearError();
  $("results").hidden = true; $("sequence").focus();
});
$("analyzeBtn").addEventListener("click", runAnalysis);
$("motifBtn").addEventListener("click", runMotif);
$("gcBtn").addEventListener("click", runGC);
$("compareBtn").addEventListener("click", runCompare);
$("orfBtn").addEventListener("click", runORFs);
$("reportBtn").addEventListener("click", downloadReport);
document.querySelectorAll("[data-copy]").forEach(button => {
  button.addEventListener("click", async () => {
    const text = $(button.dataset.copy).textContent;
    try {
      await navigator.clipboard.writeText(text);
      showToast("Copied to clipboard.");
    } catch {
      const area = document.createElement("textarea");
      area.value = text; document.body.append(area); area.select();
      document.execCommand("copy"); area.remove(); showToast("Copied to clipboard.");
    }
  });
});
updateCount();
