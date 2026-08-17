/* GToTree interactive dashboard — dependency-free (tree, alignment, sunburst).
 * Plotly is used only for the SCG presence heatmap and is bundled inline by the
 * Python generator. All data is injected as window.GT_DATA.
 *
 * window.GT_DATA = {
 *   title, gene_set, generated,
 *   newick: "<newick string or ''>",
 *   alignment: { ids:[...], seqs:[...] } | null,   // aligned .faa
 *   summary: { columns:[...], rows:[[...]] } | null,
 *   hits: { genomes:[...], genes:[...], matrix:[[...]] } | null
 * }
 */
(function () {
  "use strict";
  const D = window.GT_DATA || {};
  const $ = (s, r) => (r || document).querySelector(s);
  const $$ = (s, r) => Array.from((r || document).querySelectorAll(s));
  const NS = "http://www.w3.org/2000/svg";
  const el = (t, a) => { const e = document.createElementNS(NS, t); if (a) for (const k in a) e.setAttribute(k, a[k]); return e; };

  /* ---------------- theme ---------------- */
  const rootEl = document.documentElement;
  function setTheme(t) { rootEl.setAttribute("data-theme", t); try { localStorage.setItem("gt-report-theme", t); } catch (e) {} const b = $("#themeBtn"); if (b) b.textContent = t === "dark" ? "\u2600\ufe0f" : "\u263e"; }
  (function initTheme() { let t = "light"; try { t = localStorage.getItem("gt-report-theme") || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light"); } catch (e) {} setTheme(t); })();

  /* ---------------- nav / tabs ---------------- */
  function showTab(id) {
    $$(".tabpane").forEach(p => p.classList.toggle("active", p.id === "tab-" + id));
    $$(".nav-link").forEach(a => a.classList.toggle("active", a.dataset.tab === id));
    if (id === "tree") renderTree();
    if (id === "alignment") renderAlignment();
    if (id === "scg") renderHeatmap();
    if (id === "taxonomy") renderSunburst();
    if (id === "genomes") renderTable();
  }

  /* =========================================================
   * Newick parser
   * =======================================================*/
  function parseNewick(s) {
    if (!s) return null;
    s = s.trim();
    let i = 0;
    function node() { return { name: "", length: null, children: [] }; }
    function readLabel() {
      let start = i;
      while (i < s.length && !"(),:;".includes(s[i])) i++;
      return s.slice(start, i).trim();
    }
    function parse() {
      const n = node();
      if (s[i] === "(") {
        i++;
        while (true) {
          n.children.push(parse());
          if (s[i] === ",") { i++; continue; }
          if (s[i] === ")") { i++; break; }
          break;
        }
      }
      n.name = readLabel().replace(/'/g, "");
      if (s[i] === ":") { i++; let start = i; while (i < s.length && !"(),:;".includes(s[i])) i++; n.length = parseFloat(s.slice(start, i)); }
      return n;
    }
    const root = parse();
    return root;
  }

  function treeStats(root) {
    let tips = 0, internal = 0, maxDepth = 0;
    (function walk(n, d) {
      if (!n.children.length) { tips++; maxDepth = Math.max(maxDepth, d); }
      else { internal++; n.children.forEach(c => walk(c, d + (c.length || 0))); }
    })(root, 0);
    return { tips, internal, maxDepth };
  }

  /* =========================================================
   * Tree renderer (rectangular phylogram/cladogram + radial)
   * =======================================================*/
  const T = { root: null, layout: "phylogram", radial: false, forcePhylo: false, showSupport: true, zoom: 1, highlight: "", rerootNode: null, collapsed: new Set(), effectiveCladogram: false };

  function assignLeafOrder(n, counter) {
    if (!n.children.length || T.collapsed.has(n)) { n._y = counter.i++; n._leaf = true; return n._y; }
    let sum = 0; n.children.forEach(c => sum += assignLeafOrder(c, counter));
    n._y = sum / n.children.length; n._leaf = false; return n._y;
  }
  function assignX(n, x, depthIdx) {
    n._depthIdx = depthIdx;
    n._xLen = x + (n.length || 0);
    if (n.children.length && !T.collapsed.has(n)) n.children.forEach(c => assignX(c, n._xLen, depthIdx + 1));
  }
  function maxXLen(n, m) { m = Math.max(m, n._xLen || 0); if (n.children.length && !T.collapsed.has(n)) n.children.forEach(c => { m = maxXLen(c, m); }); return m; }
  function maxDepthIdx(n, m) { m = Math.max(m, n._depthIdx || 0); if (n.children.length && !T.collapsed.has(n)) n.children.forEach(c => { m = maxDepthIdx(c, m); }); return m; }

  function collectLeaves(n, arr) { if (!n.children.length || T.collapsed.has(n)) arr.push(n); else n.children.forEach(c => collectLeaves(c, arr)); return arr; }

  function reroot(root, target) {
    // simple midpoint-free reroot on the branch leading to target
    if (target === root) return root;
    const path = [];
    (function find(n, acc) { acc.push(n); if (n === target) { path.push(...acc); return true; } for (const c of n.children) { if (find(c, acc)) return true; } acc.pop(); return false; })(root, []);
    if (!path.length) return root;
    // rebuild: make a new root between target and its parent
    for (let k = path.length - 1; k > 0; k--) {
      const child = path[k], parent = path[k - 1];
      parent.children = parent.children.filter(c => c !== child);
    }
    const newRoot = { name: "", length: null, children: [target] };
    // re-attach the old chain as sibling
    let attach = newRoot, prevLen = target.length;
    for (let k = path.length - 2; k >= 0; k--) {
      const node = path[k];
      const nn = { name: node.name, length: prevLen, children: node.children.slice() };
      prevLen = node.length;
      attach.children.push(nn);
      attach = nn;
    }
    target.length = 0;
    return newRoot;
  }

  function renderTree() {
    const host = $("#treeCanvas");
    if (!host) return;
    if (!T.root) { host.innerHTML = '<p class="empty">No tree file (run was alignment-only, or -N was set).</p>'; return; }
    host.innerHTML = "";
    const leaves = [];
    const counter = { i: 0 };
    assignLeafOrder(T.root, counter);
    assignX(T.root, 0, 0);
    collectLeaves(T.root, leaves);
    const nLeaves = leaves.length;

    const rowH = 20 * T.zoom;
    const labelPad = 8;
    const rawMaxLen = maxXLen(T.root, 0);
    const maxLen = rawMaxLen || 1;
    const maxIdx = maxDepthIdx(T.root, 0) || 1;
    const margin = { top: 24, right: 220, bottom: 24, left: 24 };
    const plotW = 620 * T.zoom;
    const height = margin.top + margin.bottom + nLeaves * rowH;
    const width = margin.left + margin.right + plotW;

    // Degenerate / low-resolution tree: when many tips have ~0-length branches
    // (near-identical genomes), a phylogram draws them with no visible branch.
    // Detect that and fall back to cladogram spacing so every tip gets a line,
    // and tell the user why. (Fully-zero trees also hit this.)
    let zeroTips = 0;
    leaves.forEach(l => { if (!l.length || Math.abs(l.length) < 1e-9) zeroTips++; });
    const lowRes = rawMaxLen < 1e-9 || (leaves.length > 0 && zeroTips / leaves.length >= 0.5);
    const useCladogram = T.layout === "cladogram" || (lowRes && !T.radial && !T.forcePhylo);
    T.effectiveCladogram = useCladogram;
    const xScale = v => margin.left + (useCladogram ? 0 : (v / maxLen) * plotW);
    const xScaleIdx = idx => margin.left + (idx / maxIdx) * plotW;
    const getX = n => useCladogram ? xScaleIdx(n._depthIdx) : xScale(n._xLen);
    const getY = n => margin.top + n._y * rowH + rowH / 2;

    if (T.radial) { renderRadial(host, leaves, maxLen, useCladogram, maxIdx); return; }

    const svg = el("svg", { width, height, viewBox: `0 0 ${width} ${height}`, class: "treeSvg" });
    const hl = (T.highlight || "").toLowerCase();

    (function draw(n) {
      const x = getX(n), y = getY(n);
      if (n.children.length && !T.collapsed.has(n)) {
        const ys = n.children.map(getY);
        const vx = x;
        svg.appendChild(el("path", { d: `M${vx},${Math.min(...ys)} L${vx},${Math.max(...ys)}`, class: "branch" }));
        n.children.forEach(c => {
          const cx = getX(c), cy = getY(c);
          svg.appendChild(el("path", { d: `M${vx},${cy} L${cx},${cy}`, class: "branch" }));
          draw(c);
        });
        // support value label on internal node
        if (T.showSupport && n.name && !isNaN(parseFloat(n.name))) {
          const t = el("text", { x: x + 3, y: y - 3, class: "support" }); t.textContent = n.name; svg.appendChild(t);
        }
        // collapse toggle marker
        const dot = el("circle", { cx: x, cy: y, r: 3.2, class: "node-dot" });
        dot.style.cursor = "pointer";
        dot.addEventListener("click", () => { T.collapsed.has(n) ? T.collapsed.delete(n) : T.collapsed.add(n); renderTree(); });
        svg.appendChild(dot);
      } else {
        // leaf (or collapsed clade)
        const collapsed = T.collapsed.has(n) && n.children.length;
        const label = collapsed ? (n.name || "clade") + " \u25b8 (" + collectLeaves(n, []).length + ")" : n.name;
        const matched = hl && label.toLowerCase().includes(hl);
        if (collapsed) {
          const tri = el("path", { d: `M${x},${y} l10,-6 l0,12 z`, class: "collapsed-tri" });
          tri.style.cursor = "pointer";
          tri.addEventListener("click", () => { T.collapsed.delete(n); renderTree(); });
          svg.appendChild(tri);
        }
        const t = el("text", { x: x + labelPad + (collapsed ? 12 : 0), y: y + 4, class: "tiplabel" + (matched ? " match" : "") });
        t.textContent = label;
        t.style.cursor = "pointer";
        t.addEventListener("click", () => { T.rerootNode = n; T.root = reroot(T.root, n); T.collapsed.clear(); renderTree(); });
        svg.appendChild(t);
        svg.appendChild(el("circle", { cx: x, cy: y, r: matched ? 3.5 : 2, class: "tip-dot" + (matched ? " match" : "") }));
      }
    })(T.root);

    // scale bar (phylogram only)
    if (!useCladogram && maxLen > 0) {
      const barLen = niceScale(maxLen * 0.25);
      const px = (barLen / maxLen) * plotW;
      const by = height - 10, bx = margin.left;
      svg.appendChild(el("path", { d: `M${bx},${by} L${bx + px},${by}`, class: "scalebar" }));
      const tl = el("text", { x: bx + px / 2, y: by - 5, class: "scalelabel", "text-anchor": "middle" }); tl.textContent = barLen; svg.appendChild(tl);
    }
    if (lowRes && useCladogram && !T.radial) {
      const note = document.createElement("div");
      note.className = "tree-note";
      note.innerHTML = `${zeroTips} of ${leaves.length} tips have ~0-length branches (near-identical genomes), which a true phylogram would draw with no visible branch. Showing <strong>cladogram</strong> spacing (topology only) so every tip is visible. Click <strong>Phylogram</strong> to force branch-length scaling anyway.`;
      host.appendChild(note);
    }
    host.appendChild(svg);
    syncTreeBtns();
  }

  function renderRadial(host, leaves, maxLen, useCladogram, maxIdx) {
    const n = leaves.length;
    const size = Math.max(520, 26 * n) * T.zoom;
    const cx = size / 2, cy = size / 2, R = size / 2 - 130;
    const svg = el("svg", { width: size, height: size, viewBox: `0 0 ${size} ${size}`, class: "treeSvg" });
    const ang = node => (node._y / Math.max(1, n - 1)) * 2 * Math.PI;
    const rad = node => useCladogram ? (node._depthIdx / (maxIdx || 1)) * R : (node._xLen / (maxLen || 1)) * R;
    const pt = (r, a) => [cx + r * Math.cos(a - Math.PI / 2), cy + r * Math.sin(a - Math.PI / 2)];
    const hl = (T.highlight || "").toLowerCase();
    (function draw(node) {
      const a = ang(node), r = rad(node);
      const [x, y] = pt(r, a);
      if (node.children.length && !T.collapsed.has(node)) {
        const angs = node.children.map(ang);
        // arc connecting children at radius r
        const a0 = Math.min(...angs), a1 = Math.max(...angs);
        const steps = 12; let d = "";
        for (let k = 0; k <= steps; k++) { const aa = a0 + (a1 - a0) * k / steps; const [px, py] = pt(r, aa); d += (k ? "L" : "M") + px + "," + py; }
        svg.appendChild(el("path", { d, class: "branch", fill: "none" }));
        node.children.forEach(c => { const [cxp, cyp] = pt(rad(c), ang(c)); const [rx, ry] = pt(r, ang(c)); svg.appendChild(el("path", { d: `M${rx},${ry} L${cxp},${cyp}`, class: "branch" })); draw(c); });
      } else {
        const label = node.name;
        const matched = hl && label && label.toLowerCase().includes(hl);
        const [lx, ly] = pt(r + 6, a);
        const deg = (a * 180 / Math.PI) - 90;
        const flip = a > Math.PI;
        const t = el("text", { x: lx, y: ly, class: "tiplabel" + (matched ? " match" : ""), transform: `rotate(${flip ? deg + 180 : deg},${lx},${ly})`, "text-anchor": flip ? "end" : "start" });
        t.textContent = label; svg.appendChild(t);
        svg.appendChild(el("circle", { cx: x, cy: y, r: matched ? 3.5 : 1.8, class: "tip-dot" + (matched ? " match" : "") }));
      }
    })(T.root);
    host.appendChild(svg);
    syncTreeBtns();
  }

  function niceScale(v) { const p = Math.pow(10, Math.floor(Math.log10(v))); const f = v / p; let nf = f >= 5 ? 5 : f >= 2 ? 2 : 1; return +(nf * p).toPrecision(2); }

  function downloadNewick() {
    function toNwk(n) { let s = ""; if (n.children.length) s = "(" + n.children.map(toNwk).join(",") + ")"; s += (n.name || ""); if (n.length != null && !isNaN(n.length)) s += ":" + n.length; return s; }
    const txt = toNwk(T.root) + ";";
    dl(txt, (D.title || "tree") + ".tre", "text/plain");
  }
  // Inline the computed colors so the exported SVG renders standalone (page CSS
  // is not carried into a downloaded .svg file).
  const SVG_STYLE = `
    .branch{stroke:#5b6b7b;stroke-width:1.3;fill:none}
    .tiplabel{font:12px 'JetBrains Mono',ui-monospace,monospace;fill:#0f1e2e}
    .tiplabel.match{fill:#d6336c;font-weight:700}
    .tip-dot{fill:#5b6b7b}
    .tip-dot.match{fill:#d6336c}
    .node-dot{fill:#ffffff;stroke:#5b6b7b;stroke-width:1.2}
    .collapsed-tri{fill:#0891b2}
    .support{font:9px 'JetBrains Mono',monospace;fill:#566878}
    .scalebar{stroke:#0f1e2e;stroke-width:1.5}
    .scalelabel{font:10px 'JetBrains Mono',monospace;fill:#5b6b7b}
    .arc-label{font:10px 'Inter',sans-serif;fill:#fff}
  `;
  // Rebuild the full alignment FASTA from the embedded data and download it,
  // so it works offline straight from the report (no bucket round-trip).
  function downloadAlignment() {
    if (!D.alignment || !D.alignment.ids) return;
    let out = "";
    for (let i = 0; i < D.alignment.ids.length; i++) {
      out += ">" + D.alignment.ids[i] + "\n";
      const seq = D.alignment.seqs[i] || "";
      for (let p = 0; p < seq.length; p += 60) out += seq.slice(p, p + 60) + "\n";
    }
    dl(out, (D.title || "gtotree") + "_Aligned_SCGs.faa", "text/plain");
  }
  // Rebuild the genome summary TSV from the embedded table and download it.
  function downloadSummary() {
    if (!D.summary || !D.summary.columns) return;
    const lines = [D.summary.columns.join("\t")];
    D.summary.rows.forEach(r => lines.push(r.join("\t")));
    dl(lines.join("\n") + "\n", (D.title || "gtotree") + "_Genomes_summary_info.tsv", "text/tab-separated-values");
  }
  function dl(text, name, type) { const b = new Blob([text], { type }); const u = URL.createObjectURL(b); const a = document.createElement("a"); a.href = u; a.download = name; a.click(); URL.revokeObjectURL(u); }

  /* =========================================================
   * Alignment viewer
   * =======================================================*/
  const AA_GROUPS = { hydrophobic: "AVLIMFWC", polar: "STNQ", positive: "KRH", negative: "DE", special: "GP", other: "YXBZ*-" };
  const AA_COLOR = {};
  (function () { const c = { hydrophobic: "#2b8a9e", polar: "#4c9a52", positive: "#3161c9", negative: "#c0392b", special: "#b7791f", other: "#8a8f98" }; for (const g in AA_GROUPS) for (const ch of AA_GROUPS[g]) AA_COLOR[ch] = c[g]; })();
  const ALN = { start: 0, cols: 120, colorOn: true, sortByTree: true };

  function renderAlignment() {
    const host = $("#alnCanvas");
    if (!host) return;
    if (!D.alignment || !D.alignment.ids || !D.alignment.ids.length) { host.innerHTML = '<p class="empty">No alignment file available.</p>'; return; }
    let ids = D.alignment.ids.slice(), seqs = D.alignment.seqs.slice();
    // order by current tree leaf order if available
    if (ALN.sortByTree && T.root) {
      const order = collectLeaves(T.root, []).map(l => l.name);
      const idx = new Map(order.map((n, i) => [n, i]));
      const pack = ids.map((id, i) => ({ id, seq: seqs[i] }));
      pack.sort((a, b) => (idx.has(a.id) ? idx.get(a.id) : 1e9) - (idx.has(b.id) ? idx.get(b.id) : 1e9));
      ids = pack.map(p => p.id); seqs = pack.map(p => p.seq);
    }
    const alnLen = Math.max(...seqs.map(s => s.length));
    ALN.start = Math.min(ALN.start, Math.max(0, alnLen - ALN.cols));
    const start = ALN.start, end = Math.min(alnLen, start + ALN.cols);
    const cw = 9, rh = 16, labelW = 180;

    const rows = [];
    // Simple range header — the per-column tick numbers never lined up cleanly
    // across proportional layouts, so we just state the visible column range.
    const rangeHdr = `<div class="aln-row aln-range"><span class="aln-label">columns</span>` +
      `<span class="aln-seq">${start + 1}\u2013${end} of ${alnLen}</span></div>`;
    rows.push(rangeHdr);

    for (let r = 0; r < ids.length; r++) {
      const seq = seqs[r];
      let cells = "";
      for (let c = start; c < end; c++) {
        const ch = (seq[c] || "-").toUpperCase();
        const col = ALN.colorOn && ch !== "-" ? AA_COLOR[ch] || "#8a8f98" : "";
        const cls = ch === "-" ? "aln-cell gap" : "aln-cell";
        cells += `<i class="${cls}" ${col ? `style="background:${col}"` : ""}>${ch}</i>`;
      }
      rows.push(`<div class="aln-row"><span class="aln-label" title="${escapeHtml(ids[r])}">${escapeHtml(ids[r])}</span><span class="aln-seq">${cells}</span></div>`);
    }
    host.innerHTML = `<div class="aln-grid" style="--cw:${cw}px;--rh:${rh}px;--labelW:${labelW}px">${rows.join("")}</div>`;
    $("#alnMeta").textContent = `${ids.length} genomes \u00b7 concatenated SCG alignment ${alnLen} aa wide \u00b7 viewing columns ${start + 1}\u2013${end} of ${alnLen}`;
    const sl = $("#alnSlider"); if (sl) { sl.max = Math.max(0, alnLen - ALN.cols); sl.value = start; }
  }

  /* =========================================================
   * SCG presence heatmap (Plotly)
   * =======================================================*/
  function renderHeatmap() {
    const host = $("#scgCanvas");
    if (!host) return;
    if (!D.hits || !D.hits.matrix || !D.hits.matrix.length) { host.innerHTML = '<p class="empty">No SCG_hit_counts.tsv available.</p>'; return; }
    if (host.dataset.rendered) return;
    const dark = rootEl.getAttribute("data-theme") === "dark";
    const z = D.hits.matrix, x = D.hits.genes, y = D.hits.genomes;
    // Hit counts are integers: 0 = absent, 1 = single copy (ideal for an SCG),
    // 2+ = multiple copies (possible redundancy/paralog). Use a fixed discrete
    // scale anchored at 0..2 so the legend is meaningful even when every value
    // is 1 (a single-value matrix would otherwise auto-scale to 0.6..1.4).
    const zmax = Math.max(2, ...z.map(row => Math.max(...row)));
    const discrete = [
      [0.0, dark ? "#3a1220" : "#f7d6d6"],   // 0 absent  (red-ish)
      [0.999 / zmax, dark ? "#3a1220" : "#f7d6d6"],
      [1 / zmax, dark ? "#164e63" : "#0891b2"], // 1 present (cyan)
      [1.999 / zmax, dark ? "#164e63" : "#0891b2"],
      [2 / zmax, dark ? "#7c3aed" : "#6d28d9"],  // 2+ multi-copy (purple)
      [1.0, dark ? "#7c3aed" : "#6d28d9"],
    ];
    const data = [{
      z, x, y, type: "heatmap", colorscale: discrete, zmin: 0, zmax,
      colorbar: {
        title: "copies", thickness: 14, len: 0.6,
        tickmode: "array", tickvals: [0, 1, 2], ticktext: ["0 (absent)", "1 (present)", "2+ (multi)"],
      },
      hovertemplate: "<b>%{y}</b><br>%{x}: %{z} copy(ies)<extra></extra>",
      xgap: 1, ygap: 1,
    }];
    const layout = {
      margin: { l: 200, r: 30, t: 10, b: 130 },
      height: Math.max(360, 20 * y.length + 160),
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { color: dark ? "#e6edf3" : "#0f1e2e", family: "Inter, system-ui, sans-serif", size: 11 },
      xaxis: { tickangle: -55, automargin: true, side: "bottom" },
      yaxis: { automargin: true },
    };
    Plotly.newPlot(host, data, layout, { displayModeBar: true, responsive: true, displaylogo: false });
    host.dataset.rendered = "1";
  }

  /* =========================================================
   * Genome QC table
   * =======================================================*/
  const TBL = { sortCol: null, asc: true, filter: "" };
  function renderTable() {
    const host = $("#genomesCanvas");
    if (!host) return;
    if (!D.summary || !D.summary.rows || !D.summary.rows.length) { host.innerHTML = '<p class="empty">No Genomes_summary_info.tsv available.</p>'; return; }
    const cols = D.summary.columns;
    let rows = D.summary.rows.slice();
    const f = TBL.filter.toLowerCase();
    if (f) rows = rows.filter(r => r.some(c => String(c).toLowerCase().includes(f)));
    if (TBL.sortCol != null) {
      const ci = TBL.sortCol;
      rows.sort((a, b) => { const av = a[ci], bv = b[ci]; const an = parseFloat(av), bn = parseFloat(bv); const both = !isNaN(an) && !isNaN(bn); const r = both ? an - bn : String(av).localeCompare(String(bv)); return TBL.asc ? r : -r; });
    }
    let html = '<div class="tbl-wrap"><table class="qc"><thead><tr>';
    cols.forEach((c, i) => { html += `<th data-ci="${i}" class="${TBL.sortCol === i ? (TBL.asc ? "sort-asc" : "sort-desc") : ""}">${escapeHtml(c)}</th>`; });
    html += "</tr></thead><tbody>";
    rows.forEach(r => { html += "<tr>" + r.map(c => `<td>${escapeHtml(String(c))}</td>`).join("") + "</tr>"; });
    html += "</tbody></table></div>";
    host.innerHTML = html;
    $$("#genomesCanvas th").forEach(th => th.addEventListener("click", () => { const ci = +th.dataset.ci; if (TBL.sortCol === ci) TBL.asc = !TBL.asc; else { TBL.sortCol = ci; TBL.asc = true; } renderTable(); }));
    $("#genomesMeta").textContent = `${rows.length} genomes \u00b7 ${cols.length} columns`;
  }

  /* =========================================================
   * Taxonomy sunburst (dependency-free SVG)
   * =======================================================*/
  function buildTaxonomyTree() {
    // Look for lineage columns in the summary table.
    if (!D.summary) return null;
    const cols = D.summary.columns.map(c => c.toLowerCase());
    const rankNames = ["domain", "phylum", "class", "order", "family", "genus", "species"];
    const rankIdx = rankNames.map(r => cols.indexOf(r)).filter(i => i >= 0);
    if (rankIdx.length < 2) return null;
    const root = { name: "root", children: {}, count: 0 };
    D.summary.rows.forEach(row => {
      let node = root; node.count++;
      for (const ci of rankIdx) {
        let v = (row[ci] || "").trim(); if (!v || v.toUpperCase() === "NA") v = "unclassified";
        node.children[v] = node.children[v] || { name: v, children: {}, count: 0 };
        node = node.children[v]; node.count++;
      }
    });
    return root;
  }
  function renderSunburst() {
    const host = $("#taxCanvas");
    if (!host) return;
    const tree = buildTaxonomyTree();
    if (!tree) { host.innerHTML = '<p class="empty">No taxonomy/lineage columns in the summary (run with -t or -D to add NCBI/GTDB taxonomy).</p>'; return; }
    host.innerHTML = "";
    const size = 560, cx = size / 2, cy = size / 2, rings = maxDepthTax(tree), ringW = (size / 2 - 20) / Math.max(1, rings);
    const svg = el("svg", { width: size, height: size, viewBox: `0 0 ${size} ${size}`, class: "sunburst" });
    const palette = ["#0891b2", "#2b8a9e", "#4c9a52", "#b7791f", "#c0392b", "#7c3aed", "#3161c9", "#d6336c", "#0e7490", "#4d908e"];
    let ci = 0;
    function arc(x, y, r0, r1, a0, a1) {
      const p = (r, a) => [x + r * Math.cos(a), y + r * Math.sin(a)];
      const large = (a1 - a0) > Math.PI ? 1 : 0;
      const [x0, y0] = p(r0, a0), [x1, y1] = p(r1, a0), [x2, y2] = p(r1, a1), [x3, y3] = p(r0, a1);
      return `M${x0},${y0} L${x1},${y1} A${r1},${r1} 0 ${large} 1 ${x2},${y2} L${x3},${y3} A${r0},${r0} 0 ${large} 0 ${x0},${y0} Z`;
    }
    function draw(node, depth, a0, a1, color) {
      const kids = Object.values(node.children);
      if (depth > 0) {
        const r0 = (depth - 1) * ringW + 14, r1 = depth * ringW + 14;
        const path = el("path", { d: arc(cx, cy, r0, r1, a0, a1), fill: color, class: "arc" });
        const title = el("title"); title.textContent = `${node.name} — ${node.count} genome(s)`; path.appendChild(title);
        path.addEventListener("mouseenter", () => { $("#taxLabel").textContent = `${node.name} · ${node.count} genome(s)`; });
        svg.appendChild(path);
        // label if wide enough
        if ((a1 - a0) > 0.16) {
          const mid = (a0 + a1) / 2, rr = (r0 + r1) / 2;
          const lx = cx + rr * Math.cos(mid), ly = cy + rr * Math.sin(mid);
          const deg = mid * 180 / Math.PI;
          const t = el("text", { x: lx, y: ly, class: "arc-label", "text-anchor": "middle", transform: `rotate(${deg > 90 && deg < 270 ? deg + 180 : deg},${lx},${ly})` });
          t.textContent = node.name.length > 14 ? node.name.slice(0, 13) + "\u2026" : node.name; svg.appendChild(t);
        }
      }
      let a = a0; const total = node.count || 1;
      kids.forEach(k => { const span = (a1 - a0) * (k.count / total); const col = depth === 0 ? palette[ci++ % palette.length] : shade(color, depth % 2 ? 0.14 : -0.1); draw(k, depth + 1, a, a + span, col); a += span; });
    }
    draw(tree, 0, 0, 2 * Math.PI, "#0891b2");
    host.appendChild(svg);
    const lbl = document.createElement("div"); lbl.id = "taxLabel"; lbl.className = "tax-label"; lbl.textContent = `${D.summary.rows.length} genomes across ${rings} ranks`; host.appendChild(lbl);
  }
  function maxDepthTax(node, d) { d = d || 0; const kids = Object.values(node.children); if (!kids.length) return d; return Math.max(...kids.map(k => maxDepthTax(k, d + 1))); }
  function shade(hex, amt) { const c = hex.replace("#", ""); const n = parseInt(c, 16); const cl = v => Math.max(0, Math.min(255, v)); let r = cl((n >> 16) + Math.round(255 * amt)), g = cl(((n >> 8) & 255) + Math.round(255 * amt)), b = cl((n & 255) + Math.round(255 * amt)); return "#" + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1); }

  function escapeHtml(s) { return String(s).replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c])); }

  /* ---------------- wire up ---------------- */
  document.addEventListener("DOMContentLoaded", function () {
    // parse tree once
    T.root = parseNewick(D.newick);
    // KPIs
    if (T.root) { const st = treeStats(T.root); setText("kpiTips", st.tips); setText("kpiInternal", st.internal); }
    else { setText("kpiTips", "\u2014"); setText("kpiInternal", "\u2014"); }
    setText("kpiAln", D.alignment ? (D.alignment.ids.length + " \u00d7 " + Math.max(...D.alignment.seqs.map(s => s.length)) + " aa") : "\u2014");
    setText("kpiGeneSet", D.gene_set || "\u2014");

    // nav
    $$(".nav-link").forEach(a => a.addEventListener("click", e => { e.preventDefault(); showTab(a.dataset.tab); }));
    const tb = $("#themeBtn"); if (tb) tb.addEventListener("click", () => setTheme(rootEl.getAttribute("data-theme") === "dark" ? "light" : "dark"));

    // tree controls
    bind("#btnPhylo", () => { T.layout = "phylogram"; T.radial = false; T.forcePhylo = true; syncTreeBtns(); renderTree(); });
    bind("#btnClado", () => { T.layout = "cladogram"; T.radial = false; T.forcePhylo = false; syncTreeBtns(); renderTree(); });
    bind("#btnRadial", () => { T.radial = true; T.forcePhylo = false; syncTreeBtns(); renderTree(); });
    bind("#btnZoomIn", () => { T.zoom = Math.min(4, T.zoom * 1.2); renderTree(); });
    bind("#btnZoomOut", () => { T.zoom = Math.max(0.3, T.zoom / 1.2); renderTree(); });
    bind("#btnSupport", () => { T.showSupport = !T.showSupport; syncTreeBtns(); renderTree(); });
    bind("#btnResetRoot", () => { T.root = parseNewick(D.newick); T.collapsed.clear(); renderTree(); });
    bind("#btnNewick", downloadNewick);
    const hi = $("#treeSearch"); if (hi) hi.addEventListener("input", () => { T.highlight = hi.value; renderTree(); });

    // alignment controls
    const sl = $("#alnSlider"); if (sl) sl.addEventListener("input", () => { ALN.start = +sl.value; renderAlignment(); });
    bind("#alnColor", () => { ALN.colorOn = !ALN.colorOn; renderAlignment(); });
    bind("#alnSort", () => { ALN.sortByTree = !ALN.sortByTree; renderAlignment(); });
    bind("#alnPrev", () => { ALN.start = Math.max(0, ALN.start - ALN.cols); renderAlignment(); });
    bind("#alnNext", () => { ALN.start = ALN.start + ALN.cols; renderAlignment(); });
    bind("#btnAlnDownload", downloadAlignment);
    bind("#btnSummaryDownload", downloadSummary);

    // table filter
    const tf = $("#tblFilter"); if (tf) tf.addEventListener("input", () => { TBL.filter = tf.value; renderTable(); });

    // re-render heatmap on theme change
    const mo = new MutationObserver(() => { const h = $("#scgCanvas"); if (h && h.dataset.rendered && $("#tab-scg").classList.contains("active")) { delete h.dataset.rendered; renderHeatmap(); } });
    mo.observe(rootEl, { attributes: true, attributeFilter: ["data-theme"] });

    // hide tabs with no data
    if (!D.hits) hideTab("scg");
    if (!D.summary) { hideTab("genomes"); hideTab("taxonomy"); }
    else if (!buildTaxonomyTree()) hideTab("taxonomy");
    if (!D.alignment) hideTab("alignment");

    syncTreeBtns();
    showTab(T.root ? "tree" : (D.alignment ? "alignment" : "genomes"));
  });

  function setText(id, v) { const e = document.getElementById(id); if (e) e.textContent = v; }
  function bind(sel, fn) { const e = $(sel); if (e) e.addEventListener("click", fn); }
  function hideTab(id) { const a = $(`.nav-link[data-tab="${id}"]`); if (a) a.style.display = "none"; }
  function syncTreeBtns() {
    // Reflect the EFFECTIVE view: an auto-cladogram (from a low-res tree) lights
    // the Cladogram button even though T.layout is still "phylogram".
    setActive("#btnPhylo", !T.radial && !T.effectiveCladogram);
    setActive("#btnClado", !T.radial && T.effectiveCladogram);
    setActive("#btnRadial", T.radial);
    setActive("#btnSupport", T.showSupport);
  }
  function setActive(sel, on) { const e = $(sel); if (e) e.classList.toggle("active", !!on); }
})();
