// Shared tool engine for spt.html and cpt.html.
// Each page calls initTool({ engine, soilColumn, sbtColumn }).

const $ = s => document.querySelector(s);
const fmt = (v, d = 1) => v.toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d });
const isClay = p => p.c > 0 && p.f === 0;

// Robertson (2010) non-normalized SBT zones — clay(brown) -> sand(yellow) ramp
const SBT = {
  0: { label: 'Unclassified',          color: '#9AA3AB' },
  2: { label: 'Organic soils',         color: '#4E3B2A' },
  3: { label: 'Clay',                  color: '#8A5A3B' },
  4: { label: 'Silt mixture',          color: '#B08E5A' },
  5: { label: 'Sand mixture',          color: '#C9A14A' },
  6: { label: 'Sand',                  color: '#DCC079' },
  7: { label: 'Gravelly / dense sand', color: '#E3D3A0' },
};
const sbtOf = p => SBT[p.sbt] || SBT[0];

export async function initTool(cfg) {
  const statusEl = $('#status');
  const runBtn = $('#run');
  let engine = null;

  const setStatus = (kind, msg) => { statusEl.className = 'status ' + kind; statusEl.textContent = msg; };

  // load the WebAssembly engine
  try {
    const mod = await import(cfg.engine);
    engine = await mod.default();
    setStatus('ok', 'Engine ready — run the analysis.');
    runBtn.disabled = false;
  } catch (e) {
    runBtn.disabled = true;
    const file = cfg.engine.replace('./', '');
    setStatus('err', `Engine (${file} / .wasm) not found. Build it with emcc first, then open this page from a local server.`);
    console.error(e);
  }

  // CSV file loader
  $('#pick').addEventListener('click', () => $('#file').click());
  $('#file').addEventListener('change', async ev => {
    const f = ev.target.files[0];
    if (f) $('#csv').value = await f.text();
  });

  // run
  runBtn.addEventListener('click', () => {
    if (!engine) return;
    const csv = $('#csv').value;
    const layersEl = $('#layers');
    const n = layersEl ? Math.max(1, parseInt(layersEl.value || '8', 10)) : 8;
    let data;
    try {
      data = JSON.parse(engine.analyzeCSV(csv, n));
    } catch (e) {
      setStatus('err', 'An error occurred during analysis. Please check the input format.');
      console.error(e); return;
    }
    if (!data.layers.length) {
      setStatus('warn', 'No valid data rows. Please check the input format.'); return;
    }
    setStatus('ok', `Done — ${data.layers.length} layers.`);
    renderTable(data.layers, cfg);
    renderChart(data.profile, cfg);
    renderLegend(cfg, data.profile);
    $('#results').classList.add('show');
  });
}

function renderLegend(cfg, profile) {
  const el = $('#legend');
  if (cfg.soilColumn) {
    el.classList.remove('hidden');
    el.innerHTML = `<span><span class="swatch" style="background:var(--clay)"></span>Clay</span>
      <span><span class="swatch" style="background:var(--sand)"></span>Sand</span>`;
  } else if (cfg.sbtColumn) {
    el.classList.remove('hidden');
    const zones = [...new Set(profile.map(p => p.sbt))].filter(z => z > 0).sort();
    el.innerHTML = zones.map(z => `<span><span class="swatch" style="background:${SBT[z].color}"></span>${SBT[z].label}</span>`).join('');
  } else {
    el.classList.add('hidden');
  }
}

function renderTable(layers, cfg) {
  let prev = 0;
  const rows = layers.map((p, i) => {
    const top = prev; prev = p.depth;
    let soilCell = '';
    if (cfg.soilColumn) soilCell = `<td class="l ${isClay(p) ? 'soil-clay' : 'soil-sand'}">${isClay(p) ? 'Clay' : 'Sand'}</td>`;
    else if (cfg.sbtColumn) { const s = sbtOf(p); soilCell = `<td class="l" style="color:${s.color};font-weight:600">${s.label}</td>`; }
    return `<tr><td>${i + 1}</td><td class="l">${fmt(top)}–${fmt(p.depth)}</td>
      <td>${fmt(p.c)}</td><td>${fmt(p.f)}</td><td>${fmt(p.g)}</td><td>${fmt(p.E, 0)}</td>${soilCell}</tr>`;
  }).join('');
  const soilHead = (cfg.soilColumn || cfg.sbtColumn) ? '<th class="l">Soil Classification</th>' : '';
  $('#table').innerHTML =
    `<thead><tr><th>#</th><th class="l">Depth (m)</th>
     <th>c (kPa)</th><th>φ (°)</th><th>γ (kN/m³)</th><th>E (kPa)</th>${soilHead}</tr></thead><tbody>${rows}</tbody>`;
}

function renderChart(profile, cfg) {
  const maxDepth = profile[profile.length - 1].depth || 1;
  const panels = [
    { key: 'f', label: 'φ (°)',     color: '#1F6F78' },
    { key: 'c', label: 'c (kPa)',   color: '#A0562D' },
    { key: 'g', label: 'γ (kN/m³)', color: '#5B6470' },
    { key: 'E', label: 'E (kPa)',   color: '#6B5BA6' },
  ];
  const padTop = 26, padBottom = 24, padLeft = 42, bandW = 16, gap = 20, panelW = 120, plotH = 300;
  const H = padTop + plotH + padBottom, W = padLeft + bandW + 8 + panels.length * (panelW + gap);
  const yOf = d => padTop + (d / maxDepth) * plotH;
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="Depth profile of soil parameters">`;

  const tickStep = maxDepth <= 6 ? 1 : (maxDepth <= 15 ? 2 : 5);
  for (let d = 0; d <= maxDepth + 1e-9; d += tickStep) {
    const y = yOf(d);
    svg += `<line class="gridline" x1="${padLeft}" y1="${y.toFixed(1)}" x2="${W - gap}" y2="${y.toFixed(1)}"/>`;
    svg += `<text class="axis-label" x="${padLeft - 6}" y="${(y + 3).toFixed(1)}" text-anchor="end">${d}</text>`;
  }
  svg += `<text class="axis-label" x="${padLeft - 30}" y="${padTop - 10}">m</text>`;

  for (let i = 0; i < profile.length - 1; i++) {
    const p = profile[i], y1 = yOf(p.depth), y2 = yOf(profile[i + 1].depth);
    if (y2 <= y1) continue;
    let col = '#9AA3AB';
    if (cfg.soilColumn) col = isClay(p) ? '#1F6F78' : '#C9A14A';
    else if (cfg.sbtColumn) col = sbtOf(p).color;
    svg += `<rect x="${padLeft}" y="${y1.toFixed(1)}" width="${bandW}" height="${(y2 - y1 + 0.6).toFixed(1)}" fill="${col}" opacity="0.82"/>`;
  }
  svg += `<rect x="${padLeft}" y="${padTop}" width="${bandW}" height="${plotH}" fill="none" stroke="var(--line)"/>`;

  panels.forEach((pan, idx) => {
    const x0 = padLeft + bandW + 8 + idx * (panelW + gap);
    const vals = profile.map(p => p[pan.key]);
    const vmax = Math.max(...vals, 1e-6);
    const xOf = v => x0 + (v / vmax) * panelW;
    svg += `<text class="panel-title" x="${x0}" y="${padTop - 10}">${pan.label}</text>`;
    svg += `<line class="gridline" x1="${x0}" y1="${padTop}" x2="${x0}" y2="${padTop + plotH}"/>`;
    svg += `<line class="gridline" x1="${x0 + panelW}" y1="${padTop}" x2="${x0 + panelW}" y2="${padTop + plotH}"/>`;
    svg += `<text class="axis-label" x="${x0}" y="${padTop + plotH + 14}">0</text>`;
    svg += `<text class="axis-label" x="${x0 + panelW}" y="${padTop + plotH + 14}" text-anchor="end">${vmax >= 100 ? Math.round(vmax) : vmax.toFixed(0)}</text>`;
    let d = '';
    profile.forEach((p, i) => { d += (i === 0 ? 'M' : 'L') + xOf(p[pan.key]).toFixed(1) + ',' + yOf(p.depth).toFixed(1); });
    svg += `<path class="curve" d="${d}" stroke="${pan.color}"/>`;
  });

  svg += `</svg>`;
  $('#chart').innerHTML = svg;
}
