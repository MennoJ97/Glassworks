// Connects the controls, the solver and the plots.

import { MATERIALS, PRESETS, PROBLEMS, buildRequests, decodeState, defaultState, encodeState, layer, problemOf, totalLength } from './model.js';
import { FieldPlot, LineChart, seriesColors } from './plot.js';
import { Solver } from './solver.js';

const $ = (id) => document.getElementById(id);
const THZ_MICRONS = 299.792458; // c, in THz times microns.
const WALLS = { pec: 'Perfect conductor', pmc: 'Magnetic wall', periodic: 'Periodic', absorbing: 'Open (no reflection)' };
const PERIOD_MS = 3000;
const VIEW_HELP = {
  spectrum: 'Reflection, transmission and absorption across a range of wavelengths, for a plane wave from the left.',
  field: 'The field at one wavelength, from a plane wave or from a current sheet inside the structure.',
  cavity: 'The wavelengths at which the structure resonates between its two ends.',
  bands: 'The waves that can travel through an endless repeat of this structure, and the gaps where none can.',
  waveguide: 'The modes that travel along the layers at a fixed wavelength, and how fast each one propagates.',
};
const PROBLEM_TEXT = {
  driven: {
    settings: 'Light source',
    field: 'This is the field the source produces.',
    learn: ['learn/07-driven.html', 'How driven problems work'],
  },
  eigen: {
    settings: 'Mode search',
    field: 'This is the shape of the selected mode. Its overall size is arbitrary, because no source sets it.',
    learn: ['learn/06-eigenproblems.html', 'How eigenproblems work'],
  },
};

const state = decodeState(location.hash.slice(1)) ?? defaultState();
const solver = new Solver();
// Switching problem returns to the view last used there.
const lastView = { driven: 'spectrum', eigen: 'cavity' };
let selectedLayer = -1;
let selectedField = 0;
let results = null;
let visible = new Set();
let playing = !matchMedia('(prefers-reduced-motion: reduce)').matches;
let timer = null;
let linkTimer = null;

const fieldPlot = new FieldPlot($('field-plot'), {
  onSelect: (i) => { selectedLayer = i; renderLayers(); fieldPlot.selected = i; fieldPlot.draw(); },
  onInterface: moveInterface,
  onSource: (x) => {
    state.sourceX = Math.min(totalLength(state), Math.max(0, x));
    fieldPlot.source = state.sourceX;
    changed();
  },
  onSplit: (i, offset) => {
    const part = state.layers[i];
    if (offset < 0.005 || part.thickness - offset < 0.005) return;
    state.layers.splice(i, 1, { ...part, thickness: offset }, { ...part, thickness: part.thickness - offset });
    selectedLayer = i + 1;
    renderLayers();
    changed();
  },
});

const chart = new LineChart($('chart'), (hit) => {
  if (state.analysis === 'spectrum') {
    state.wavelength = hit.x;
    $('wavelength').value = hit.x.toFixed(3);
  } else {
    state.bandPick = { k: hit.x, band: hit.seriesIndex };
  }
  changed();
});

function moveInterface(index, x) {
  const [a, b] = [state.layers[index], state.layers[index + 1]];
  const start = state.layers.slice(0, index).reduce((s, l) => s + l.thickness, 0);
  const sum = a.thickness + b.thickness;
  const wanted = Math.round((x - start) * 1000) / 1000;
  a.thickness = Math.min(sum - 0.005, Math.max(0.005, wanted));
  b.thickness = +(sum - a.thickness).toFixed(6);
  renderLayers();
  changed();
}

// Layer editor -----------------------------------------------------------------

function numberInput(value, step, min, onInput, label) {
  const input = Object.assign(document.createElement('input'), { type: 'number', value, step, min });
  input.setAttribute('aria-label', label);
  input.addEventListener('input', () => {
    const v = parseFloat(input.value);
    if (Number.isFinite(v) && v >= min) onInput(v);
  });
  return input;
}

function button(text, title, onClick) {
  const b = Object.assign(document.createElement('button'), { type: 'button', textContent: text, title });
  b.setAttribute('aria-label', title);
  b.addEventListener('click', onClick);
  return b;
}

function restructure(fn) {
  fn();
  renderLayers();
  changed();
}

function renderLayers() {
  const focused = document.activeElement?.closest?.('#layers') ? document.activeElement.getAttribute('aria-label') : null;
  const body = $('layers');
  body.replaceChildren();
  state.layers.forEach((l, i) => {
    const row = document.createElement('tr');
    row.classList.toggle('selected', i === selectedLayer);
    row.addEventListener('focusin', () => { selectedLayer = i; fieldPlot.selected = i; fieldPlot.draw(); });

    const select = document.createElement('select');
    select.setAttribute('aria-label', `Layer ${i + 1} material`);
    for (const m of MATERIALS) select.add(new Option(m.name, m.name, false, m.name === l.name));
    select.add(new Option('Custom', 'Custom', false, !MATERIALS.some((m) => m.name === l.name)));
    select.addEventListener('change', () => restructure(() => {
      const m = MATERIALS.find((m) => m.name === select.value);
      Object.assign(l, m ?? { name: 'Custom' });
    }));
    const custom = (key) => (v) => { l[key] = v; l.name = 'Custom'; select.value = 'Custom'; changed(); };

    const cells = [
      select,
      numberInput(+l.thickness.toFixed(4), 0.01, 0.001, (v) => { l.thickness = v; changed(); }, `Layer ${i + 1} thickness`),
      numberInput(l.n, 0.01, 0.01, custom('n'), `Layer ${i + 1} refractive index`),
      numberInput(l.k, 0.01, 0, custom('k'), `Layer ${i + 1} absorption`),
    ];
    const actions = document.createElement('div');
    actions.className = 'row-actions';
    actions.append(
      button('↑', `Move layer ${i + 1} left`, () => restructure(() => {
        if (i > 0) state.layers.splice(i - 1, 2, state.layers[i], state.layers[i - 1]);
      })),
      button('↓', `Move layer ${i + 1} right`, () => restructure(() => {
        if (i + 1 < state.layers.length) state.layers.splice(i, 2, state.layers[i + 1], state.layers[i]);
      })),
      button('Copy', `Duplicate layer ${i + 1}`, () => restructure(() => state.layers.splice(i + 1, 0, { ...l }))),
      button('Remove', `Remove layer ${i + 1}`, () => restructure(() => {
        if (state.layers.length > 1) state.layers.splice(i, 1);
      })),
    );
    for (const c of [...cells, actions]) {
      const td = document.createElement('td');
      td.append(c);
      row.append(td);
    }
    body.append(row);
  });
  if (focused) body.querySelector(`[aria-label="${CSS.escape(focused).replace(/\\ /g, ' ')}"]`)?.focus();
  $('repeat-from').max = $('repeat-to').max = state.layers.length;
  fieldPlot.layers = state.layers;
}

$('add-layer').addEventListener('click', () => restructure(() => state.layers.push(layer('Silica glass', 0.2))));
$('repeat').addEventListener('click', () => restructure(() => {
  const from = Math.max(1, +$('repeat-from').value) - 1;
  const to = Math.min(state.layers.length, +$('repeat-to').value);
  const times = Math.min(50, Math.max(1, +$('repeat-times').value));
  if (from >= to) return;
  const block = state.layers.slice(from, to);
  const copies = Array.from({ length: times }, () => block.map((l) => ({ ...l }))).flat();
  state.layers.splice(to, 0, ...copies);
}));

// Settings ---------------------------------------------------------------------

function renderSettings() {
  const a = state.analysis;
  const problem = problemOf(a);
  const text = PROBLEM_TEXT[problem];
  lastView[problem] = a;
  document.querySelectorAll('[data-problem]').forEach((t) => {
    t.setAttribute('aria-selected', t.dataset.problem === problem);
  });
  document.querySelectorAll('[data-views]').forEach((nav) => { nav.hidden = nav.dataset.views !== problem; });
  document.querySelectorAll('[data-analysis]').forEach((t) => {
    t.setAttribute('aria-selected', t.dataset.analysis === a);
  });
  $('view-text').textContent = VIEW_HELP[a];
  $('view-learn').href = text.learn[0];
  $('view-learn').textContent = text.learn[1];
  $('settings-title').textContent = text.settings;
  $('field-note').textContent = text.field;
  const show = (id, on) => { $(id).hidden = !on; };
  show('wavelength-field', a === 'waveguide' || a === 'spectrum' || a === 'field');
  show('angle-field', a === 'spectrum' || (a === 'field' && state.source === 'plane_wave'));
  show('modes-field', a === 'cavity' || a === 'bands' || a === 'waveguide');
  show('source-field', a === 'field');
  show('range-field', a === 'spectrum');
  $('polarization').value = state.polarization;
  $('wavelength').value = state.wavelength;
  $('angle').value = state.angle;
  $('n-modes').value = state.nModes;
  $('order').value = state.order;
  $('cells').value = state.cellsPerWavelength;
  $('source').value = state.source;
  $('range-from').value = state.spectrum.from;
  $('range-to').value = state.spectrum.to;
  renderBoundaries();
}

function renderBoundaries() {
  const a = state.analysis;
  const closed = ['pec', 'pmc'];
  const options = { cavity: [...closed, 'periodic'], waveguide: closed, spectrum: ['absorbing', ...closed], field: ['absorbing', ...closed] }[a];
  const planeWave = a === 'spectrum' || (a === 'field' && state.source === 'plane_wave');
  for (const side of ['left', 'right']) {
    const select = $(`${side}-wall`);
    const fixed = !options || (side === 'left' && planeWave);
    select.disabled = fixed;
    select.replaceChildren();
    if (fixed) {
      select.add(new Option(a === 'bands' ? 'Periodic' : 'Incoming plane wave'));
    } else {
      for (const o of options) select.add(new Option(WALLS[o], o, false, state.bc[a][side] === o));
    }
  }
}

function onWall(side) {
  const value = $(`${side}-wall`).value;
  const bc = state.bc[state.analysis];
  const other = side === 'left' ? 'right' : 'left';
  bc[side] = value;
  // Periodicity ties the two ends together, so they change as one.
  if (value === 'periodic') bc[other] = 'periodic';
  else if (bc[other] === 'periodic') bc[other] = value;
  renderBoundaries();
  changed();
}
$('left-wall').addEventListener('change', () => onWall('left'));
$('right-wall').addEventListener('change', () => onWall('right'));

function bindNumber(id, apply, min = 0) {
  $(id).addEventListener('input', () => {
    const v = parseFloat($(id).value);
    if (Number.isFinite(v) && v >= min) { apply(v); changed(); }
  });
}
bindNumber('wavelength', (v) => { state.wavelength = v; }, 0.01);
bindNumber('angle', (v) => { state.angle = Math.min(89, v); });
bindNumber('n-modes', (v) => { state.nModes = Math.min(30, Math.max(1, Math.round(v))); }, 1);
bindNumber('order', (v) => { state.order = Math.min(5, Math.max(1, Math.round(v))); }, 1);
bindNumber('cells', (v) => { state.cellsPerWavelength = Math.min(200, v); }, 1);
bindNumber('range-from', (v) => { state.spectrum.from = v; }, 0.01);
bindNumber('range-to', (v) => { state.spectrum.to = v; }, 0.01);
$('polarization').addEventListener('change', () => { state.polarization = $('polarization').value; visible = new Set(); changed(); });
$('source').addEventListener('change', () => { state.source = $('source').value; renderSettings(); changed(); });

function selectView(analysis) {
  state.analysis = analysis;
  selectedField = 0;
  renderSettings();
  changed();
}
document.querySelectorAll('[data-analysis]').forEach((t) => t.addEventListener('click', () => selectView(t.dataset.analysis)));
document.querySelectorAll('[data-problem]').forEach((t) => t.addEventListener('click', () => {
  if (t.dataset.problem !== problemOf(state.analysis)) selectView(lastView[t.dataset.problem]);
}));

for (const [problem, { name: groupName }] of Object.entries(PROBLEMS)) {
  const group = Object.assign(document.createElement('optgroup'), { label: groupName });
  for (const [name, preset] of Object.entries(PRESETS)) {
    if (problemOf(preset.analysis) === problem) group.append(new Option(name, name));
  }
  $('preset').append(group);
}
$('preset').addEventListener('change', () => {
  const preset = PRESETS[$('preset').value];
  if (!preset) return;
  Object.assign(state, structuredClone(preset));
  selectedLayer = -1;
  selectedField = 0;
  renderLayers();
  renderSettings();
  changed();
});

$('play').addEventListener('click', () => {
  playing = !playing;
  $('play').textContent = playing ? 'Pause' : 'Play';
  if (playing) requestAnimationFrame(animate);
});

// Sharing ----------------------------------------------------------------------

$('copy-link').addEventListener('click', async () => {
  await navigator.clipboard.writeText(location.href);
  setStatus('Copied a link to this structure.');
});

function download(name, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  Object.assign(document.createElement('a'), { href: url, download: name }).click();
  URL.revokeObjectURL(url);
}

$('save-json').addEventListener('click', () => download('structure.json', JSON.stringify(state, null, 2), 'application/json'));
$('load-json').addEventListener('change', async (e) => {
  const file = e.target.files[0];
  if (!file) return;
  const loaded = decodeState(encodeURIComponent(await file.text()));
  if (!loaded) return setError('That file is not a saved structure. Choose a JSON file saved from this page.');
  Object.assign(state, loaded);
  renderLayers();
  renderSettings();
  changed();
});

$('save-csv').addEventListener('click', () => {
  const source = results?.[0];
  if (!source) return;
  const names = Object.keys(source.channels).filter((n) => {
    return !n.startsWith('field') || n.startsWith(`field${selectedField}/`);
  });
  const columns = names.map((n) => source.channels[n]);
  const rows = Array.from({ length: Math.max(...columns.map((c) => c.length)) }, (_, i) => columns.map((c) => c[i] ?? '').join(','));
  download('results.csv', [names.join(','), ...rows].join('\n'), 'text/csv');
});

// Solving and showing results ----------------------------------------------------

function changed() {
  fieldPlot.draw();
  clearTimeout(timer);
  timer = setTimeout(solveNow, 60);
  // Browsers limit URL updates (Safari throws past about 100 in 30 s), so write
  // the link only once editing pauses, not on every step of a drag.
  clearTimeout(linkTimer);
  linkTimer = setTimeout(() => history.replaceState(null, '', `#${encodeState(state)}`), 500);
}

async function solveNow() {
  const reply = await solver.run(buildRequests(state));
  if (!reply) return;
  if (reply.error) {
    setError(`${reply.error[0].toUpperCase()}${reply.error.slice(1)}.`);
    clearResults();
    return;
  }
  results = reply.results;
  $('error').hidden = true;
  $('save-csv').hidden = false;
  const { n_dofs: dofs } = results[results.length - 1].meta;
  setStatus(`Solved ${dofs.toLocaleString()} unknowns in ${reply.ms < 10 ? reply.ms.toFixed(1) : Math.round(reply.ms)} ms, in your browser.`);
  showResults();
}

// Round-off can leave tiny negative powers, which should read as zero.
const percent = (p) => `${Math.abs(100 * p).toFixed(1)}%`;
const setStatus = (text) => { $('status').textContent = text; };
function setError(text) {
  $('error').textContent = text;
  $('error').hidden = false;
}

// A failed solve must not leave one problem's results showing under another's.
function clearResults() {
  results = null;
  setStatus('No results for these settings.');
  $('results-title').textContent = 'Results';
  $('results').replaceChildren();
  $('empty').hidden = true;
  $('save-csv').hidden = true;
  $('chart-wrap').hidden = true;
  chart.spec = null;
  fieldPlot.field = null;
  fieldPlot.source = null;
  renderComponents(null);
  fieldPlot.draw();
}

function fieldFrom(result, i) {
  const meta = result.meta.fields[i];
  if (!meta) return null;
  const channel = (name) => result.channels[`field${i}/${name}`];
  const components = {};
  for (const name of meta.components) components[name] = { re: channel(`${name}/re`), im: channel(`${name}/im`) };
  return { x: channel('x'), components, names: meta.components, meta };
}

function showResults() {
  const a = state.analysis;
  const fieldResult = results[results.length - 1];
  const count = fieldResult.meta.fields.length;
  if (a === 'bands') selectedField = Math.min(state.bandPick.band, count - 1);
  selectedField = Math.max(0, Math.min(selectedField, count - 1));

  const field = fieldFrom(fieldResult, selectedField);
  fieldPlot.field = field;
  fieldPlot.source = a === 'field' && state.source === 'line' ? state.sourceX : null;
  if (field && !field.names.some((n) => visible.has(n))) visible = new Set([field.names[0]]);
  fieldPlot.visible = visible;
  renderComponents(field);
  renderTable(fieldResult);
  renderChart();
  fieldPlot.draw();
}

function renderComponents(field) {
  const box = $('components');
  box.replaceChildren();
  if (!field) return;
  const colors = seriesColors();
  field.names.forEach((name, n) => {
    const wrap = document.createElement('label');
    wrap.className = 'chip';
    const check = Object.assign(document.createElement('input'), { type: 'checkbox', checked: visible.has(name) });
    check.addEventListener('change', () => {
      if (check.checked) visible.add(name);
      else visible.delete(name);
      fieldPlot.draw();
    });
    const swatch = document.createElement('span');
    swatch.className = 'swatch';
    swatch.style.background = colors[n];
    const sub = document.createElement('sub');
    sub.textContent = name[1];
    const text = document.createElement('span');
    text.append(name[0] === 'H' ? 'η₀H' : 'E', sub);
    wrap.append(check, swatch, text);
    box.append(wrap);
  });
}

function resultsTitle() {
  switch (state.analysis) {
    case 'spectrum':
    case 'field': return `Response at ${+state.wavelength.toFixed(3)} µm`;
    case 'waveguide': return `Guided modes at ${+state.wavelength.toFixed(3)} µm`;
    case 'bands': return `Modes at kL/π = ${state.bandPick.k.toFixed(2)}`;
    default: return 'Resonant modes';
  }
}

function renderTable(result) {
  const a = state.analysis;
  const table = $('results');
  $('results-title').textContent = resultsTitle();
  table.replaceChildren();
  const fields = result.meta.fields;
  const addRow = (cells, tag = 'td') => {
    const row = document.createElement('tr');
    for (const c of cells) row.append(Object.assign(document.createElement(tag), { textContent: c }));
    table.append(row);
    return row;
  };
  $('empty').hidden = fields.length > 0;
  if (fields.length === 0) {
    $('empty').textContent = a === 'waveguide'
      ? 'No guided modes at this wavelength. Give an inner layer a higher index than both outer layers, or make it thicker.'
      : 'No modes found. Add a layer or widen the structure.';
    return;
  }
  if (a === 'spectrum' || a === 'field') {
    const power = fields[0].power;
    if (power) {
      addRow(['Reflected', 'Transmitted', 'Absorbed'], 'th');
      addRow(power.map(percent));
    } else {
      $('empty').hidden = false;
      $('empty').textContent = 'The source radiates in both directions. Drag it to move it.';
    }
    return;
  }
  const heading = a === 'waveguide' ? ['Mode', 'Effective index', 'Propagation β (1/µm)'] : ['Mode', 'Wavelength (µm)', 'Frequency (THz)'];
  addRow(heading, 'th');
  fields.forEach((f, i) => {
    const wavelength = f.k0 > 1e-3 ? (2 * Math.PI) / f.k0 : Infinity;
    const cells = a === 'waveguide'
      ? [i, f.n_eff.toFixed(5), (f.n_eff * f.k0).toFixed(4)]
      : [i + 1, Number.isFinite(wavelength) ? wavelength.toFixed(4) : 'Static', Number.isFinite(wavelength) ? (THZ_MICRONS / wavelength).toFixed(1) : '0'];
    const row = addRow(cells);
    row.tabIndex = 0;
    row.classList.add('pickable');
    row.classList.toggle('selected', i === selectedField);
    const pick = () => {
      selectedField = i;
      if (a === 'bands') state.bandPick.band = i;
      showResults();
    };
    row.addEventListener('click', pick);
    row.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); pick(); }
    });
  });
}

function renderChart() {
  const a = state.analysis;
  const wrap = $('chart-wrap');
  wrap.hidden = a !== 'spectrum' && a !== 'bands';
  if (wrap.hidden) { chart.spec = null; return; }
  const channels = results[0].channels;
  const colors = seriesColors();
  if (a === 'spectrum') {
    const x = channels.wavelength;
    $('chart-title').textContent = 'Click the chart to see the field at a wavelength.';
    chart.spec = {
      xMin: x[0], xMax: x[x.length - 1], yMin: 0, yMax: 1, xLabel: 'Wavelength (µm)', endLabels: true,
      marker: state.wavelength,
      series: ['R', 'T', 'A'].map((n, i) => ({
        name: ['Reflected', 'Transmitted', 'Absorbed'][i], x, y: channels[n].map((v) => Math.min(1, Math.max(0, v))), color: colors[i],
      })),
      format: (index) => [`${x[index].toFixed(3)} µm`, ...['R', 'T', 'A'].map((n, i) => `${['Reflected', 'Transmitted', 'Absorbed'][i]} ${percent(channels[n][index])}`)],
    };
  } else {
    const L = results[0].meta.length;
    const x = Array.from(channels.k, (k) => (k * L) / Math.PI);
    const scale = L / (2 * Math.PI);
    const bands = Object.keys(channels).filter((n) => n.startsWith('band')).map((n) => channels[n].map((k0) => k0 * scale));
    const gaps = [];
    for (let i = 0; i + 1 < bands.length; i++) {
      const [lo, hi] = [Math.max(...bands[i]), Math.min(...bands[i + 1])];
      if (hi - lo > 1e-3) gaps.push([lo, hi]);
    }
    $('chart-title').textContent = 'Band diagram. Shaded ranges are band gaps, where no wave can travel. Click a point to see its mode.';
    chart.spec = {
      xMin: 0, xMax: 1, yMin: 0, yMax: Math.max(...bands.flat()) * 1.05 || 1, xLabel: 'Bloch wavevector kL/π', gaps, marker: state.bandPick.k,
      series: bands.map((y, i) => ({ name: `Band ${i + 1}`, x, y, color: colors[0] })),
      format: (index, s) => [`Band ${s + 1}`, `kL/π ${x[index].toFixed(3)}`, `L/λ ${bands[s][index].toFixed(4)}`],
    };
  }
  chart.draw();
}

// Animation and start-up -----------------------------------------------------------

let lastFrame = 0;
function animate(now) {
  if (!playing) return;
  fieldPlot.phase = (fieldPlot.phase + (2 * Math.PI * Math.min(50, now - lastFrame)) / PERIOD_MS) % (2 * Math.PI);
  lastFrame = now;
  if (fieldPlot.field && !document.hidden) fieldPlot.draw();
  requestAnimationFrame(animate);
}

new ResizeObserver(() => { fieldPlot.draw(); chart.draw(); }).observe($('field-plot'));
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { if (results) showResults(); });

$('play').textContent = playing ? 'Pause' : 'Play';
renderLayers();
renderSettings();
changed();
requestAnimationFrame(animate);
