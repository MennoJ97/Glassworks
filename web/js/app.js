// Connects the controls, the solver and the plots.

import { MATERIALS, PRESETS, PROBLEMS, blochAt, buildRequests, decodeState, defaultState, encodeState, layer, problemOf, totalLength } from './model.js';
import { EvolutionMap, FieldMap, FieldPlot, LineChart, seriesColors } from './plot.js';
import { Solver } from './solver.js';
import { STAGES, TOTAL_ITERATIONS, finishedLayers, gradientRequests, meanReflection, optimize, pixelCount, plotLayers, references, sharpnessAt } from './design.js';

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
  cavity2d: 'The modes of a 2D structure inside walls, uniform along the third direction.',
  bands2d: 'The waves that can travel through an endless square or rectangular lattice of this cell, along the edge of its symmetry zone.',
  mirror: 'Find the layers that reflect the most light at one wavelength, within a fixed thickness.',
  antireflection: 'Find the coating that reflects the least light, averaged across a band of wavelengths.',
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
  design: {
    settings: 'Goal',
    field: 'This is the field at the goal wavelength. The design starts as a grey mix of both materials and sharpens into layers.',
  },
};

const state = decodeState(location.hash.slice(1)) ?? defaultState();
const solver = new Solver();
// Switching problem returns to the view last used there.
const lastView = { driven: 'spectrum', eigen: 'cavity', design: 'mirror' };
// The design in progress. Changing any setting starts it over.
const design = { rho: null, frames: [], history: [], shown: 0, references: [], finished: false, running: false, run: 0 };
let selectedLayer = -1;
let selectedField = 0;
let results = null;
// Results outlive a change of view until the next solve, and fit no other view.
let resultsView = null;
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

const fieldMap = new FieldMap($('field-map'));
const planar = () => Boolean(state.plane[state.analysis]);

const chart = new LineChart($('chart'), (hit) => {
  if (problemOf(state.analysis) === 'design') { showFrame(hit.index); return; }
  if (state.analysis === 'bands2d') {
    state.plane.bands2d.pick = { s: hit.x, band: hit.seriesIndex };
    changed();
    return;
  }
  if (state.analysis === 'spectrum') {
    state.wavelength = hit.x;
    $('wavelength').value = hit.x.toFixed(3);
  } else {
    state.bandPick = { k: hit.x, band: hit.seriesIndex };
  }
  changed();
});

const evolution = new EvolutionMap($('evolution'), (frame) => showFrame(frame));

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
  $('view-learn').hidden = !text.learn;
  if (text.learn) [$('view-learn').href, $('view-learn').textContent] = text.learn;
  $('settings-title').textContent = text.settings;
  $('field-note').textContent = text.field;
  const show = (id, on) => { $(id).hidden = !on; };
  const designing = problem === 'design';
  show('wavelength-field', a === 'waveguide' || a === 'spectrum' || a === 'field' || a === 'mirror');
  show('angle-field', a === 'spectrum' || (a === 'field' && state.source === 'plane_wave') || designing);
  show('design-fields', designing);
  show('design-band', a === 'antireflection');
  show('design-run', designing);
  show('left-field', !designing && !planar());
  show('right-field', !designing && !planar());
  show('accuracy', !designing && !planar());
  const flat = planar();
  show('layer-builder', !designing && !flat);
  show('plane-builder', flat);
  show('plane-fields', flat);
  show('field-plot', !flat);
  show('field-map', flat);
  if (flat) {
    renderPlane();
    show('plane-walls-field', a === 'cavity2d');
    $('field-note').textContent = 'Colour shows the field at this moment, with red and blue as opposite signs. Its overall size is arbitrary.';
  }
  // TE and TM keep their meaning for the coefficients, but read differently in a plane.
  $('polarization').options[0].text = flat ? 'Ez (electric field out of the plane, often called TM)' : 'TE (electric field along the layers)';
  $('polarization').options[1].text = flat ? 'Hz (magnetic field out of the plane, often called TE)' : 'TM (magnetic field along the layers)';
  document.querySelector('main').classList.toggle('design', designing);
  $('polarization-title').textContent = designing || flat ? 'Polarization' : 'Polarization and ends';
  $('field-edit').hidden = designing || flat;
  $('field-layers').hidden = flat;
  fieldPlot.locked = designing;
  if (designing) renderDesignSettings();
  else fieldPlot.layers = state.layers;
  show('modes-field', a === 'cavity' || a === 'bands' || a === 'waveguide' || flat);
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

// The 2D domain and its shapes.
function renderPlane() {
  const p = state.plane[state.analysis];
  $('plane-width').value = p.width;
  $('plane-height').value = p.height;
  $('plane-cell').value = p.cell;
  $('plane-background').value = p.background;
  $('plane-walls').replaceChildren(...['pec', 'pmc', 'periodic'].map((w) => new Option(WALLS[w], w, false, w === p.walls)));
  const body = $('shapes');
  const focused = document.activeElement?.closest?.('#shapes') ? document.activeElement.getAttribute('aria-label') : null;
  body.replaceChildren();
  p.shapes.forEach((shape, i) => {
    const row = document.createElement('tr');
    const kind = document.createElement('select');
    kind.setAttribute('aria-label', `Shape ${i + 1} type`);
    kind.add(new Option('Circle', 'circle', false, shape.type === 'circle'));
    kind.add(new Option('Rectangle', 'rect', false, shape.type === 'rect'));
    kind.addEventListener('change', () => { shape.type = kind.value; renderPlane(); changed(); });
    const edit = (key, min, label) => numberInput(shape[key], 0.01, min, (v) => { shape[key] = v; changed(); }, `Shape ${i + 1} ${label}`);
    const height = edit('h', 0.001, 'height');
    height.disabled = shape.type === 'circle';
    const remove = button('Remove', `Remove shape ${i + 1}`, () => { p.shapes.splice(i, 1); renderPlane(); changed(); });
    for (const c of [kind, edit('x', -100, 'centre x'), edit('y', -100, 'centre y'), edit('w', 0.001, 'width or diameter'), height, edit('n', 0.01, 'index'), remove]) {
      const td = document.createElement('td');
      td.append(c);
      row.append(td);
    }
    body.append(row);
  });
  if (focused) body.querySelector(`[aria-label="${CSS.escape(focused).replace(/\\ /g, ' ')}"]`)?.focus();
}

const planeSettings = () => state.plane[state.analysis];
bindNumber('plane-width', (v) => { planeSettings().width = v; }, 0.05);
bindNumber('plane-height', (v) => { planeSettings().height = v; }, 0.05);
bindNumber('plane-cell', (v) => { planeSettings().cell = v; }, 0.001);
bindNumber('plane-background', (v) => { planeSettings().background = v; }, 0.01);
$('plane-walls').addEventListener('change', () => { planeSettings().walls = $('plane-walls').value; changed(); });
$('add-shape').addEventListener('click', () => {
  const p = planeSettings();
  p.shapes.push({ type: 'circle', x: p.width / 2, y: p.height / 2, w: Math.min(p.width, p.height) / 4, h: Math.min(p.width, p.height) / 4, n: 2 });
  renderPlane();
  changed();
});

function renderDesignSettings() {
  const d = state.design[state.analysis];
  for (const key of ['thickness', 'radius', 'seed']) $(`design-${key}`).value = d[key];
  for (const key of ['low', 'high', 'left', 'right']) {
    const select = $(`design-${key}`);
    select.replaceChildren(...MATERIALS.filter((m) => m.k === 0).map((m) => new Option(m.name, m.name, false, m.name === d[key])));
  }
  if (d.from !== undefined) [$('design-from').value, $('design-to').value] = [d.from, d.to];
}

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
const designSettings = () => state.design[state.analysis];
bindNumber('design-thickness', (v) => { designSettings().thickness = Math.min(5, v); }, 0.05);
bindNumber('design-radius', (v) => { designSettings().radius = Math.min(100, v); }, 5);
bindNumber('design-seed', (v) => { designSettings().seed = Math.round(v); }, 1);
bindNumber('design-from', (v) => { designSettings().from = v; }, 0.05);
bindNumber('design-to', (v) => { designSettings().to = v; }, 0.05);
for (const key of ['low', 'high', 'left', 'right']) {
  $(`design-${key}`).addEventListener('change', () => { designSettings()[key] = $(`design-${key}`).value; changed(); });
}
$('polarization').addEventListener('change', () => { state.polarization = $('polarization').value; visible = new Set(); changed(); });
$('source').addEventListener('change', () => { state.source = $('source').value; renderSettings(); changed(); });

// Nothing on screen fits a new view until its first solve returns.
function clearForNewView() {
  clearResults();
  $('error').hidden = true;
  setStatus('Solving.');
}

function selectView(analysis) {
  state.analysis = analysis;
  clearForNewView();
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
  clearForNewView();
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
  resetDesign();
  fieldPlot.draw();
  clearTimeout(timer);
  timer = setTimeout(solveNow, 60);
  // Browsers limit URL updates (Safari throws past about 100 in 30 s), so write
  // the link only once editing pauses, not on every step of a drag.
  clearTimeout(linkTimer);
  linkTimer = setTimeout(() => history.replaceState(null, '', `#${encodeState(state)}`), 500);
}

async function solveNow() {
  if (problemOf(state.analysis) === 'design') return previewDesign();
  const view = state.analysis;
  const reply = await solver.run(buildRequests(state));
  if (!reply || view !== state.analysis) return;
  if (reply.error) {
    setError(`${reply.error[0].toUpperCase()}${reply.error.slice(1)}.`);
    clearResults();
    return;
  }
  results = reply.results;
  resultsView = view;
  $('error').hidden = true;
  $('save-csv').hidden = false;
  const { n_dofs: dofs } = results[results.length - 1].meta;
  setStatus(`Solved ${dofs.toLocaleString()} unknowns in ${reply.ms < 10 ? reply.ms.toFixed(1) : Math.round(reply.ms)} ms, in your browser.`);
  showResults();
}

// Inverse design -----------------------------------------------------------------

function resetDesign() {
  design.run++;
  design.running = false;
  design.finished = false;
  design.history = [];
  design.frames = [];
  design.rho = null;
  $('design-replay').textContent = 'Replay';
  $('design-start').textContent = 'Start';
  $('design-use').hidden = true;
}

// Repaints the map from the stored frames, as when the colour scheme changes.
function restartEvolution() {
  const d = state.design[state.analysis];
  evolution.reset(pixelCount(d), TOTAL_ITERATIONS + 1, d.thickness, STAGES);
  design.frames.forEach((rho) => evolution.add(rho));
}

function showDesign(rho, reply) {
  design.rho = rho;
  results = reply.results;
  resultsView = state.analysis;
  $('save-csv').hidden = false;
  $('error').hidden = true;
  fieldPlot.layers = plotLayers(state.design[state.analysis], rho);
  showResults();
}

// Solves the hand-made references and the current design, without optimizing.
async function previewDesign() {
  const run = design.run;
  const found = [];
  for (const { name, requests } of references(state)) {
    const reply = await solver.run(requests);
    if (!reply || run !== design.run) return;
    if (reply.error) return setError(reply.error);
    found.push({ name, value: meanReflection(reply.results) });
  }
  design.references = found;
  const rho = design.rho ?? new Float64Array(pixelCount(state.design[state.analysis])).fill(0.5);
  const reply = await solver.run(gradientRequests(state, rho));
  if (!reply || run !== design.run) return;
  if (reply.error) return setError(reply.error);
  setStatus('Press Start to optimize this design.');
  showDesign(rho, reply);
}

async function startDesign() {
  const run = ++design.run;
  design.running = true;
  design.finished = false;
  design.history = [];
  design.frames = [];
  restartEvolution();
  $('design-start').textContent = 'Stop';
  $('design-replay').textContent = 'Replay';
  $('design-use').hidden = true;
  const started = performance.now();
  let solving = 0;
  const failure = await optimize(state, (requests) => solver.run(requests), ({ rho, reply, iteration, beta }) => {
    if (run !== design.run) return false;
    solving += reply.ms;
    design.history.push(meanReflection(reply.results));
    design.frames.push(rho);
    design.shown = design.frames.length - 1;
    evolution.add(rho);
    design.finished = beta === Infinity;
    const seconds = ((performance.now() - started) / 1000).toFixed(1);
    setStatus(design.finished
      ? `Finished in ${seconds} s, of which the solver took ${(solving / 1000).toFixed(1)} s.`
      : `Iteration ${iteration} of ${TOTAL_ITERATIONS}, sharpness ${beta}.`);
    showDesign(rho, reply);
    return true;
  });
  if (run !== design.run) return;
  design.running = false;
  $('design-start').textContent = 'Start';
  $('design-use').hidden = !design.finished;
  if (failure?.error) setError(failure.error);
  renderChart();
}

// Shows one stored iteration: its design, field and place on the chart and map.
async function solveFrame(frame, run) {
  const rho = design.frames[frame];
  const reply = await solver.run(gradientRequests(state, rho));
  if (!reply || reply.error || run !== design.run) return;
  design.shown = frame;
  const sharp = sharpnessAt(frame);
  setStatus(`Showing iteration ${frame + 1}, ${sharp === Infinity ? 'fully sharpened' : `sharpness ${sharp}`}.`);
  showDesign(rho, reply);
}

function showFrame(frame) {
  if (design.running || !design.frames[frame]) return;
  $('design-replay').textContent = 'Replay';
  return solveFrame(frame, ++design.run);
}

$('design-frame').addEventListener('input', () => showFrame(+$('design-frame').value));

$('design-replay').addEventListener('click', async () => {
  const run = ++design.run;
  if ($('design-replay').textContent === 'Stop replay') {
    $('design-replay').textContent = 'Replay';
    return;
  }
  $('design-replay').textContent = 'Stop replay';
  for (let frame = 0; frame < design.frames.length && run === design.run; frame++) {
    await solveFrame(frame, run);
    await new Promise((resolve) => setTimeout(resolve, 30));
  }
  if (run === design.run) $('design-replay').textContent = 'Replay';
});

$('design-start').addEventListener('click', () => {
  if (!design.running) return startDesign();
  design.run++;
  design.running = false;
  $('design-start').textContent = 'Start';
  setStatus('Stopped. Press Start to begin again from the random start.');
});

$('design-use').addEventListener('click', () => {
  const d = state.design[state.analysis];
  state.layers = finishedLayers(d, design.frames.at(-1));
  state.spectrum = state.analysis === 'mirror'
    ? { from: +(0.5 * state.wavelength).toFixed(3), to: +(1.8 * state.wavelength).toFixed(3) }
    : { from: d.from, to: d.to };
  state.wavelength = state.analysis === 'mirror' ? state.wavelength : (d.from + d.to) / 2;
  state.bc.spectrum.right = 'absorbing';
  selectedLayer = -1;
  renderLayers();
  selectView('spectrum');
});

function renderDesignTable(table, addRow) {
  addRow(['Structure', 'Reflected'], 'th');
  const last = design.finished && design.shown === design.frames.length - 1;
  const count = last && finishedLayers(state.design[state.analysis], design.rho).length - 2;
  const name = count ? `This design, ${count} layers` : design.frames.length ? `Iteration ${design.shown + 1}` : 'This design';
  addRow([name, percent(meanReflection(results))]);
  for (const { name, value: reference } of design.references) addRow([name, percent(reference)]);
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
  $('evolution-wrap').hidden = true;
  chart.spec = null;
  fieldPlot.field = null;
  fieldPlot.source = null;
  fieldMap.clear();
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
  if (!results || resultsView !== a) return;
  const fieldResult = results[results.length - 1];
  const count = fieldResult.meta.fields.length;
  if (a === 'bands') selectedField = Math.min(state.bandPick.band, count - 1);
  if (a === 'bands2d') selectedField = Math.min(state.plane.bands2d.pick.band, count - 1);
  selectedField = Math.max(0, Math.min(selectedField, count - 1));

  if (planar()) {
    const p = state.plane[a];
    const { channels } = fieldResult;
    fieldMap.setMesh({ x: channels.x, y: channels.y, triangles: channels.triangles, width: p.width, height: p.height }, p.shapes);
    const name = fieldResult.meta.fields[selectedField]?.components[0];
    if (name) fieldMap.setField(channels[`field${selectedField}/${name}/re`], channels[`field${selectedField}/${name}/im`]);
    fieldPlot.field = null;
    $('components').replaceChildren(name ?? '');
    renderTable(fieldResult);
    renderChart();
    fieldMap.phase = fieldPlot.phase;
    fieldMap.draw();
    return;
  }
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
    case 'bands2d': {
      const [x, y] = blochAt(state.plane.bands2d.pick.s).map((v) => (v / Math.PI).toFixed(2));
      return `Modes at k = (${x}, ${y}) π/a`;
    }
    case 'mirror': return `Reflection at ${+state.wavelength.toFixed(3)} µm`;
    case 'antireflection': return `Mean reflection, ${designSettings().from} to ${designSettings().to} µm`;
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
  if (problemOf(a) === 'design') return renderDesignTable(table, addRow);
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
      if (a === 'bands2d') state.plane.bands2d.pick.band = i;
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
  const designing = problemOf(a) === 'design';
  wrap.hidden = designing ? design.history.length < 2 : a !== 'spectrum' && a !== 'bands' && a !== 'bands2d';
  $('evolution-wrap').hidden = !designing || design.frames.length < 2;
  if (wrap.hidden) { chart.spec = null; return; }
  const channels = results[0].channels;
  const colors = seriesColors();
  if (designing) {
    // Both goals become light lost, on a log scale, where a mirror's last 1% shows.
    const lost = (r) => Math.log10(Math.max(1e-5, a === 'mirror' ? 1 - r : r));
    const x = design.history.map((_, i) => i + 1);
    const y = design.history.map(lost);
    const levels = design.references.map(({ name, value }) => ({ name, y: lost(value) }));
    const yMin = Math.floor(Math.min(...y, ...levels.map((l) => l.y)));
    const yMax = Math.max(Math.ceil(Math.max(...y)), yMin + 1);
    $('chart-title').textContent = a === 'mirror'
      ? 'Light that gets through the mirror at each iteration, lower is better. It jumps where the pixels get sharper.'
      : 'Light reflected by the coating at each iteration, lower is better. It jumps where the pixels get sharper.';
    chart.spec = {
      xMin: 1, xMax: TOTAL_ITERATIONS + 1, yMin, yMax, xLabel: 'Iteration', levels, marker: design.shown + 1,
      yTicks: Array.from({ length: yMax - yMin + 1 }, (_, i) => yMin + i),
      yFormat: (v) => `${+(100 * 10 ** v).toPrecision(1)}%`,
      series: [{ name: 'Lost', x, y, color: colors[0] }],
      format: (index) => [`Iteration ${x[index]}`, `Reflected ${percent(design.history[index])}`],
    };
    evolution.marker = design.shown;
    evolution.draw();
    const slider = $('design-frame');
    slider.max = design.frames.length - 1;
    slider.value = design.shown;
    slider.disabled = $('design-replay').disabled = design.running;
  } else if (a === 'spectrum') {
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
    const x = Array.from(a === 'bands2d' ? channels.s : channels.k, (k) => (a === 'bands2d' ? k : (k * L) / Math.PI));
    const scale = L / (2 * Math.PI);
    const bands = Object.keys(channels).filter((n) => n.startsWith('band')).map((n) => channels[n].map((k0) => k0 * scale));
    const gaps = [];
    for (let i = 0; i + 1 < bands.length; i++) {
      const [lo, hi] = [Math.max(...bands[i]), Math.min(...bands[i + 1])];
      if (hi - lo > 1e-3) gaps.push([lo, hi]);
    }
    $('chart-title').textContent = 'Band diagram, with frequency as period over wavelength. Shaded ranges are band gaps, where no wave can travel. Click a point to see its mode.';
    const flat = a === 'bands2d';
    chart.spec = {
      xMin: 0, xMax: flat ? 3 : 1, yMin: 0, yMax: Math.max(...bands.flat()) * 1.05 || 1, gaps,
      xLabel: flat ? 'Bloch wavevector around the zone' : 'Bloch wavevector kL/π',
      marker: flat ? state.plane.bands2d.pick.s : state.bandPick.k,
      ...(flat ? { xTicks: [0, 1, 2, 3], xFormat: (v) => ['Γ', 'X', 'M', 'Γ'][v] } : {}),
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
  if (planar() && !document.hidden) {
    fieldMap.phase = fieldPlot.phase;
    fieldMap.draw();
  }
  requestAnimationFrame(animate);
}

const resized = new ResizeObserver(() => { fieldPlot.draw(); fieldMap.draw(); chart.draw(); evolution.draw(); });
resized.observe($('field-plot'));
resized.observe($('field-map'));
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
  if (design.frames.length) restartEvolution();
  if (results) showResults();
});

$('play').textContent = playing ? 'Pause' : 'Play';
renderLayers();
renderSettings();
changed();
requestAnimationFrame(animate);
