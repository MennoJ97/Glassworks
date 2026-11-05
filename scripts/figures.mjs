// Draws the tutorial's figures into docs/figures as SVG: `node scripts/figures.mjs`.
//
// Every number in them comes from the solver. Fields, error rates and the design run come
// from the WebAssembly build in web/pkg (run ./scripts/build-web.sh first), and meshes and
// orderings from the Rust example `figure_data`. The diagrams in chapters 3, 4 and 6 are
// drawn from their formulas.
//
// Colours are the site's own variables, with light-mode values as fallbacks, so a figure
// follows the reader's theme once build-docs.mjs inlines it, and still reads as a plain
// image on GitHub. Pictures made of pixels are greyscale masks over a coloured rectangle,
// for the same reason.

import { execFileSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { crc32, deflateSync } from 'node:zlib';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const out = join(root, 'docs', 'figures');
const at = (path) => pathToFileURL(join(root, path)).href;

const { default: init, solve, solve_2d } = await import(at('web/pkg/glassworks_wasm.js'));
await init({ module_or_path: readFileSync(join(root, 'web/pkg/glassworks_wasm_bg.wasm')) });
const { defaultState, buildRequests } = await import(at('web/js/model.js'));
const design = await import(at('web/js/design.js'));

function run(request) {
  const solution = (request.width === undefined ? solve : solve_2d)(JSON.stringify(request));
  const channels = {};
  solution.channel_names().forEach((name, i) => (channels[name] = solution.channel(i)));
  const meta = JSON.parse(solution.meta());
  solution.free();
  return { meta, channels };
}

const scratch = mkdtempSync(join(tmpdir(), 'glassworks-figures-'));
execFileSync('cargo', ['run', '--release', '-q', '-p', 'glassworks-core', '--example', 'figure_data', '--', join(scratch, 'data.json')], { cwd: root, stdio: 'inherit' });
const data = JSON.parse(readFileSync(join(scratch, 'data.json'), 'utf8'));
rmSync(scratch, { recursive: true });

// ---- Drawing ----

const W = 640;
const FONT = 'var(--font, "Helvetica Neue", Helvetica, Arial, sans-serif)';
const MATH = '"STIX Two Math", "Cambria Math", "Latin Modern Math", serif';
const colour = {
  bg: 'var(--bg, #F5F5F3)', ink: 'var(--ink, #111111)', muted: 'var(--muted, #666663)', rule: 'var(--rule, #D6D6D2)',
  grid: 'var(--grid, #E2E2DE)', glass: 'var(--glass, #0F6E78)', red: 'var(--series-1, #b5173f)', blue: 'var(--series-2, #2458d3)',
  amber: 'var(--series-3, #9a6a00)',
};
const STYLE = `<style>
.fig text { font: 13px ${FONT}; fill: ${colour.muted}; }
.fig .title { fill: ${colour.ink}; font-weight: 700; }
.fig .ink-text { fill: ${colour.ink}; }
.fig .math { font-family: ${MATH}; font-style: italic; font-size: 15px; fill: ${colour.ink}; }
.fig .upright { font-style: normal; }
.fig .halo { paint-order: stroke; stroke: ${colour.bg}; stroke-width: 5px; stroke-linejoin: round; }
${Object.entries(colour).map(([name, value]) => `.fig .s-${name} { stroke: ${value}; fill: none; }\n.fig .f-${name} { fill: ${value}; }`).join('\n')}
</style>`;

const n = (v) => +v.toFixed(2);
const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;');
const text = (x, y, s, cls = '', anchor = 'start') => `<text x="${n(x)}" y="${n(y)}"${cls ? ` class="${cls}"` : ''}${anchor !== 'start' ? ` text-anchor="${anchor}"` : ''}>${s}</text>`;
const line = (x1, y1, x2, y2, cls, extra = '') => `<line x1="${n(x1)}" y1="${n(y1)}" x2="${n(x2)}" y2="${n(y2)}" class="${cls}"${extra}/>`;
const poly = (points, cls, extra = '') => `<polyline points="${points.map(([x, y]) => `${n(x)},${n(y)}`).join(' ')}" class="${cls}"${extra}/>`;
const rect = (x, y, w, h, cls, extra = '') => `<rect x="${n(x)}" y="${n(y)}" width="${n(w)}" height="${n(h)}" class="${cls}"${extra}/>`;
const circle = (x, y, r, cls, extra = '') => `<circle cx="${n(x)}" cy="${n(y)}" r="${r}" class="${cls}"${extra}/>`;
const group = (x, y, body) => `<g transform="translate(${n(x)} ${n(y)})">${body}</g>`;
const thousands = (v) => v.toLocaleString('en-GB');

function figure(name, height, body) {
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" class="fig" viewBox="0 0 ${W} ${height}" width="${W}" height="${height}">
${STYLE}
<rect width="${W}" height="${height}" class="f-bg"/>
${body}
</svg>
`;
  writeFileSync(join(out, `${name}.svg`), svg);
  console.log(`${name}.svg ${(svg.length / 1024).toFixed(0)} kB`);
}

// A greyscale PNG, as a data URI. 255 is full colour, 0 is none.
function png(width, height, grey) {
  const raw = Buffer.alloc((width + 1) * height);
  for (let y = 0; y < height; y++) raw.set(grey.subarray(y * width, (y + 1) * width), y * (width + 1) + 1);
  const chunk = (type, body) => {
    const head = Buffer.alloc(4);
    head.writeUInt32BE(body.length);
    const typed = Buffer.concat([Buffer.from(type), body]);
    const crc = Buffer.alloc(4);
    crc.writeUInt32BE(crc32(typed));
    return Buffer.concat([head, typed, crc]);
  };
  const header = Buffer.alloc(13);
  header.writeUInt32BE(width, 0);
  header.writeUInt32BE(height, 4);
  header[8] = 8; // bits per sample
  const file = Buffer.concat([Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]), chunk('IHDR', header), chunk('IDAT', deflateSync(raw, { level: 9 })), chunk('IEND', Buffer.alloc(0))]);
  return `data:image/png;base64,${file.toString('base64')}`;
}

// Paints `fill` wherever the greyscale picture is light. Ids must be unique on a page.
function masked(id, x, y, w, h, width, height, grey, fill, sharp = false) {
  const style = sharp ? ' style="image-rendering: pixelated"' : '';
  return `<mask id="${id}" maskUnits="userSpaceOnUse" x="${n(x)}" y="${n(y)}" width="${n(w)}" height="${n(h)}"><image href="${png(width, height, grey)}" x="${n(x)}" y="${n(y)}" width="${n(w)}" height="${n(h)}" preserveAspectRatio="none"${style}/></mask>${rect(x, y, w, h, `f-${fill}`, ` mask="url(#${id})"`)}`;
}

// The solver page's tint for a layer of index n: deeper for a higher index.
const tint = (index) => Math.min(0.55, (0.3 * Math.max(0, index - 1)) / 1.2);

mkdirSync(out, { recursive: true });

// ---- Chapter 1: u kinks where p jumps ----
{
  const layers = [0.5, 0.5, 0.5];
  const eps = [1, 4, 1];
  const r = run({
    layers: layers.map((thickness, i) => ({ thickness, eps: [eps[i], 0] })),
    polarization: 'tm', left: 'port', right: 'absorbing', mesh: { order: 4, max_cell: 0.05 },
    analysis: { type: 'driven', wavelength: 1, angle_deg: 0, source: { type: 'plane_wave' } },
  });
  const c = (name) => [r.channels[`field0/${name}/re`], r.channels[`field0/${name}/im`]];
  const x = r.channels['field0/x'];
  // Any phase of the wave is a valid snapshot. Show the one whose slopes at the interfaces
  // are steepest, where a corner is easiest to see. p u' = -i k0 Ez, so E_z's matching
  // part is its imaginary one.
  const [[ur, ui], [zr, zi]] = [c('Hy'), c('Ez')];
  const at = [0.5, 1].map((v) => x.findIndex((p) => Math.abs(p - v) < 1e-9));
  const part = (re, im, phase) => re.map((v, i) => Math.cos(phase) * v - Math.sin(phase) * im[i]);
  const steepness = (phase) => {
    const f = part(ur, ui, phase);
    return at.reduce((sum, i) => sum + Math.abs(f[i] - f[i - 1]) + Math.abs(f[i + 1] - f[i]), 0);
  };
  const phase = Array.from({ length: 180 }, (_, k) => (k * Math.PI) / 180).reduce((best, p) => (steepness(p) > steepness(best) ? p : best));
  const u = part(ur, ui, phase);
  const ez = part(zi, zr.map((v) => -v), phase);
  const [left, right, top, bottom] = [40, 600, 34, 200];
  const sx = (v) => left + (v / 1.5) * (right - left);
  const mid = (top + bottom) / 2;
  const scale = (values) => {
    const peak = Math.max(...values.map(Math.abs));
    return (v) => mid - (v / peak) * (bottom - top) * 0.44;
  };
  const [su, sz] = [scale(u), scale(ez)];
  let body = rect(sx(0.5), top, sx(1) - sx(0.5), bottom - top, 'f-glass', ` fill-opacity="${tint(2)}"`);
  body += line(left, mid, right, mid, 's-grid');
  body += line(left, bottom, right, bottom, 's-rule');
  for (const v of [0, 0.5, 1, 1.5]) body += line(sx(v), bottom, sx(v), bottom + 4, 's-rule') + text(sx(v), bottom + 17, v, '', 'middle');
  body += text(right, bottom + 34, 'Position (µm)', '', 'end');
  body += poly([...x].map((v, i) => [sx(v), sz(ez[i])]), 's-blue', ' stroke-width="1.6"');
  body += poly([...x].map((v, i) => [sx(v), su(u[i])]), 's-red', ' stroke-width="2"');
  for (const i of at) body += circle(sx(x[i]), su(u[i]), 7, 's-ink', ' stroke-width="1.2"');
  body += text(sx(0.75), top - 10, 'glass, ', '', 'middle').replace('</text>', `<tspan class="math">ε</tspan> = 4</text>`);
  body += text(sx(0.25), top - 10, 'air', '', 'middle') + text(sx(1.25), top - 10, 'air', '', 'middle');
  body += `<text x="${right + 6}" y="${n(su(u.at(-1)) + 4)}" class="math">u</text>`;
  body += `<text x="${right + 6}" y="${n(sz(ez.at(-1)) + 4)}" class="math">E<tspan font-size="11" dy="3">z</tspan></text>`;
  figure('interface-kink', 244, body);
}

// ---- Chapter 3: Lagrange basis functions ----
{
  const lagrange = (order, i, xi) => {
    let value = 1;
    for (let j = 0; j <= order; j++) if (j !== i) value *= (xi - (-1 + (2 * j) / order)) / ((2 * i) / order - (2 * j) / order);
    return value;
  };
  const [left, right] = [110, 610];
  const sx = (v) => left + v * (right - left);
  const row = (y0, order, cells, highlight, label) => {
    const h = 1 / cells;
    const sy = (v) => y0 + 90 - v * 72;
    let body = text(20, y0 + 56, label, 'ink-text');
    body += line(sx(0), sy(0), sx(1), sy(0), 's-rule');
    const functions = new Map();
    for (let c = 0; c < cells; c++) {
      for (let i = 0; i <= order; i++) {
        const g = c * order + i;
        const points = Array.from({ length: 61 }, (_, k) => {
          const xi = -1 + (2 * k) / 60;
          return [sx(c * h + ((xi + 1) / 2) * h), sy(lagrange(order, i, xi))];
        });
        functions.set(g, [...(functions.get(g) ?? []), points]);
      }
    }
    for (const [g, pieces] of functions) {
      if (g === highlight) continue;
      for (const p of pieces) body += poly(p, 's-muted', ' stroke-width="1" stroke-opacity="0.55"');
    }
    const lit = functions.get(highlight);
    const outline = lit.flat();
    body += `<polygon points="${[[outline[0][0], sy(0)], ...outline, [outline.at(-1)[0], sy(0)]].map(([x, y]) => `${n(x)},${n(y)}`).join(' ')}" class="f-glass" fill-opacity="0.14"/>`;
    for (const p of lit) body += poly(p, 's-glass', ' stroke-width="2.4"');
    for (let g = 0; g <= cells * order; g++) {
      const vertex = g % order === 0;
      body += circle(sx(g / (cells * order)), sy(0), vertex ? 3.5 : 3, vertex ? 'f-ink' : 's-ink', vertex ? '' : ' stroke-width="1.2" fill="none"');
    }
    for (let c = 0; c <= cells; c++) body += line(sx(c * h), sy(0) + 6, sx(c * h), sy(0) + 12, 's-muted');
    const peak = highlight / (cells * order);
    body += `<text x="${n(sx(peak) + 8)}" y="${n(sy(1) + 2)}" class="math">φ<tspan font-size="11" dy="4" class="upright">${highlight}</tspan></text>`;
    return body;
  };
  let body = row(4, 1, 4, 2, 'Order 1');
  body += row(124, 2, 2, 2, 'Order 2');
  figure('basis-functions', 232, body);
}

// ---- Chapter 4: two cells' matrices add into one ----
{
  const size = 40;
  const grid = (x0, y0, rows, cols, entries, labels) => {
    let body = '';
    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        const e = entries(i, j);
        const [x, y] = [x0 + j * size, y0 + i * size];
        if (e.split) {
          body += `<polygon points="${x},${y} ${x + size},${y} ${x},${y + size}" class="f-blue" fill-opacity="0.16"/>`;
          body += `<polygon points="${x + size},${y} ${x + size},${y + size} ${x},${y + size}" class="f-amber" fill-opacity="0.2"/>`;
        } else if (e.fill) body += rect(x, y, size, size, `f-${e.fill}`, ` fill-opacity="${e.fill === 'amber' ? 0.2 : 0.16}"`);
        body += rect(x, y, size, size, 's-rule');
        if (e.text) body += text(x + size / 2, y + size / 2 + 5, e.text, 'ink-text', 'middle');
      }
    }
    labels.forEach((l, k) => {
      body += text(x0 - 8, y0 + k * size + size / 2 + 4, l, '', 'end');
      body += text(x0 + k * size + size / 2, y0 - 8, l, '', 'middle');
    });
    return body;
  };
  // The mesh: three nodes, two cells.
  let body = line(40, 40, 260, 40, 's-ink', ' stroke-width="1.5"');
  body += line(40, 40, 150, 40, 's-blue', ' stroke-width="4" stroke-opacity="0.5"') + line(150, 40, 260, 40, 's-amber', ' stroke-width="4" stroke-opacity="0.55"');
  [40, 150, 260].forEach((x, k) => (body += circle(x, 40, 4, 'f-ink') + text(x, 62, `node ${k}`, '', 'middle')));
  body += text(95, 28, 'cell 1', '', 'middle') + text(205, 28, 'cell 2', '', 'middle');
  const local = (fill) => (i, j) => ({ fill, text: i === j ? '1' : '−1' });
  body += grid(54, 112, 2, 2, local('blue'), ['0', '1']);
  body += grid(184, 112, 2, 2, local('amber'), ['1', '2']);
  body += text(94, 214, 'cell 1', '', 'middle') + text(224, 214, 'cell 2', '', 'middle');
  body += `<path d="M300 152 H352" class="s-muted" stroke-width="1.5"/><path d="M346 146 L353 152 L346 158" class="s-muted" stroke-width="1.5"/>`;
  const global = (i, j) => {
    if (i === 1 && j === 1) return { split: true, text: '1 + 1' };
    if (Math.abs(i - j) > 1) return {};
    return { fill: i + j < 2 || (i + j === 2 && i !== j) ? 'blue' : 'amber', text: i === j ? '1' : '−1' };
  };
  body += grid(420, 92, 3, 3, global, ['0', '1', '2']);
  body += text(480, 230, 'whole line', '', 'middle');
  figure('assembly', 244, body);
}

// ---- Chapter 5: reverse Cuthill-McKee folds a ring ----
{
  const { n: size, natural, rcm, order } = data.ring;
  const newOf = Array(size);
  order.forEach((old, k) => (newOf[old] = k));
  const half = (x0, title, number, entries) => {
    const [cx, cy, r] = [x0 + 78, 118, 58];
    const at = (k) => [cx + r * Math.sin((2 * Math.PI * k) / size), cy - r * Math.cos((2 * Math.PI * k) / size)];
    let body = text(x0 + 4, 24, title, 'title');
    for (let k = 0; k < size; k++) {
      const [a, b] = [at(k), at((k + 1) % size)];
      const far = Math.abs(number(k) - number((k + 1) % size)) > 2;
      body += line(...a, ...b, far ? 's-red' : 's-rule', ` stroke-width="${far ? 2 : 1.5}"`);
    }
    for (let k = 0; k < size; k++) {
      const [x, y] = at(k);
      body += circle(x, y, 10, 'f-bg') + circle(x, y, 10, 's-muted');
      body += text(x, y + 4.5, number(k), 'ink-text', 'middle');
    }
    const cell = 8;
    const [mx, my] = [x0 + 168, 70];
    body += rect(mx, my, size * cell, size * cell, 's-rule');
    for (const [i, j] of entries) {
      const far = Math.abs(i - j) > 2;
      body += rect(mx + j * cell + 0.5, my + i * cell + 0.5, cell - 1, cell - 1, far ? 'f-red' : 'f-ink');
    }
    return body;
  };
  let body = half(0, 'Numbered along the line', (k) => k, natural);
  body += half(330, 'Reverse Cuthill–McKee', (k) => newOf[k], rcm);
  figure('ring-ordering', 196, body);
}

// ---- Chapter 6: shift and invert ----
{
  const sigma = 130;
  const lambdas = [1, 2, 3, 4, 5, 6].map((m) => (m * Math.PI) ** 2);
  const [left, right] = [60, 600];
  const sl = (v) => left + (v / 400) * (right - left);
  const st = (v) => left + ((v + 0.03) / 0.07) * (right - left);
  const [y1, y2] = [56, 150];
  let body = line(left, y1, right, y1, 's-ink') + line(left, y2, right, y2, 's-ink');
  for (const v of [0, 100, 200, 300, 400]) body += line(sl(v), y1, sl(v), y1 + 4, 's-ink');
  for (const v of [-0.02, 0, 0.02, 0.04]) body += line(st(v), y2, st(v), y2 + 4, 's-ink');
  body += `<text x="${left - 14}" y="${y1 + 5}" class="math" text-anchor="end">λ</text>`;
  body += `<text x="${left - 14}" y="${y2 + 5}" class="math" text-anchor="end">θ</text>`;
  body += line(sl(sigma), y1 - 22, sl(sigma), y1, 's-ink', ' stroke-dasharray="3 3"') + `<text x="${n(sl(sigma))}" y="${y1 - 28}" class="math" text-anchor="middle">σ</text>`;
  const theta = lambdas.map((l) => 1 / (l - sigma));
  const nearest = theta.reduce((best, t, i) => (Math.abs(t) > Math.abs(theta[best]) ? i : best), 0);
  lambdas.forEach((l, i) => {
    const lit = i === nearest;
    body += line(sl(l), y1, st(theta[i]), y2, lit ? 's-glass' : 's-muted', lit ? ' stroke-width="2"' : ' stroke-opacity="0.5"');
    body += circle(sl(l), y1, 4.5, lit ? 'f-glass' : 'f-muted') + circle(st(theta[i]), y2, 4.5, lit ? 'f-glass' : 'f-muted');
    body += `<text x="${n(sl(l))}" y="${y1 - 10}" class="math" text-anchor="middle" font-size="13">λ<tspan font-size="10" dy="3" class="upright">${i + 1}</tspan></text>`;
  });
  // Numbers last, haloed, so the lines crossing them don't strike them out.
  for (const v of [0, 100, 200, 300, 400]) body += text(sl(v), y1 + 18, v, 'halo', 'middle');
  for (const v of [-0.02, 0, 0.02, 0.04]) body += text(st(v), y2 + 18, v.toString().replace('-', '−'), 'halo', 'middle');
  figure('shift-invert', 184, body);
}

// ---- Chapter 8: error against cells, for each order ----
{
  const cells = [2, 4, 8, 16, 32];
  const errors = [1, 2, 3, 4, 5].map((order) =>
    cells.map((c) => {
      const k = run({ layers: [{ thickness: 1, eps: [1, 0] }], polarization: 'te', left: 'pec', right: 'pec', mesh: { order, max_cell: 1 / c }, analysis: { type: 'cavity', n_modes: 1 } }).meta.fields[0].k0;
      return Math.abs(k * k - Math.PI * Math.PI);
    }),
  );
  console.log('  errors', errors.map((row) => row.map((e) => e.toExponential(2)).join(' ')).join(' | '));
  const [left, right, top, bottom] = [70, 520, 34, 256];
  const sx = (c) => left + (Math.log2(c / 2) / 4) * (right - left);
  const sy = (e) => top + ((1 - Math.log10(Math.max(e, 1e-15))) / 16) * (bottom - top);
  let body = rect(left, sy(1e-12), right - left, bottom - sy(1e-12), 'f-grid');
  body += text(right - 8, bottom - 10, 'round-off', '', 'end');
  for (const p of [0, -4, -8, -12]) body += line(left, sy(10 ** p), right, sy(10 ** p), 's-grid') + `<text x="${left - 8}" y="${n(sy(10 ** p) + 4)}" text-anchor="end">10<tspan font-size="9" dy="-6">${String(p).replace('-', '−')}</tspan></text>`;
  for (const c of cells) body += line(sx(c), bottom, sx(c), bottom + 4, 's-rule') + text(sx(c), bottom + 18, c, '', 'middle');
  body += line(left, bottom, right, bottom, 's-rule') + line(left, top, left, bottom, 's-rule');
  body += text(right, bottom + 36, 'Cells', '', 'end');
  body += text(left, top - 14, 'Error in the lowest eigenvalue');
  errors.forEach((row, k) => {
    const opacity = 0.35 + (0.65 * k) / 4;
    const points = row.map((e, i) => [sx(cells[i]), sy(e)]);
    body += poly(points, 's-glass', ` stroke-width="2" stroke-opacity="${n(opacity)}"`);
    for (const [x, y] of points) body += circle(x, y, 3.2, 'f-glass', ` fill-opacity="${n(opacity)}"`);
    const [x, y] = points[1];
    body += text(x + 6, y - 8, `order ${k + 1}`, 'ink-text halo');
  });
  figure('convergence', 302, body);
}

// ---- Chapter 9: a staircase and a fitted mesh ----
{
  const panel = (x0, title, mesh) => {
    const [size, y0] = [256, 34];
    const p = (v) => [x0 + mesh.coords[2 * v] * size, y0 + (1 - mesh.coords[2 * v + 1]) * size];
    let body = text(x0, 20, title, 'title');
    for (let c = 0; c < mesh.region.length; c++) {
      const vs = mesh.cells.slice(3 * c, 3 * c + 3).map(p);
      body += `<polygon points="${vs.map(([x, y]) => `${n(x)},${n(y)}`).join(' ')}" class="${mesh.region[c] ? 'f-glass s-muted' : 's-rule'}"${mesh.region[c] ? ` fill-opacity="${tint(2)}" stroke-opacity="0.6"` : ''} stroke-width="0.8"/>`;
    }
    body += circle(x0 + 0.5 * size, y0 + 0.5 * size, 0.3 * size, 's-ink', ' stroke-width="1.4" stroke-dasharray="4 3"');
    return body;
  };
  figure('fitted-mesh', 300, panel(40, 'Staircase', data.staircase) + panel(344, 'Fitted', data.fitted));
}

// ---- Chapter 10: nested dissection ----
{
  const { n: size, positions, triangles, dissected } = data.mesh;
  const box = 176;
  let body = text(10, 20, 'Order of elimination', 'title');
  const p = (v) => [10 + positions[v][0] * box, 34 + (1 - positions[v][1]) * box];
  for (const t of triangles) body += `<polygon points="${t.map(p).map(([x, y]) => `${n(x)},${n(y)}`).join(' ')}" class="s-rule" stroke-width="0.6"/>`;
  const byOrder = [...positions.keys()].sort((a, b) => dissected[a] - dissected[b]);
  for (const v of byOrder) {
    const late = (dissected[v] / (size - 1)) ** 6;
    body += circle(...p(v), 2.1, 'f-glass', ` fill-opacity="${n(0.12 + 0.88 * late)}"`);
  }
  body += text(10, 34 + box + 20, 'pale first, dark last');
  const spy = (x0, title, entries, id) => {
    const pixels = 175;
    const grey = new Uint8Array(pixels * pixels);
    const cover = new Float32Array(pixels * pixels);
    const per = (size / pixels) ** 2;
    for (const [i, j] of entries) {
      for (const [r, c] of i === j ? [[i, j]] : [[i, j], [j, i]]) cover[Math.floor((r * pixels) / size) * pixels + Math.floor((c * pixels) / size)] += 1 / per;
    }
    cover.forEach((v, k) => (grey[k] = Math.round(255 * Math.min(1, Math.sqrt(v * 3)))));
    let svg = text(x0, 20, title, 'title');
    svg += rect(x0, 34, box, box, 's-rule');
    svg += masked(id, x0, 34, box, box, pixels, pixels, grey, 'ink');
    svg += text(x0, 34 + box + 20, `${thousands(entries.length)} entries in L`);
    return svg;
  };
  body += spy(224, 'Banded, after RCM', data.mesh.banded_factor, 'nd-banded');
  body += spy(438, 'Nested dissection', data.mesh.dissected_factor, 'nd-dissected');
  figure('nested-dissection', 236, body);
}

// ---- Chapter 11: a line source radiating into an absorbing layer ----
{
  const s = defaultState();
  s.analysis = 'scatter2d';
  s.source = 'line';
  Object.assign(s.plane.scatter2d, { shapes: [], sourceX: 1.5, sourceY: 1.5, cell: 0.05 });
  const [request] = buildRequests(s);
  const { channels } = run(request);
  const { x, y, triangles } = channels;
  const field = channels['field0/Ez/re'];
  const [lo, hi] = [Math.min(...x), Math.max(...x)];
  const pml = request.analysis.pml;
  const pixels = 240;
  const value = new Float32Array(pixels * pixels).fill(NaN);
  const px = (v) => ((v - lo) / (hi - lo)) * pixels - 0.5;
  for (let t = 0; t < triangles.length; t += 3) {
    const [a, b, c] = [triangles[t], triangles[t + 1], triangles[t + 2]];
    const [ax, ay, bx, by, cx, cy] = [px(x[a]), px(y[a]), px(x[b]), px(y[b]), px(x[c]), px(y[c])];
    const det = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy);
    for (let j = Math.max(0, Math.ceil(Math.min(ay, by, cy))); j <= Math.min(pixels - 1, Math.floor(Math.max(ay, by, cy))); j++) {
      for (let i = Math.max(0, Math.ceil(Math.min(ax, bx, cx))); i <= Math.min(pixels - 1, Math.floor(Math.max(ax, bx, cx))); i++) {
        const l1 = ((by - cy) * (i - cx) + (cx - bx) * (j - cy)) / det;
        const l2 = ((cy - ay) * (i - cx) + (ax - cx) * (j - cy)) / det;
        const l3 = 1 - l1 - l2;
        if (l1 < -1e-9 || l2 < -1e-9 || l3 < -1e-9) continue;
        value[(pixels - 1 - j) * pixels + i] = l1 * field[a] + l2 * field[b] + l3 * field[c];
      }
    }
  }
  // The source's own peak would wash everything else out.
  const sorted = [...value].map(Math.abs).sort((p, q) => p - q);
  const peak = sorted[Math.floor(sorted.length * 0.985)];
  const part = (sign) => Uint8Array.from(value, (v) => Math.round(255 * Math.min(1, Math.max(0, (sign * v) / peak)) ** 0.85));
  const [mx, my, size] = [10, 30, 250];
  const domain = (v) => mx + ((v - lo) / (hi - lo)) * size;
  let body = text(mx, 18, 'The field, in the domain and the layer', 'title');
  body += masked('pml-plus', mx, my, size, size, pixels, pixels, part(1), 'red');
  body += masked('pml-minus', mx, my, size, size, pixels, pixels, part(-1), 'blue');
  body += rect(mx, my, size, size, 's-rule');
  body += rect(domain(0), my + size - (domain(3) - mx), domain(3) - domain(0), domain(3) - domain(0), 's-ink', ' stroke-dasharray="4 3" stroke-width="1.2"');
  // The cut along the middle row of nodes, through the source.
  const row = [...x.keys()].filter((k) => Math.abs(y[k] - 1.5) < 1e-6).sort((p, q) => x[p] - x[q]);
  const [left, right, top, bottom] = [320, 620, 40, 250];
  const sx = (v) => left + ((v - lo) / (hi - lo)) * (right - left);
  const cutPeak = Math.max(...row.map((k) => Math.abs(field[k])));
  const sy = (v) => (top + bottom) / 2 - (v / cutPeak) * ((bottom - top) / 2 - 22);
  body += text(left, 18, 'Along the middle', 'title');
  for (const [a, b] of [[lo, 0], [3, hi]]) body += rect(sx(a), top, sx(b) - sx(a), bottom - top, 'f-grid');
  body += text(sx(lo + pml / 2), top + 16, 'layer', '', 'middle') + text(sx(3 + pml / 2), top + 16, 'layer', '', 'middle');
  body += line(left, (top + bottom) / 2, right, (top + bottom) / 2, 's-rule');
  body += poly(row.map((k) => [sx(x[k]), sy(field[k])]), 's-red', ' stroke-width="1.8"');
  body += line(left, bottom, right, bottom, 's-rule');
  for (const v of [0, 1.5, 3]) body += line(sx(v), bottom, sx(v), bottom + 4, 's-rule') + text(sx(v), bottom + 18, v, '', 'middle');
  body += text(right, bottom + 36, 'Position (µm)', '', 'end');
  figure('absorbing-layer', 296, body);
}

// ---- Chapter 13: a bend's density, sharpening ----
{
  const s = defaultState();
  s.analysis = 'bend2d';
  const solveAll = async (requests) => ({ results: requests.map(run) });
  const stage = design.stageLength(s);
  const keep = new Map([[stage, 'β = 1'], [3 * stage, 'β = 4'], [5 * stage, 'β = 16'], [7 * stage, 'β = 64']]);
  const shots = [];
  const started = performance.now();
  await design.optimize(s, solveAll, ({ rho, iteration, beta }) => {
    if (beta === Infinity) shots.push({ rho, label: 'Final', iteration });
    else if (keep.has(iteration)) shots.push({ rho, label: keep.get(iteration), iteration });
    return true;
  });
  console.log(`  bend design: ${shots.length} snapshots in ${((performance.now() - started) / 1000).toFixed(0)} s`);
  const side = Math.round(Math.sqrt(shots[0].rho.length));
  const [box, gap] = [104, 22];
  let body = '';
  shots.forEach(({ rho, label, iteration }, k) => {
    const x0 = 10 + k * (box + gap);
    // Pixel rows run up the page in the solver, down in an image.
    const grey = new Uint8Array(side * side);
    for (let j = 0; j < side; j++) for (let i = 0; i < side; i++) grey[(side - 1 - j) * side + i] = Math.round(255 * rho[j * side + i]);
    body += masked(`density-${k}`, x0, 10, box, box, side, side, grey, 'glass', true);
    body += rect(x0, 10, box, box, 's-rule');
    body += text(x0, box + 30, label, 'ink-text');
    body += text(x0, box + 48, `iteration ${iteration}`);
  });
  figure('design-sharpening', 172, body);
}
