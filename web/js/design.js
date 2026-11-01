// Density topology optimization of a layered structure. The design region is a row
// of thin pixels, each with a density that mixes two materials. The solver supplies
// the adjoint gradient; everything else here is optimization, not physics.

import { MATERIALS, scaled, toShape } from './model.js';

const PIXEL = 0.005;
const AMBIENT = 0.3;
const SHARPNESS = [1, 2, 4, 8, 16, 32, 64, 128];
const ITERATIONS_PER_STAGE = 40;
// A 2D solve costs far more than a 1D one, and these converge sooner: 25 did no better
// over three starts, while 7 sometimes fell short, and so did skipping every other sharpness.
export const stageLength = (s) => (s.analysis === 'cavityDesign' ? CAVITY.iterations : s.analysis === 'bend2d' ? 10 : ITERATIONS_PER_STAGE);
const sharpnessOf = (s) => (s.analysis === 'cavityDesign' ? CAVITY.stages.map(([, beta]) => beta) : SHARPNESS);
export const totalIterations = (s) => sharpnessOf(s).length * stageLength(s);
export const STAGES = SHARPNESS.map((beta, i) => ({ frame: i * ITERATIONS_PER_STAGE, name: `β ${beta}` }));
export const sharpnessAt = (s, frame) => sharpnessOf(s)[Math.floor(frame / stageLength(s))] ?? Infinity;

const material = (name) => MATERIALS.find((m) => m.name === name);
const eps = (name) => material(name).n ** 2;
const solverLayer = (thickness, e) => ({ thickness, eps: [e, 0], mu: [1, 0] });

export const pixelCount = (d) => (d.cavity ? CAVITY.pixels ** 2 : d.thickness === undefined ? BEND.pixels ** 2 : Math.max(1, Math.round(d.thickness / PIXEL)));

// A slab of index 3.5 turning a corner through a square of pixels, in units of the band's
// centre wavelength, so a band anywhere gets the same design at its own size.
// Everything lies on a grid of 1/24, so each pixel is a whole number of pairs of triangles.
const BEND = { pixels: 24, x0: 1, y0: 0.5, size: 1, width: 0.125, domain: 2.5, n: 3.5, pml: 0.5, span: 0.35 };
const bendCentre = { x: BEND.x0 + BEND.size / 2, y: BEND.y0 + BEND.size / 2 };
const rectangle = (x0, y0, x1, y1) => ({ type: 'rect', x: (x0 + x1) / 2, y: (y0 + y1) / 2, w: x1 - x0, h: y1 - y0, n: BEND.n });

export const bendScale = (s) => (s.design.bend2d.from + s.design.bend2d.to) / 2;

// Cells per pixel side, and cells per wavelength in the slab at the band's short end. Three
// per pixel is as many as the solver's limit on unknowns allows.
function bendRefinement(s) {
  const shortest = s.design.bend2d.from / bendScale(s);
  const cell = BEND.size / BEND.pixels;
  const refine = Math.min(3, Math.ceil((4.5 * BEND.n * cell) / shortest));
  return { refine, cells: (shortest * refine) / (BEND.n * cell) };
}

// Why a band's percentages can't be trusted, if they can't. A slab this wide carries a
// second mode below 0.84 of the centre, and the port launches only the first.
export function bendWarning(s) {
  if (bendRefinement(s).cells < 4) return 'This band is too wide for the mesh to resolve its short end, so these percentages aren\'t reliable.';
  if (s.design.bend2d.from / bendScale(s) < 0.84) return 'The guide carries a second mode at the short end of this band, so these percentages are less reliable there. A narrower band avoids that.';
  return null;
}

// What the plot needs, and each solve shares.
export function bendGeometry(straight, scale) {
  const { x0, y0, size, width, domain, span } = BEND;
  const { x, y } = bendCentre;
  const far = domain + 3;
  const shapes = straight
    ? [rectangle(-3, y - width / 2, far, y + width / 2)]
    : [rectangle(-3, y - width / 2, x0, y + width / 2), rectangle(x - width / 2, y0 + size, x + width / 2, far)];
  const end = domain - 0.4;
  return scaled({
    size: domain, shapes,
    source: { type: 'port', x: 0.4, y0: y - 0.7, y1: y + 0.7 },
    mark: [0.4, y],
    monitor: straight ? { x0: end, y0: y - span, x1: end, y1: y + span } : { x0: x + span, y0: end, x1: x - span, y1: end },
    region: { x0, y0, width: size, height: size, nx: BEND.pixels, ny: BEND.pixels },
  }, scale);
}

function bendRequest(s, wavelength, straight, rho) {
  const scale = bendScale(s);
  const g = bendGeometry(straight, scale);
  const n2 = BEND.n ** 2;
  return {
    width: g.size,
    height: g.size,
    cell: (BEND.size * scale) / BEND.pixels / bendRefinement(s).refine,
    background: { eps: 1, mu: 1 },
    shapes: g.shapes.map(toShape),
    polarization: 'te',
    walls: 'absorbing',
    monitors: [g.monitor],
    analysis: { type: 'driven', wavelength, source: g.source, pml: BEND.pml * scale },
    ...(rho ? { density: { ...g.region, eps_low: 1, eps_high: n2, rho: Array.from(rho) } } : {}),
  };
}

// A straight guide's power at each wavelength, by its `k0`, which turns power into transmission.
const straightPower = new Map();

// Solves the straight guide once for these wavelengths. Returns a reply holding any error.
export async function prepareDesign(s, solve) {
  if (s.analysis !== 'bend2d') return {};
  const wavelengths = goalWavelengths(s);
  const reply = await solve(wavelengths.map((wavelength) => bendRequest(s, wavelength, true)));
  if (!reply || reply.error) return reply;
  reply.results.forEach((r) => straightPower.set(r.meta.fields[0].k0, r.channels.monitors[0]));
  return reply;
}

// A square of pixels around a line source, in units of the wavelength, as big as 7 x 7 rods.
// `stages` raises the bandwidth's Q to a level a sharpening design can follow, then sharpens.
export const CAVITY_Q = 1000;
const CAVITY = {
  pixels: 64, cellsPerPixel: 1, margin: 0.25, pml: 0.5, eps: 12, radius: 2, iterations: 20,
  stages: [[10, 1], [20, 2], [50, 4], [100, 8], [200, 16], [300, 32], [500, 64], [500, 128], [500, 256], [500, 512], [500, 1024]],
  // Measured at this Q once finished, and in the table.
  Q: CAVITY_Q,
};
const cavitySide = CAVITY.pixels / 24;
// How much a finished design was resized to put its resonance on the target, by its densities.
const tuned = new WeakMap();

export function cavityGeometry(wavelength) {
  const { margin } = CAVITY;
  const size = cavitySide + 2 * margin;
  return scaled({
    size, shapes: [],
    source: { type: 'line', x: size / 2, y: size / 2 },
    mark: [size / 2, size / 2],
    region: { x0: margin, y0: margin, width: cavitySide, height: cavitySide, nx: CAVITY.pixels, ny: CAVITY.pixels },
  }, wavelength);
}

// What the plot needs for a 2D design.
export const designGeometry = (s, rho) => (s.analysis === 'bend2d'
  ? bendGeometry(false, bendScale(s))
  : cavityGeometry(s.wavelength * (tuned.get(rho) ?? 1)));

function cavityRequest(s, rho, Q, wavelength = s.wavelength) {
  const scale = s.wavelength * (tuned.get(rho) ?? 1);
  const g = cavityGeometry(scale);
  return {
    width: g.size,
    height: g.size,
    cell: (cavitySide * scale) / CAVITY.pixels,
    background: { eps: 1, mu: 1 },
    shapes: [],
    polarization: 'te',
    walls: 'absorbing',
    monitors: [],
    analysis: { type: 'driven', wavelength, source: g.source, pml: CAVITY.pml * scale, average_q: Q },
    density: { ...g.region, eps_low: 1, eps_high: CAVITY.eps, rho: Array.from(rho) },
  };
}

const times = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const over = (a, b) => {
  const d = b[0] * b[0] + b[1] * b[1];
  return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
};
const complexK = (wavelength, Q) => [(2 * Math.PI) / wavelength, Math.PI / (wavelength * Q)];

// Where a design's resonance is, and how each density moves it. Close to a pole, 1 / u_0 is
// nearly linear in k, so two solves place it. The solver's two gradient channels are the slope
// of `2 u_0`, which through `1 / u_0` is also the slope of the pole's own k.
export async function findResonance(s, rho, Q, solve, wavelength = s.wavelength) {
  const beside = wavelength * (1 + 1 / (4 * Q) + 1e-3);
  const reply = await solve([cavityRequest(s, rho, Q, wavelength), cavityRequest(s, rho, Q, beside)]);
  if (!reply || reply.error) return { reply };
  const [here, there] = reply.results;
  const response = here.meta.fields[0].response;
  const [k, kBeside] = [complexK(wavelength, Q), complexK(beside, Q)];
  const [g, gBeside] = [over([1, 0], response), over([1, 0], there.meta.fields[0].response)];
  const slope = over([kBeside[0] - k[0], kBeside[1] - k[1]], [gBeside[0] - g[0], gBeside[1] - g[1]]);
  const factor = times(slope, over([1, 0], times(response, response)));
  const { gradient, gradient_im: imaginary } = here.channels;
  return {
    reply: { ...reply, results: [here] },
    poleK: k[0] - times(g, slope)[0],
    moves: Float64Array.from(gradient, (re, i) => (factor[0] * re - factor[1] * imaginary[i]) / 2),
  };
}

// Pixels as ordinary rectangles, which need nothing of densities to solve.
function cavityShapes(rho, scale) {
  const { pixels, margin } = CAVITY;
  const h = cavitySide / pixels;
  const shapes = [];
  for (let j = 0; j < pixels; j++) {
    for (let i = 0; i < pixels;) {
      if (rho[i + j * pixels] <= 0.5) { i++; continue; }
      let end = i;
      while (end < pixels && rho[end + j * pixels] > 0.5) end++;
      shapes.push({ material: { eps: CAVITY.eps, mu: 1 }, type: 'rect', x0: (margin + i * h) * scale, y0: (margin + j * h) * scale, x1: (margin + end * h) * scale, y1: (margin + (j + 1) * h) * scale });
      i = end;
    }
  }
  return shapes;
}

// Eight-fold symmetry around the source, which is its own adjoint.
function mirrored(n) {
  const images = (i, j) => [[i, j], [n - 1 - i, j], [i, n - 1 - j], [n - 1 - i, n - 1 - j], [j, i], [n - 1 - j, i], [j, n - 1 - i], [n - 1 - j, n - 1 - i]];
  const lists = Array.from({ length: n * n }, (_, k) => images(k % n, Math.floor(k / n)).map(([x, y]) => x + y * n));
  return (v) => Float64Array.from(lists, (list) => list.reduce((sum, m) => sum + v[m], 0) / 8);
}

// A cavity's goal is LDOS, as in Liang & Johnson 2013, with the frequency held on target. A
// high Q's peak is narrow, and a free step moves it further than its width, so each step here
// has its slope along the frequency removed. A correction along that slope then keeps the
// resonance on target, and scaling the finished design takes out what sharpening leaves.
async function optimizeCavity(s, solve, onStep) {
  const n = CAVITY.pixels ** 2;
  const filter = coneFilter2d(CAVITY.pixels, CAVITY.pixels, CAVITY.radius);
  const symmetric = mirrored(CAVITY.pixels);
  const through = (g, slope) => symmetric(filter.adjoint(g.map((v, k) => v * slope[k])));
  const target = (2 * Math.PI) / s.wavelength;
  let x = randomStart(n, s.design.cavityDesign.seed);
  let iteration = 0;
  for (const [Q, beta] of CAVITY.stages) {
    const reach = Math.min(0.1, 0.2 / beta);
    const mma = new Mma(reach);
    for (let i = 0; i < CAVITY.iterations; i++) {
      const { rho, slope } = project(filter.apply(symmetric(x)), beta);
      const { reply, poleK, moves } = await findResonance(s, rho, Q, solve);
      if (!reply || reply.error) return reply;
      if (!onStep({ rho, reply, iteration: ++iteration, beta, bandwidth: Q })) return null;
      const goal = through(reply.results[0].channels.gradient.map((v) => (-2 * v) / target), slope);
      const frequency = through(moves, slope);
      const square = frequency.reduce((sum, v) => sum + v * v, 0);
      const along = goal.reduce((sum, v, k) => sum + v * frequency[k], 0) / square;
      const next = mma.step(x, goal.map((v, k) => v - along * frequency[k]));
      const moved = next.reduce((sum, v, k) => sum + (v - x[k]) * frequency[k], 0);
      const largest = frequency.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
      const correction = Math.max(-reach / largest, Math.min(reach / largest, (0.5 * (target - poleK) - moved) / square));
      x = next.map((v, k) => Math.min(1, Math.max(0, v + correction * frequency[k])));
    }
  }
  const rho = project(filter.apply(symmetric(x)), 1000).rho.map((v) => (v > 0.5 ? 1 : 0));
  // Find the peak, and place it precisely, then resize so that it's on the target.
  const scan = Array.from({ length: 21 }, (_, i) => s.wavelength * (0.96 + 0.004 * i));
  const swept = await solve(scan.map((wavelength) => cavityRequest(s, rho, 100, wavelength)));
  if (!swept || swept.error) return swept;
  const strengths = swept.results.map((r, i) => r.meta.fields[0].ldos * scan[i]);
  let peak = scan[strengths.indexOf(Math.max(...strengths))];
  for (const Q of [100, 300, 1000, 1000]) {
    const { reply, poleK } = await findResonance(s, rho, Q, solve, peak);
    if (!reply || reply.error) return reply;
    const next = (2 * Math.PI) / poleK;
    if (Math.abs(next - peak) < 0.02 * peak) peak = next;
  }
  tuned.set(rho, s.wavelength / peak);
  const reply = await solve(gradientRequests(s, rho));
  if (!reply || reply.error) return reply;
  // Checks it as a drawn structure, with the search for leaky modes, which shares nothing of this.
  const scale = tuned.get(rho) * s.wavelength;
  const g = cavityGeometry(scale);
  const found = await solve([{ width: g.size, height: g.size, cell: g.region.width / CAVITY.pixels, background: { eps: 1, mu: 1 }, shapes: cavityShapes(rho, scale), polarization: 'te', walls: 'absorbing', monitors: [], analysis: { type: 'resonance', target_k0: target, x: g.source.x, y: g.source.y, pml: CAVITY.pml * scale } }]);
  const field = found?.results?.[0]?.meta.fields[0];
  const mode = field && { wavelength: (2 * Math.PI) / field.k0, q: field.q, volume: field.volume / s.wavelength ** 2 };
  onStep({ rho, reply, iteration, beta: Infinity, mode, bandwidth: CAVITY.Q });
  return null;
}

// Goal wavelengths, with the centre last so the plotted field comes from it.
function goalWavelengths(s) {
  if (s.analysis === 'mirror') return [s.wavelength];
  const { from, to } = s.design[s.analysis];
  const count = s.analysis === 'bend2d' ? 3 : 11;
  const wavelengths = Array.from({ length: count }, (_, i) => from + ((to - from) * i) / (count - 1));
  wavelengths.push(wavelengths.splice((count - 1) / 2, 1)[0]);
  return wavelengths;
}

function requests(s, layers, makeAnalysis) {
  const wavelengths = goalWavelengths(s);
  const max_cell = Math.min(...wavelengths) / 6;
  return wavelengths.map((wavelength) => ({
    layers,
    polarization: s.polarization,
    left: 'port',
    right: 'absorbing',
    mesh: { order: 2, max_cell },
    analysis: makeAnalysis(wavelength),
  }));
}

function pixelLayers(d, rho) {
  const [low, high] = [eps(d.low), eps(d.high)];
  return [
    solverLayer(AMBIENT, eps(d.left)),
    ...Array.from(rho, (r) => solverLayer(PIXEL, low + r * (high - low))),
    solverLayer(AMBIENT, eps(d.right)),
  ];
}

export const gradientRequests = (s, rho) => (s.analysis === 'cavityDesign' ? [cavityRequest(s, rho, CAVITY.Q)] : s.analysis === 'bend2d'
  ? goalWavelengths(s).map((wavelength) => bendRequest(s, wavelength, false, rho))
  : requests(s, pixelLayers(s.design[s.analysis], rho),
    (wavelength) => ({ type: 'driven_gradient', wavelength, angle_deg: s.angle })));

// What a design is judged by: reflection for layers, transmission for a bend, and for a
// cavity, LDOS relative to vacuum, where a line source's `F` is k / 2.
const measure = ({ meta, channels }) => (meta.fields[0].ldos != null
  ? (2 * meta.fields[0].ldos) / meta.fields[0].k0
  : channels.monitors
  ? channels.monitors[0] / straightPower.get(meta.fields[0].k0)
  : meta.fields[0].power[0]);
export const meanReflection = (results) => results.reduce((sum, r) => sum + measure(r), 0) / results.length;

// The structures a designer would try by hand, solved at the same wavelengths.
export function references(s) {
  if (s.analysis === 'cavityDesign') {
    // Rods with one removed. A period of 0.3781 puts this crystal's mode on the wavelength.
    const a = 0.3781 * s.wavelength;
    const { margin, pml } = CAVITY;
    const edge = margin * s.wavelength;
    const rods = [];
    for (let i = 0; i < 49; i++) {
      if (i !== 24) rods.push({ type: 'circle', x: edge + (i % 7 + 0.5) * a, y: edge + (Math.floor(i / 7) + 0.5) * a, w: 0.4 * a, h: 0.4 * a, n: Math.sqrt(CAVITY.eps) });
    }
    const size = 7 * a + 2 * edge;
    return [{
      name: 'Crystal of 7 x 7 rods',
      requests: [{
        width: size, height: size, cell: a / 12, background: { eps: 1, mu: 1 }, shapes: rods.map(toShape), polarization: 'te', walls: 'absorbing', monitors: [],
        analysis: { type: 'driven', wavelength: s.wavelength, source: { type: 'line', x: size / 2, y: size / 2 }, pml: pml * s.wavelength, average_q: CAVITY.Q },
      }],
    }];
  }
  if (s.analysis === 'bend2d') {
    return [{ name: 'Plain corner', requests: gradientRequests(s, new Float64Array(BEND.pixels ** 2)) }];
  }
  const d = s.design[s.analysis];
  const quarter = (name, wavelength) => wavelength / (4 * material(name).n);
  const wrap = (inner) => [
    solverLayer(AMBIENT, eps(d.left)),
    ...inner.map(([name, t]) => solverLayer(t, eps(name))),
    solverLayer(AMBIENT, eps(d.right)),
  ];
  const plain = (name, inner) => ({
    name,
    requests: requests(s, wrap(inner), (wavelength) => ({ type: 'driven', wavelength, angle_deg: s.angle, source: { type: 'plane_wave' } })),
  });
  if (s.analysis === 'mirror') {
    const [hi, lo] = [quarter(d.high, s.wavelength), quarter(d.low, s.wavelength)];
    const pairs = Math.max(0, Math.floor((d.thickness - hi) / (hi + lo)));
    const stack = [[d.high, hi], ...Array.from({ length: pairs }, () => [[d.low, lo], [d.high, hi]]).flat()];
    return [plain(`Quarter-wave stack, ${stack.length} layers`, stack)];
  }
  return [
    plain('No coating', []),
    plain('One quarter-wave layer', [[d.low, quarter(d.low, (d.from + d.to) / 2)]]),
  ];
}

function randomStart(n, seed) {
  // Mulberry32: the same seed always gives the same start.
  let a = Math.imul(seed, 0x9e3779b1) >>> 0;
  return Float64Array.from({ length: n }, () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = Math.imul(a ^ (a >>> 15), a | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return 0.5 + 0.4 * (((t ^ (t >>> 14)) >>> 0) / 2 ** 32 - 0.5);
  });
}

// A cone filter averages each pixel with its neighbours, which removes single-pixel
// noise. `adjoint` carries gradients back through the same weights.
export function coneFilter(n, radius) {
  const reach = Math.ceil(radius) - 1;
  const weight = (i, j) => radius - Math.abs(i - j);
  const window = (i) => [Math.max(0, i - reach), Math.min(n - 1, i + reach)];
  const total = Float64Array.from({ length: n }, (_, i) => {
    const [a, b] = window(i);
    let sum = 0;
    for (let j = a; j <= b; j++) sum += weight(i, j);
    return sum;
  });
  return {
    apply(x) {
      return Float64Array.from({ length: n }, (_, i) => {
        const [a, b] = window(i);
        let sum = 0;
        for (let j = a; j <= b; j++) sum += weight(i, j) * x[j];
        return sum / total[i];
      });
    },
    adjoint(g) {
      const out = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        const [a, b] = window(i);
        for (let j = a; j <= b; j++) out[j] += (weight(i, j) * g[i]) / total[i];
      }
      return out;
    },
  };
}

// The same filter over a rectangle of pixels, stored by rows.
export function coneFilter2d(nx, ny, radius) {
  const reach = Math.ceil(radius) - 1;
  const neighbours = Array.from({ length: nx * ny }, (_, i) => {
    const [cx, cy] = [i % nx, Math.floor(i / nx)];
    const list = [];
    for (let y = Math.max(0, cy - reach); y <= Math.min(ny - 1, cy + reach); y++) {
      for (let x = Math.max(0, cx - reach); x <= Math.min(nx - 1, cx + reach); x++) {
        const weight = radius - Math.hypot(x - cx, y - cy);
        if (weight > 0) list.push([x + y * nx, weight]);
      }
    }
    const total = list.reduce((sum, [, w]) => sum + w, 0);
    return list.map(([j, w]) => [j, w / total]);
  });
  return {
    apply: (x) => Float64Array.from(neighbours, (list) => list.reduce((sum, [j, w]) => sum + w * x[j], 0)),
    adjoint(g) {
      const out = new Float64Array(nx * ny);
      neighbours.forEach((list, i) => list.forEach(([j, w]) => { out[j] += w * g[i]; }));
      return out;
    },
  };
}

// Pushes densities toward 0 or 1. Larger `beta` is sharper.
export function project(x, beta) {
  const half = Math.tanh(beta / 2);
  const rho = x.map((v) => (half + Math.tanh(beta * (v - 0.5))) / (2 * half));
  const slope = x.map((v) => beta / (2 * half * Math.cosh(beta * (v - 0.5)) ** 2));
  return { rho, slope };
}

// Svanberg's method of moving asymptotes, for bounds only. Each step minimizes a
// convex approximation, which has a closed form when there are no other constraints.
export class Mma {
  // `move` is the largest change to any variable in one step.
  constructor(move) {
    this.move = move;
    this.history = [];
  }

  // A new `x` in [0, 1] that should reduce an objective with gradient `g` at `x`.
  step(x, g) {
    const n = x.length;
    const [low, high, next] = [new Float64Array(n), new Float64Array(n), new Float64Array(n)];
    const [prev, prev2] = [this.history.at(-1), this.history.at(-2)];
    for (let i = 0; i < n; i++) {
      if (prev2) {
        // Oscillating variables get tighter asymptotes, steady ones looser.
        const gamma = (x[i] - prev[i]) * (prev[i] - prev2[i]) < 0 ? 0.7 : 1.2;
        low[i] = Math.min(x[i] - 0.01, Math.max(x[i] - 10, x[i] - gamma * (prev[i] - this.low[i])));
        high[i] = Math.max(x[i] + 0.01, Math.min(x[i] + 10, x[i] + gamma * (this.high[i] - prev[i])));
      } else {
        low[i] = x[i] - 2 * this.move;
        high[i] = x[i] + 2 * this.move;
      }
      const p = (high[i] - x[i]) * Math.sqrt(Math.max(g[i], 0)) + 1e-7;
      const q = (x[i] - low[i]) * Math.sqrt(Math.max(-g[i], 0)) + 1e-7;
      const best = (p * low[i] + q * high[i]) / (p + q);
      const lower = Math.max(0, x[i] - this.move, low[i] + 0.1 * (x[i] - low[i]));
      const upper = Math.min(1, x[i] + this.move, high[i] - 0.1 * (high[i] - x[i]));
      next[i] = Math.min(upper, Math.max(lower, best));
    }
    [this.low, this.high] = [low, high];
    this.history.push(x);
    return next;
  }
}

// The slope of what the optimizer minimizes, for each pixel's density before projection,
// where `slope` is the projection's own. A mirror and a bend want more of what they measure.
export function goalGradient(s, results, slope) {
  const d = s.design[s.analysis];
  const gradient = new Float64Array(slope.length);
  for (const { channels, meta } of results) {
    const [slopes, scale] = s.analysis === 'bend2d'
      ? [channels.gradient, -1 / straightPower.get(meta.fields[0].k0)]
      : [channels['dR/deps_re'], (s.analysis === 'mirror' ? -1 : 1) * (eps(d.high) - eps(d.low))];
    for (let k = 0; k < gradient.length; k++) gradient[k] += (scale * slopes[k] * slope[k]) / results.length;
  }
  return gradient;
}

// Runs until the schedule ends or `onStep` returns false. `solve` returns a solver
// reply, or null if a newer request replaced it. Returns a reply holding any error.
export async function optimize(s, solve, onStep) {
  if (s.analysis === 'cavityDesign') return optimizeCavity(s, solve, onStep);
  const d = s.design[s.analysis];
  const bend = s.analysis === 'bend2d';
  const n = pixelCount(d);
  const filter = bend ? coneFilter2d(BEND.pixels, BEND.pixels, d.radius) : coneFilter(n, Math.max(d.radius / (PIXEL * 1000), 1));
  const prepared = await prepareDesign(s, solve);
  if (!prepared || prepared.error) return prepared;
  let x = randomStart(n, d.seed);
  let iteration = 0;
  for (const beta of SHARPNESS) {
    // Projection amplifies a change in x by about beta, so steps shrink to match.
    const mma = new Mma(Math.min(0.1, 0.2 / beta));
    for (let i = 0; i < stageLength(s); i++) {
      const { rho, slope } = project(filter.apply(x), beta);
      const reply = await solve(gradientRequests(s, rho));
      if (!reply || reply.error) return reply;
      if (!onStep({ rho, reply, iteration: ++iteration, beta })) return null;
      x = mma.step(x, filter.adjoint(goalGradient(s, reply.results, slope)));
    }
  }
  const { rho } = project(filter.apply(x), 1000);
  const reply = await solve(gradientRequests(s, rho));
  if (!reply || reply.error) return reply;
  onStep({ rho, reply, iteration, beta: Infinity });
  return null;
}

// Merges neighbouring pixels of nearly equal density into layers for plotting.
export function plotLayers(d, rho) {
  const [low, high] = [eps(d.low), eps(d.high)];
  const out = [{ ...material(d.left), thickness: AMBIENT }];
  let start = 0;
  for (let i = 1; i <= rho.length; i++) {
    if (i === rho.length || Math.abs(rho[i] - rho[start]) > 0.02) {
      out.push({ n: Math.sqrt(low + rho[start] * (high - low)), k: 0, thickness: (i - start) * PIXEL });
      start = i;
    }
  }
  out.push({ ...material(d.right), thickness: AMBIENT });
  return out;
}

// A finished design as ordinary editable layers. Densities must already be sharp.
export function finishedLayers(d, rho) {
  const out = [{ ...material(d.left), thickness: 0.5 }];
  let start = 0;
  for (let i = 1; i <= rho.length; i++) {
    if (i === rho.length || rho[i] > 0.5 !== rho[start] > 0.5) {
      out.push({ ...material(rho[start] > 0.5 ? d.high : d.low), thickness: +((i - start) * PIXEL).toFixed(6) });
      start = i;
    }
  }
  out.push({ ...material(d.right), thickness: 0.5 });
  return out;
}
