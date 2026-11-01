// Density topology optimization of a layered structure. The design region is a row
// of thin pixels, each with a density that mixes two materials. The solver supplies
// the adjoint gradient; everything else here is optimization, not physics.

import { MATERIALS, toShape } from './model.js';

const PIXEL = 0.005;
const AMBIENT = 0.3;
const SHARPNESS = [1, 2, 4, 8, 16, 32, 64, 128];
const ITERATIONS_PER_STAGE = 40;
// A 2D solve costs far more than a 1D one, and these converge sooner.
const stageLength = (s) => (s.analysis === 'bend2d' ? 25 : ITERATIONS_PER_STAGE);
export const totalIterations = (s) => SHARPNESS.length * stageLength(s);
export const STAGES = SHARPNESS.map((beta, i) => ({ frame: i * ITERATIONS_PER_STAGE, name: `β ${beta}` }));
export const sharpnessAt = (s, frame) => SHARPNESS[Math.floor(frame / stageLength(s))] ?? Infinity;

const material = (name) => MATERIALS.find((m) => m.name === name);
const eps = (name) => material(name).n ** 2;
const solverLayer = (thickness, e) => ({ thickness, eps: [e, 0], mu: [1, 0] });

export const pixelCount = (d) => (d.thickness === undefined ? BEND.pixels ** 2 : Math.max(1, Math.round(d.thickness / PIXEL)));

// A slab of index 3.5 turning a corner through a square of pixels, in wavelengths of about 1.
// Everything lies on a grid of 1/24, so each pixel is exactly a pair of triangles.
const BEND = { pixels: 24, x0: 1, y0: 0.5, size: 1, width: 0.125, domain: 2.5, n: 3.5, pml: 0.5, span: 0.35 };
const bendCentre = { x: BEND.x0 + BEND.size / 2, y: BEND.y0 + BEND.size / 2 };
const rectangle = (x0, y0, x1, y1) => ({ type: 'rect', x: (x0 + x1) / 2, y: (y0 + y1) / 2, w: x1 - x0, h: y1 - y0, n: BEND.n });

// What the plot needs, and each solve shares.
export function bendGeometry(straight = false) {
  const { x0, y0, size, width, domain, span } = BEND;
  const { x, y } = bendCentre;
  const far = domain + 3;
  const shapes = straight
    ? [rectangle(-3, y - width / 2, far, y + width / 2)]
    : [rectangle(-3, y - width / 2, x0, y + width / 2), rectangle(x - width / 2, y0 + size, x + width / 2, far)];
  const end = domain - 0.4;
  return {
    size: domain, shapes,
    source: { type: 'port', x: 0.4, y0: y - 0.7, y1: y + 0.7 },
    mark: [0.4, y],
    monitor: straight ? { x0: end, y0: y - span, x1: end, y1: y + span } : { x0: x + span, y0: end, x1: x - span, y1: end },
    region: { x0, y0, width: size, height: size, nx: BEND.pixels, ny: BEND.pixels },
  };
}

function bendRequest(wavelength, straight, rho) {
  const g = bendGeometry(straight);
  const n2 = BEND.n ** 2;
  return {
    width: g.size,
    height: g.size,
    cell: BEND.size / BEND.pixels,
    background: { eps: 1, mu: 1 },
    shapes: g.shapes.map(toShape),
    polarization: 'te',
    walls: 'absorbing',
    monitors: [g.monitor],
    analysis: { type: 'driven', wavelength, source: g.source, pml: BEND.pml },
    ...(rho ? { density: { ...g.region, eps_low: 1, eps_high: n2, rho: Array.from(rho) } } : {}),
  };
}

// A straight guide's power at each wavelength, by its `k0`, which turns power into transmission.
const straightPower = new Map();

// Solves the straight guide once for these wavelengths. Returns a reply holding any error.
export async function prepareDesign(s, solve) {
  if (s.analysis !== 'bend2d') return {};
  const wavelengths = goalWavelengths(s);
  const reply = await solve(wavelengths.map((wavelength) => bendRequest(wavelength, true)));
  if (!reply || reply.error) return reply;
  reply.results.forEach((r) => straightPower.set(r.meta.fields[0].k0, r.channels.monitors[0]));
  return reply;
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

export const gradientRequests = (s, rho) => (s.analysis === 'bend2d'
  ? goalWavelengths(s).map((wavelength) => bendRequest(wavelength, false, rho))
  : requests(s, pixelLayers(s.design[s.analysis], rho),
    (wavelength) => ({ type: 'driven_gradient', wavelength, angle_deg: s.angle })));

// What a design is judged by: reflection for layers, and transmission for a bend.
const measure = ({ meta, channels }) => (channels.monitors
  ? channels.monitors[0] / straightPower.get(meta.fields[0].k0)
  : meta.fields[0].power[0]);
export const meanReflection = (results) => results.reduce((sum, r) => sum + measure(r), 0) / results.length;

// The structures a designer would try by hand, solved at the same wavelengths.
export function references(s) {
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
function coneFilter(n, radius) {
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

// Runs until the schedule ends or `onStep` returns false. `solve` returns a solver
// reply, or null if a newer request replaced it. Returns a reply holding any error.
export async function optimize(s, solve, onStep) {
  const d = s.design[s.analysis];
  const bend = s.analysis === 'bend2d';
  const n = pixelCount(d);
  const filter = bend ? coneFilter2d(BEND.pixels, BEND.pixels, d.radius) : coneFilter(n, Math.max(d.radius / (PIXEL * 1000), 1));
  const contrast = bend ? 0 : eps(d.high) - eps(d.low);
  const direction = s.analysis === 'mirror' ? -1 : 1;
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
      const gradient = new Float64Array(n);
      for (const { channels, meta } of reply.results) {
        // The bend's goal is more power through, so minimize its negative.
        const [slopes, scale] = bend
          ? [channels.gradient, -1 / straightPower.get(meta.fields[0].k0)]
          : [channels['dR/deps_re'], direction * contrast];
        for (let k = 0; k < n; k++) gradient[k] += (scale * slopes[k] * slope[k]) / reply.results.length;
      }
      x = mma.step(x, filter.adjoint(gradient));
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
