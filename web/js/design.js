// Density topology optimization of a layered structure. The design region is a row
// of thin pixels, each with a density that mixes two materials. The solver supplies
// the adjoint gradient; everything else here is optimization, not physics.

import { MATERIALS } from './model.js';

const PIXEL = 0.005;
const AMBIENT = 0.3;
const SHARPNESS = [1, 2, 4, 8, 16, 32, 64, 128];
const ITERATIONS_PER_STAGE = 40;
export const TOTAL_ITERATIONS = SHARPNESS.length * ITERATIONS_PER_STAGE;
export const STAGES = SHARPNESS.map((beta, i) => ({ frame: i * ITERATIONS_PER_STAGE, name: `β ${beta}` }));
export const sharpnessAt = (frame) => SHARPNESS[Math.floor(frame / ITERATIONS_PER_STAGE)] ?? Infinity;

const material = (name) => MATERIALS.find((m) => m.name === name);
const eps = (name) => material(name).n ** 2;
const solverLayer = (thickness, e) => ({ thickness, eps: [e, 0], mu: [1, 0] });

export const pixelCount = (d) => Math.max(1, Math.round(d.thickness / PIXEL));

// Goal wavelengths, with the centre last so the plotted field comes from it.
function goalWavelengths(s) {
  if (s.analysis === 'mirror') return [s.wavelength];
  const { from, to } = s.design.antireflection;
  const wavelengths = Array.from({ length: 11 }, (_, i) => from + ((to - from) * i) / 10);
  wavelengths.push(wavelengths.splice(5, 1)[0]);
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

export const gradientRequests = (s, rho) => requests(s, pixelLayers(s.design[s.analysis], rho),
  (wavelength) => ({ type: 'driven_gradient', wavelength, angle_deg: s.angle }));

export const meanReflection = (results) => results.reduce((sum, r) => sum + r.meta.fields[0].power[0], 0) / results.length;

// The structures a designer would try by hand, solved at the same wavelengths.
export function references(s) {
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

// Pushes densities toward 0 or 1. Larger `beta` is sharper.
function project(x, beta) {
  const half = Math.tanh(beta / 2);
  const rho = x.map((v) => (half + Math.tanh(beta * (v - 0.5))) / (2 * half));
  const slope = x.map((v) => beta / (2 * half * Math.cosh(beta * (v - 0.5)) ** 2));
  return { rho, slope };
}

// Svanberg's method of moving asymptotes, for bounds only. Each step minimizes a
// convex approximation, which has a closed form when there are no other constraints.
class Mma {
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
  const n = pixelCount(d);
  const filter = coneFilter(n, Math.max(d.radius / (PIXEL * 1000), 1));
  const contrast = eps(d.high) - eps(d.low);
  const direction = s.analysis === 'mirror' ? -1 : 1;
  let x = randomStart(n, d.seed);
  let iteration = 0;
  for (const beta of SHARPNESS) {
    // Projection amplifies a change in x by about beta, so steps shrink to match.
    const mma = new Mma(Math.min(0.1, 0.2 / beta));
    for (let i = 0; i < ITERATIONS_PER_STAGE; i++) {
      const { rho, slope } = project(filter.apply(x), beta);
      const reply = await solve(gradientRequests(s, rho));
      if (!reply || reply.error) return reply;
      if (!onStep({ rho, reply, iteration: ++iteration, beta })) return null;
      const gradient = new Float64Array(n);
      for (const { channels } of reply.results) {
        const dR = channels['dR/deps_re'];
        for (let k = 0; k < n; k++) gradient[k] += (direction * dR[k] * contrast * slope[k]) / reply.results.length;
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
