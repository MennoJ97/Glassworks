// Designing a cavity: LDOS as the goal, with its frequency held on target. Chapter 17 of the
// tutorial explains this, and it builds on the optimizer that chapter 13 explains.

import { scaled } from './model.js';
import { Mma, coneFilter2d, project, randomStart } from './optimizer.js';

// A square of pixels around an emitter, in units of the wavelength, as big as 7 x 7 rods.
// `stages` raises the bandwidth's Q to a level a sharpening design can follow, then sharpens.
export const CAVITY_Q = 1000;
export const CAVITY = {
  pixels: 64, margin: 0.25, pml: 0.25, eps: 12, radius: 2, iterations: 20,
  stages: [[10, 1], [20, 2], [50, 4], [100, 8], [200, 16], [300, 32], [500, 64], [500, 128]],
  // Measured at this Q once finished, and in the table.
  Q: CAVITY_Q,
};
const cavitySide = CAVITY.pixels / 24;
// How much a finished design was resized to put its resonance on the target, by its densities.
export const tuned = new WeakMap();

export function cavityGeometry(wavelength) {
  const { margin } = CAVITY;
  const size = cavitySide + 2 * margin;
  return scaled({
    size, shapes: [],
    mark: [size / 2, size / 2],
    region: { x0: margin, y0: margin, width: cavitySide, height: cavitySide, nx: CAVITY.pixels, ny: CAVITY.pixels },
  }, wavelength);
}

// An emitter's electric field: a line source for E_z, and for H_z, a dipole in the plane, since
// a line source there is a magnetic current.
const emitter = (s, g) => (s.polarization === 'tm'
  ? { type: 'dipole', x: g.mark[0], y: g.mark[1], angle_deg: 0 }
  : { type: 'line', x: g.mark[0], y: g.mark[1] });

export function cavityRequest(s, rho, Q, wavelength = s.wavelength) {
  const scale = s.wavelength * (tuned.get(rho) ?? 1);
  const g = cavityGeometry(scale);
  return {
    width: g.size,
    height: g.size,
    cell: (cavitySide * scale) / CAVITY.pixels,
    background: { eps: 1, mu: 1 },
    shapes: [],
    polarization: s.polarization,
    walls: 'absorbing',
    monitors: [],
    analysis: { type: 'driven', wavelength, source: emitter(s, g), pml: CAVITY.pml * scale, average_q: Q },
    density: { ...g.region, eps_low: 1, eps_high: CAVITY.eps, rho: Array.from(rho) },
  };
}

const times = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const over = (a, b) => {
  const d = b[0] * b[0] + b[1] * b[1];
  return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
};
const complexK = (wavelength, Q) => [(2 * Math.PI) / wavelength, Math.PI / (wavelength * Q)];

const plus = (a, b) => [a[0] + b[0], a[1] + b[1]];
const minus = (a, b) => [a[0] - b[0], a[1] - b[1]];

// `r = B + R / (k - p)` through three points, as `B` and `R`, and `p`. It rearranges to a straight
// line in the unknowns `(p, B, C)`, with `r k = r p + B k - C`, which Cramer's rule solves.
function fitPole(points) {
  const det3 = (m) => [[0, 1, 2], [1, 2, 0], [2, 0, 1]].reduce((sum, [i, j, l]) => plus(sum, times(m[0][i], minus(times(m[1][j], m[2][l]), times(m[1][l], m[2][j])))), [0, 0]);
  const rows = points.map(({ r, k }) => [r, k, [-1, 0]]);
  const rhs = points.map(({ r, k }) => times(r, k));
  const swap = (col) => rows.map((row, i) => row.map((v, j) => (j === col ? rhs[i] : v)));
  const d = det3(rows);
  const [p, B, C] = [0, 1, 2].map((col) => over(det3(swap(col)), d));
  return { pole: p, background: B, residue: minus(times(B, p), C) };
}

// Where a design's resonance is, and how each density moves it. Close to a pole, 1 / u_0 is
// nearly linear in k, so two solves place it. The solver's two gradient channels are the slope
// of `2 u_0`, which through `1 / u_0` is also the slope of the pole's own k. A dipole's response
// also holds a large background from its near field, so it takes three solves to fit that too.
export async function findResonance(s, rho, Q, solve, wavelength = s.wavelength) {
  const offset = 1 / (4 * Q) + 1e-3;
  const dipole = s.polarization === 'tm';
  const wavelengths = [wavelength, wavelength * (1 + offset), ...(dipole ? [wavelength * (1 - offset)] : [])];
  const reply = await solve(wavelengths.map((w) => cavityRequest(s, rho, Q, w)));
  if (!reply || reply.error) return { reply };
  const here = reply.results[0];
  const points = reply.results.map((result, i) => ({ r: result.meta.fields[0].response, k: complexK(wavelengths[i], Q) }));
  const { gradient, gradient_im: imaginary } = here.channels;
  const resonant = dipole ? fitPole(points) : null;
  // The pole's own slope is `R / r^2` times the response's, with `r` less any background.
  const r = dipole ? minus(points[0].r, resonant.background) : points[0].r;
  const g = over([1, 0], r);
  const slope = dipole
    ? resonant.residue
    : over(minus(points[1].k, points[0].k), minus(over([1, 0], points[1].r), g));
  const factor = times(slope, over([1, 0], times(r, r)));
  return {
    reply: { ...reply, results: [here] },
    poleK: dipole ? resonant.pole[0] : points[0].k[0] - times(g, slope)[0],
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
// LDOS is read at the mode's own wavelength: past one linewidth, raising it would lower Q.
export async function optimizeCavity(s, solve, onStep) {
  const n = CAVITY.pixels ** 2;
  const filter = coneFilter2d(CAVITY.pixels, CAVITY.pixels, CAVITY.radius);
  const symmetric = mirrored(CAVITY.pixels);
  const through = (g, slope) => symmetric(filter.adjoint(g.map((v, k) => v * slope[k])));
  const target = (2 * Math.PI) / s.wavelength;
  let x = randomStart(n, s.design.cavityDesign.seed);
  let iteration = 0;
  let own = s.wavelength;
  for (const [Q, beta] of CAVITY.stages) {
    const reach = Math.min(0.1, 0.2 / beta);
    const mma = new Mma(reach);
    for (let i = 0; i < CAVITY.iterations; i++) {
      const { rho, slope } = project(filter.apply(symmetric(x)), beta);
      const { reply, poleK, moves } = await findResonance(s, rho, Q, solve, own);
      if (!reply || reply.error) return reply;
      if (!onStep({ rho, reply, iteration: ++iteration, beta, bandwidth: Q })) return null;
      if (Math.abs((2 * Math.PI) / poleK - s.wavelength) < 0.04 * s.wavelength) own = (2 * Math.PI) / poleK;
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
  // Grey pixels matter more to H_z, whose field crosses their edges, so removing them moves its
  // mode further: nine starts were resized by up to 10.4%, where E_z needs 1.2%.
  const scan = s.polarization === 'tm'
    ? Array.from({ length: 61 }, (_, i) => s.wavelength * (0.85 + 0.005 * i))
    : Array.from({ length: 21 }, (_, i) => s.wavelength * (0.96 + 0.004 * i));
  const swept = await solve(scan.map((wavelength) => cavityRequest(s, rho, 100, wavelength)));
  if (!swept || swept.error) return swept;
  const strengths = swept.results.map((r, i) => r.meta.fields[0].ldos * scan[i]);
  let peak = scan[strengths.indexOf(Math.max(...strengths))];
  if (s.polarization === 'tm') {
    // A dipole's background makes the fit below less sure of a broad mode, so narrow in by LDOS.
    for (const Q of [300, 1000, 3000]) {
      const near = Array.from({ length: 9 }, (_, i) => peak * (1 + ((i - 4) * 0.5) / Q));
      const local = await solve(near.map((wavelength) => cavityRequest(s, rho, Q, wavelength)));
      if (!local || local.error) return local;
      const values = local.results.map((r, i) => r.meta.fields[0].ldos * near[i]);
      peak = near[values.indexOf(Math.max(...values))];
    }
  } else {
    for (const Q of [100, 300, 1000, 1000]) {
      const { reply, poleK } = await findResonance(s, rho, Q, solve, peak);
      if (!reply || reply.error) return reply;
      const next = (2 * Math.PI) / poleK;
      if (Math.abs(next - peak) < 0.02 * peak) peak = next;
    }
  }
  tuned.set(rho, s.wavelength / peak);
  const reply = await solve([cavityRequest(s, rho, CAVITY.Q)]);
  if (!reply || reply.error) return reply;
  // Checks it as a drawn structure, with the search for leaky modes, which shares nothing of this.
  const scale = tuned.get(rho) * s.wavelength;
  const g = cavityGeometry(scale);
  const found = await solve([{ width: g.size, height: g.size, cell: g.region.width / CAVITY.pixels, background: { eps: 1, mu: 1 }, shapes: cavityShapes(rho, scale), polarization: s.polarization, walls: 'absorbing', monitors: [], analysis: { type: 'resonance', target_k0: target, x: g.mark[0], y: g.mark[1], pml: CAVITY.pml * scale, dipole_deg: s.polarization === 'tm' ? 0 : null } }]);
  const field = found?.results?.[0]?.meta.fields[0];
  const mode = field && { wavelength: (2 * Math.PI) / field.k0, q: field.q, volume: field.volume / s.wavelength ** 2 };
  onStep({ rho, reply, iteration, beta: Infinity, mode, bandwidth: CAVITY.Q });
  return null;
}
