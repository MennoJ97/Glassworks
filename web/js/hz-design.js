// What chapter 19 adds to a cavity's optimizer for H_z, kept apart so that the optimizer of
// chapter 13 and the cavity of chapter 17 read without it: a smoothed projection, a step that
// holds the frequency inside MMA, and edges moved at random while optimizing.

import { Mma, randomNumbers } from './optimizer.js';

// Subpixel-smoothed projection, from Hammond et al. 2025, in its second-order form from Romano
// et al. 2026. Each pixel of a filtered field `x`, of `n` by `n`, finds from its value and slope
// how far an edge at 0.5 is, and mixes the projections either side of it by how much of the pixel
// it covers. So an edge stays about a pixel wide at any sharpness, and its gradient stays finite.
// The curvature keeps that width from vanishing where two edges meet and the slope is zero.
// `back` carries a gradient through it, as the transpose of its derivative.
const SMOOTHING = 0.55;
const covered = (d) => {
  if (d >= SMOOTHING) return 0;
  if (d <= -SMOOTHING) return 1;
  const t = d / SMOOTHING;
  return 0.5 - (15 / 16) * t + (5 / 8) * t ** 3 - (3 / 16) * t ** 5;
};
function smoothed(v, size, beta) {
  const sharp = (u) => (beta === Infinity
    ? (u > 0.5 ? 1 : u < 0.5 ? 0 : 0.5)
    : (Math.tanh(beta / 2) + Math.tanh(beta * (u - 0.5))) / (2 * Math.tanh(beta / 2)));
  const d = (0.5 - v) / Math.max(size, 1e-9);
  if (Math.abs(d) >= SMOOTHING) return sharp(v);
  const F = covered(d);
  return (1 - F) * sharp(v - SMOOTHING * size * F) + F * sharp(v + SMOOTHING * size * covered(-d));
}

export function smoothProject(x, n, beta) {
  const R = SMOOTHING, h = 1e-6, N = n * n;
  const index = (i, j) => Math.max(0, Math.min(n - 1, i)) + Math.max(0, Math.min(n - 1, j)) * n;
  const at = (i, j) => x[index(i, j)];
  const rho = new Float64Array(N), byValue = new Float64Array(N), bySize = new Float64Array(N);
  // Each pixel's slope and curvature, scaled by how they enter its size.
  const parts = Array.from({ length: 5 }, () => new Float64Array(N));
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < n; i++) {
      const k = i + j * n, v = x[k];
      const gx = (at(i + 1, j) - at(i - 1, j)) / 2, gy = (at(i, j + 1) - at(i, j - 1)) / 2;
      const hxx = at(i + 1, j) - 2 * v + at(i - 1, j), hyy = at(i, j + 1) - 2 * v + at(i, j - 1);
      const hxy = (at(i + 1, j + 1) - at(i + 1, j - 1) - at(i - 1, j + 1) + at(i - 1, j - 1)) / 4;
      const size = Math.sqrt(gx * gx + gy * gy + R * R * (hxx * hxx + hyy * hyy + 2 * hxy * hxy));
      rho[k] = smoothed(v, size, beta);
      byValue[k] = (smoothed(v + h, size, beta) - smoothed(v - h, size, beta)) / (2 * h);
      const low = Math.max(0, size - h);
      // Divided by the size, since each part's own slope of it is that part over the size.
      bySize[k] = size > 1e-8 ? (smoothed(v, size + h, beta) - smoothed(v, low, beta)) / (size + h - low) / size : 0;
      [gx, gy, R * R * hxx, R * R * hyy, 2 * R * R * hxy].forEach((p, m) => { parts[m][k] = p; });
    }
  }
  const back = (u) => {
    const out = new Float64Array(N);
    const add = (i, j, value) => { out[index(i, j)] += value; };
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        const k = i + j * n, a = u[k] * bySize[k];
        const [gx, gy, xx, yy, xy] = parts.map((p) => a * p[k]);
        // In this order, a run repeats the measured designs exactly: round-off changes them.
        out[k] += u[k] * byValue[k];
        add(i + 1, j, gx / 2); add(i - 1, j, -gx / 2);
        add(i, j + 1, gy / 2); add(i, j - 1, -gy / 2);
        add(i + 1, j, xx); add(i - 1, j, xx); out[k] -= 2 * xx;
        add(i, j + 1, yy); add(i, j - 1, yy); out[k] -= 2 * yy;
        add(i + 1, j + 1, xy / 4); add(i + 1, j - 1, -xy / 4); add(i - 1, j + 1, -xy / 4); add(i - 1, j - 1, xy / 4);
      }
    }
    return out;
  };
  return { rho, back };
}

// MMA's step from optimizer.js, with its subproblem opened up: `choose`, if given, receives the
// step for any gradient on the same asymptotes, and returns the one it wants. Without `choose` it
// steps exactly as `Mma` does, which a test checks, so keep the two alike.
export class HeldMma extends Mma {
  step(x, g, choose) {
    const n = x.length;
    const [low, high] = [new Float64Array(n), new Float64Array(n)];
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
    }
    const stepFor = (gradient) => Float64Array.from(x, (v, i) => {
      const p = (high[i] - v) * Math.sqrt(Math.max(gradient[i], 0)) + 1e-7;
      const q = (v - low[i]) * Math.sqrt(Math.max(-gradient[i], 0)) + 1e-7;
      const best = (p * low[i] + q * high[i]) / (p + q);
      const lower = Math.max(0, v - this.move, low[i] + 0.1 * (v - low[i]));
      const upper = Math.min(1, v + this.move, high[i] - 0.1 * (high[i] - v));
      return Math.min(upper, Math.max(lower, best));
    });
    const next = choose ? choose(stepFor) : stepFor(g);
    [this.low, this.high] = [low, high];
    this.history.push(x);
    return next;
  }
}

// MMA's step for `goal`, holding one constraint: `shift(step)` should be `want`. Adding `lambda`
// times the constraint's gradient `along` to the goal moves the step along the constraint, so
// a bisection on `lambda` finds it. This is what MMA's dual does for a constraint, and each try is
// one closed-form step. Where the bounds keep it out of reach, it takes the closest end.
export function holdingStep(stepFor, goal, along, shift, want) {
  const largest = (v) => v.reduce((m, a) => Math.max(m, Math.abs(a)), 0);
  const at = (lambda) => stepFor(goal.map((v, k) => v + lambda * along[k]));
  const error = (lambda) => shift(at(lambda)) - want;
  let span = largest(goal) / largest(along);
  let [a, b] = [-span, span];
  let [ea, eb] = [error(a), error(b)];
  for (let t = 0; t < 16 && ea * eb > 0; t++) {
    span *= 4;
    [a, b] = [-span, span];
    [ea, eb] = [error(a), error(b)];
  }
  if (ea * eb > 0) return at(Math.abs(ea) < Math.abs(eb) ? a : b);
  for (let t = 0; t < 30; t++) {
    const m = (a + b) / 2, em = error(m);
    if (ea * em <= 0) b = m;
    else [a, ea] = [m, em];
  }
  return at((a + b) / 2);
}

// What each H_z iteration samples: edges moved by this standard deviation, correlated over this
// distance, both in µm. A fraction of a pixel moves H_z's mode about twice as far as a drawn edge
// does, so `grey` scales it. Chapter 19 gives the measurements.
export const DISORDER = { samples: 3, sigma: 0.002, correlation: 0.1, grey: 0.5, waves: 100 };

// A smooth random field, in pixels of a design with pixels of `h`, from a sum of random waves.
export function edgeField(random, h) {
  const gaussian = () => Math.sqrt(-2 * Math.log(random() + 1e-300)) * Math.cos(2 * Math.PI * random());
  const { sigma, correlation, grey, waves: count } = DISORDER;
  const waves = Array.from({ length: count }, () => [gaussian() / correlation, gaussian() / correlation, 2 * Math.PI * random()]);
  const size = ((grey * sigma) / h) * Math.sqrt(2 / count);
  return (x, y) => size * waves.reduce((sum, [a, b, phase]) => sum + Math.cos(a * x * h + b * y * h + phase), 0);
}

// Moves each edge between pixels by `move(x, y)` of a pixel, with positive growing material. It
// follows each edge's contrast, so a uniform region stays as it is.
export function movedEdges(rho, n, move) {
  const out = Float64Array.from(rho);
  const edge = (a, b, d) => {
    const contrast = Math.abs(rho[a] - rho[b]);
    if (contrast < 1e-6) return;
    const [high, low] = rho[a] > rho[b] ? [a, b] : [b, a];
    if (d > 0) out[low] += contrast * d;
    else out[high] += contrast * d;
  };
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < n; i++) {
      const a = i + j * n;
      if (i + 1 < n) edge(a, a + 1, move(i + 1, j + 0.5));
      if (j + 1 < n) edge(a, a + n, move(i + 0.5, j + 1));
    }
  }
  return out.map((v) => Math.min(1, Math.max(0, v)));
}

// The slope of the LDOS averaged over disordered copies, each read at the wavelength that `moves`
// predicts for it, since a copy read on its peak's wing would trade Q for LDOS. `request` builds a
// cavity's solve for a design, and `h` is a pixel's width, of `n` by `n`.
export async function disorderedSlope(s, rho, { Q, poleK, moves, own, iteration }, solve, { request, h, n }) {
  const random = randomNumbers(iteration * 7907 + s.design.cavityDesign.seed);
  const copies = Array.from({ length: DISORDER.samples }, () => movedEdges(rho, n, edgeField(random, h)));
  const wavelengthOf = (copy) => {
    const k = copy.reduce((sum, v, q) => sum + moves[q] * (v - rho[q]), poleK);
    const wavelength = (2 * Math.PI) / k;
    return Math.abs(wavelength - s.wavelength) < 0.05 * s.wavelength ? wavelength : own;
  };
  const reply = await solve(copies.map((copy) => request(s, copy, Q, wavelengthOf(copy))));
  if (!reply || reply.error) return { reply };
  const mean = new Float64Array(rho.length);
  reply.results.forEach((r) => r.channels.gradient.forEach((v, k) => { mean[k] += v / reply.results.length; }));
  return { mean };
}
