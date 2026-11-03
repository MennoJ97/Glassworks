// Getting from grey to sharp, with no physics in it: a filter, a projection and the method of
// moving asymptotes. Chapter 13 of the tutorial explains each one.

// Mulberry32: the same state always gives the same numbers, from 0 up to 1.
export function randomNumbers(state) {
  let a = state;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = Math.imul(a ^ (a >>> 15), a | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 2 ** 32;
  };
}

export function randomStart(n, seed) {
  const random = randomNumbers(Math.imul(seed, 0x9e3779b1) >>> 0);
  return Float64Array.from({ length: n }, () => 0.5 + 0.4 * (random() - 0.5));
}

// A cone filter averages each pixel with its neighbours, which removes single-pixel
// noise. `adjoint` carries gradients back through the same weights.
export const coneFilter = (n, radius) => coneFilter2d(n, 1, radius);

// Over a rectangle of pixels, stored by rows.
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

// Svanberg's method of moving asymptotes, for bounds only. Each step minimizes a
// convex approximation, which has a closed form when there are no other constraints.
export class Mma {
  // `move` is the largest change to any variable in one step.
  constructor(move) {
    this.move = move;
    this.history = [];
  }

  // A new `x` in [0, 1] that should reduce an objective with gradient `g` at `x`. `choose`, if
  // given, receives the step for any gradient on the same asymptotes, and returns the one it wants.
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
    // Variant MIX: a constraint's part `extra` mixes into each side's coefficient, as MMA's dual does.
    const stepFor = (gradient, extra) => Float64Array.from(x, (v, i) => {
      const e = extra ? extra[i] : 0;
      const p = (high[i] - v) * Math.sqrt(Math.max(gradient[i], 0) + Math.max(e, 0)) + 1e-7;
      const q = (v - low[i]) * Math.sqrt(Math.max(-gradient[i], 0) + Math.max(-e, 0)) + 1e-7;
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
  const at = process.env.MIX ? (lambda) => stepFor(goal, along.map((v) => lambda * v)) : (lambda) => stepFor(goal.map((v, k) => v + lambda * along[k]));
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
