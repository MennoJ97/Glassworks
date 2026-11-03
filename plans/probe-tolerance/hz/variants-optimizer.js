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
// Subpixel-smoothed projection, Hammond et al. 2025: uses the filtered field and its gradient, so
// that an edge stays about a pixel wide at any sharpness, with a gradient that stays finite.
const fill = (d, R) => { if (d >= R) return 0; if (d <= -R) return 1; const t = d / R; return 0.5 - (15 / 16) * t + (5 / 8) * t ** 3 - (3 / 16) * t ** 5; };
export function sspValue(v, g, beta) {
  const R = 0.55, gm = Math.max(g, 1e-9), d = (0.5 - v) / gm;
  const P = (u) => (beta === Infinity ? (u > 0.5 ? 1 : u < 0.5 ? 0 : 0.5) : (Math.tanh(beta / 2) + Math.tanh(beta * (u - 0.5))) / (2 * Math.tanh(beta / 2)));
  if (Math.abs(d) >= R) return P(v);
  const F = fill(d, R);
  return (1 - F) * P(v - R * gm * fill(d, R)) + F * P(v + R * gm * fill(-d, R));
}
export function ssp(x, n, beta, second = false) {
  const at = (i, j) => x[Math.max(0, Math.min(n - 1, i)) + Math.max(0, Math.min(n - 1, j)) * n];
  const R = 0.55, h = 1e-6, N = n * n;
  const rho = new Float64Array(N), dv = new Float64Array(N), dG = new Float64Array(N);
  const parts = Array.from({ length: 5 }, () => new Float64Array(N));
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
    const k = i + j * n, v = x[k];
    const gx = (at(i + 1, j) - at(i - 1, j)) / 2, gy = (at(i, j + 1) - at(i, j - 1)) / 2;
    const hxx = second ? at(i + 1, j) - 2 * v + at(i - 1, j) : 0, hyy = second ? at(i, j + 1) - 2 * v + at(i, j - 1) : 0;
    const hxy = second ? (at(i + 1, j + 1) - at(i + 1, j - 1) - at(i - 1, j + 1) + at(i - 1, j - 1)) / 4 : 0;
    // Romano et al. 2026, Eqs. 11 and 14: the Hessian keeps this size from vanishing where edges merge.
    const G = Math.sqrt(gx * gx + gy * gy + R * R * (hxx * hxx + hyy * hyy + 2 * hxy * hxy));
    rho[k] = sspValue(v, G, beta);
    dv[k] = (sspValue(v + h, G, beta) - sspValue(v - h, G, beta)) / (2 * h);
    const low = Math.max(0, G - h);
    dG[k] = G > 1e-8 ? (sspValue(v, G + h, beta) - sspValue(v, low, beta)) / (G + h - low) / G : 0;
    parts[0][k] = gx; parts[1][k] = gy; parts[2][k] = R * R * hxx; parts[3][k] = R * R * hyy; parts[4][k] = 2 * R * R * hxy;
  }
  const back = (u) => {
    const out = new Float64Array(N);
    const add = (i, j, value) => { out[Math.max(0, Math.min(n - 1, i)) + Math.max(0, Math.min(n - 1, j)) * n] += value; };
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
      const k = i + j * n, a = u[k] * dG[k];
      out[k] += u[k] * dv[k];
      add(i + 1, j, (a * parts[0][k]) / 2); add(i - 1, j, -(a * parts[0][k]) / 2);
      add(i, j + 1, (a * parts[1][k]) / 2); add(i, j - 1, -(a * parts[1][k]) / 2);
      if (second) {
        add(i + 1, j, a * parts[2][k]); add(i - 1, j, a * parts[2][k]); out[k] -= 2 * a * parts[2][k];
        add(i, j + 1, a * parts[3][k]); add(i, j - 1, a * parts[3][k]); out[k] -= 2 * a * parts[3][k];
        add(i + 1, j + 1, (a * parts[4][k]) / 4); add(i + 1, j - 1, -(a * parts[4][k]) / 4); add(i - 1, j + 1, -(a * parts[4][k]) / 4); add(i - 1, j - 1, (a * parts[4][k]) / 4);
      }
    }
    return out;
  };
  return { rho, back };
}

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
  // `choose`, if given, receives a function from a gradient to its step, and returns the step it wants.
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
    const solveFor = (gradient) => {
      const next = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        const p = (high[i] - x[i]) * Math.sqrt(Math.max(gradient[i], 0)) + 1e-7;
        const q = (x[i] - low[i]) * Math.sqrt(Math.max(-gradient[i], 0)) + 1e-7;
        const best = (p * low[i] + q * high[i]) / (p + q);
        const lower = Math.max(0, x[i] - this.move, low[i] + 0.1 * (x[i] - low[i]));
        const upper = Math.min(1, x[i] + this.move, high[i] - 0.1 * (high[i] - x[i]));
        next[i] = Math.min(upper, Math.max(lower, best));
      }
      return next;
    };
    const next = choose ? choose(solveFor) : solveFor(g);
    [this.low, this.high] = [low, high];
    this.history.push(x);
    return next;
  }
}
