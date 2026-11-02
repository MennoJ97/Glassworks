// Disorder correlated over `L`, so that nearby sides move together, as they do when a process drifts. Same measures as tol.mjs.
import { load, rng, gauss } from './pix.mjs';
import { F, pole3 } from './hz.mjs';
const tag = process.argv[2]; const num = process.argv[3]; const id = tag === 'plain' ? num : `${tag}-${num}`; const L = +process.argv[4]; const draws = +process.argv[5]; const sub = 2;
const D = load(id); const { n, h, m, mat } = D; const flat = (f) => { const o = []; const over = 1e-6 * h; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { if (!mat(i, j)) continue; const [xl, xr, yb, yt] = [m + i * h, m + (i + 1) * h, m + j * h, m + (j + 1) * h]; const e = (a, b, x, y) => (mat(a, b) ? over : f(x, y)); o.push({ type: 'rect', x0: xl - e(i - 1, j, xl, (yb + yt) / 2), x1: xr + e(i + 1, j, xr, (yb + yt) / 2), y0: yb - e(i, j - 1, (xl + xr) / 2, yb), y1: yt + e(i, j + 1, (xl + xr) / 2, yt), material: { eps: 12, mu: 1 } }); } return o; };
const b = pole3(D, flat(() => 0), sub, 1.0045, 45, [100, 300, 1000, 3000, 1e4, 3e4]);
const b0 = F(D, flat(() => 0), sub, b.w, 0), b90 = F(D, flat(() => 0), sub, b.w, 90);
console.log(JSON.stringify({ tag, seed: +num, base: true, w: b.w, Q: b.Q, F0: b0, F90: b90 }));
const atF = (s, w, d) => F(D, s, sub, w, d);
const sweep = (s, c, half, d) => { const ws = Array.from({ length: 21 }, (_, i) => c * (1 + ((i - 10) / 10) * half)); const vs = ws.map((w) => atF(s, w, d)); const best = ws[vs.indexOf(Math.max(...vs))]; const near = Array.from({ length: 9 }, (_, i) => best * (1 + ((i - 4) / 4) * (half / 10))); const vn = near.map((w) => atF(s, w, d)); return Math.max(...vn) >= Math.max(...vs) ? { w: near[vn.indexOf(Math.max(...vn))], F: Math.max(...vn) } : { w: best, F: Math.max(...vs) }; };
for (const sigma of [2, 5]) for (let k = 1; k <= draws; k++) { const u = rng(104729 * k + 37 * +num + sigma); const K = 200; const waves = Array.from({ length: K }, () => [gauss(u) / L, gauss(u) / L, 2 * Math.PI * u()]); const field = (x, y) => (sigma / 1000) * Math.sqrt(2 / K) * waves.reduce((sum, [a, c, p]) => sum + Math.cos(a * x + c * y + p), 0);
  const s = flat(field); const half = 0.004 + 0.006 * sigma; const x = sweep(s, b.w, half, 0), y = { w: 0, F: 0 }; let px = {}; try { px = pole3(D, s, sub, x.w, 0, [1000, 3000, 1e4]); } catch (e) { /* no pole found */ }
  console.log(JSON.stringify({ tag, seed: +num, sigma, draw: k, L, Fx: x.F, wx: x.w, Fy: y.F, wy: y.w, Fx_orig: atF(s, b.w, 0), Fy_orig: atF(s, b.w, 90), Qx: px.Q, wpx: px.w })); }
console.log(JSON.stringify({ tag, seed: +num, finished: true }));
