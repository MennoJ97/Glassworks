import { one, PML } from './lib.mjs';
import { load, pixels } from './pix.mjs';
const name = process.argv[2]; const D = load(name); const sub = 2; const { size, h, m, n, mat, scale } = D;
const r = one({ width: size, height: size, cell: h / sub, background: { eps: 1, mu: 1 }, shapes: pixels(D, () => 0), polarization: 'te', walls: 'absorbing', monitors: [], analysis: { type: 'resonance', target_k0: 2 * Math.PI, x: size / 2, y: size / 2, pml: PML * scale } });
const f = r.meta.fields[0], c = r.channels, re = c['field0/Ez/re'], im = c['field0/Ez/im'], x = c.x, y = c.y, tri = c.triangles;
const inMat = (X, Y) => mat(Math.floor((X - m) / h), Math.floor((Y - m) / h));
let air = 0, die = 0; for (let t = 0; t < tri.length; t += 3) { const [a, b, q] = [tri[t], tri[t + 1], tri[t + 2]]; const gx = (x[a] + x[b] + x[q]) / 3, gy = (y[a] + y[b] + y[q]) / 3; if (gx < 0 || gy < 0 || gx > size || gy > size) continue; const A = Math.abs((x[b] - x[a]) * (y[q] - y[a]) - (x[q] - x[a]) * (y[b] - y[a])) / 2; const e2 = [a, b, q].reduce((s, v) => s + re[v] * re[v] + im[v] * im[v], 0) / 3; if (inMat(gx, gy)) die += 12 * e2 * A; else air += e2 * A; }
const cx = size / 2; let peak = 0, pi = 0, src = 0, bestSrc = 1e9; for (let i = 0; i < x.length; i++) { const e = re[i] * re[i] + im[i] * im[i]; if (e > peak) { peak = e; pi = i; } const d = Math.hypot(x[i] - cx, y[i] - cx); if (d < bestSrc) { bestSrc = d; src = e; } }
// Second-brightest lobe, away from the centre, and the field at the nearest material or air boundary.
let far = 0, farR = 0; for (let i = 0; i < x.length; i++) { const d = Math.hypot(x[i] - cx, y[i] - cx); if (d > 0.12 * scale) { const e = re[i] * re[i] + im[i] * im[i]; if (e > far) { far = e; farR = d; } } }
let bound = 0; const pts = []; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { if (!mat(i, j)) continue; for (const [a, b, sx, sy] of [[i + 1, j, 1, .5], [i - 1, j, 0, .5], [i, j + 1, .5, 1], [i, j - 1, .5, 0]]) if (!mat(a, b)) pts.push([m + (i + sx) * h, m + (j + sy) * h]); }
const near = Math.min(...pts.map((p) => Math.hypot(p[0] - cx, p[1] - cx)));
const grid = new Map(); for (let i = 0; i < x.length; i++) grid.set(`${Math.round(x[i] / (h / sub / 2))},${Math.round(y[i] / (h / sub / 2))}`, i);
for (const p of pts) { const i = grid.get(`${Math.round(p[0] / (h / sub / 2))},${Math.round(p[1] / (h / sub / 2))}`); if (i !== undefined) bound += (re[i] * re[i] + im[i] * im[i]) * h; }
const total = air + die; const sh = 16; const rows = []; for (let j = n / 2 + sh - 1; j >= n / 2 - sh; j--) { let s = ''; for (let i = n / 2 - sh; i < n / 2 + sh; i++) s += mat(i, j) ? '#' : '.'; rows.push(s); }
console.log(JSON.stringify({ name, Q: f.q, V: f.volume, srcE2: src, peakE2: peak, peakR: Math.hypot(x[pi] - cx, y[pi] - cx), peakInMat: inMat(x[pi] + 1e-6, y[pi] + 1e-6), srcPixels: [[n / 2 - 1, n / 2 - 1], [n / 2, n / 2 - 1], [n / 2 - 1, n / 2], [n / 2, n / 2]].map(([i, j]) => mat(i, j) ? 1 : 0).join(''), energyInMat: die / total, E2perEnergy: src / total, nearestEdge: near / scale * 1000, farLobe: far / peak, farR, boundary: bound / total, rows }));
