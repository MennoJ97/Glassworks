// More of the crystal's disorder draws, along x, to firm up its medians. One draw per process, with `k` from 7 on.
import { rng, gauss } from './pix.mjs';
import { F, crystal, crystalD } from './hz.mjs';
const [sigma, k] = process.argv.slice(2).map(Number); const D = crystalD, r0 = 0.4 * 0.2759, sub = 2, w0 = 0.9945529287;
const u = rng(7919 * k + 31 * 1000 + Math.round(sigma * 100)); const v = Array.from({ length: 200 }, () => r0 - (sigma / 1000) * gauss(u)); const s = crystal((i) => v[i]);
const half = 0.003 + 0.005 * sigma; const ws = Array.from({ length: 21 }, (_, i) => w0 * (1 + ((i - 10) / 10) * half)); const vs = ws.map((w) => F(D, s, sub, w, 0)); const best = ws[vs.indexOf(Math.max(...vs))]; const near = Array.from({ length: 9 }, (_, i) => best * (1 + ((i - 4) / 4) * (half / 10))); const vn = near.map((w) => F(D, s, sub, w, 0));
console.log(JSON.stringify({ id: 'crystal', sigma, draw: k, Fx: Math.max(...vs, ...vn) }));
