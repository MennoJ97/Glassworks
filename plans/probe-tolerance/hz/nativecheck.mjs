// Native against WebAssembly on one request, channel by channel.
import { readFileSync } from 'node:fs';
import { one as native } from './native.mjs';
const g = await import(`${process.cwd()}/pkg/eigenmenno_wasm.js`); g.initSync({ module: readFileSync(`${process.cwd()}/pkg/eigenmenno_wasm_bg.wasm`) });
const N = 64, side = 64 / 24, margin = 0.25, size = side + 2 * margin; let r = 7;
const rho = Array.from({ length: N * N }, () => (r = (r * 16807) % 2147483647) / 2147483647);
const req = { width: size, height: size, cell: side / N, background: { eps: 1, mu: 1 }, shapes: [], polarization: 'tm', walls: 'absorbing', monitors: [], analysis: { type: 'driven', wavelength: 1, source: { type: 'dipole', x: size / 2, y: size / 2, angle_deg: 0 }, pml: 0.25, average_q: 1000 }, density: { x0: margin, y0: margin, width: side, height: side, nx: N, ny: N, eps_low: 1, eps_high: 12, rho } };
const w = g.solve_2d(JSON.stringify(req)); const wm = JSON.parse(w.meta()); const wg = w.channel(w.channel_names().indexOf('gradient')); w.free();
const t = Date.now(); const n = native(req); const sec = (Date.now() - t) / 1000;
const rel = Math.max(...wg.map((v, i) => Math.abs(v - n.channels.gradient[i]))) / Math.max(...wg.map(Math.abs));
console.log(JSON.stringify({ ldosWasm: wm.fields[0].ldos, ldosNative: n.meta.fields[0].ldos, gradientRelDiff: rel, nativeSec: sec }));
