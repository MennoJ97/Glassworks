import { one } from './lib.mjs';
import { rng, gauss } from './pix.mjs';
// The crystal: rods with one removed, whose radius changes uniformly, or randomly rod by rod.
const a = 0.3781, e = 0.25; const size = 7 * a + 2 * e; const cell = a / +(process.argv[2] ?? 24);
const rods = (r) => { const o = []; for (let k = 0; k < 49; k++) if (k !== 24) o.push({ type: 'circle', x: e + (k % 7 + 0.5) * a, y: e + (Math.floor(k / 7) + 0.5) * a, r: r(k), material: { eps: 12, mu: 1 } }); return o; };
const req = (shapes, analysis) => ({ width: size, height: size, cell, background: { eps: 1, mu: 1 }, shapes, polarization: 'te', walls: 'absorbing', monitors: [], analysis });
const reso = (shapes, t) => { const f = one(req(shapes, { type: 'resonance', target_k0: (2 * Math.PI) / t, x: size / 2, y: size / 2, pml: 0.25 })).meta.fields[0]; return { w: (2 * Math.PI) / f.k0, Q: f.q, V: f.volume }; };
const F = (shapes, w) => { const f = one(req(shapes, { type: 'driven', wavelength: w, source: { type: 'line', x: size / 2, y: size / 2 }, pml: 0.25, average_q: 1000 })).meta.fields[0]; return (2 * f.ldos) / f.k0; };
const r0 = 0.2 * a; const b = reso(rods(() => r0), 1); console.log(JSON.stringify({ crystal: true, nm: 0, cell: a / cell, ...b, Fown: F(rods(() => r0), b.w), F1: F(rods(() => r0), 1) }));
for (const nm of [-10, -5, -2, -1, 1, 2, 5, 10]) { const s = rods(() => r0 + nm / 1000); const m = reso(s, b.w * (1 + 0.003 * nm)); console.log(JSON.stringify({ crystal: true, nm, ...m, Fown: F(s, m.w), F1: F(s, b.w) })); }
for (const sigma of [1, 2, 5]) for (let k = 1; k <= 6; k++) { const u = rng(101 * k + sigma); const vals = Array.from({ length: 49 }, () => r0 + (sigma / 1000) * gauss(u)); const s = rods((i) => vals[i]); const m = reso(s, b.w); console.log(JSON.stringify({ crystal: true, sigma, draw: k, ...m, Fown: F(s, m.w), F1: F(s, b.w) })); }
