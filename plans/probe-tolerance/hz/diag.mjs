// Same design, same cell, four ways of drawing it: where does the LDOS difference come from?
import { readFileSync } from 'node:fs';
import { load, pixels } from './pix.mjs';
import { one, base, PML, MARGIN, SIDE } from './hz.mjs';
const seed = process.argv[2]; const D = load(seed); const w = D.des.w; const n = D.n, h = D.h, m = D.m;
const src = { type: 'dipole', x: D.size / 2, y: D.size / 2, angle_deg: 0 };
const drv = (r) => ({ ...base(D.size, h, r.shapes ?? [], { type: 'driven', wavelength: w, source: src, pml: PML * D.scale, average_q: 1000 }), ...(r.density ? { density: r.density } : {}) });
const rd = (r) => { const f = one(drv(r)).meta.fields[0]; return { F: 2 * f.ldos / f.k0, resp: f.response }; };
const merged = []; for (let j = 0; j < n; j++) for (let i = 0; i < n;) { if (!D.mat(i, j)) { i++; continue; } let e = i; while (e < n && D.mat(e, j)) e++; merged.push({ material: { eps: 12, mu: 1 }, type: 'rect', x0: m + i * h, y0: m + j * h, x1: m + e * h, y1: m + (j + 1) * h }); i = e; }
console.log('density', rd({ density: { x0: m, y0: m, width: n * h, height: n * h, nx: n, ny: n, eps_low: 1, eps_high: 12, rho: D.des.rho } }));
console.log('merged rects', rd({ shapes: merged }));
console.log('pixels ov=1e-6', rd({ shapes: pixels(D, () => 0) }));
console.log('around source:'); for (let j = n / 2 - 4; j < n / 2 + 4; j++) { let s = ''; for (let i = n / 2 - 4; i < n / 2 + 4; i++) s += D.mat(i, j) ? '#' : '.'; console.log(s); }
