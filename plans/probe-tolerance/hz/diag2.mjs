// Pixel drawing against merged rectangles: which nodes move, and how much does the overlap matter?
import { load } from './pix.mjs';
import { one, base, PML } from './hz.mjs';
const seed = process.argv[2]; const D = load(seed); const w = D.des.w; const { n, h, m, mat } = D;
const src = { type: 'dipole', x: D.size / 2, y: D.size / 2, angle_deg: 0 };
const run = (shapes, pol = 'tm') => one({ ...base(D.size, h, shapes, pol === 'tm' ? { type: 'driven', wavelength: w, source: src, pml: PML * D.scale, average_q: 1000 } : { type: 'driven', wavelength: w, source: { type: 'line', x: D.size / 2, y: D.size / 2 }, pml: PML * D.scale, average_q: 1000 }), polarization: pol });
const merged = []; for (let j = 0; j < n; j++) for (let i = 0; i < n;) { if (!mat(i, j)) { i++; continue; } let e = i; while (e < n && mat(e, j)) e++; merged.push({ material: { eps: 12, mu: 1 }, type: 'rect', x0: m + i * h, y0: m + j * h, x1: m + e * h, y1: m + (j + 1) * h }); i = e; }
const pix = (ov) => { const o = []; const over = ov * h; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { if (!mat(i, j)) continue; const e = (a, b) => (mat(a, b) ? over : 0); o.push({ type: 'rect', x0: m + i * h - e(i - 1, j), x1: m + (i + 1) * h + e(i + 1, j), y0: m + j * h - e(i, j - 1), y1: m + (j + 1) * h + e(i, j + 1), material: { eps: 12, mu: 1 } }); } return o; };
const a = run(merged); const Fa = 2 * a.meta.fields[0].ldos / a.meta.fields[0].k0;
console.log('merged', Fa, 'nodes', a.channels.x.length);
for (const ov of [0, 1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]) { const b = run(pix(ov)); const f = b.meta.fields[0]; let moved = 0, big = 0; for (let i = 0; i < a.channels.x.length; i++) { const d = Math.hypot(a.channels.x[i] - b.channels.x[i], a.channels.y[i] - b.channels.y[i]); if (d > 1e-12) moved++; big = Math.max(big, d); } console.log('ov', ov, 'F', (2 * f.ldos / f.k0).toFixed(3), 'moved', moved, 'largest/h', (big / h).toExponential(2), 'tris', b.channels.triangles?.length); }
const ea = run(merged, 'te'), eb = run(pix(1e-6), 'te'); console.log('E_z, same geometry: merged', 2 * ea.meta.fields[0].ldos / ea.meta.fields[0].k0, 'pixels 1e-6', 2 * eb.meta.fields[0].ldos / eb.meta.fields[0].k0);
