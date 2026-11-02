// One finished design: shift per nm on half cells, where its source sits, and its reading on finer cells if asked.
import { load, pixels } from './pix.mjs';
import { F, pole3 } from './hz.mjs';
const [id, finer] = process.argv.slice(2); const D = load(id); const s0 = pixels(D, () => 0); const c = D.n / 2;
const fill = D.des.rho.reduce((a, v) => a + (v > 0.5), 0) / D.des.rho.length;
const a = pole3(D, s0, 2, 1.0045, 45); const up = pole3(D, pixels(D, () => 0.001), 2, a.w * 1.009, 45), dn = pole3(D, pixels(D, () => -0.001), 2, a.w * 0.991, 45);
const out = { id, source: D.mat(c, c) ? 'material' : 'air', fill, w: a.w, Q: a.Q, F0: F(D, s0, 2, a.w, 0), F45: F(D, s0, 2, a.w, 45), shift: (up.w - dn.w) / 2 / a.w * 100, Qpage: D.des.Q, Fpage: D.des.F, Vpage: D.des.V };
if (finer) { const b = pole3(D, s0, 3, a.w * 1.001, 45); const b4 = pole3(D, s0, 4, b.w * 1.0005, 45); Object.assign(out, { w3: b.w, Q3: b.Q, F3: F(D, s0, 3, b.w, 0), w4: b4.w, Q4: b4.Q, F4: F(D, s0, 4, b4.w, 0) }); }
console.log(JSON.stringify(out));
