// Is a fraction of a pixel a fair model of a moved edge in H_z? Compares a grey edge's shift, the gradient's prediction and a drawn edge.
import { readFileSync } from 'node:fs';
import { load, pixels } from './pix.mjs';
import { one, pole3 } from './hz.mjs';
const root = `${process.cwd()}/src`;
const { defaultState } = await import(`${root}/web/js/model.js`);
const { tuned, findResonance, CAVITY } = await import(`${root}/web/js/cavity-design.js`);
const seed = process.argv[2]; const D = load(seed); const n = 64;
const s = { ...defaultState(), analysis: 'cavityDesign', wavelength: 1, polarization: 'tm' };
const solve = async (rs) => ({ results: rs.map(one), ms: 0 });
const edgeMoved = (rho, d) => { const out = Float64Array.from(rho); const edge = (a, b) => { const diff = rho[a] - rho[b]; if (Math.abs(diff) < 1e-6) return; const [hi, lo] = diff > 0 ? [a, b] : [b, a]; const amount = Math.abs(diff) * Math.abs(d); if (d > 0) out[lo] += amount; else out[hi] -= amount; }; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const a = i + j * n; if (i + 1 < n) edge(a, a + 1); if (j + 1 < n) edge(a, a + n); } return out.map((v) => Math.min(1, Math.max(0, v))); };
const rho0 = Float64Array.from(D.des.rho); tuned.set(rho0, D.scale);
const fit = async (rho, w0) => { let w = w0, r; for (const q of [300, 1000, 3000, 10000]) { r = await findResonance(s, rho, q, solve, w); const nx = 2 * Math.PI / r.poleK; if (Math.abs(nx - w) < 0.04) w = nx; } return { w, moves: r.moves }; };
const b = await fit(rho0, D.des.w);
for (const nm of [1, 2]) { const px = (nm / 1000) / D.h; const up = edgeMoved(rho0, px), dn = edgeMoved(rho0, -px); tuned.set(up, D.scale); tuned.set(dn, D.scale);
  const pred = (r) => { let dk = 0; for (let q = 0; q < r.length; q++) dk += b.moves[q] * (r[q] - rho0[q]); return (2 * Math.PI) / (2 * Math.PI / b.w + dk); };
  const u = await fit(up, pred(up)), v = await fit(dn, pred(dn));
  console.log(JSON.stringify({ seed: +seed, nm, grey_shift: (u.w - v.w) / 2 / nm / b.w * 100, grad_shift: (pred(up) - pred(dn)) / 2 / nm / b.w * 100 })); }
for (const sub of [1, 2]) { const a = pole3(D, pixels(D, () => 0), sub, sub === 1 ? 1 : 1.0045, 45); const u = pole3(D, pixels(D, () => 0.001), sub, a.w * 1.009, 45), v = pole3(D, pixels(D, () => -0.001), sub, a.w * 0.991, 45); console.log(JSON.stringify({ seed: +seed, sub, drawn_shift: (u.w - v.w) / 2 / a.w * 100 })); }
