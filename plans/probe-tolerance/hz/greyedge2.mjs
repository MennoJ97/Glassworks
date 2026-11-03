// greyedge.mjs again, with the anisotropic-averaging probe: a grey edge's shift per nm, with
// today's scalar and with the tensor, to set against the drawn shift measured before.
// Usage: node greyedge2.mjs <seed>
import { readFileSync } from 'node:fs';
import { load } from './pix.mjs';
const A = process.env.PKG ?? `${process.cwd()}/../aniso/pkg`;
const g = await import(`${A}/glassworks_wasm.js`);
g.initSync({ module: readFileSync(`${A}/glassworks_wasm_bg.wasm`) });
const one = (r) => { const x = g.solve_2d(JSON.stringify(r)); const channels = {}; x.channel_names().forEach((n, i) => { channels[n] = x.channel(i); }); const meta = JSON.parse(x.meta()); x.free(); return { meta, channels }; };
const root = `${process.cwd()}/src`;
const { defaultState } = await import(`${root}/web/js/model.js`);
const { tuned, findResonance } = await import(`${root}/web/js/cavity-design.js`);
const seed = process.argv[2]; const D = load(seed); const n = 64;
const s = { ...defaultState(), analysis: 'cavityDesign', wavelength: 1, polarization: 'tm' };
const rho0 = Float64Array.from(D.des.rho); tuned.set(rho0, D.scale);
// Each pixel beside an edge gets the edge's normal. Opposite edges cancel, leaving it isotropic.
const sum = Array.from({ length: n * n }, () => [0, 0]);
for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
  const a = i + j * n;
  if (i + 1 < n && rho0[a] !== rho0[a + 1]) { sum[a][0] += 1; sum[a + 1][0] += 1; }
  if (j + 1 < n && rho0[a] !== rho0[a + n]) { sum[a][1] += 1; sum[a + n][1] += 1; }
}
// Two edges at a corner give the diagonal; nn^T only needs the axis, so signs don't matter here.
const normals = sum.map(([x, y]) => (x || y ? Math.atan2(y, x) : NaN));
const edgeMoved = (rho, d) => { const out = Float64Array.from(rho); const edge = (a, b) => { const diff = rho[a] - rho[b]; if (Math.abs(diff) < 1e-6) return; const [hi, lo] = diff > 0 ? [a, b] : [b, a]; const amount = Math.abs(diff) * Math.abs(d); if (d > 0) out[lo] += amount; else out[hi] -= amount; }; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const a = i + j * n; if (i + 1 < n) edge(a, a + 1); if (j + 1 < n) edge(a, a + n); } return out.map((v) => Math.min(1, Math.max(0, v))); };
for (const tensor of [false, true]) {
  const solve = async (rs) => ({ results: rs.map((r) => one(tensor ? { ...r, density: { ...r.density, normal: normals.map((v) => (Number.isFinite(v) ? v : 99)) } } : r)), ms: 0 });
  const fit = async (rho, w0) => { let w = w0, r; for (const q of [300, 1000, 3000, 10000]) { r = await findResonance(s, rho, q, solve, w); const nx = (2 * Math.PI) / r.poleK; if (Math.abs(nx - w) < 0.04) w = nx; } return { w }; };
  const b = await fit(rho0, D.des.w);
  for (const nm of [1, 2]) {
    const px = nm / 1000 / D.h; const up = edgeMoved(rho0, px), dn = edgeMoved(rho0, -px); tuned.set(up, D.scale); tuned.set(dn, D.scale);
    const u = await fit(up, b.w * (1 + 0.008 * nm)), v = await fit(dn, b.w * (1 - 0.008 * nm));
    console.log(JSON.stringify({ seed: +seed, tensor, nm, base: b.w, grey_shift: ((u.w - v.w) / 2 / nm / b.w) * 100 }));
  }
}
