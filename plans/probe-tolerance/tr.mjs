import { readFileSync, writeFileSync } from 'node:fs';
const [tag, seedS] = process.argv.slice(2); const seed = +seedS; const root = `${process.cwd()}/trade`;
const g = await import(`${root}/web/pkg/glassworks_wasm.js`); g.initSync({ module: readFileSync(`${root}/web/pkg/glassworks_wasm_bg.wasm`) });
const { defaultState } = await import(`${root}/web/js/model.js`); const { optimize, designGeometry } = await import(`${root}/web/js/design.js`);
const solveOne = (r) => { const x = g.solve_2d(JSON.stringify(r)); const channels = {}; x.channel_names().forEach((n, i) => { channels[n] = x.channel(i); }); const meta = JSON.parse(x.meta()); x.free(); return { meta, channels }; };
const s = { ...defaultState(), analysis: 'cavityDesign', wavelength: 1 }; s.design = { ...s.design, cavityDesign: { ...s.design.cavityDesign, seed } };
let last; await optimize(s, async (rs) => ({ results: rs.map(solveOne), ms: 0 }), (st) => { last = st; return true; });
const name = `${tag}-${seed}`; writeFileSync(`des/${name}.json`, JSON.stringify({ seed, rho: Array.from(last.rho), scale: designGeometry(s, last.rho).size / designGeometry(s).size, Q: last.mode?.q, V: last.mode?.volume }));
const { load, pixels, reso } = await import('./pix.mjs'); const { one, PML } = await import('./lib.mjs');
const D = load(name); const sub = 2;
const F = (shapes, w) => { const f = one({ width: D.size, height: D.size, cell: D.h / sub, background: { eps: 1, mu: 1 }, shapes, polarization: 'te', walls: 'absorbing', monitors: [], analysis: { type: 'driven', wavelength: w, source: { type: 'line', x: D.size / 2, y: D.size / 2 }, pml: PML * D.scale, average_q: 1000 } }).meta.fields[0]; return (2 * f.ldos) / f.k0; };
try { const b = reso(D, pixels(D, () => 0), sub, 1.0005); const up = reso(D, pixels(D, () => 0.0005), sub, b.w * 1.006), dn = reso(D, pixels(D, () => -0.0005), sub, b.w * 0.994); const fill = last.rho.reduce((a, v) => a + v, 0) / last.rho.length;
  console.log(JSON.stringify({ tag, seed, w: b.w, Q: b.Q, V: b.V, F: F(pixels(D, () => 0), b.w), slope: (up.w - dn.w) / b.w * 100, Vup: up.V, Vdn: dn.V, fill, scale: D.scale })); } catch (e) { console.log(JSON.stringify({ tag, seed, fail: String(e).slice(0, 60) })); }
