import { readFileSync, writeFileSync } from 'node:fs';
const [tag, seedS] = process.argv.slice(2); const seed = +seedS; const root = `${process.cwd()}/rob`;
const g = await import(`${root}/web/pkg/glassworks_wasm.js`); g.initSync({ module: readFileSync(`${root}/web/pkg/glassworks_wasm_bg.wasm`) });
const { defaultState } = await import(`${root}/web/js/model.js`); const { optimize, designGeometry } = await import(`${root}/web/js/design.js`);
const solveOne = (r) => { const x = g.solve_2d(JSON.stringify(r)); const channels = {}; x.channel_names().forEach((n, i) => { channels[n] = x.channel(i); }); const meta = JSON.parse(x.meta()); x.free(); return { meta, channels }; };
const s = { ...defaultState(), analysis: 'cavityDesign', wavelength: 1 }; s.design = { ...s.design, cavityDesign: { ...s.design.cavityDesign, seed } };
let last; const t0 = Date.now(); await optimize(s, async (rs) => ({ results: rs.map(solveOne), ms: 0 }), (st) => { last = st; return true; });
const name = `${tag}-${seed}`; writeFileSync(`des/${name}.json`, JSON.stringify({ seed, rho: Array.from(last.rho), scale: designGeometry(s, last.rho).size / designGeometry(s).size, Q: last.mode?.q, V: last.mode?.volume }));
const { load, pixels, reso, exposed, rng, gauss } = await import('./pix.mjs'); const { one, PML } = await import('./lib.mjs');
const D = load(name), sub = 2; const emit = (o) => console.log(JSON.stringify({ tag, seed, ...o }));
const F = (shapes, w) => { const f = one({ width: D.size, height: D.size, cell: D.h / sub, background: { eps: 1, mu: 1 }, shapes, polarization: 'te', walls: 'absorbing', monitors: [], analysis: { type: 'driven', wavelength: w, source: { type: 'line', x: D.size / 2, y: D.size / 2 }, pml: PML * D.scale, average_q: 1000 } }).meta.fields[0]; return (2 * f.ldos) / f.k0; };
try {
  const b = reso(D, pixels(D, () => 0), sub, 1.0005); const up = reso(D, pixels(D, () => 0.0005), sub, b.w * 1.003), dn = reso(D, pixels(D, () => -0.0005), sub, b.w * 0.997);
  emit({ base: true, w: b.w, Q: b.Q, V: b.V, F: F(pixels(D, () => 0), b.w), slope: (up.w - dn.w) / b.w * 100, fill: D.des.rho.reduce((a, v) => a + v, 0) / D.des.rho.length, secs: (Date.now() - t0) / 1000 });
  for (const nm of [1, -1, 2, -2]) { const sh = pixels(D, () => nm / 1000); emit({ nm, F1: F(sh, b.w) }); }
  for (const sigma of [1, 2, 5]) for (let k = 1; k <= 6; k++) { const u = rng(7919 * k + 31 * seed + Math.round(sigma * 100)); const list = Array.from({ length: exposed(D) }, () => (sigma / 1000) * gauss(u)); const sh = pixels(D, (i) => list[i]); try { const r = reso(D, sh, sub, b.w); emit({ sigma, draw: k, w: r.w, Q: r.Q, Fown: F(sh, r.w), F1: F(sh, b.w) }); } catch (e) { emit({ sigma, draw: k, fail: 1, F1: F(sh, b.w) }); } }
} catch (e) { emit({ fail: String(e).slice(0, 60) }); }
