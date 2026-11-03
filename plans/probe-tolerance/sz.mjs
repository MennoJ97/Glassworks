import { readFileSync } from 'node:fs';
const [P, seed] = process.argv.slice(2); const root = `${process.cwd()}/size${P}`;
const g = await import(`${root}/web/pkg/glassworks_wasm.js`); g.initSync({ module: readFileSync(`${root}/web/pkg/glassworks_wasm_bg.wasm`) });
const { defaultState } = await import(`${root}/web/js/model.js`); const { optimize, gradientRequests } = await import(`${root}/web/js/design.js`);
const one = (r) => { const x = g.solve_2d(JSON.stringify(r)); const channels = {}; x.channel_names().forEach((n, i) => { channels[n] = x.channel(i); }); const meta = JSON.parse(x.meta()); x.free(); return { meta, channels }; };
const s = { ...defaultState(), analysis: 'cavityDesign', wavelength: 1 }; s.design = { ...s.design, cavityDesign: { ...s.design.cavityDesign, seed: +seed } };
let last; const t = Date.now(); await optimize(s, async (rs) => ({ results: rs.map(one), ms: 0 }), (st) => { last = st; return true; });
const q = gradientRequests(s, last.rho)[0]; const a = one(q);
console.log(JSON.stringify({ P: +P, seed: +seed, Q: last.mode?.q ?? 0, V: last.mode?.volume ?? 0, F: 2 * a.meta.fields[0].ldos / a.meta.fields[0].k0, s: (Date.now() - t) / 1000 }));
