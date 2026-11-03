// Optimizes an H_z cavity with the page's own code, and saves it.
import { readFileSync, writeFileSync } from 'node:fs';
const S = process.cwd(); const root = `${S}/rob3`; const seed = +process.argv[2]; const tag = process.argv[3];
const { one: native } = await import(`${S}/native.mjs`);

const { defaultState } = await import(`${root}/web/js/model.js`);
const { optimize } = await import(`${root}/web/js/design.js`);
const { tuned, cavityRequest, CAVITY } = await import(`${root}/web/js/cavity-design.js`);
const one = native;
const s = { ...defaultState(), analysis: 'cavityDesign', wavelength: 1, polarization: 'tm' }; s.design = { ...s.design, cavityDesign: { ...s.design.cavityDesign, seed } };
const t0 = Date.now(); let last, loop;
await optimize(s, async (rs) => ({ results: rs.map(one), ms: 0 }), (st) => { last = st; if (st.beta !== Infinity) loop = st; return true; });
const a = one(cavityRequest(s, last.rho, CAVITY.Q));
const rec = { loop: process.env.KEEP ? Array.from(loop.rho) : undefined, loopX: process.env.KEEP ? Array.from(loop.x) : undefined, loopBeta: loop.beta, seed, rho: Array.from(last.rho), scale: tuned.get(last.rho), Q: last.mode?.q, V: last.mode?.volume, w: last.mode?.wavelength, F: 2 * a.meta.fields[0].ldos / a.meta.fields[0].k0 };
writeFileSync(`des3/${tag}-${seed}.json`, JSON.stringify(rec));
console.log(JSON.stringify({ seed, Q: rec.Q, V: rec.V, w: rec.w, F: rec.F, scale: rec.scale, sec: (Date.now() - t0) / 1000 }));
