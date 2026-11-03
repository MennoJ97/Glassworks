// How much the finish costs: the last optimized design, grey at its edges, read at its own peak
// at Q 1000, against the finished design. Usage: node diag3.mjs <tag> <seed>
import { readFileSync } from 'node:fs';
const S = process.cwd();
const g = await import(`${S}/pkg/eigenmenno_wasm.js`);
g.initSync({ module: readFileSync(`${S}/pkg/eigenmenno_wasm_bg.wasm`) });
const { defaultState } = await import(`${S}/rob3/web/js/model.js`);
const { cavityRequest } = await import(`${S}/rob3/web/js/cavity-design.js`);
const one = (r) => { const x = g.solve_2d(JSON.stringify(r)); const meta = JSON.parse(x.meta()); x.free(); return meta.fields[0]; };
const [tag, seed] = [process.argv[2], +process.argv[3]];
const d = JSON.parse(readFileSync(`des3/${tag}-${seed}.json`));
const s = { ...defaultState(), analysis: 'cavityDesign', wavelength: 1, polarization: 'tm' };
const rho = Float64Array.from(d.loop);
const F = (w, Q) => { const f = one(cavityRequest(s, rho, Q, w)); return (2 * f.ldos) / f.k0; };
let peak = 1, best = -1;
for (let i = 0; i < 61; i++) { const w = 0.85 + 0.005 * i; const v = F(w, 100) * w; if (v > best) { best = v; peak = w; } }
for (const Q of [300, 1000, 3000]) {
  const near = Array.from({ length: 9 }, (_, i) => peak * (1 + ((i - 4) * 0.5) / Q));
  const vals = near.map((w) => F(w, Q) * w);
  peak = near[vals.indexOf(Math.max(...vals))];
}
const grey = rho.filter((v) => v > 0.02 && v < 0.98).length;
const mid = rho.filter((v) => v > 0.2 && v < 0.8).length;
const flips = d.rho.reduce((n, v, k) => n + ((v > 0.5) !== (rho[k] > 0.5) ? 1 : 0), 0);
console.log(JSON.stringify({ tag, seed, loopBeta: d.loopBeta, greyF: F(peak, 1000), greyPeak: peak, finalF: d.F, finalQ: d.Q, grey, mid, flips }));
