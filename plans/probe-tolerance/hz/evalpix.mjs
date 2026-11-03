// Reads a finished design of any pixel size at its own peak, on cells of the area over `cells`,
// so designs with different pixels are read on one mesh. Usage: node evalpix.mjs <tag> <seed> [cells]
import { readFileSync } from 'node:fs';
const S = process.cwd();
const g = await import(`${S}/pkg/glassworks_wasm.js`);
g.initSync({ module: readFileSync(`${S}/pkg/glassworks_wasm_bg.wasm`) });
const { defaultState } = await import(`${S}/rob3/web/js/model.js`);
const { cavityRequest, CAVITY } = await import(`${S}/rob3/web/js/cavity-design.js`);
const [tag, seed, cells] = [process.argv[2], +process.argv[3], +(process.argv[4] ?? 256)];
const d = JSON.parse(readFileSync(`des3/${tag}-${seed}.json`));
const n = Math.round(Math.sqrt(d.rho.length));
const s = { ...defaultState(), analysis: 'cavityDesign', wavelength: d.scale, polarization: 'tm' };
const rho = Float64Array.from(d.rho);
const F = (w, Q) => {
  const req = cavityRequest(s, rho, Q, w);
  req.cell = (req.cell * n) / cells;
  const x = g.solve_2d(JSON.stringify(req)); const f = JSON.parse(x.meta()).fields[0]; x.free();
  return (2 * f.ldos) / f.k0;
};
let peak = 1, best = -1;
for (let i = 0; i < 31; i++) { const w = 0.97 + 0.002 * i; const v = F(w, 300) * w; if (v > best) { best = v; peak = w; } }
for (const Q of [1000, 3000, 10000]) {
  const near = Array.from({ length: 9 }, (_, i) => peak * (1 + ((i - 4) * 0.5) / Q));
  const vals = near.map((w) => F(w, Q) * w);
  peak = near[vals.indexOf(Math.max(...vals))];
}
console.log(JSON.stringify({ tag, seed, n, cells, ownF: d.F, ownQ: d.Q, F: F(peak, CAVITY.Q), peak }));
