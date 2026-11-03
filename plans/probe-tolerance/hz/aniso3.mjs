// Does anisotropic averaging make an optimized grey design read as its drawn edges do? Reads one
// kept design four ways at its own peak, all on cells half a pixel wide, with the probe solver:
// grey with today's scalar, grey with the tensor, its 0.5 contour drawn twice as finely, and the
// page's finished pixels. Usage: node aniso3.mjs <tag> <seed> [readings]
import { readFileSync } from 'node:fs';
const S = process.cwd(), A = `${S}/../aniso/pkg`;
const g = await import(`${A}/glassworks_wasm.js`);
g.initSync({ module: readFileSync(`${A}/glassworks_wasm_bg.wasm`) });
const { defaultState } = await import(`${S}/rob3/web/js/model.js`);
const { cavityRequest, CAVITY, mirrored } = await import(`${S}/rob3/web/js/cavity-design.js`);
const { coneFilter2d, smoothProject } = await import(`${S}/rob3/web/js/optimizer.js`);
const [tag, seed] = [process.argv[2], +process.argv[3]];
const which = (process.argv[4] ?? 'greyS,greyT,drawn2,snap').split(',');
const d = JSON.parse(readFileSync(`des3/${tag}-${seed}.json`));
const n = CAVITY.pixels;
const filtered = coneFilter2d(n, n, CAVITY.radius).apply(mirrored(n)(Float64Array.from(d.loopX)));
const ssp = smoothProject(filtered, n, d.loopBeta).rho;
const at = (i, j) => filtered[Math.max(0, Math.min(n - 1, i)) + Math.max(0, Math.min(n - 1, j)) * n];
// The edge's normal from the filtered field's slope, as SSP itself reads it.
const normal = Array.from({ length: n * n }, (_, k) => {
  const [i, j] = [k % n, Math.floor(k / n)];
  const [gx, gy] = [(at(i + 1, j) - at(i - 1, j)) / 2, (at(i, j + 1) - at(i, j - 1)) / 2];
  return Math.hypot(gx, gy) > 1e-6 ? Math.atan2(gy, gx) : 99;
});
const upsampled = (m) => {
  const N = n * m;
  return Float64Array.from({ length: N * N }, (_, k) => {
    const [px, py] = [((k % N) + 0.5) / m - 0.5, (Math.floor(k / N) + 0.5) / m - 0.5];
    const [i, j] = [Math.floor(px), Math.floor(py)];
    const [s, t] = [px - i, py - j];
    const v = (1 - s) * (1 - t) * at(i, j) + s * (1 - t) * at(i + 1, j) + (1 - s) * t * at(i, j + 1) + s * t * at(i + 1, j + 1);
    return v > 0.5 ? 1 : 0;
  });
};
const s = { ...defaultState(), analysis: 'cavityDesign', wavelength: 1, polarization: 'tm' };
// `divide` is how many cells each pixel of `rho` gets, so that every reading's cells are half a pixel wide.
const read = (rho, divide, normals) => {
  const F = (w, Q) => {
    const req = cavityRequest(s, rho, Q, w);
    req.cell /= divide;
    if (normals) req.density.normal = normals;
    const x = g.solve_2d(JSON.stringify(req)); const f = JSON.parse(x.meta()).fields[0]; x.free();
    return (2 * f.ldos) / f.k0;
  };
  let peak = 1, best = -1;
  for (let i = 0; i < 61; i++) { const w = 0.85 + 0.005 * i; const v = F(w, 100) * w; if (v > best) { best = v; peak = w; } }
  for (const Q of [300, 1000, 3000]) {
    const near = Array.from({ length: 9 }, (_, i) => peak * (1 + ((i - 4) * 0.5) / Q));
    const vals = near.map((w) => F(w, Q) * w);
    peak = near[vals.indexOf(Math.max(...vals))];
  }
  return { F: F(peak, CAVITY.Q), peak };
};
const out = { tag, seed, pageF: d.F, sameAsLoop: Math.max(...ssp.map((v, k) => Math.abs(v - d.loop[k]))) };
for (const r of which) {
  if (r === 'greyS') out.greyS = read(ssp, 2);
  if (r === 'greyT') out.greyT = read(ssp, 2, normal);
  if (r === 'drawn2') out.drawn2 = read(upsampled(2), 1);
  // Its pixels are a quarter wide, and the mesh must line up with them, so its cells are too.
  if (r === 'drawn4') out.drawn4 = read(upsampled(4), 1);
  if (r === 'snap') out.snap = read(Float64Array.from(d.rho), 2);
}
console.log(JSON.stringify(out));
