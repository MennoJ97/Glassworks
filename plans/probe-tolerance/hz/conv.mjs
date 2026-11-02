// Mesh convergence: wavelength, Q, LDOS and the shift per nm against mesh, at pixel corners.
import { load, pixels } from './pix.mjs';
import { reso, F, peak } from './hz.mjs';
const seed = process.argv[2]; const subs = process.argv[3].split(',').map(Number);
const D = load(seed); const shapes = pixels(D, () => 0);
const c = D.n / 2; const atSource = D.mat(c, c) ? 'material' : 'air';
const out = { seed: +seed, atSource, scale: D.scale, Qpage: D.des.Q, Fpage: D.des.F, wpage: D.des.w };
for (const sub of subs) { const t0 = Date.now();
  try { const pk = peak(D, shapes, sub, 1, 0, 0.03); const b = reso(D, shapes, sub, pk.w); const b90 = reso(D, shapes, sub, pk.w, 90);
    const up = reso(D, pixels(D, () => 0.001), sub, b.w * 1.003), dn = reso(D, pixels(D, () => -0.001), sub, b.w * 0.997);
    console.log(JSON.stringify({ ...out, sub, pk, w: b.w, Q: b.Q, V: b.V, w90: b90.w, Q90: b90.Q, Fown: F(D, shapes, sub, b.w), F90: F(D, shapes, sub, b.w, 90), wup: up.w, wdn: dn.w, slope: (up.w - dn.w) / 2 / b.w * 100, sec: (Date.now() - t0) / 1000 })); }
  catch (e) { console.log(JSON.stringify({ ...out, sub, fail: String(e).slice(0, 80) })); } }
