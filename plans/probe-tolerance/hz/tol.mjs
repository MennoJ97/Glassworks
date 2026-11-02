// One structure at half-cell mesh: a uniform edge move, then random disorder, by driven sweeps and poles at both dipole directions.
import { load, pixels, exposed, rng, gauss } from './pix.mjs';
import { F, pole3, peak, merged, crystal, crystalD, one } from './hz.mjs';
const id = process.argv[2]; const sub = 2; const isCrystal = id === 'crystal'; const draws = +(process.argv[3] ?? 6);
const r0 = 0.4 * 0.2759;
const D = isCrystal ? crystalD : load(id);
const flat = isCrystal ? crystal() : pixels(D, () => 0);
const uniform = (nm) => (isCrystal ? crystal(() => r0 - nm / 1000) : pixels(D, () => nm / 1000));
const emit = (o) => console.log(JSON.stringify({ id, ...o }));
const deg0 = isCrystal ? 0 : 45;
const g = isCrystal ? peak(D, flat, sub, 1, 0, 0.04).w : 1.0045;
const b = pole3(D, flat, sub, g, deg0, [100, 300, 1000, 3000, 1e4, 3e4]);
const b90 = pole3(D, flat, sub, isCrystal ? peak(D, flat, sub, b.w, 90, 0.04).w : b.w, isCrystal ? 90 : 135, [100, 300, 1000, 3000, 1e4, 3e4]);
const atF = (s, w, d) => F(D, s, sub, w, d);
emit({ base: true, w: b.w, Q: b.Q, w90: b90.w, Q90: b90.Q, F0: atF(flat, b.w, 0), F90: atF(flat, b.w, 90), Fmode: atF(flat, b.w, deg0) });
// A uniform move: the shift, and what's left at the wavelength that was intended.
let last = [[0, b.w]];
for (const sign of [1, -1]) { last = [[0, b.w]]; for (const m of [0.5, 1, 2, 5]) { const nm = sign * m; const [p, q] = [last.at(-1), last.at(-2)]; const guess = q ? p[1] + ((p[1] - q[1]) / (p[0] - q[0])) * (nm - p[0]) : b.w * (1 + 0.009 * nm); const s = uniform(nm); try { const r = pole3(D, s, sub, guess, deg0, [300, 1000, 3000, 1e4, 3e4]); last.push([nm, r.w]); emit({ nm, w: r.w, Q: r.Q, Fown: atF(s, r.w, deg0), Forig: atF(s, b.w, deg0) }); } catch (e) { emit({ nm, fail: String(e).slice(0, 60) }); } } }
// Disorder: a dipole's LDOS at each direction, at the intended wavelength and at its own best one.
const sweep = (s, c, half, d) => { const ws = Array.from({ length: 21 }, (_, i) => c * (1 + ((i - 10) / 10) * half)); const vs = ws.map((w) => atF(s, w, d)); let best = ws[vs.indexOf(Math.max(...vs))]; const step = half / 10; const near = Array.from({ length: 9 }, (_, i) => best * (1 + ((i - 4) / 4) * step)); const vn = near.map((w) => atF(s, w, d)); const k = vn.indexOf(Math.max(...vn)); return vn[k] >= Math.max(...vs) ? { w: near[k], F: vn[k] } : { w: best, F: Math.max(...vs) }; };
const nExposed = isCrystal ? 0 : exposed(D);
for (const sigma of [1, 2, 5]) for (let k = 1; k <= draws; k++) { const u = rng(7919 * k + 31 * (isCrystal ? 1000 : +id) + Math.round(sigma * 100));
  let s; if (isCrystal) { const v = Array.from({ length: 200 }, () => r0 - (sigma / 1000) * gauss(u)); s = crystal((i) => v[i]); } else { const list = Array.from({ length: nExposed }, () => (sigma / 1000) * gauss(u)); s = pixels(D, (i) => list[i]); }
  const half = 0.003 + 0.005 * sigma; const x = sweep(s, b.w, half, 0), y = sweep(s, b.w, half, 90);
  let px, py; try { px = pole3(D, s, sub, x.w, 0, [1000, 3000, 1e4]); } catch (e) { px = {}; } try { py = pole3(D, s, sub, y.w, 90, [1000, 3000, 1e4]); } catch (e) { py = {}; }
  emit({ sigma, draw: k, Fx: x.F, wx: x.w, Fy: y.F, wy: y.w, Fx_orig: atF(s, b.w, 0), Fy_orig: atF(s, b.w, 90), wpx: px.w, Qx: px.Q, wpy: py.w, Qy: py.Q }); }
emit({ finished: true });
