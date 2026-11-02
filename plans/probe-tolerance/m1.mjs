import { one, PML } from './lib.mjs';
import { load, pixels, reso, exposed, rng, gauss } from './pix.mjs';
const [seed, job, a1, a2, a3] = process.argv.slice(2);
const emit = (o) => console.log(JSON.stringify({ seed: +seed, ...o }));
const D = load(seed);
const drive = (shapes, sub, w) => { const { size, h, scale } = D; const f = one({ width: size, height: size, cell: h / sub, background: { eps: 1, mu: 1 }, shapes, polarization: 'te', walls: 'absorbing', monitors: [], analysis: { type: 'driven', wavelength: w, source: { type: 'line', x: size / 2, y: size / 2 }, pml: PML * scale, average_q: 1000 } }).meta.fields[0]; return (2 * f.ldos) / f.k0; };
const base = (sub) => reso(D, pixels(D, () => 0), sub, 1.0005);
const measure = (shapes, sub, target, tag) => { try { const r = reso(D, shapes, sub, target); return { ...tag, sub, w: r.w, Q: r.Q, V: r.V, Fown: drive(shapes, sub, r.w), F1: drive(shapes, sub, base0.w) }; } catch (e) { return { ...tag, sub, fail: String(e).slice(0, 50) }; } };
const sub = +(a1 ?? 2);
const base0 = base(sub);
if (job === 'uniform') { emit(measure(pixels(D, () => 0), sub, base0.w, { nm: 0 })); for (const sign of [1, -1]) { let hist = [[0, base0.w]]; for (const m of [0.5, 1, 2, 5, 10]) { const nm = sign * m; const last = hist.at(-1), prev = hist.at(-2); const guess = prev ? last[1] + ((last[1] - prev[1]) / (last[0] - prev[0])) * (nm - last[0]) : base0.w * (1 + 0.012 * nm); const r = measure(pixels(D, () => nm / 1000), sub, guess, { nm, guess }); emit(r); if (r.w) hist.push([nm, r.w]); } } }
if (job === 'random') { const sigma = +a2, draws = +a3; emit({ base: true, sub, w: base0.w, Q: base0.Q, V: base0.V, Fown: drive(pixels(D, () => 0), sub, base0.w) }); for (let k = 1; k <= draws; k++) { const u = rng(7919 * k + 31 * +seed + Math.round(sigma * 100)); const n = exposed(D); const list = Array.from({ length: n }, () => (sigma / 1000) * gauss(u)); emit(measure(pixels(D, (i) => list[i]), sub, base0.w, { sigma, draw: k })); } }
