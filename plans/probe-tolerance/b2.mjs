import { one } from './lib.mjs';
import { rng, gauss } from './pix.mjs';
const R = '/Users/mennojansen/source/eigenmenno';
const { defaultState } = await import(`${R}/web/js/model.js`);
const { optimize, gradientRequests, meanReflection, prepareDesign, bendScale } = await import(`${R}/web/js/design.js`);
const seed = +process.argv[2];
const solve = async (rs) => ({ results: rs.map(one), ms: 0 });
const s = { ...defaultState(), analysis: 'bend2d' }; s.design = { ...s.design, bend2d: { ...s.design.bend2d, seed } };
await prepareDesign(s, solve);
let last; await optimize(s, solve, (st) => { last = st; return true; });
const rho = last.rho; const reqs = gradientRequests(s, rho); const r0 = reqs[0]; const d0 = r0.density; const n = d0.nx, h = d0.width / n; const scale = bendScale(s);
const mat = (i, j) => i >= 0 && j >= 0 && i < n && j < n && rho[i + j * n] > 0.5;
const eps = d0.eps_high; const arms = r0.shapes;
// A side outside the region facing an arm stays joined to it.
const touchesArm = (i, j) => { const x = d0.x0 + (i + 0.5) * h, y = d0.y0 + (j + 0.5) * h; return arms.some((a) => a.outline ? false : false) || armAt(x, y); };
const armBox = arms.map((a) => { const o = JSON.parse(JSON.stringify(a)); return o; });
const bbox = armBox.map((a) => [a.x0, a.y0, a.x1, a.y1]);
function armAt(x, y) { return bbox.some(([p, q, u, v]) => x > p && x < u && y > q && y < v); }
const pixels = (offset) => { const over = 1e-6 * h; const out = []; let side = 0; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { if (!mat(i, j)) continue; const e = (a, b) => (mat(a, b) || touchesArm(a, b) ? over : offset(side++)); const [L, RR, B, T] = [e(i - 1, j), e(i + 1, j), e(i, j - 1), e(i, j + 1)]; out.push({ type: 'rect', x0: d0.x0 + i * h - L, x1: d0.x0 + (i + 1) * h + RR, y0: d0.y0 + j * h - B, y1: d0.y0 + (j + 1) * h + T, material: { eps, mu: 1 } }); } return out; };
const count = () => { let c = 0; pixels(() => { c++; return 0; }); return c; };
const T = (shapes, armGrow, finer = 1) => meanReflection(reqs.map((r) => { const { density, ...rest } = r; const a = armBox.map((b) => { const horiz = b.x1 - b.x0 > b.y1 - b.y0; return { ...b, x0: b.x0 - (horiz ? 0 : armGrow), x1: b.x1 + (horiz ? 0 : armGrow), y0: b.y0 - (horiz ? armGrow : 0), y1: b.y1 + (horiz ? armGrow : 0) }; }); return one({ ...rest, cell: rest.cell / finer, shapes: [...a, ...shapes] }); }));
const pct = (v) => (100 * v).toFixed(2);
console.log(`bend seed ${seed}: density ${pct(meanReflection(last.reply.results))}% | drawn ${pct(T(pixels(() => 0), 0))}%, cell ${r0.cell.toFixed(4)} µm, pixel ${(1000 * h).toFixed(1)} nm, scale ${scale}, sides ${count()}, wavelengths ${reqs.map((r) => r.analysis.wavelength.toFixed(3)).join(' ')}`);
for (const nm of [-20, -10, -5, -2, -1, 1, 2, 5, 10, 20]) console.log(`  uniform ${String(nm).padStart(3)} nm, arms fixed: ${pct(T(pixels(() => nm / 1000), 0))}%  | 1/2 cell: ${pct(T(pixels(() => nm / 1000), 0, 1.5))}%`);
