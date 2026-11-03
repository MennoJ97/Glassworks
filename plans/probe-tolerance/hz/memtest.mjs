// One cavity LDOS solve on a random design of N x N pixels over the page's area, with cells a pixel wide.
import { readFileSync } from 'node:fs';
const g = await import(`${process.cwd()}/pkg/eigenmenno_wasm.js`);
const wasm = g.initSync({ module: readFileSync(`${process.cwd()}/pkg/eigenmenno_wasm_bg.wasm`) });
const N = +process.argv[2], side = 64 / 24, margin = 0.25, size = side + 2 * margin;
let r = 7; const rho = Array.from({ length: N * N }, () => ((r = (r * 16807) % 2147483647) / 2147483647 > 0.5 ? 1 : 0));
const req = { width: size, height: size, cell: side / N, background: { eps: 1, mu: 1 }, shapes: [], polarization: 'tm', walls: 'absorbing', monitors: [],
  analysis: { type: 'driven', wavelength: 1, source: { type: 'dipole', x: size / 2, y: size / 2, angle_deg: 0 }, pml: 0.25, average_q: 1000 },
  density: { x0: margin, y0: margin, width: side, height: side, nx: N, ny: N, eps_low: 1, eps_high: 12, rho } };
for (let k = 0; k < 2; k++) { const t = Date.now(); const x = g.solve_2d(JSON.stringify(req)); const m = JSON.parse(x.meta()); x.free();
  console.log(JSON.stringify({ N, dofs: m.n_dofs, sec: (Date.now() - t) / 1000, wasmMB: Math.round(wasm.memory.buffer.byteLength / 2 ** 20), rssMB: Math.round(process.memoryUsage().rss / 2 ** 20) })); }
