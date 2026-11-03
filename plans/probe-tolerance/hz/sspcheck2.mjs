// SSP's and SSP2's transposes against finite differences, and how each behaves where two edges merge.
const { ssp, randomStart, coneFilter2d } = await import(`${process.cwd()}/rob2/web/js/optimizer.js`);
const n = 64; const f = coneFilter2d(n, n, 2); const x = f.apply(randomStart(n * n, 3).map((v) => 0.5 + 3 * (v - 0.5)));
for (const second of [false, true]) for (const beta of [8, 128, Infinity]) {
  const { rho, back } = ssp(x, n, beta, second); const u = rho.map((_, k) => Math.sin(k * 0.37)); const g = back(u);
  let worst = 0, checked = 0; const J = (a) => ssp(a, n, beta, second).rho.reduce((s, v, k) => s + v * u[k], 0);
  for (let t = 0; t < 400; t++) { const k = 100 + t * 9; const d = 1e-5; const p = x.slice(), m = x.slice(); p[k] += d; m[k] -= d; const fd = (J(p) - J(m)) / (2 * d); if (Math.abs(fd) > 1e-6 || Math.abs(g[k]) > 1e-6) { checked++; worst = Math.max(worst, Math.abs(fd - g[k]) / Math.max(Math.abs(fd), 1e-3)); } }
  console.log(JSON.stringify({ second, beta: String(beta), checked, worst, grey: rho.filter((v) => v > 0.01 && v < 0.99).length / rho.length }));
}
// A parabola's minimum crossing the level set: the middle pixel's density as two edges merge.
for (const second of [false, true]) { const row = []; for (const a of [0.47, 0.49, 0.499, 0.501, 0.51, 0.53]) { const y = new Float64Array(n * n); for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) y[i + j * n] = a + 0.02 * (i - 32) ** 2; row.push(+ssp(y, n, Infinity, second).rho[32 + 32 * n].toFixed(3)); } console.log(second ? 'SSP2' : 'SSP ', row.join(' ')); }
