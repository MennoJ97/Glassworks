// SSP's transpose against finite differences, on a random smooth field.
const { ssp, randomStart, coneFilter2d } = await import(`${process.cwd()}/rob2/web/js/optimizer.js`);
const n = 64; const f = coneFilter2d(n, n, 2); const x = f.apply(randomStart(n * n, 3).map((v) => 0.5 + 3 * (v - 0.5)));
for (const beta of [8, 64, 128, Infinity]) {
  const { rho, back } = ssp(x, n, beta); const u = rho.map((_, k) => Math.sin(k * 0.37)); const g = back(u);
  let worst = 0, checked = 0; const J = (a) => ssp(a, n, beta).rho.reduce((s, v, k) => s + v * u[k], 0);
  for (let t = 0; t < 400; t++) { const k = 100 + t * 9; const h = 1e-5; const p = x.slice(), m = x.slice(); p[k] += h; m[k] -= h; const fd = (J(p) - J(m)) / (2 * h); if (Math.abs(fd) > 1e-6 || Math.abs(g[k]) > 1e-6) { checked++; worst = Math.max(worst, Math.abs(fd - g[k]) / Math.max(Math.abs(fd), 1e-3)); } }
  const grey = rho.filter((v) => v > 0.01 && v < 0.99).length / rho.length; const norm = Math.sqrt(g.reduce((s, v) => s + v * v, 0));
  console.log(JSON.stringify({ beta: String(beta), checked, worstRelativeError: worst, greyFraction: grey, gradientNorm: norm }));
}
