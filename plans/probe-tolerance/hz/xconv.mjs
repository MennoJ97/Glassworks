// The crystal's mesh convergence, for comparison with the designs'.
import { F, pole3, peak, crystal, crystalD } from './hz.mjs';
const D = crystalD, r0 = 0.4 * 0.2759; const g = { 1: 0.99, 2: 0.9946, 3: 0.996 };
for (const sub of [1, 2, 3]) { const pk = peak(D, crystal(), sub, g[sub], 0, 0.02).w; const b = pole3(D, crystal(), sub, pk, 0, [100, 300, 1000, 3000, 1e4, 3e4]); const up = pole3(D, crystal(() => r0 - 0.001), sub, b.w * 1.006, 0, [300, 1000, 3000, 1e4]), dn = pole3(D, crystal(() => r0 + 0.001), sub, b.w * 0.994, 0, [300, 1000, 3000, 1e4]); console.log(JSON.stringify({ sub, w: b.w, Q: b.Q, F: F(D, crystal(), sub, b.w, 0), shift: (up.w - dn.w) / 2 / b.w * 100 })); }
