// The pixel drawing against merged rectangles as the mesh refines, and how far a nanometre moves the mode.
import { load, pixels } from './pix.mjs';
import { F, pole3, merged } from './hz.mjs';
const seed = process.argv[2]; const subs = process.argv[3].split(',').map(Number); const D = load(seed); const m = merged(D), p0 = pixels(D, () => 0);
const guess = { 1: 1.0000, 2: 1.0045, 3: 1.0055, 4: 1.006 };
for (const sub of subs) { const a = pole3(D, m, sub, guess[sub], 45); const b = pole3(D, p0, sub, a.w, 45); const up = pole3(D, pixels(D, () => 0.001), sub, b.w * 1.005, 45), dn = pole3(D, pixels(D, () => -0.001), sub, b.w * 0.995, 45);
  console.log(JSON.stringify({ seed: +seed, sub, merged: { w: a.w, Q: a.Q, F45: F(D, m, sub, a.w, 45), F0: F(D, m, sub, a.w, 0) }, pixels: { w: b.w, Q: b.Q, F45: F(D, p0, sub, b.w, 45), F0: F(D, p0, sub, b.w, 0) }, up: up.w, dn: dn.w, shift_pct_per_nm: (up.w - dn.w) / 2 / b.w * 100, Qup: up.Q, Qdn: dn.Q })); }
