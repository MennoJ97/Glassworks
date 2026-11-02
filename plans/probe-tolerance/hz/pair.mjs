// The mode a dipole at the centre excites belongs to a pair. Along the mesh's own mirror lines, 45 and 135 degrees, each member is excited alone.
import { load } from './pix.mjs';
import { F, peak, pole3, merged } from './hz.mjs';
const seed = process.argv[2]; const subs = process.argv[3].split(',').map(Number); const D = load(seed); const s = merged(D);
for (const sub of subs) { const t = Date.now(); const pk = peak(D, s, sub, 1, 0, 0.03); const a = pole3(D, s, sub, pk.w, 45), b = pole3(D, s, sub, pk.w, 135);
  console.log(JSON.stringify({ seed: +seed, sub, peak0: pk, w45: a.w, Q45: a.Q, w135: b.w, Q135: b.Q, split_pct: (a.w - b.w) / pk.w * 100, F0: F(D, s, sub, pk.w, 0), F90: F(D, s, sub, pk.w, 90), F45_at45: F(D, s, sub, a.w, 45), F135_at135: F(D, s, sub, b.w, 135), sec: (Date.now() - t) / 1000 })); }
