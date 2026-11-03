// The pixels around the source, and the smallest gap and line along the rows through it.
import { readFileSync } from 'node:fs';
const [id] = process.argv.slice(2); const d = JSON.parse(readFileSync(`des/${id}.json`)); const n = 64; const m = (i, j) => d.rho[i + j * n] > 0.5;
let s = ''; for (let j = 22; j < 42; j++) { for (let i = 22; i < 42; i++) s += m(i, j) ? '#' : '.'; s += '\n'; }
const px = (64 / 24) * d.scale / n * 1000; let gap = 0; while (gap < 32 && !m(32 + gap, 32)) gap++;
console.log(`${id}: pixel ${px.toFixed(1)} nm, air along x from source ${gap} px (${(gap * px).toFixed(0)} nm), F ${d.F.toFixed(0)}, Q ${d.Q.toFixed(0)}, V ${d.V.toFixed(3)}\n${s}`);
