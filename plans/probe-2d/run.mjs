import { readFileSync } from 'node:fs';
const bytes = readFileSync('sparse.wasm');
const fresh = async () => (await WebAssembly.instantiate(bytes, {})).instance.exports;
const attempt = async (name, nx) => {
  const m = await fresh();
  const t = performance.now();
  try { const x = m[name](nx); return `${(performance.now() - t).toFixed(0)} ms (${x.toFixed(3)})`; }
  catch (e) { return `failed: ${String(e.message ?? e).slice(0, 60)}`; }
};
for (const nx of [25, 50, 100, 200, 400, 800]) {
  const m = await fresh();
  const t = performance.now(); const n = m.wasm_assemble(nx); const asm = performance.now() - t;
  console.log(`nx ${nx}: dofs ${n}, assemble ${asm.toFixed(0)} ms, faer incl. assembly ${await attempt('wasm_faer', nx)}, banded incl. assembly ${await attempt('wasm_banded', nx)}`);
}
