// Times the WebAssembly build, which holds no clock of its own. Each measurement gets a fresh
// instance, so one solver's memory growth and JIT warm-up don't reach the other's.
import { readFileSync } from 'node:fs';
const bytes = readFileSync(process.argv[2] ?? 'sparse.wasm');
const fresh = async () => (await WebAssembly.instantiate(bytes, {})).instance.exports;
const time = (f) => { const t = performance.now(); const v = f(); return [performance.now() - t, v]; };
const count = (await fresh()).cases();
for (let i = 0; i < count; i++) {
  try {
    const a = await fresh(); const n = a.build(i);
    const [analyse] = time(() => a.analyse());
    const [m1, r1] = time(() => a.run_mf());
    const [m2] = time(() => a.run_mf());
    const b = await fresh(); b.build(i);
    const [f1, r2] = time(() => b.run_faer());
    const [f2] = time(() => b.run_faer());
    console.log(`${i} n ${n}: ours analyse ${analyse.toFixed(0)}, first ${m1.toFixed(0)}, again ${m2.toFixed(0)} ms, residual ${r1.toExponential(1)}; faer first ${f1.toFixed(0)}, again ${f2.toFixed(0)} ms, residual ${r2.toExponential(1)}`);
  } catch (e) { console.log(`${i}: failed: ${String(e.message ?? e).slice(0, 80)}`); }
}
