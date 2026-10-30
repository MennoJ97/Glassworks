// Runs the WebAssembly solver off the main thread, so the page stays responsive.
import init, { solve } from './pkg/eigenmenno_wasm.js';

const ready = init();

onmessage = async ({ data: { id, requests } }) => {
  await ready;
  const start = performance.now();
  try {
    const transfer = [];
    const results = requests.map((request) => {
      const solution = solve(JSON.stringify(request));
      const channels = {};
      solution.channel_names().forEach((name, i) => {
        channels[name] = solution.channel(i);
        transfer.push(channels[name].buffer);
      });
      const meta = JSON.parse(solution.meta());
      solution.free();
      return { meta, channels };
    });
    postMessage({ id, results, ms: performance.now() - start }, transfer);
  } catch (e) {
    // A trap means the module's memory is unusable, so the worker must restart.
    postMessage({ id, error: e.message ?? String(e), fatal: e instanceof WebAssembly.RuntimeError });
  }
};
