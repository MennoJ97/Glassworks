// The probe's native solver as a drop-in for the WebAssembly `one`: same request, same
// `{ meta, channels }`, with no 4 GB limit.
import { spawnSync } from 'node:child_process';
const BIN = `${process.cwd()}/../aniso/target-native/release/examples/serve`;
export const one = (req) => {
  const r = spawnSync(BIN, [], { input: JSON.stringify(req), maxBuffer: 2 ** 31 });
  if (r.status !== 0) throw new Error(`native solver exited with ${r.status}: ${r.stderr}`);
  const out = r.stdout, end = out.indexOf(10);
  const head = JSON.parse(out.subarray(0, end).toString());
  if (head.error) throw new Error(head.error);
  const channels = {};
  let at = end + 1;
  for (const [name, n] of head.channels) {
    channels[name] = new Float64Array(out.buffer.slice(out.byteOffset + at, out.byteOffset + at + 8 * n));
    at += 8 * n;
  }
  return { meta: head.meta, channels };
};
