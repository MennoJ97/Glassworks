// The optimizer in web/js/design.js, against the WebAssembly solver the page uses. The physics
// has its own tests in Rust, and these check that the optimizer gets the most out of it.
// Run `scripts/build-web.sh` first, then `npm test`.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { defaultState } from '../web/js/model.js';
import { Mma, coneFilter, coneFilter2d, findResonance, goalGradient, gradientRequests, meanReflection, optimize, pixelCount, prepareDesign, project, references, totalIterations } from '../web/js/design.js';
import { initSync, solve as solve1d, solve_2d } from '../web/pkg/eigenmenno_wasm.js';

initSync({ module: readFileSync(new URL('../web/pkg/eigenmenno_wasm_bg.wasm', import.meta.url)) });

const solveOne = (request) => {
  const solution = (request.width === undefined ? solve1d : solve_2d)(JSON.stringify(request));
  const channels = {};
  solution.channel_names().forEach((name, i) => { channels[name] = solution.channel(i); });
  const meta = JSON.parse(solution.meta());
  solution.free();
  return { meta, channels };
};
const solve = async (requests) => ({ results: requests.map(solveOne), ms: 0 });
const dot = (a, b) => a.reduce((sum, v, i) => sum + v * b[i], 0);
const noise = (n, seed) => Float64Array.from({ length: n }, (_, i) => Math.sin(12.9898 * (i + 1) * seed) * 0.5 + 0.5);
const state = (analysis) => ({ ...defaultState(), analysis });

test('each filter\'s adjoint is its true transpose, and it keeps a constant', () => {
  for (const [filter, n] of [[coneFilter(40, 3.5), 40], [coneFilter2d(7, 6, 2.5), 42]]) {
    const [x, y] = [noise(n, 1), noise(n, 2)];
    assert.ok(Math.abs(dot(filter.apply(x), y) - dot(x, filter.adjoint(y))) < 1e-12);
    filter.apply(new Float64Array(n).fill(0.3)).forEach((v) => assert.ok(Math.abs(v - 0.3) < 1e-12));
  }
});

test('projection keeps its ends and centre, and reports its true slope', () => {
  for (const beta of [1, 8, 64]) {
    const { rho } = project(Float64Array.of(0, 0.5, 1), beta);
    [0, 0.5, 1].forEach((want, i) => assert.ok(Math.abs(rho[i] - want) < 1e-12));
    const at = Float64Array.of(0.2, 0.45, 0.7);
    const { slope } = project(at, beta);
    const h = 1e-6;
    const [up, down] = [project(at.map((v) => v + h), beta).rho, project(at.map((v) => v - h), beta).rho];
    slope.forEach((v, i) => assert.ok(Math.abs(v - (up[i] - down[i]) / (2 * h)) < 1e-5 * Math.max(1, v)));
  }
});

test('MMA finds a bowl\'s minimum, stays in bounds, and respects its move limit', () => {
  const target = Float64Array.of(0.1, 0.5, 0.9, 1.4, -0.3);
  const best = target.map((v) => Math.min(1, Math.max(0, v)));
  const move = 0.1;
  const mma = new Mma(move);
  let x = new Float64Array(5).fill(0.5);
  for (let i = 0; i < 60; i++) {
    const next = mma.step(x, x.map((v, k) => 2 * (v - target[k])));
    next.forEach((v, k) => {
      assert.ok(v >= 0 && v <= 1);
      assert.ok(Math.abs(v - x[k]) <= move + 1e-12);
    });
    x = next;
  }
  x.forEach((v, k) => assert.ok(Math.abs(v - best[k]) < 1e-2, `${v} against ${best[k]}`));
});

// The whole chain the optimizer follows, from `x` through filter, projection and solver to its
// goal, against finite differences. A wrong sign anywhere would make it walk the wrong way.
for (const [analysis, goal, radius] of [['mirror', (t) => -t, 4], ['antireflection', (t) => t, 4], ['bend2d', (t) => -t, 3]]) {
  test(`the ${analysis} gradient matches finite differences through the whole chain`, async () => {
    const s = state(analysis);
    await prepareDesign(s, solve);
    const n = pixelCount(s.design[analysis]);
    const filter = analysis === 'bend2d' ? coneFilter2d(24, 24, radius) : coneFilter(n, radius);
    const beta = 4;
    const x = noise(n, 3).map((v) => 0.3 + 0.4 * v);
    const value = async (at) => goal(meanReflection((await solve(gradientRequests(s, project(filter.apply(at), beta).rho))).results));
    const { rho, slope } = project(filter.apply(x), beta);
    const exact = filter.adjoint(goalGradient(s, (await solve(gradientRequests(s, rho))).results, slope));
    const h = 1e-4;
    const picks = [Math.floor(n * 0.3), Math.floor(n * 0.52), Math.floor(n * 0.77)];
    for (const k of picks) {
      const shift = (d) => { const y = x.slice(); y[k] += d; return y; };
      const numeric = ((await value(shift(h))) - (await value(shift(-h)))) / (2 * h);
      assert.ok(Math.abs(numeric - exact[k]) < 2e-3 * Math.abs(exact[k]) + 1e-7, `pixel ${k}: ${numeric} against ${exact[k]}`);
    }
  });
}

// A run from start to finish, returning the last frame and each iteration's value.
async function run(s) {
  const seen = [];
  let last = null;
  const failure = await optimize(s, solve, ({ rho, reply, iteration, beta }) => {
    seen.push(meanReflection(reply.results));
    last = { rho, iteration, beta };
    return true;
  });
  assert.equal(failure?.error, undefined);
  return { ...last, seen, grey: last.rho.reduce((sum, r) => sum + (r > 0.02 && r < 0.98), 0) / last.rho.length };
}
const referenceValues = async (s) => Promise.all(references(s).map(async ({ requests }) => meanReflection((await solve(requests)).results)));

test('a mirror ends sharp, and reflects about as well as the quarter-wave stack that fits', async () => {
  const s = state('mirror');
  const { seen, grey, iteration, beta } = await run(s);
  assert.equal(beta, Infinity);
  assert.equal(iteration, totalIterations(s));
  assert.ok(grey < 0.05, `grey ${grey}`);
  const [stack] = await referenceValues(s);
  assert.ok(seen.at(-1) > 0.98 * stack, `${seen.at(-1)} against ${stack}`);
});

test('a coating ends sharp, and beats both no coating and one quarter-wave layer', async () => {
  const s = state('antireflection');
  const { seen, grey } = await run(s);
  assert.ok(grey < 0.05, `grey ${grey}`);
  const [none, single] = await referenceValues(s);
  assert.ok(seen.at(-1) < single && single < none, `${seen.at(-1)}, ${single}, ${none}`);
});

test('a bend ends sharp, and carries nearly all the light that a plain corner loses', async () => {
  const s = state('bend2d');
  const { seen, grey } = await run(s);
  assert.ok(grey < 0.05, `grey ${grey}`);
  const [plain] = await referenceValues(s);
  assert.ok(plain < 0.1, `plain ${plain}`);
  // Three starts in plans/2d-solver.md reach 98.1 to 98.4%.
  assert.ok(seen.at(-1) > 0.97, `T ${seen.at(-1)}`);
});

test('the same start makes the same design, which a band at twice the size also gets', async () => {
  const at = async (from, to) => {
    const s = state('bend2d');
    s.design = { ...s.design, bend2d: { ...s.design.bend2d, from, to } };
    await prepareDesign(s, solve);
    const rho = noise(576, 4);
    return { rho, s, results: (await solve(gradientRequests(s, rho))).results };
  };
  const [a, again, big] = [await at(0.95, 1.05), await at(0.95, 1.05), await at(1.9, 2.1)];
  assert.equal(meanReflection(a.results), meanReflection(again.results));
  const slope = new Float64Array(576).fill(1);
  const relative = (x) => { const g = goalGradient(x.s, x.results, slope); const peak = Math.max(...g.map(Math.abs)); return g.map((v) => v / peak); };
  const [ga, gb] = [relative(a), relative(big)];
  assert.ok(Math.abs(meanReflection(a.results) - meanReflection(big.results)) < 1e-10);
  ga.forEach((v, i) => assert.ok(Math.abs(v - gb[i]) < 1e-9));
});

const cavityState = () => ({ ...state('cavityDesign'), wavelength: 1 });

test('a cavity\'s frequency slope predicts how far a density moves its resonance', async () => {
  // A true resonance: the crystal's rods as pixels, close to their own pole. A grey structure
  // at a low Q has no single pole for this to describe.
  const s = cavityState();
  const [n, h, a] = [64, 1 / 24, 0.3781];
  const inset = (64 / 24 - 7 * a) / 2;
  const rho = Float64Array.from({ length: n * n }, (_, k) => {
    const [x, y] = [(k % n + 0.5) * h - inset, (Math.floor(k / n) + 0.5) * h - inset];
    const [i, j] = [Math.floor(x / a), Math.floor(y / a)];
    return i >= 0 && j >= 0 && i < 7 && j < 7 && !(i === 3 && j === 3) && Math.hypot(x - (i + 0.5) * a, y - (j + 0.5) * a) < 0.2 * a ? 1 : 0;
  });
  const Q = 2000;
  let wavelength = 1;
  for (let i = 0; i < 4; i++) wavelength = (2 * Math.PI) / (await findResonance(s, rho, Q, solve, wavelength)).poleK;
  const { poleK, moves } = await findResonance(s, rho, Q, solve, wavelength);
  const strongest = moves.reduce((best, v, k) => (Math.abs(v) > Math.abs(moves[best]) ? k : best), 0);
  const step = 1e-3;
  const nudged = rho.slice();
  nudged[strongest] += rho[strongest] > 0.5 ? -step : step;
  const after = (await findResonance(s, nudged, Q, solve, wavelength)).poleK;
  const [measured, predicted] = [after - poleK, moves[strongest] * (nudged[strongest] - rho[strongest])];
  assert.ok(Math.abs(measured - predicted) < 0.02 * Math.abs(predicted), `${measured} against ${predicted}`);
});

test('a cavity beats the crystal in its footprint, on target, and as drawn rectangles', async () => {
  const s = cavityState();
  let last = null;
  let early = 0;
  const failure = await optimize(s, solve, (step) => {
    last = step;
    if (step.iteration === 80) early = meanReflection(step.reply.results);
    return true;
  });
  assert.equal(failure?.error, undefined);
  assert.equal(last.beta, Infinity);
  // Minimizing can finish nearly as well as maximizing (see plans/2d-solver.md). So the goal's sign
  // shows best here, at Q 100, where 18 starts reach 15 to 98 maximizing and 0.15 to 0.4 minimizing.
  assert.ok(early > 20, `early ${early}`);
  last.rho.forEach((v) => assert.ok(v === 0 || v === 1));
  const ldos = meanReflection(last.reply.results);
  const [crystal] = await referenceValues(s);
  // In plans/2d-solver.md, the crystal reads about 300, and 18 starts 886 to 1,075.
  assert.ok(crystal > 250 && crystal < 350, `crystal ${crystal}`);
  assert.ok(ldos > 1.5 * crystal, `${ldos} against ${crystal}`);
  assert.ok(Math.abs(last.mode.wavelength - 1) < 2e-3, `mode at ${last.mode.wavelength}`);
  assert.ok(last.mode.q > 1000, `Q ${last.mode.q}`);
});
