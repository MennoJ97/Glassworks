// Density topology optimization: a row of thin pixels for layers, or a square of them for a
// bend, each with a density that mixes two materials. The solver supplies the adjoint gradient,
// and `optimizer.js` does the rest. Cavities have their own file, `cavity-design.js`.

import { CAVITY, cavityGeometry, cavityReferences, cavityRequest, optimizeCavity, tuned } from './cavity-design.js';
import { MATERIALS, scaled, toShape } from './model.js';
import { Mma, coneFilter, coneFilter2d, project, randomStart } from './optimizer.js';

const PIXEL = 0.005;
const AMBIENT = 0.3;
const SHARPNESS = [1, 2, 4, 8, 16, 32, 64, 128];
const ITERATIONS_PER_STAGE = 40;
// A 2D solve costs far more than a 1D one, and these converge sooner: for E_z, 25 did no better
// over three starts, while 5 sometimes fell short. H_z converges more slowly: its mean over
// three starts rose from 97.2% at 10 to 98.0% at 20 and 98.4% at 40.
const bendIterations = (s) => (s.polarization === 'tm' ? 40 : 10);
export const stageLength = (s) => (s.analysis === 'cavityDesign' ? CAVITY.iterations : s.analysis === 'bend2d' ? bendIterations(s) : ITERATIONS_PER_STAGE);
const sharpnessOf = (s) => (s.analysis === 'cavityDesign' ? CAVITY.stages.map(([, beta]) => beta) : SHARPNESS);
export const totalIterations = (s) => sharpnessOf(s).length * stageLength(s);
export const STAGES = SHARPNESS.map((beta, i) => ({ frame: i * ITERATIONS_PER_STAGE, name: `β ${beta}` }));
export const sharpnessAt = (s, frame) => sharpnessOf(s)[Math.floor(frame / stageLength(s))] ?? Infinity;

const material = (name) => MATERIALS.find((m) => m.name === name);
const eps = (name) => material(name).n ** 2;
const solverLayer = (thickness, e) => ({ thickness, eps: [e, 0], mu: [1, 0] });

export const pixelCount = (d) => (d.cavity ? CAVITY.pixels ** 2 : d.thickness === undefined ? BEND.pixels ** 2 : Math.max(1, Math.round(d.thickness / PIXEL)));

// A slab of index 3.5 turning a corner through a square of pixels, in units of the band's
// centre wavelength, so a band anywhere gets the same design at its own size. The guides meet
// the square `turn` from its inner corner, which keeps a smooth bend tight, and the rest of the
// square reaches past the corner for the design to use.
// Everything lies on a grid of 1/24, so each pixel is a whole number of pairs of triangles.
const BEND = { pixels: 18, x0: 1, y0: 0.75, size: 0.75, turn: 0.25, domain: 2.5, n: 3.5, pml: 0.5, span: 0.35 };
const bendCentre = { x: BEND.x0 + BEND.turn, y: BEND.y0 + BEND.size - BEND.turn };
const rectangle = (x0, y0, x1, y1) => ({ type: 'rect', x: (x0 + x1) / 2, y: (y0 + y1) / 2, w: x1 - x0, h: y1 - y0, n: BEND.n });

export const bendScale = (s) => (s.design.bend2d.from + s.design.bend2d.to) / 2;

// Cells per pixel side, and cells per wavelength in the slab at the band's short end. Three
// per pixel is as many as the solver's limit on unknowns allows.
function bendRefinement(s) {
  const shortest = s.design.bend2d.from / bendScale(s);
  const cell = BEND.size / BEND.pixels;
  const refine = Math.min(3, Math.ceil((4.5 * BEND.n * cell) / shortest));
  return { refine, cells: (shortest * refine) / (BEND.n * cell) };
}

// A thin slab barely holds H_z, whose electric field crosses its walls: at 0.125 its effective
// index is 1.32, against 2.73 for E_z. Twice as wide gives H_z a similar 2.92, so that the two
// polarizations ask comparable things of a design.
const bendWidth = (polarization) => (polarization === 'tm' ? 0.25 : 0.125);

// Why a band's percentages can't be trusted, if they can't. Below 0.84 of the centre, the slab
// carries a mode more than it does at the centre, and the port launches only the first.
export function bendWarning(s) {
  if (bendRefinement(s).cells < 4) return 'This band is too wide for the mesh to resolve its short end, so these percentages aren\'t reliable.';
  if (s.design.bend2d.from / bendScale(s) < 0.84) return 'The guide carries another mode at the short end of this band, so these percentages are less reliable there. A narrower band avoids that.';
  return null;
}

// A quarter circle joining the guides through the square, as wide as fits: a ring about the
// square's inner corner, with air laid over the three quarters not wanted, under the straight guides.
function curvedBend(width) {
  const [cx, cy, r] = [BEND.x0, BEND.y0 + BEND.size, BEND.turn];
  const [outer, inner] = [r + width / 2, r - width / 2];
  const disc = (radius, n) => ({ type: 'circle', x: cx, y: cy, w: 2 * radius, h: 2 * radius, n });
  const air = (x0, y0, x1, y1) => ({ ...rectangle(x0, y0, x1, y1), n: 1 });
  return [disc(outer, BEND.n), disc(inner, 1), air(cx - outer, cy - outer, cx, cy + outer), air(cx, cy, cx + outer, cy + outer)];
}

// What the plot needs, and each solve shares. `corner`, if given, fills the square with a
// 'sharp' corner or a 'curved' bend for comparison.
export function bendGeometry(straight, scale, polarization, corner) {
  const { x0, y0, size, domain, span } = BEND;
  const width = bendWidth(polarization);
  const { x, y } = bendCentre;
  const far = domain + 3;
  const shapes = straight
    ? [rectangle(-3, y - width / 2, far, y + width / 2)]
    : corner === 'sharp'
    ? [rectangle(-3, y - width / 2, x + width / 2, y + width / 2), rectangle(x - width / 2, y - width / 2, x + width / 2, far)]
    : [
      ...(corner === 'curved' ? curvedBend(width) : []),
      rectangle(-3, y - width / 2, x0, y + width / 2), rectangle(x - width / 2, y0 + size, x + width / 2, far),
    ];
  const end = domain - 0.4;
  return scaled({
    size: domain, shapes,
    source: { type: 'port', x: 0.4, y0: y - 0.7, y1: y + 0.7 },
    mark: [0.4, y],
    monitor: straight ? { x0: end, y0: y - span, x1: end, y1: y + span } : { x0: x + span, y0: end, x1: x - span, y1: end },
    region: { x0, y0, width: size, height: size, nx: BEND.pixels, ny: BEND.pixels },
  }, scale);
}

function bendRequest(s, wavelength, straight, rho, corner) {
  const scale = bendScale(s);
  const g = bendGeometry(straight, scale, s.polarization, corner);
  const n2 = BEND.n ** 2;
  return {
    width: g.size,
    height: g.size,
    cell: (BEND.size * scale) / BEND.pixels / bendRefinement(s).refine,
    background: { eps: 1, mu: 1 },
    shapes: g.shapes.map(toShape),
    polarization: s.polarization,
    walls: 'absorbing',
    monitors: [g.monitor],
    analysis: { type: 'driven', wavelength, source: g.source, pml: BEND.pml * scale },
    ...(rho ? { density: { ...g.region, eps_low: 1, eps_high: n2, rho: Array.from(rho) } } : {}),
  };
}

// A straight guide's power at each wavelength, by its `k0`, which turns power into transmission.
const straightPower = new Map();

// Solves the straight guide once for these wavelengths. Returns a reply holding any error.
export async function prepareDesign(s, solve) {
  if (s.analysis !== 'bend2d') return {};
  const wavelengths = goalWavelengths(s);
  const reply = await solve(wavelengths.map((wavelength) => bendRequest(s, wavelength, true)));
  if (!reply || reply.error) return reply;
  reply.results.forEach((r) => straightPower.set(r.meta.fields[0].k0, r.channels.monitors[0]));
  return reply;
}

export const designGeometry = (s, rho) => (s.analysis === 'bend2d'
  ? bendGeometry(false, bendScale(s), s.polarization)
  : cavityGeometry(s.wavelength * (tuned.get(rho) ?? 1)));

// Goal wavelengths, with the centre last so the plotted field comes from it.
function goalWavelengths(s) {
  if (s.analysis === 'mirror') return [s.wavelength];
  const { from, to } = s.design[s.analysis];
  const count = s.analysis === 'bend2d' ? 3 : 11;
  const wavelengths = Array.from({ length: count }, (_, i) => from + ((to - from) * i) / (count - 1));
  wavelengths.push(wavelengths.splice((count - 1) / 2, 1)[0]);
  return wavelengths;
}

function requests(s, layers, makeAnalysis) {
  const wavelengths = goalWavelengths(s);
  const max_cell = Math.min(...wavelengths) / 6;
  return wavelengths.map((wavelength) => ({
    layers,
    polarization: s.polarization,
    left: 'port',
    right: 'absorbing',
    mesh: { order: 2, max_cell },
    analysis: makeAnalysis(wavelength),
  }));
}

function pixelLayers(d, rho) {
  const [low, high] = [eps(d.low), eps(d.high)];
  return [
    solverLayer(AMBIENT, eps(d.left)),
    ...Array.from(rho, (r) => solverLayer(PIXEL, low + r * (high - low))),
    solverLayer(AMBIENT, eps(d.right)),
  ];
}

export const gradientRequests = (s, rho) => (s.analysis === 'cavityDesign' ? [cavityRequest(s, rho, CAVITY.Q)] : s.analysis === 'bend2d'
  ? goalWavelengths(s).map((wavelength) => bendRequest(s, wavelength, false, rho))
  : requests(s, pixelLayers(s.design[s.analysis], rho),
    (wavelength) => ({ type: 'driven_gradient', wavelength, angle_deg: s.angle })));

// What a design is judged by: reflection for layers, transmission for a bend, and for a
// cavity, LDOS relative to vacuum, where a line source's `F` is k / 2.
const measure = ({ meta, channels }) => (meta.fields[0].ldos != null
  ? (2 * meta.fields[0].ldos) / meta.fields[0].k0
  : channels.monitors
  ? channels.monitors[0] / straightPower.get(meta.fields[0].k0)
  : meta.fields[0].power[0]);
export const meanReflection = (results) => results.reduce((sum, r) => sum + measure(r), 0) / results.length;

// The structures a designer would try by hand, solved at the same wavelengths.
export function references(s) {
  if (s.analysis === 'cavityDesign') return cavityReferences(s);
  if (s.analysis === 'bend2d') {
    return [
      { name: 'Empty region', requests: gradientRequests(s, new Float64Array(BEND.pixels ** 2)) },
      { name: 'Sharp corner', requests: goalWavelengths(s).map((wavelength) => bendRequest(s, wavelength, false, undefined, 'sharp')) },
      { name: 'Curved bend', requests: goalWavelengths(s).map((wavelength) => bendRequest(s, wavelength, false, undefined, 'curved')) },
    ];
  }
  const d = s.design[s.analysis];
  const quarter = (name, wavelength) => wavelength / (4 * material(name).n);
  const wrap = (inner) => [
    solverLayer(AMBIENT, eps(d.left)),
    ...inner.map(([name, t]) => solverLayer(t, eps(name))),
    solverLayer(AMBIENT, eps(d.right)),
  ];
  const plain = (name, inner) => ({
    name,
    requests: requests(s, wrap(inner), (wavelength) => ({ type: 'driven', wavelength, angle_deg: s.angle, source: { type: 'plane_wave' } })),
  });
  if (s.analysis === 'mirror') {
    const [hi, lo] = [quarter(d.high, s.wavelength), quarter(d.low, s.wavelength)];
    const pairs = Math.max(0, Math.floor((d.thickness - hi) / (hi + lo)));
    const stack = [[d.high, hi], ...Array.from({ length: pairs }, () => [[d.low, lo], [d.high, hi]]).flat()];
    return [plain(`Quarter-wave stack, ${stack.length} layers`, stack)];
  }
  return [
    plain('No coating', []),
    plain('One quarter-wave layer', [[d.low, quarter(d.low, (d.from + d.to) / 2)]]),
  ];
}

// The slope of what the optimizer minimizes, for each pixel's density before projection,
// where `slope` is the projection's own. A mirror and a bend want more of what they measure.
export function goalGradient(s, results, slope) {
  const d = s.design[s.analysis];
  const gradient = new Float64Array(slope.length);
  for (const { channels, meta } of results) {
    const [slopes, scale] = s.analysis === 'bend2d'
      ? [channels.gradient, -1 / straightPower.get(meta.fields[0].k0)]
      : [channels['dR/deps_re'], (s.analysis === 'mirror' ? -1 : 1) * (eps(d.high) - eps(d.low))];
    for (let k = 0; k < gradient.length; k++) gradient[k] += (scale * slopes[k] * slope[k]) / results.length;
  }
  return gradient;
}

// Runs until the schedule ends or `onStep` returns false. `solve` returns a solver
// reply, or null if a newer request replaced it. Returns a reply holding any error.
export async function optimize(s, solve, onStep) {
  if (s.analysis === 'cavityDesign') return optimizeCavity(s, solve, onStep);
  const d = s.design[s.analysis];
  const bend = s.analysis === 'bend2d';
  const n = pixelCount(d);
  const filter = bend ? coneFilter2d(BEND.pixels, BEND.pixels, d.radius) : coneFilter(n, Math.max(d.radius / (PIXEL * 1000), 1));
  const prepared = await prepareDesign(s, solve);
  if (!prepared || prepared.error) return prepared;
  let x = randomStart(n, d.seed);
  let iteration = 0;
  for (const beta of SHARPNESS) {
    // Projection amplifies a change in x by up to about beta / 2, so steps shrink to match.
    const mma = new Mma(Math.min(0.1, 0.2 / beta));
    for (let i = 0; i < stageLength(s); i++) {
      const { rho, slope } = project(filter.apply(x), beta);
      const reply = await solve(gradientRequests(s, rho));
      if (!reply || reply.error) return reply;
      if (!onStep({ rho, reply, iteration: ++iteration, beta })) return null;
      x = mma.step(x, filter.adjoint(goalGradient(s, reply.results, slope)));
    }
  }
  const { rho } = project(filter.apply(x), 1000);
  const reply = await solve(gradientRequests(s, rho));
  if (!reply || reply.error) return reply;
  onStep({ rho, reply, iteration, beta: Infinity });
  return null;
}

// Merges neighbouring pixels of nearly equal density into layers for plotting.
export function plotLayers(d, rho) {
  const [low, high] = [eps(d.low), eps(d.high)];
  const out = [{ ...material(d.left), thickness: AMBIENT }];
  let start = 0;
  for (let i = 1; i <= rho.length; i++) {
    if (i === rho.length || Math.abs(rho[i] - rho[start]) > 0.02) {
      out.push({ n: Math.sqrt(low + rho[start] * (high - low)), k: 0, thickness: (i - start) * PIXEL });
      start = i;
    }
  }
  out.push({ ...material(d.right), thickness: AMBIENT });
  return out;
}

// A finished design as ordinary editable layers. Densities must already be sharp.
export function finishedLayers(d, rho) {
  const out = [{ ...material(d.left), thickness: 0.5 }];
  let start = 0;
  for (let i = 1; i <= rho.length; i++) {
    if (i === rho.length || rho[i] > 0.5 !== rho[start] > 0.5) {
      out.push({ ...material(rho[start] > 0.5 ? d.high : d.low), thickness: +((i - start) * PIXEL).toFixed(6) });
      start = i;
    }
  }
  out.push({ ...material(d.right), thickness: 0.5 });
  return out;
}
