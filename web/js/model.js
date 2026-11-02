// The page's state: the geometry the visitor builds, and how it becomes solver
// requests. A 2D version would add a new geometry here, emitting regions instead
// of layers, while the materials and analyses stay the same.

export const MATERIALS = [
  { name: 'Air', n: 1, k: 0 },
  { name: 'Magnesium fluoride', n: 1.38, k: 0 },
  { name: 'Silica glass', n: 1.45, k: 0 },
  { name: 'Silicon nitride', n: 2, k: 0 },
  { name: 'Titanium dioxide', n: 2.4, k: 0 },
  { name: 'Silicon', n: 3.48, k: 0 },
  { name: 'Absorbing film', n: 2, k: 0.2 },
  { name: 'Gold-like metal', n: 0.25, k: 6 },
];

export const layer = (name, thickness) => ({ ...MATERIALS.find((m) => m.name === name), thickness });

const H = () => layer('Titanium dioxide', 0.104);
const L = () => layer('Silica glass', 0.172);
const air = (t = 0.5) => layer('Air', t);

// A square lattice of rods with the middle one missing, in open space. Periodic walls
// turn the same shapes into a supercell. The target sits inside the lattice's band gap.
function defectCavity(n) {
  const rods = Array.from({ length: n * n }, (_, i) => ({ type: 'circle', x: i % n + 0.5, y: Math.floor(i / n) + 0.5, w: 0.4, h: 0.4, n: 3.4641 }));
  rods.splice((n * n - 1) / 2, 1);
  return { width: n, height: n, cell: 0.1, background: 1, walls: 'absorbing', target: 2.65, sourceX: n / 2, sourceY: n / 2, shapes: rods };
}

// The waveguide bends a view can show, with what each needs. Bends in rods are from Mekis et
// al. 1996, whose guided band runs from wavelengths of 3.2a down to 2.26a. Their layer is
// thick, because rods running through it reflect more than free space does.
export const GUIDES = {
  sharp: { name: 'Sharp bend in a photonic crystal', wavelength: 2.833, cell: 0.16667, n: 3.4, band: [2.4, 3.1] },
  stepped: { name: 'Stepped bend in a photonic crystal', wavelength: 2.833, cell: 0.16667, n: 3.4, band: [2.4, 3.1] },
  slab: { name: 'Sharp bend in an ordinary waveguide', wavelength: 1, cell: 0.05, n: 3.5, band: [1, 1.5] },
};
const guideSetup = (bend) => ({ bend, scale: 1, cell: GUIDES[bend].cell, background: 1, walls: 'absorbing', shapes: [] });
const LATTICE = { arm: 5, clad: 3, pml: 5 };
const SLAB = { size: 4, x: 2.5, y: 1.5, width: 0.143, n: 3.5, pml: 1 };

// Multiplies every length in a geometry. Maxwell's equations have no length of their own,
// so a structure scaled with its wavelength behaves exactly as it did.
export function scaled(value, factor, key) {
  if (typeof value === 'number') return ['n', 'nx', 'ny'].includes(key) ? value : value * factor;
  if (Array.isArray(value)) return value.map((v) => scaled(v, factor));
  if (value && typeof value === 'object') return Object.fromEntries(Object.entries(value).map(([k, v]) => [k, scaled(v, factor, k)]));
  return value;
}

// Links from before scaling existed have no scale.
export const guideScale = (s) => s.plane.guide2d.scale ?? 1;

// Cells per wavelength inside the densest material: the page aims for its target, and a
// result below `GUIDE_RESOLVED` isn't converged. Second-order cells need only a few.
const GUIDE_TARGET = 4.5;
const GUIDE_RESOLVED = 4;

// The cell a guide's mesh uses: finer than asked if this wavelength needs it, as far as
// the solver's limit on unknowns allows, which is about 310 cells across.
export function guideCell(s) {
  const p = s.plane.guide2d;
  const { size, pml } = guideGeometry(p.bend, false, guideScale(s));
  return Math.min(p.cell, Math.max(s.wavelength / (GUIDES[p.bend].n * GUIDE_TARGET), (size + 2 * pml) / 310));
}
export const guideResolved = (s) => s.wavelength / (GUIDES[s.plane.guide2d.bend].n * guideCell(s)) >= GUIDE_RESOLVED;

// A guide's shapes, source and output monitor, or those of a straight one to compare with.
export function guideGeometry(bend, straight, scale = 1) {
  return scaled(unitGuide(bend, straight), scale);
}

function unitGuide(bend, straight) {
  if (bend === 'slab') {
    const { size, x, y, width, n, pml } = SLAB;
    const far = size + pml + 1;
    const end = straight ? far : x + width / 2;
    const rect = (x0, y0, x1, y1) => ({ type: 'rect', x: (x0 + x1) / 2, y: (y0 + y1) / 2, w: x1 - x0, h: y1 - y0, n });
    const shapes = [rect(-pml - 1, y - width / 2, end, y + width / 2)];
    if (!straight) shapes.push(rect(x - width / 2, y - width / 2, x + width / 2, far));
    return {
      size, pml, shapes,
      source: { type: 'port', x: 0.5, y0: y - 1.2, y1: y + 1.2 },
      mark: [0.5, y],
      monitor: straight ? { x0: 3.5, y0: y - 1, x1: 3.5, y1: y + 1 } : { x0: x + 1, y0: 3.5, x1: x - 1, y1: 3.5 },
    };
  }
  const { arm, clad, pml } = LATTICE;
  const size = arm + 1 + clad;
  const [i0, j0] = [arm, clad];
  const removed = (i, j) => {
    if (straight) return j === j0;
    if (bend === 'sharp') return (j === j0 && i <= i0) || (i === i0 && j >= j0);
    return (j === j0 && i < i0) || (i === i0 - 1 && j === j0 + 1) || (i === i0 && j > j0);
  };
  const shapes = [];
  for (let i = -pml - 1; i < size + pml + 1; i++) {
    for (let j = -pml - 1; j < size + pml + 1; j++) {
      if (!removed(i, j)) shapes.push({ type: 'circle', x: i + 0.5, y: j + 0.5, w: 0.36, h: 0.36, n: 3.4 });
    }
  }
  const [x, y, span, far] = [i0 + 0.5, j0 + 0.5, 3, size - 1.03];
  return {
    size, pml, shapes,
    source: { type: 'line', x: 1.5, y },
    mark: [1.5, y],
    monitor: straight ? { x0: far, y0: y - span, x1: far, y1: y + span } : { x0: x + span, y0: far, x1: x - span, y1: far },
  };
}

// Each preset starts on the view that shows it best. A `planeSetup` replaces that view's domain.
export const PRESETS = {
  'Mirror with a resonant gap': {
    analysis: 'spectrum', wavelength: 1,
    layers: [air(), H(), L(), H(), L(), H(), layer('Silica glass', 0.345), H(), L(), H(), L(), H(), air()],
  },
  'Quarter-wave mirror': {
    analysis: 'spectrum', wavelength: 1,
    layers: [air(), H(), L(), H(), L(), H(), L(), H(), L(), H(), layer('Silica glass', 0.5)],
  },
  'Anti-reflection coating': {
    analysis: 'spectrum', wavelength: 0.55,
    layers: [air(), layer('Magnesium fluoride', 0.1), layer('Silica glass', 0.6)],
  },
  'Metal film': {
    analysis: 'field', wavelength: 1,
    layers: [air(), layer('Gold-like metal', 0.04), air()],
  },
  'Layered resonator': {
    analysis: 'cavity',
    layers: [air(0.6), layer('Silicon', 0.3), air(0.6)],
  },
  'Photonic crystal': {
    analysis: 'bands',
    layers: [air(0.5), layer('Silicon', 0.25)],
  },
  'Slab waveguide': {
    analysis: 'waveguide', wavelength: 1,
    layers: [layer('Silica glass', 3), layer('Silicon nitride', 0.8), layer('Silica glass', 3)],
  },
  'Scattering by a rod': { analysis: 'scatter2d', source: 'plane_wave', wavelength: 1, angle: 0 },
  'Diffraction grating': { analysis: 'grating2d', wavelength: 1, angle: 20 },
  'Sharp bend in a photonic crystal': { analysis: 'guide2d', wavelength: GUIDES.sharp.wavelength, planeSetup: guideSetup('sharp') },
  'Bend in an ordinary waveguide': { analysis: 'guide2d', wavelength: GUIDES.slab.wavelength, planeSetup: guideSetup('slab') },
  'Empty square cavity': { analysis: 'cavity2d', planeSetup: { width: 1, height: 1, cell: 0.04, background: 1, walls: 'pec', target: null, sourceX: 0.5, sourceY: 0.5, shapes: [] } },
  'Point-defect cavity': { analysis: 'cavity2d', planeSetup: defectCavity(5) },
  'Lattice of rods': { analysis: 'bands2d' },
  'Design a mirror': { analysis: 'mirror', wavelength: 1 },
  'Design a broadband coating': { analysis: 'antireflection' },
  'Design a waveguide bend': { analysis: 'bend2d' },
  'Design a cavity': { analysis: 'cavityDesign', wavelength: 1 },
};

// The kinds of problem. Driven views solve one linear system for a source the visitor
// picks; eigenmode views have no source and solve for the modes themselves; design
// views repeat driven solves, using their gradients to change the structure.
export const PROBLEMS = {
  driven: { name: 'Driven response', views: ['spectrum', 'field', 'scatter2d', 'guide2d', 'grating2d'] },
  eigen: { name: 'Eigenmodes', views: ['cavity', 'bands', 'waveguide', 'cavity2d', 'bands2d'] },
  design: { name: 'Inverse design', views: ['mirror', 'antireflection', 'bend2d', 'cavityDesign'] },
};
export const problemOf = (analysis) => Object.keys(PROBLEMS).find((p) => PROBLEMS[p].views.includes(analysis));

export function defaultState() {
  return {
    layers: [],
    analysis: 'spectrum',
    polarization: 'te',
    wavelength: 1,
    angle: 0,
    nModes: 6,
    order: 4,
    cellsPerWavelength: 6,
    spectrum: { from: 0.5, to: 1.8 },
    bc: {
      cavity: { left: 'pec', right: 'pec' },
      waveguide: { left: 'pec', right: 'pec' },
      spectrum: { right: 'absorbing' },
      field: { left: 'absorbing', right: 'absorbing' },
    },
    source: 'plane_wave',
    sourceX: 0.5,
    bandPick: { k: 0.5, band: 0 },
    plane: {
      scatter2d: {
        width: 3, height: 3, cell: 0.06, background: 1, walls: 'absorbing', sourceX: 0.75, sourceY: 1.5, averageQ: null,
        shapes: [{ type: 'circle', x: 1.5, y: 1.5, w: 1, h: 1, n: 2 }],
      },
      guide2d: guideSetup('sharp'),
      // The substrate runs out through the layer below, as anything reaching it must.
      grating2d: {
        width: 1.5, height: 2, cell: 0.025, background: 1, walls: 'periodic',
        shapes: [{ type: 'rect', x: 0.75, y: -1.25, w: 3.5, h: 3.5, n: 1.5 }, { type: 'rect', x: 0.75, y: 0.75, w: 0.75, h: 0.5, n: 1.5 }],
      },
      bend2d: { background: 1, walls: 'absorbing', shapes: [] },
      cavityDesign: { background: 1, walls: 'absorbing', shapes: [] },
      cavity2d: { width: 1, height: 1, cell: 0.04, background: 1, walls: 'pec', target: null, sourceX: 0.5, sourceY: 0.5, shapes: [] },
      bands2d: {
        width: 1, height: 1, cell: 0.04, background: 1, walls: 'periodic', pick: { s: 1, band: 0 },
        shapes: [{ type: 'circle', x: 0.5, y: 0.5, w: 0.4, h: 0.4, n: 3.4641 }],
      },
    },
    design: {
      mirror: { thickness: 2, low: 'Silica glass', high: 'Titanium dioxide', left: 'Air', right: 'Air', radius: 20, seed: 1 },
      antireflection: { thickness: 0.6, low: 'Magnesium fluoride', high: 'Titanium dioxide', left: 'Air', right: 'Silica glass', radius: 15, seed: 1, from: 0.4, to: 0.8 },
      bend2d: { radius: 3, seed: 1, from: 0.95, to: 1.05 },
      cavityDesign: { cavity: true, seed: 1 },
    },
    ...structuredClone(PRESETS['Mirror with a resonant gap']),
  };
}

export const totalLength = (s) => s.layers.reduce((sum, l) => sum + l.thickness, 0);

function toLayer(l) {
  return { thickness: l.thickness, eps: [l.n * l.n - l.k * l.k, 2 * l.n * l.k], mu: [1, 0] };
}

const material = (n) => ({ eps: n * n, mu: 1 });

// The Bloch phases per period at a place `s` on the path from Gamma to X to M to Gamma.
export function blochAt(s) {
  if (s <= 1) return [Math.PI * s, 0];
  if (s <= 2) return [Math.PI, Math.PI * (s - 1)];
  return [Math.PI * (3 - s), Math.PI * (3 - s)];
}

export const toShape = (c) => ({
  material: material(c.n),
  ...(c.type === 'circle'
    ? { type: 'circle', x: c.x, y: c.y, r: c.w / 2 }
    : { type: 'rect', x0: c.x - c.w / 2, y0: c.y - c.h / 2, x1: c.x + c.w / 2, y1: c.y + c.h / 2 }),
});

// What a pointer at (x, y) would grab, with later shapes on top: an outline within `reach`
// resizes, and inside moves. A small shape keeps a middle to hold it by.
export function grabShape(shapes, x, y, reach) {
  for (let index = shapes.length - 1; index >= 0; index--) {
    const s = shapes[index];
    const [dx, dy] = [x - s.x, y - s.y];
    if (s.type === 'circle') {
      const rim = Math.hypot(dx, dy) - s.w / 2;
      if (rim < reach && rim > -Math.min(reach, s.w / 4)) return { index, part: 'rim' };
      if (rim < 0) return { index, part: 'body', dx, dy };
    } else {
      const [ex, ey] = [Math.abs(dx) - s.w / 2, Math.abs(dy) - s.h / 2];
      if (ex >= reach || ey >= reach) continue;
      const sideX = ex > -Math.min(reach, s.w / 4) ? Math.sign(dx) : 0;
      const sideY = ey > -Math.min(reach, s.h / 4) ? Math.sign(dy) : 0;
      return sideX || sideY ? { index, part: 'edge', sideX, sideY } : { index, part: 'body', dx, dy };
    }
  }
  return null;
}

// Applies a grab to its shape. Dragging a rectangle's edge leaves the opposite edge in place.
export function dragShape(s, grab, x, y) {
  const [snap, smallest] = [(v) => Math.round(v * 1000) / 1000, 0.005];
  if (grab.part === 'body') {
    [s.x, s.y] = [snap(x - grab.dx), snap(y - grab.dy)];
  } else if (grab.part === 'rim') {
    s.w = Math.max(smallest, snap(2 * Math.hypot(x - s.x, y - s.y)));
  } else {
    const resize = (centre, size, side, at) => {
      const fixed = centre - (side * size) / 2;
      const next = Math.max(smallest, snap(side * (at - fixed)));
      return [+(fixed + (side * next) / 2).toFixed(6), next];
    };
    if (grab.sideX) [s.x, s.w] = resize(s.x, s.w, grab.sideX, x);
    if (grab.sideY) [s.y, s.h] = resize(s.y, s.h, grab.sideY, y);
  }
}

// A straight guide, then the bend. The same source in each gives the bend's transmission.
function guideRequests(s) {
  const p = s.plane.guide2d;
  return [true, false].map((straight) => {
    const g = guideGeometry(p.bend, straight, guideScale(s));
    return {
      width: g.size,
      height: g.size,
      cell: guideCell(s),
      background: material(1),
      polarization: 'te',
      walls: 'absorbing',
      shapes: g.shapes.map(toShape),
      monitors: [g.monitor],
      analysis: { type: 'driven', wavelength: s.wavelength, source: g.source, pml: g.pml },
    };
  });
}

// Where a grating's lines go: reflection just above the structure with the sheet above that,
// and transmission halfway down the space below it.
export function gratingLines(s) {
  const p = s.plane.grating2d;
  const inside = p.shapes.map((c) => [c.y - c.h / 2, c.y + c.h / 2]);
  const top = Math.max(0, ...inside.map(([, high]) => high));
  const lows = inside.map(([low]) => low).filter((low) => low > 0);
  const wavelength = s.wavelength / p.background;
  const reflect = top + Math.max(2 * p.cell, 0.05 * wavelength);
  return { reflect, source: reflect + 0.5 * wavelength, transmit: Math.min(p.height, ...lows) / 2 };
}

function planeRequests(s) {
  if (s.analysis === 'guide2d') return guideRequests(s);
  const p = s.plane[s.analysis];
  const base = {
    width: p.width,
    height: p.height,
    cell: p.cell,
    background: material(p.background),
    polarization: s.polarization,
    walls: p.walls,
    shapes: p.shapes.map(toShape),
  };
  const targetK0 = p.target ? (2 * Math.PI) / p.target : null;
  const modes = (bloch) => ({ ...base, analysis: { type: 'cavity', n_modes: s.nModes, bloch, target_k0: targetK0 } });
  if (s.analysis === 'cavity2d') {
    if (p.walls !== 'absorbing') return [modes([0, 0])];
    // A leaky mode needs somewhere to look, so start from the domain's fundamental. A quarter
    // wavelength of layer changed the rod crystal's Q by under 0.2% against a whole one.
    const wavelength = p.target || 2 * Math.max(p.width, p.height);
    return [{ ...base, analysis: { type: 'resonance', target_k0: (2 * Math.PI) / wavelength, x: p.sourceX, y: p.sourceY, pml: wavelength / 4 } }];
  }
  if (s.analysis === 'grating2d') {
    const { reflect, source, transmit } = gratingLines(s);
    const pml = s.wavelength / p.background;
    return [{ ...base, analysis: { type: 'grating', wavelength: s.wavelength, angle_deg: s.angle, pml, source_y: source, reflect_y: reflect, transmit_y: transmit } }];
  }
  if (s.analysis === 'scatter2d') {
    const source = s.source === 'line' ? { type: 'line', x: p.sourceX, y: p.sourceY } : { type: 'plane_wave', angle_deg: s.angle };
    // In tests, this much layer on the mesh's cells left errors well under 1%.
    const averageQ = s.source === 'line' && p.averageQ ? p.averageQ : null;
    return [{ ...base, analysis: { type: 'driven', wavelength: s.wavelength, source, pml: 0.6 * s.wavelength, average_q: averageQ } }];
  }
  return [{ ...base, analysis: { type: 'bands', n_bands: s.nModes, n_k: 5 } }, modes(blochAt(p.pick.s))];
}

// One or two solves: a chart's data, then the field of the point picked on it.
export function buildRequests(s) {
  if (s.plane[s.analysis]) return planeRequests(s);
  const optical = s.layers.reduce((sum, l) => sum + l.thickness * Math.max(l.n, 1), 0);
  const make = (analysis, left, right, referenceWavelength) => ({
    layers: s.layers.map(toLayer),
    polarization: s.polarization,
    left,
    right,
    mesh: { order: s.order, max_cell: referenceWavelength / s.cellsPerWavelength },
    analysis,
  });
  const nearestModes = 2 * optical / s.nModes;
  const bc = s.bc[s.analysis];
  const planeWave = { type: 'plane_wave' };
  switch (s.analysis) {
    case 'cavity':
      return [make({ type: 'cavity', n_modes: s.nModes }, bc.left, bc.right, nearestModes)];
    case 'bands':
      return [
        make({ type: 'bands', n_bands: s.nModes, n_k: 41 }, 'periodic', 'periodic', nearestModes),
        make({ type: 'cavity', n_modes: s.nModes, bloch_phase: Math.PI * s.bandPick.k }, 'periodic', 'periodic', nearestModes),
      ];
    case 'waveguide':
      return [make({ type: 'waveguide', wavelength: s.wavelength, n_modes: s.nModes }, bc.left, bc.right, s.wavelength)];
    case 'spectrum': {
      const { from, to } = s.spectrum;
      const wavelengths = Array.from({ length: 301 }, (_, i) => from + ((to - from) * i) / 300);
      return [
        make({ type: 'sweep', wavelengths, angle_deg: s.angle }, 'port', bc.right, Math.min(from, s.wavelength)),
        make({ type: 'driven', wavelength: s.wavelength, angle_deg: s.angle, source: planeWave }, 'port', bc.right, s.wavelength),
      ];
    }
    default: {
      const sheet = s.source === 'line';
      const source = sheet ? { type: 'line', x: s.sourceX } : planeWave;
      return [make({ type: 'driven', wavelength: s.wavelength, angle_deg: s.angle, source }, sheet ? bc.left : 'port', bc.right, s.wavelength)];
    }
  }
}

// Geometry lives in the URL, so a link shares a structure with no server involved.
export const encodeState = (s) => encodeURIComponent(JSON.stringify(s));

export function decodeState(text) {
  try {
    const parsed = JSON.parse(decodeURIComponent(text));
    const valid = Array.isArray(parsed.layers) && parsed.layers.length > 0 &&
      parsed.layers.every((l) => [l.n, l.k, l.thickness].every(Number.isFinite));
    if (!valid) return null;
    const fresh = defaultState();
    // Links from before a view existed have no settings for it.
    return { ...fresh, ...parsed, plane: { ...fresh.plane, ...parsed.plane }, design: { ...fresh.design, ...parsed.design } };
  } catch {
    return null;
  }
}
