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
  'Empty square cavity': { analysis: 'cavity2d', planeSetup: { width: 1, height: 1, cell: 0.04, background: 1, walls: 'pec', target: null, sourceX: 0.5, sourceY: 0.5, shapes: [] } },
  'Point-defect cavity': { analysis: 'cavity2d', planeSetup: defectCavity(5) },
  'Lattice of rods': { analysis: 'bands2d' },
  'Design a mirror': { analysis: 'mirror', wavelength: 1 },
  'Design a broadband coating': { analysis: 'antireflection' },
};

// The kinds of problem. Driven views solve one linear system for a source the visitor
// picks; eigenmode views have no source and solve for the modes themselves; design
// views repeat driven solves, using their gradients to change the structure.
export const PROBLEMS = {
  driven: { name: 'Driven response', views: ['spectrum', 'field', 'scatter2d'] },
  eigen: { name: 'Eigenmodes', views: ['cavity', 'bands', 'waveguide', 'cavity2d', 'bands2d'] },
  design: { name: 'Inverse design', views: ['mirror', 'antireflection'] },
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
      cavity2d: { width: 1, height: 1, cell: 0.04, background: 1, walls: 'pec', target: null, sourceX: 0.5, sourceY: 0.5, shapes: [] },
      bands2d: {
        width: 1, height: 1, cell: 0.04, background: 1, walls: 'periodic', pick: { s: 1, band: 0 },
        shapes: [{ type: 'circle', x: 0.5, y: 0.5, w: 0.4, h: 0.4, n: 3.4641 }],
      },
    },
    design: {
      mirror: { thickness: 2, low: 'Silica glass', high: 'Titanium dioxide', left: 'Air', right: 'Air', radius: 20, seed: 1 },
      antireflection: { thickness: 0.6, low: 'Magnesium fluoride', high: 'Titanium dioxide', left: 'Air', right: 'Silica glass', radius: 15, seed: 1, from: 0.4, to: 0.8 },
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

function planeRequests(s) {
  const p = s.plane[s.analysis];
  const base = {
    width: p.width,
    height: p.height,
    cell: p.cell,
    background: material(p.background),
    polarization: s.polarization,
    walls: p.walls,
    shapes: p.shapes.map((c) => ({
      material: material(c.n),
      ...(c.type === 'circle'
        ? { type: 'circle', x: c.x, y: c.y, r: c.w / 2 }
        : { type: 'rect', x0: c.x - c.w / 2, y0: c.y - c.h / 2, x1: c.x + c.w / 2, y1: c.y + c.h / 2 }),
    })),
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
    return valid ? { ...defaultState(), ...parsed } : null;
  } catch {
    return null;
  }
}
