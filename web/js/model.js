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

// Each preset starts on the view that shows it best.
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
};

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
    ...structuredClone(PRESETS['Mirror with a resonant gap']),
  };
}

export const totalLength = (s) => s.layers.reduce((sum, l) => sum + l.thickness, 0);

function toLayer(l) {
  return { thickness: l.thickness, eps: [l.n * l.n - l.k * l.k, 2 * l.n * l.k], mu: [1, 0] };
}

// One or two solves: a chart's data, then the field of the point picked on it.
export function buildRequests(s) {
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
