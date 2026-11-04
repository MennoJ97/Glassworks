//! Geometry only: no knowledge of elements or physics.

use serde::Deserialize;
use std::collections::HashMap;

/// `cells` holds `verts_per_cell` vertex ids per cell, and `coords` holds x, or x
/// and y, per vertex. Lines and triangles share this layout.
pub struct Mesh {
    pub verts_per_cell: usize,
    pub coords: Vec<f64>,
    pub cells: Vec<usize>,
    /// Material region of each cell.
    pub region: Vec<usize>,
}

/// A closed shape in the plane, positioned in the first quadrant's unit cell.
#[derive(Deserialize, Clone, Copy)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Outline {
    Circle { x: f64, y: f64, r: f64 },
    Rect { x0: f64, y0: f64, x1: f64, y1: f64 },
}

impl Outline {
    /// The same shape, moved.
    pub fn shifted(self, dx: f64, dy: f64) -> Self {
        match self {
            Outline::Circle { x, y, r } => Outline::Circle { x: x + dx, y: y + dy, r },
            Outline::Rect { x0, y0, x1, y1 } => {
                Outline::Rect { x0: x0 + dx, y0: y0 + dy, x1: x1 + dx, y1: y1 + dy }
            }
        }
    }

    /// Whether the shape stays within `[0, lx] x [0, ly]`.
    pub fn inside(&self, lx: f64, ly: f64) -> bool {
        match *self {
            Outline::Circle { x, y, r } => x - r >= 0.0 && x + r <= lx && y - r >= 0.0 && y + r <= ly,
            Outline::Rect { x0, y0, x1, y1 } => x0 >= 0.0 && x1 <= lx && y0 >= 0.0 && y1 <= ly,
        }
    }

    /// Negative inside and positive outside. Only its sign and rough size are used. This is
    /// the mesh's hottest call, and `hypot` is about 30 times slower than this in WebAssembly.
    pub(crate) fn level(&self, px: f64, py: f64) -> f64 {
        match *self {
            Outline::Circle { x, y, r } => ((px - x) * (px - x) + (py - y) * (py - y)).sqrt() - r,
            Outline::Rect { x0, y0, x1, y1 } => (x0 - px).max(px - x1).max(y0 - py).max(py - y1),
        }
    }

    /// The nearest point on the outline.
    fn project(&self, px: f64, py: f64) -> (f64, f64) {
        match *self {
            Outline::Circle { x, y, r } => {
                let distance = (px - x).hypot(py - y).max(1e-300);
                (x + (px - x) * r / distance, y + (py - y) * r / distance)
            }
            Outline::Rect { x0, y0, x1, y1 } => {
                let clamped = (px.clamp(x0, x1), py.clamp(y0, y1));
                if clamped != (px, py) {
                    return clamped;
                }
                let gaps = [px - x0, x1 - px, py - y0, y1 - py];
                let least = gaps.iter().copied().fold(f64::MAX, f64::min);
                if least == gaps[0] {
                    (x0, py)
                } else if least == gaps[1] {
                    (x1, py)
                } else if least == gaps[2] {
                    (px, y0)
                } else {
                    (px, y1)
                }
            }
        }
    }

    /// `[x0, y0, x1, y1]` around the shape.
    fn bounds(&self) -> [f64; 4] {
        match *self {
            Outline::Circle { x, y, r } => [x - r, y - r, x + r, y + r],
            Outline::Rect { x0, y0, x1, y1 } => [x0, y0, x1, y1],
        }
    }

    fn corners(&self) -> Vec<(f64, f64)> {
        match *self {
            Outline::Circle { .. } => vec![],
            Outline::Rect { x0, y0, x1, y1 } => vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)],
        }
    }
}

impl Mesh {
    /// A line through consecutive layers, with a vertex on every interface so
    /// material jumps fall exactly on cell boundaries.
    pub fn layered_line(thickness: &[f64], max_cell: &[f64]) -> Self {
        let mut mesh = Mesh { verts_per_cell: 2, coords: vec![0.0], cells: vec![], region: vec![] };
        let mut start = 0.0;
        for (layer, (&width, &max_h)) in thickness.iter().zip(max_cell).enumerate() {
            let n_cells = (width / max_h).ceil().max(1.0) as usize;
            for i in 1..=n_cells {
                let vertex = mesh.coords.len();
                mesh.coords.push(start + width * i as f64 / n_cells as f64);
                mesh.cells.extend([vertex - 1, vertex]);
                mesh.region.push(layer);
            }
            start += width;
        }
        mesh
    }

    /// Triangles on an `nx` by `ny` grid over `[0, lx] x [0, ly]`, fitted to the
    /// outlines. Region 0 is the background, and region `i + 1` is outline `i`, with
    /// later outlines on top. Outer-edge nodes never move, so opposite edges stay
    /// matched for periodic ties.
    pub fn fitted_grid(nx: usize, ny: usize, lx: f64, ly: f64, outlines: &[Outline]) -> Self {
        let mut mesh = Mesh { verts_per_cell: 3, coords: vec![], cells: vec![], region: vec![] };
        for j in 0..=ny {
            for i in 0..=nx {
                mesh.coords.extend([lx * i as f64 / nx as f64, ly * j as f64 / ny as f64]);
            }
        }
        mesh.snap(nx, ny, lx, ly, outlines);
        let crosses = |mesh: &Mesh, from: usize, to: usize| mesh.crosses(from, to, outlines);
        let vertex_at = |i: usize, j: usize| i + j * (nx + 1);
        for j in 0..ny {
            for i in 0..nx {
                // The square's corners, counter-clockwise from its lower left.
                let (ll, lr, ur, ul) =
                    (vertex_at(i, j), vertex_at(i + 1, j), vertex_at(i + 1, j + 1), vertex_at(i, j + 1));
                // Alternating diagonals give the mesh no preferred direction. With no
                // grid line crossing an outline, one of the two diagonals is always clear.
                let rising = if (i + j) % 2 == 0 { !crosses(&mesh, ll, ur) } else { crosses(&mesh, lr, ul) };
                if rising {
                    mesh.cells.extend([ll, lr, ur, ll, ur, ul]);
                } else {
                    mesh.cells.extend([ll, lr, ul, lr, ur, ul]);
                }
            }
        }
        mesh.flip_to_delaunay(outlines);
        mesh.region = (0..mesh.cells.len() / 3)
            .map(|cell| {
                let (x, y) = mesh.centroid(cell);
                // Where an outline meets the domain's edge, whose nodes can't move, a cell's
                // centroid can lie exactly on it. Noise mustn't choose its material.
                let tolerance = 1e-9 * lx / nx as f64;
                outlines.iter().rposition(|outline| outline.level(x, y) < -tolerance).map_or(0, |i| i + 1)
            })
            .collect();
        mesh
    }

    /// `(x, y)` of a vertex of a 2D mesh.
    pub fn vertex(&self, v: usize) -> (f64, f64) {
        (self.coords[2 * v], self.coords[2 * v + 1])
    }

    pub fn centroid(&self, cell: usize) -> (f64, f64) {
        let corners = &self.cells[3 * cell..3 * cell + 3];
        let (x, y) =
            corners.iter().fold((0.0, 0.0), |sum, &v| (sum.0 + self.vertex(v).0, sum.1 + self.vertex(v).1));
        (x / 3.0, y / 3.0)
    }

    /// Whether an edge between two vertices would cross an outline.
    fn crosses(&self, from: usize, to: usize, outlines: &[Outline]) -> bool {
        let ((x0, y0), (x1, y1)) = (self.vertex(from), self.vertex(to));
        outlines.iter().any(|outline| {
            let (level0, level1) = (outline.level(x0, y0), outline.level(x1, y1));
            level0 * level1 < 0.0 && level0.abs().min(level1.abs()) > 1e-9
        })
    }

    /// Moves grid nodes onto each outline until no grid line crosses it. Corners go
    /// first, then nodes in order of distance, skipping any whose crossings an earlier
    /// move already resolved. Nodes very close to an outline join it regardless, or
    /// they'd form slivers with the snapped nodes beside them.
    fn snap(&mut self, nx: usize, ny: usize, lx: f64, ly: f64, outlines: &[Outline]) {
        let n_vertices = self.coords.len() / 2;
        // Only interior nodes move, so opposite outer edges stay matched.
        let free = |v: usize| {
            let (i, j) = (v % (nx + 1), v / (nx + 1));
            i != 0 && j != 0 && i != nx && j != ny
        };
        // Outer nodes never need testing as a mover, but they still count as neighbours.
        let neighbours = |v: usize| {
            let (i, j) = (v % (nx + 1), v / (nx + 1));
            [
                (i > 0).then(|| v - 1),
                (i < nx).then(|| v + 1),
                (j > 0).then(|| v - nx - 1),
                (j < ny).then(|| v + nx + 1),
            ]
        };
        let mut placed = vec![false; n_vertices];
        let (hx, hy) = (lx / nx as f64, ly / ny as f64);
        let close = 0.3 * hx.min(hy);
        // Grid lines one step past the range, the furthest a crossing or a close node can be.
        let lines_near = |low: f64, high: f64, h: f64, last: usize| {
            ((low / h).floor() - 1.0).max(0.0) as usize
                ..=(((high / h).ceil() + 1.0).max(0.0) as usize).min(last)
        };
        for outline in outlines {
            for (x, y) in outline.corners() {
                // Halfway between two nodes goes to the lower one, whatever rounding noise says.
                let (i, j) = (
                    (x / lx * nx as f64 - 1e-9).round() as usize,
                    (y / ly * ny as f64 - 1e-9).round() as usize,
                );
                let v = i.min(nx) + j.min(ny) * (nx + 1);
                if free(v) && !placed[v] {
                    (self.coords[2 * v], self.coords[2 * v + 1]) = (x, y);
                    placed[v] = true;
                }
            }
            // Nodes not yet placed are still on the grid, so only those near the shape matter.
            let [x0, y0, x1, y1] = outline.bounds();
            let nearby = lines_near(y0, y1, hy, ny)
                .flat_map(|j| lines_near(x0, x1, hx, nx).map(move |i| i + j * (nx + 1)));
            let mut order: Vec<usize> = nearby.filter(|&v| free(v) && !placed[v]).collect();
            // A shape on a grid has many nodes at equal distances. Rounding keeps noise from
            // choosing among them, so the mesh is the same at any scale.
            let level = |mesh: &Mesh, v: usize| {
                (outline.level(mesh.coords[2 * v], mesh.coords[2 * v + 1]) / hx * 1e9).round() * hx / 1e9
            };
            order.sort_by(|&a, &b| level(self, a).abs().total_cmp(&level(self, b).abs()));
            for v in order {
                let own = level(self, v);
                let crossing = neighbours(v)
                    .into_iter()
                    .flatten()
                    .any(|neighbour| !placed[neighbour] && own * level(self, neighbour) < 0.0);
                if crossing || own.abs() < close {
                    (self.coords[2 * v], self.coords[2 * v + 1]) =
                        outline.project(self.coords[2 * v], self.coords[2 * v + 1]);
                    placed[v] = true;
                }
            }
        }
    }

    /// Lawson's flips: swaps the diagonal of any pair of triangles whose opposite
    /// angles sum past 180 degrees. Snapping leaves slivers along curves, with three
    /// nearly collinear points, and this removes them. A flip that would put an edge
    /// across an outline is refused, because cells must stay on one side.
    fn flip_to_delaunay(&mut self, outlines: &[Outline]) {
        for _ in 0..100 {
            // Each edge seen so far, with the cell it came from and that cell's opposite vertex.
            let mut first_side: HashMap<[usize; 2], (usize, usize)> = HashMap::new();
            let mut flipped = vec![false; self.cells.len() / 3];
            let mut flips = 0;
            for cell in 0..self.cells.len() / 3 {
                for k in 0..3 {
                    let tri = [self.cells[3 * cell], self.cells[3 * cell + 1], self.cells[3 * cell + 2]];
                    let (a, b, opp) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
                    let key = [a.min(b), a.max(b)];
                    let Some(&(other, far)) = first_side.get(&key) else {
                        first_side.insert(key, (cell, opp));
                        continue;
                    };
                    if flipped[cell] || flipped[other] {
                        continue;
                    }
                    // The cotangent of the angle at `apex` facing the shared edge a-b.
                    let cot = |apex: usize| {
                        let (o, u, v) = (self.vertex(apex), self.vertex(a), self.vertex(b));
                        let (ux, uy, vx, vy) = (u.0 - o.0, u.1 - o.1, v.0 - o.0, v.1 - o.1);
                        (ux * vx + uy * vy) / (ux * vy - uy * vx).abs().max(1e-300)
                    };
                    // A negative sum means the angles at opp and far add up to more than 180°.
                    if cot(opp) + cot(far) < -1e-9 && !self.crosses(opp, far, outlines) {
                        // Keeps the counter-clockwise orientation: a, far, opp and far, b, opp.
                        self.cells[3 * cell..3 * cell + 3].copy_from_slice(&[a, far, opp]);
                        self.cells[3 * other..3 * other + 3].copy_from_slice(&[far, b, opp]);
                        flipped[cell] = true;
                        flipped[other] = true;
                        flips += 1;
                    }
                }
            }
            if flips == 0 {
                break;
            }
        }
    }

    /// Each edge once, as `[low id, high id]`, in the order cells first reach them.
    pub fn edges(&self) -> Vec<[usize; 2]> {
        let mut seen = std::collections::HashSet::new();
        let mut edges = vec![];
        for corners in self.cells.chunks(3) {
            for (a, b) in [(corners[0], corners[1]), (corners[1], corners[2]), (corners[2], corners[0])] {
                let edge = [a.min(b), a.max(b)];
                if seen.insert(edge) {
                    edges.push(edge);
                }
            }
        }
        edges
    }

    pub fn n_cells(&self) -> usize {
        self.region.len()
    }
}
