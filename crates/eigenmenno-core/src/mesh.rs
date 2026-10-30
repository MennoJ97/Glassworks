//! Geometry only: no knowledge of elements or physics.

/// `cells` holds `verts_per_cell` vertex ids per cell, so triangles can share
/// this layout when 2D arrives.
pub struct Mesh {
    pub verts_per_cell: usize,
    pub coords: Vec<f64>,
    pub cells: Vec<usize>,
    /// Material region of each cell.
    pub region: Vec<usize>,
}

impl Mesh {
    /// A line through consecutive layers, with a vertex on every interface so
    /// material jumps fall exactly on cell boundaries.
    pub fn layered_line(thickness: &[f64], max_cell: &[f64]) -> Self {
        let mut m = Mesh { verts_per_cell: 2, coords: vec![0.0], cells: vec![], region: vec![] };
        let mut start = 0.0;
        for (r, (&t, &h)) in thickness.iter().zip(max_cell).enumerate() {
            let n = (t / h).ceil().max(1.0) as usize;
            for i in 1..=n {
                let v = m.coords.len();
                m.coords.push(start + t * i as f64 / n as f64);
                m.cells.extend([v - 1, v]);
                m.region.push(r);
            }
            start += t;
        }
        m
    }

    pub fn n_cells(&self) -> usize {
        self.region.len()
    }
}
