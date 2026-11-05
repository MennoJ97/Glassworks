//! The meshes and orderings drawn in the tutorial's chapters 5, 9 and 10, from the solver's
//! own `rcm`, `Ordering` and `Mesh::fitted_grid`. `scripts/figures.mjs` runs this and turns
//! its JSON into pictures: `cargo run --release -p glassworks-core --example figure_data -- out.json`.

use glassworks_core::assemble::{assemble, Csr, DofMap};
use glassworks_core::banded::rcm;
use glassworks_core::mesh::{Mesh, Outline};
use glassworks_core::space::{H1Line, H1Triangle, Space};
use glassworks_core::sparse::Ordering;
use glassworks_core::ONE;
use serde_json::json;
use std::collections::BTreeSet;

/// The matrix's nonzeros, `(row, column)`, renumbered by `new[old]`.
fn pattern(a: &Csr, new: &[usize]) -> Vec<(usize, usize)> {
    (0..a.n).flat_map(|i| a.col[a.ptr[i]..a.ptr[i + 1]].iter().map(move |&j| (new[i], new[j]))).collect()
}

/// The nonzeros of L, `(row, column)` with row >= column, when the unknowns are eliminated in
/// the order `new`. Eliminating column k connects everything below it in that column, which
/// column k's parent in the elimination tree inherits.
fn factor_pattern(a: &Csr, new: &[usize]) -> Vec<(usize, usize)> {
    let n = a.n;
    let mut below: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    for (i, j) in pattern(a, new) {
        if i > j {
            below[j].insert(i);
        }
    }
    let mut entries = vec![];
    for k in 0..n {
        let column = std::mem::take(&mut below[k]);
        if let Some(&parent) = column.iter().next() {
            below[parent].extend(column.iter().copied().filter(|&i| i != parent));
        }
        entries.push((k, k));
        entries.extend(column.iter().map(|&i| (i, k)));
    }
    entries
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: figure_data <out.json>");

    // Chapter 5: a periodic line of 12 first-order cells, whose last node is tied to its first.
    let cells = 12;
    let mesh = Mesh::layered_line(&[1.0], &[1.0 / cells as f64]);
    let space = H1Line::new(&mesh, 1);
    let dofs = DofMap::new(space.n_dofs(), &[], &[(space.n_dofs() - 1, 0, ONE)]);
    let ring = Csr::from_triplets(dofs.n_free, assemble(&space, &dofs, |_| (ONE, ONE)));
    let natural: Vec<usize> = (0..ring.n).collect();
    let perm = rcm(&ring);
    let mut after = vec![0; ring.n];
    perm.iter().enumerate().for_each(|(k, &old)| after[old] = k);

    // Chapter 10: second-order triangles on a 12 by 12 grid, ordered as banded LU would
    // (reverse Cuthill-McKee) and by nested dissection.
    let mesh = Mesh::fitted_grid(12, 12, 1.0, 1.0, &[]);
    let space = H1Triangle::new(&mesh);
    let dofs = DofMap::new(space.n_dofs(), &[], &[]);
    let a = Csr::from_triplets(dofs.n_free, assemble(&space, &dofs, |_| (ONE, ONE)));
    let positions = space.dof_positions();
    let mut banded = vec![0; a.n];
    rcm(&a).iter().enumerate().for_each(|(k, &old)| banded[old] = k);
    let dissected = Ordering::new(&a, &positions).new_index().to_vec();

    // Chapter 9: a rod on a 16 by 16 grid, as a staircase (each cell takes the material at its
    // centre) and fitted to the outline.
    let rod = Outline::Circle { x: 0.5, y: 0.5, r: 0.3 };
    let mut staircase = Mesh::fitted_grid(16, 16, 1.0, 1.0, &[]);
    staircase.region = (0..staircase.n_cells())
        .map(|c| {
            let (x, y) = staircase.centroid(c);
            usize::from((x - 0.5).powi(2) + (y - 0.5).powi(2) < 0.09)
        })
        .collect();
    let fitted = Mesh::fitted_grid(16, 16, 1.0, 1.0, &[rod]);
    let triangles = |m: &Mesh| json!({ "coords": m.coords, "cells": m.cells, "region": m.region });

    let out = json!({
        "staircase": triangles(&staircase),
        "fitted": triangles(&fitted),
        "ring": { "n": ring.n, "natural": pattern(&ring, &natural), "rcm": pattern(&ring, &after), "order": perm },
        "mesh": {
            "n": a.n,
            "positions": positions,
            "triangles": (0..mesh.n_cells()).map(|c| {
                let mut d = [0; 6];
                space.cell_dofs(c, &mut d);
                [d[0], d[1], d[2]]
            }).collect::<Vec<_>>(),
            "dissected": dissected,
            "banded_factor": factor_pattern(&a, &banded),
            "dissected_factor": factor_pattern(&a, &dissected),
            "nonzeros": a.col.len(),
        },
    });
    std::fs::write(&path, out.to_string()).unwrap();
    eprintln!("wrote {path}");
}
