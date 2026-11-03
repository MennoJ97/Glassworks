//! Thin browser interface: JSON request in, JSON metadata plus typed arrays out.

use glassworks_core::output::Output;
use glassworks_core::problem::{Request, Request2d};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct Solution {
    meta: String,
    channels: Vec<(String, Vec<f64>)>,
}

#[wasm_bindgen]
impl Solution {
    pub fn meta(&self) -> String {
        self.meta.clone()
    }
    pub fn channel_names(&self) -> Vec<String> {
        self.channels.iter().map(|(n, _)| n.clone()).collect()
    }
    /// Copies a channel into a `Float64Array`, which the worker then transfers.
    pub fn channel(&self, index: usize) -> Vec<f64> {
        self.channels[index].1.clone()
    }
}

fn solution(out: Result<Output, String>) -> Result<Solution, JsError> {
    let out = out.map_err(|e| JsError::new(&e))?;
    Ok(Solution { meta: serde_json::to_string(&out.meta)?, channels: out.channels })
}

#[wasm_bindgen]
pub fn solve(request_json: &str) -> Result<Solution, JsError> {
    let req: Request = serde_json::from_str(request_json)?;
    solution(glassworks_core::solve1d::solve(&req))
}

#[wasm_bindgen]
pub fn solve_2d(request_json: &str) -> Result<Solution, JsError> {
    let req: Request2d = serde_json::from_str(request_json)?;
    solution(glassworks_core::solve2d::solve(&req))
}
