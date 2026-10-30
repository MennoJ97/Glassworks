//! Thin browser interface: JSON request in, JSON metadata plus typed arrays out.

use eigenmenno_core::problem::Request;
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

#[wasm_bindgen]
pub fn solve(request_json: &str) -> Result<Solution, JsError> {
    let req: Request = serde_json::from_str(request_json)?;
    let out = eigenmenno_core::solve::solve(&req).map_err(|e| JsError::new(&e))?;
    Ok(Solution { meta: serde_json::to_string(&out.meta)?, channels: out.channels })
}
