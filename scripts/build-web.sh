#!/usr/bin/env bash
# Builds the WebAssembly solver into web/pkg. The output of web/ is fully static.
set -euo pipefail
cd "$(dirname "$0")/.."

# Homebrew's Rust has no wasm32 target, so prefer a rustup toolchain if present.
for bin in "$HOME/.rustup/toolchains/stable-"*/bin; do
  [ -d "$bin" ] && export PATH="$bin:$PATH"
done

BINDGEN_VERSION=$(grep -A1 '^name = "wasm-bindgen"$' Cargo.lock | sed -n 's/version = "\(.*\)"/\1/p')
if [ "$(.tools/bin/wasm-bindgen --version 2>/dev/null | cut -d' ' -f2)" != "$BINDGEN_VERSION" ]; then
  cargo install wasm-bindgen-cli --version "$BINDGEN_VERSION" --root .tools --locked
fi

RUSTFLAGS="-C target-feature=+simd128" \
  cargo build --release --target wasm32-unknown-unknown -p glassworks-wasm
.tools/bin/wasm-bindgen target/wasm32-unknown-unknown/release/glassworks_wasm.wasm \
  --out-dir web/pkg --target web --no-typescript

# Generate the tutorial pages from docs/.
[ -d node_modules ] || npm ci
node scripts/build-docs.mjs
ls -l web/pkg
