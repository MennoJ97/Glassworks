# Third-party notices

The solver in `pkg/` is compiled from Eigenmenno's own Rust code together with the libraries
below. Each is used under the MIT License, which all of them offer, alone or as one of
their choices. Its text follows the list.

Only libraries whose code ends up in `pkg/` are listed. Tools that run while building, such as
`serde_derive`, `syn` and the tutorial's `marked` and `katex`, ship nothing to visitors.

| Library | Version | Licence offered | Copyright |
|---|---|---|---|
| Rust standard library (`core`, `alloc`, `std`) | 1.99 | MIT or Apache-2.0 | The Rust Project Contributors |
| [wasm-bindgen](https://github.com/wasm-bindgen/wasm-bindgen), with the JavaScript it generates in `pkg/` | 0.2.130 | MIT or Apache-2.0 | 2014 Alex Crichton |
| wasm-bindgen-shared | 0.2.130 | MIT or Apache-2.0 | 2014 Alex Crichton |
| [cfg-if](https://github.com/rust-lang/cfg-if) | 1.0.5 | MIT or Apache-2.0 | 2014 Alex Crichton |
| [num-complex](https://github.com/rust-num/num-complex) | 0.4.6 | MIT or Apache-2.0 | 2014 The Rust Project Developers |
| [num-traits](https://github.com/rust-num/num-traits) | 0.2.19 | MIT or Apache-2.0 | 2014 The Rust Project Developers |
| [serde](https://github.com/serde-rs/serde) and serde_core | 1.0.229 | MIT or Apache-2.0 | Erick Tryzelaar and David Tolnay |
| [serde_json](https://github.com/serde-rs/json) | 1.0.154 | MIT or Apache-2.0 | Erick Tryzelaar and David Tolnay |
| [itoa](https://github.com/dtolnay/itoa) | 1.0.18 | MIT or Apache-2.0 | David Tolnay |
| [zmij](https://github.com/dtolnay/zmij) | 1.0.23 | MIT | David Tolnay |
| [memchr](https://github.com/BurntSushi/memchr) | 2.8.3 | MIT or Unlicense | 2015 Andrew Gallant |

## The MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
