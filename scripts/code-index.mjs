// Finds the functions, types and constants in a Rust or JavaScript file, and the lines
// each one spans, so the tutorial can point at them by name. It reads just enough syntax
// to skip comments and strings when matching braces, and nothing more.

// Marks each character as code (0), comment (1) or string (2).
export function classify(src, lang) {
  const kind = new Uint8Array(src.length);
  let i = 0;
  const mark = (to, k) => {
    kind.fill(k, i, to);
    i = to;
  };
  const until = (from, end) => {
    const j = src.indexOf(end, from);
    return j < 0 ? src.length : j + end.length;
  };
  const quoted = (q) => {
    let j = i + 1;
    while (j < src.length && src[j] !== q) j += src[j] === '\\' ? 2 : 1;
    return j + 1;
  };
  while (i < src.length) {
    const c = src[i];
    const next = src[i + 1];
    if (c === '/' && next === '/') mark(until(i, '\n') - 1, 1);
    else if (c === '/' && next === '*') mark(until(i + 2, '*/'), 1);
    else if (c === '"' || (lang === 'js' && (c === "'" || c === '`'))) mark(Math.min(quoted(c), src.length), 2);
    else if (lang === 'rs' && c === 'r' && /^r#*"/.test(src.slice(i, i + 8)) && !/\w/.test(src[i - 1] ?? '')) {
      const hashes = /^r(#*)"/.exec(src.slice(i))[1];
      mark(until(i + hashes.length + 2, `"${hashes}`), 2);
    } else if (lang === 'rs' && c === "'" && (next === '\\' || src[i + 2] === "'")) mark(quoted("'"), 2);
    else i++;
  }
  return kind;
}

const RUST_ITEM = /^(\s*)(?:pub(?:\([^)]*\))?\s+)?(?:(?:const|async|unsafe)\s+)*(fn|struct|enum|trait|mod|type|const|static)\s+([A-Za-z_]\w*)/;
const RUST_IMPL = /^(\s*)impl\b(?:<[^>]*>)?\s+(?:[\w:<>, &']+?\s+for\s+)?([A-Za-z_]\w*)/;
const JS_ITEM = /^(\s*)(?:export\s+)?(?:(?:async\s+)?function\*?\s+(\w+)|class\s+(\w+)|(?:const|let)\s+(\w+)\s*=)/;

// The line an item ends on: its matching closing brace, or the first `;` outside brackets.
function itemEnd(src, kind, from) {
  let depth = 0;
  for (let i = from; i < src.length; i++) {
    if (kind[i]) continue;
    const c = src[i];
    if (c === '(' || c === '[') depth++;
    else if (c === ')' || c === ']') depth--;
    else if (c === ';' && depth === 0) return i;
    else if (c === '{' && depth === 0) {
      let braces = 0;
      for (let j = i; j < src.length; j++) {
        if (kind[j]) continue;
        if (src[j] === '{') braces++;
        else if (src[j] === '}' && --braces === 0) return j;
      }
      return src.length - 1;
    }
  }
  return src.length - 1;
}

// Items as { names, start, end } with 1-based inclusive lines. A method is named both
// `method` and `Type::method`. The span starts at the comments and attributes above it.
export function definitions(src, lang) {
  const kind = classify(src, lang);
  const lines = src.split('\n');
  const offsets = [];
  let at = 0;
  for (const l of lines) {
    offsets.push(at);
    at += l.length + 1;
  }
  const lineOf = (pos) => {
    let lo = 0;
    let hi = offsets.length - 1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (offsets[mid] <= pos) lo = mid;
      else hi = mid - 1;
    }
    return lo;
  };
  const items = [];
  const owners = [];
  lines.forEach((line, n) => {
    if (kind[offsets[n] + line.search(/\S/)]) return;
    let name;
    let owner = false;
    if (lang === 'rs') {
      const impl = RUST_IMPL.exec(line);
      const item = RUST_ITEM.exec(line);
      if (impl) [name, owner] = [impl[2], true];
      else if (item) [name, owner] = [item[3], item[2] === 'trait'];
    } else {
      const item = JS_ITEM.exec(line);
      if (item && (item[2] || item[3] || (item[4] && item[1] === ''))) [name, owner] = [item[2] ?? item[3] ?? item[4], Boolean(item[3])];
    }
    if (!name) return;
    const end = lineOf(itemEnd(src, kind, offsets[n] + line.search(/\S/)));
    let start = n;
    while (start > 0 && /^\s*(\/\/|#\[|\/\*|\*)/.test(lines[start - 1])) start--;
    const parent = owners.findLast((o) => o.start < n && o.end >= end);
    const names = owner && lang === 'rs' && RUST_IMPL.test(line) ? [] : [name];
    if (parent && names.length) names.push(`${parent.name}::${name}`);
    items.push({ names, start: start + 1, end: end + 1 });
    if (owner) owners.push({ name, start: n, end });
  });
  return { items, kind, lines, offsets };
}
