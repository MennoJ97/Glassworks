// Finds the names in the source that the tutorial's code index defines, so a reader can
// click a function, type or constant and land on its definition. Like highlight.mjs it
// reads shapes, not meaning, so it only links a name that can mean one thing: `.solve(`
// goes to a method and `solve(` to a free function, a name defined in the file itself wins
// over one defined elsewhere, and anything still defined twice stays plain text.

const WORD = /[A-Za-z_$][\w$]*/g;
const LINKED = new Set(['fn', 'ty', 'co']);

// Sets `links` on each file of the index: a map from a name's offset in the source to the
// definition it means, as { length, file, start, end }.
export function linkDefinitions(index) {
  const byName = new Map();
  for (const file of index) {
    for (const item of file.items) {
      for (const name of item.names) {
        const key = `${file.lang} ${name}`;
        byName.set(key, [...(byName.get(key) ?? []), { file, ...item }]);
      }
    }
  }
  // Another file's definition is in reach when it's `pub` or exported, or when both files
  // are the solver's own modules, where a child module may call its parent's private items.
  // Each test file is a crate of its own.
  const shared = (d) => /^\s*(pub|export)\b/.test(d.file.lines[d.line - 1]);
  const inSolver = (file) => file.lang === 'rs' && !file.path.includes('/tests/');
  const reaches = (file, d) => d.file === file || shared(d) || (inSolver(file) && inSolver(d.file));
  const find = (file, name, method) => {
    const all = (byName.get(`${file.lang} ${name}`) ?? []).filter((d) => (name.includes('::') || method === d.names.length > 1) && reaches(file, d));
    const own = all.filter((d) => d.file === file);
    return own.length === 1 ? own[0] : all.length === 1 ? all[0] : null;
  };
  for (const file of index) {
    file.links = new Map();
    const src = file.lines.join('\n');
    // Local names, such as a closure `let rect = |...|`, which a bare `rect(` means instead
    // from there to the end of the function binding it. Not a pattern: `let Some(x) = ...`.
    const bound = file.lang === 'rs' ? /\blet\s+(?:mut\s+)?([a-z_]\w*)\b(?!\s*(?:::|\(|\{))/g : /^[ \t]+(?:const|let)\s+([A-Za-z_$][\w$]*)/gm;
    const lineAt = (pos) => file.offsets.findLastIndex((o) => o <= pos) + 1;
    const locals = [...src.matchAll(bound)].map((m) => {
      const from = lineAt(m.index);
      const scope = file.items.filter((it) => it.line <= from && from <= it.end).sort((a, b) => a.end - a.start - (b.end - b.start))[0];
      return { name: m[1], from, to: scope?.end ?? file.lines.length };
    });
    const shadowed = (name, n) => locals.some((l) => l.name === name && l.from <= n && n <= l.to);
    let line = 0;
    for (const { 0: name, index: at } of src.matchAll(WORD)) {
      if (!LINKED.has(file.colours[at]) || /[\w$]/.test(src[at - 1] ?? '')) continue;
      while (file.offsets[line + 1] <= at) line++;
      const path = /(\w+)::$/.exec(src.slice(Math.max(0, at - 64), at))?.[1];
      // `Type::name` is that type's, and nothing else: `Vec::new` isn't the solver's `new`.
      // `Self::name` and `x.name` are some type's method; `module::name` and `name` are free.
      // A method on a bracket, as in `(0.0..1.0).contains(&x)`, is on a value whose type
      // the text doesn't say, and is usually the standard library's.
      if (src[at - 1] === '.' && /[)\]]/.test(src[at - 2])) continue;
      if (!path && src[at - 1] !== '.' && shadowed(name, line + 1)) continue;
      const typed = path && path !== 'Self' && /^[A-Z]/.test(path);
      const def = typed ? find(file, `${path}::${name}`) : find(file, name, path === 'Self' || src[at - 1] === '.');
      if (def && !(def.file === file && def.line === line + 1)) file.links.set(at, { length: name.length, ...def });
    }
  }
}
