// Turns docs/*.md into static pages in web/learn, so docs/ stays the only copy of
// the tutorial. Equations become MathML here, at build time, so visitors download
// no math library or fonts. A broken equation or link fails the build.
//
// Each chapter also gets the code it talks about in a panel beside it. A name in `code`
// type that the source defines points at its lines, and so does a comment like
// <!-- code: local_matrices "2.0 / h" .. "h / 2.0" --> after a paragraph or equation:
// the lines of `local_matrices` from the first one containing "2.0 / h" to the next one
// containing "h / 2.0". GitHub hides the comments. One that matches nothing fails the build.
// A link to a source file can name a definition in its title, for link text that is prose.

import { mkdirSync, readFileSync, readdirSync, rmSync, writeFileSync } from 'node:fs';
import { basename, dirname, join, posix } from 'node:path';
import { fileURLToPath } from 'node:url';
import katex from 'katex';
import { Marked } from 'marked';
import { linkDefinitions } from './code-links.mjs';
import { definitions } from './code-index.mjs';
import { highlight } from './highlight.mjs';

// With --check, it renders everything, so a broken link or anchor still fails, but writes nothing.
const check = process.argv.includes('--check');
const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const docs = join(root, 'docs');
const out = join(root, 'web', 'learn');

const escapeHtml = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

// Browsers implementing MathML Core ignore mathvariant="bold", which would make the
// vector E look like the scalar E. Unicode's bold math letters work everywhere.
const boldLetter = (c) => String.fromCodePoint(c.codePointAt(0) + (c >= 'a' ? 0x1d41a - 97 : 0x1d400 - 65));
const math = (tex, displayMode) =>
  katex
    .renderToString(tex, { displayMode, output: 'mathml', throwOnError: true })
    .replace(/<mi mathvariant="bold">([A-Za-z]+)<\/mi>/g, (_, letters) => `<mi>${[...letters].map(boldLetter).join('')}</mi>`);

const sources = [
  ...['mesh', 'space', 'assemble', 'banded', 'eigen', 'problem', 'output', 'solve1d', 'solve1d/gradient', 'sparse', 'sparse/blocked', 'sparse/symmetric', 'sparse/transposed', 'solve2d', 'solve2d/closed', 'solve2d/open', 'solve2d/driven', 'solve2d/leaky', 'solve2d/adjoint', 'solve2d/grating', 'solve2d/dipole', 'solve2d/validate'].map((n) => `crates/glassworks-core/src/${n}.rs`),
  ...['analytic', 'analytic2d', 'gratings', 'hz_design', 'cross_check'].map((n) => `crates/glassworks-core/tests/${n}.rs`),
  ...['optimizer', 'design', 'cavity-design', 'hz-design'].map((n) => `web/js/${n}.js`),
];

const index = sources.map((path) => {
  const src = readFileSync(join(root, path), 'utf8').trimEnd();
  const lang = path.endsWith('.js') ? 'js' : 'rs';
  const defs = definitions(src, lang);
  return { path, name: basename(path), lang, ...defs, colours: highlight(src, defs.kind, lang) };
});
linkDefinitions(index);

// The definitions called `name`, preferring the files the chapter links to. When several
// remain, such as a method each element type has, it keeps those of a type the chapter
// names elsewhere: chapter 3's `local_matrices` is H1Line's, because it never mentions
// H1Triangle.
function lookup(name, files, text = '') {
  const hits = index.flatMap((file) => file.items.filter((it) => it.names.includes(name)).map((it) => ({ file, start: it.start, end: it.end, owner: it.names[1]?.split('::')[0] })));
  const local = hits.filter((h) => files.has(h.file.path));
  const found = local.length ? local : hits;
  const named = found.filter((h) => h.owner && !name.includes('::') && new RegExp(`\\b${h.owner}\\b`).test(text));
  return named.length === 1 ? named : found;
}

// A whole file, named by its path or the end of it.
const wholeFile = (name) => {
  const file = index.find((f) => f.path === name || f.path.endsWith(`/${name}`));
  return file && { file, start: 1, end: file.lines.length, whole: true };
};

// What an inline name refers to, or nothing when the source doesn't define it exactly once.
function resolveName(text, files, page = '', markdown = '') {
  const name = text.replace(/\(\)$/, '');
  if (!/^[\w:./-]+$/.test(name)) return null;
  const hits = lookup(name, files, markdown);
  if (hits.length > 1) console.warn(`${page}: \`${name}\` is defined ${hits.length} times (${hits.map((h) => `${h.file.name}:${h.start}`).join(', ')}), so it isn't linked; link the file you mean from the chapter, or write Type::${name}`);
  return hits.length === 1 ? hits[0] : hits.length ? null : wholeFile(name);
}

// The lines a `<!-- code: ... -->` comment names. Unlike a name in the text, it has to match.
function resolveAnchor(spec, files, page, markdown) {
  const m = /^(\S+)(?:\s+"([^"]+)"(?:\s*\.\.\s*"([^"]+)")?)?$/.exec(spec.trim());
  if (!m) throw new Error(`${page}: can't read <!-- code: ${spec} -->`);
  const [, target, from, to] = m;
  const hits = lookup(target, files, markdown);
  if (hits.length > 1) throw new Error(`${page}: ${target} is ambiguous (${hits.map((h) => `${h.file.name}:${h.start}`).join(', ')}), so link the file you mean from the chapter, or name it Type::${target}`);
  const hit = hits[0] ?? wholeFile(target);
  if (!hit) throw new Error(`${page}: no code named ${target}`);
  if (!from) return hit;
  const find = (text, after) => {
    for (let n = after; n <= hit.end; n++) if (hit.file.lines[n - 1].includes(text)) return n;
    throw new Error(`${page}: "${text}" isn't in ${target} after line ${after}`);
  };
  const start = find(from, hit.start);
  return { file: hit.file, start, end: to ? find(to, start) : start };
}

const refAttr = (refs) => refs.map((r) => (r.whole ? r.file.name : `${r.file.name}:${r.start}-${r.end}`)).join(' ');
const refHref = (r, prefix) => `${prefix}source/${r.file.name}.html${r.whole ? '' : r.start === r.end ? `#L${r.start}` : `#L${r.start}-${r.end}`}`;

// A link to a range of lines, #L12-30, lights all of them and scrolls to the first.
const LIGHT_LINES = `<script>
const lightLines = () => {
  for (const line of document.querySelectorAll('.line.hl')) line.classList.remove('hl');
  const range = /^#L(\\d+)-(\\d+)$/.exec(location.hash);
  if (!range) return;
  for (let n = +range[1]; n <= +range[2]; n++) document.getElementById(\`L\${n}\`)?.classList.add('hl');
  document.getElementById(\`L\${range[1]}\`)?.scrollIntoView();
};
addEventListener('hashchange', lightLines);
lightLines();
</script>`;

// Source lines coloured as an editor would, from highlight.mjs, with each name that
// code-links.mjs resolves linking to its definition. `prefix` leads to the source pages.
function codeLines(file, ids, prefix) {
  return file.lines
    .map((line, n) => {
      const colour = (i) => file.colours[file.offsets[n] + i];
      const linkAt = (i) => file.links.get(file.offsets[n] + i);
      let html = '';
      for (let i = 0; i < line.length; ) {
        const link = linkAt(i);
        if (link) {
          const page = link.file === file && ids ? '' : `${prefix}${link.file.name}.html`;
          html += `<a class="def" href="${page}#L${link.start}-${link.end}"><span class="${colour(i)}">${escapeHtml(line.slice(i, i + link.length))}</span></a>`;
          i += link.length;
          continue;
        }
        const c = colour(i);
        let j = i + 1;
        while (j < line.length && colour(j) === c && !linkAt(j)) j++;
        const text = escapeHtml(line.slice(i, j));
        html += c && text.trim() ? `<span class="${c}">${text}</span>` : text;
        i = j;
      }
      return `<span class="line"${ids ? ` id="L${n + 1}"` : ''}>${html || ' '}</span>`;
    })
    .join('');
}

function codePanel(files) {
  const tabs = files.map((f, i) => `<button type="button" role="tab" data-file="${f.name}" aria-selected="${i === 0}">${f.name}</button>`).join('');
  const code = files.map((f, i) => `<pre data-file="${f.name}"${i ? ' hidden' : ''}><code>${codeLines(f, false, 'source/')}</code></pre>`).join('\n');
  return `<aside class="code-panel" aria-label="Source code for this chapter">
    <div class="code-head"><div class="code-tabs" role="tablist">${tabs}</div><a class="code-full" href="source/${files[0].name}.html">Full page</a></div>
    <div class="code-body">${code}</div>
  </aside>`;
}

// A citation pointing at a missing reference has to fail the build too.
function checkAnchor(file, id, page) {
  if (!readFileSync(join(docs, file), 'utf8').includes(`id="${id}"`)) throw new Error(`${page}: no anchor ${file}#${id}`);
}

// Maps a link in the Markdown to its page on the site.
function rewriteLink(href, page) {
  if (/^[a-z]+:/.test(href)) return href;
  if (href.startsWith('#')) {
    checkAnchor(page, href.slice(1), page);
    return href;
  }
  const [path, fragment] = href.split('#');
  if (path.endsWith('.md')) {
    if (fragment) checkAnchor(path, fragment, page);
    const target = path === 'README.md' ? 'index.html' : path.replace(/\.md$/, '.html');
    return fragment ? `${target}#${fragment}` : target;
  }
  const target = posix.normalize(posix.join('docs', href));
  if (sources.includes(target)) return `source/${basename(target)}.html`;
  throw new Error(`${page}: no page for link ${href}`);
}

function renderer(page, files, used, markdown) {
  let inLink = false;
  const use = (ref) => {
    if (!used.includes(ref.file)) used.push(ref.file);
    return ref;
  };
  const marked = new Marked({
    extensions: [
      {
        name: 'codeAnchor',
        level: 'block',
        start: (src) => src.indexOf('<!-- code:'),
        tokenizer(src) {
          const m = /^<!-- code: (.+?) -->[ \t]*(?:\n|$)/.exec(src);
          return m && { type: 'codeAnchor', raw: m[0], ref: use(resolveAnchor(m[1], files, page, markdown)) };
        },
        renderer: () => '',
      },
      {
        name: 'displayMath',
        level: 'block',
        start: (src) => src.indexOf('$$'),
        tokenizer(src) {
          const m = /^\$\$([\s\S]+?)\$\$[ \t]*(?:\n|$)/.exec(src);
          return m && { type: 'displayMath', raw: m[0], tex: m[1] };
        },
        renderer: (token) => `<div class="display-math">${math(token.tex, true)}</div>\n`,
      },
      {
        name: 'inlineMath',
        level: 'inline',
        start: (src) => src.indexOf('$'),
        tokenizer(src) {
          const m = /^\$([^$\n]+?)\$/.exec(src);
          return m && { type: 'inlineMath', raw: m[0], tex: m[1] };
        },
        renderer: (token) => math(token.tex, false),
      },
    ],
    renderer: {
      // A link to a source file shows that file, or the definition it names, beside the text.
      // A title names the definition when the link text is prose: [the test](x.rs "test_name").
      link({ href, title, tokens }) {
        const target = rewriteLink(href, page);
        inLink = true;
        const inner = this.parser.parseInline(tokens);
        inLink = false;
        const file = target.startsWith('source/') && wholeFile(target.slice('source/'.length, -'.html'.length));
        if (!file) return `<a href="${escapeHtml(target)}">${inner}</a>`;
        const name = title || (tokens.length === 1 && tokens[0].type === 'codespan' && tokens[0].text);
        const named = name && resolveName(name, new Set([file.file.path]), page, markdown);
        if (title && named?.file !== file.file) throw new Error(`${page}: ${title} isn't defined in ${file.file.name}`);
        const ref = use(named?.file === file.file ? named : file);
        return `<a class="code-ref" href="${refHref(ref, '')}" data-code="${refAttr([ref])}">${inner}</a>`;
      },
      codespan({ text }) {
        const ref = !inLink && resolveName(text, files, page, markdown);
        if (!ref) return false;
        use(ref);
        return `<a class="code-ref" href="${refHref(ref, '')}" data-code="${refAttr([ref])}"><code>${escapeHtml(text)}</code></a>`;
      },
      // A paragraph that is only a figure, ![description](figures/name.svg "caption"), puts the
      // drawing in the page itself, so it takes the reader's light or dark colours.
      paragraph({ tokens }) {
        const image = tokens.length === 1 && tokens[0].type === 'image' && tokens[0];
        if (!image || !image.href.startsWith('figures/')) return false;
        const svg = readFileSync(join(docs, image.href), 'utf8').trim().replace('<svg ', `<svg role="img" aria-label="${escapeHtml(image.text)}" `);
        const caption = image.title ? `<figcaption>${marked.parseInline(image.title)}</figcaption>` : '';
        return `<figure class="figure"><div class="figure-art">${svg}</div>${caption}</figure>\n`;
      },
      // Wide tables scroll on their own, rather than widening the whole page.
      table(token) {
        return `<div class="table-scroll">${tableHtml(this, token)}</div>\n`;
      },
    },
  });
  return marked;
}

function tableHtml(self, token) {
  const cell = (c, tag) => `<${tag}>${self.parser.parseInline(c.tokens)}</${tag}>`;
  const head = `<tr>${token.header.map((c) => cell(c, 'th')).join('')}</tr>`;
  const rows = token.rows.map((r) => `<tr>${r.map((c) => cell(c, 'td')).join('')}</tr>`).join('');
  return `<table><thead>${head}</thead><tbody>${rows}</tbody></table>`;
}

// Chapters before this one are the essentials, and each one after stands on its own.
const FIRST_OPTIONAL = '14-';

const chapters = readdirSync(docs)
  .filter((f) => /^\d+-.*\.md$/.test(f))
  .sort()
  .map((file) => ({
    file,
    href: file.replace(/\.md$/, '.html'),
    title: /^# (.+)$/m.exec(readFileSync(join(docs, file), 'utf8'))[1],
  }));

// What each kind of page may be reused under: the tutorial's text, or the solver's code.
const TEXT_LICENSE = 'Text by the Glassworks authors, under <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a>. The code it shows is under the MIT License.';
const CODE_LICENSE = 'Code by the Glassworks authors, under the MIT License.';

// The site's mark, as the top bar of the solver page draws it.
const MARK = '<svg class="mark" aria-hidden="true" focusable="false" viewBox="0 0 32 24"><rect x="12" y="2" width="3" height="20" style="fill:var(--glassworks)"/><rect x="17" y="2" width="3" height="20" style="fill:var(--muted)"/><rect x="22" y="2" width="3" height="20" style="fill:var(--glassworks)"/><rect x="27" y="2" width="3" height="20" style="fill:var(--muted)"/><path d="M1 12C2.5 5 4.5 5 6 12S9.5 19 11 12S13.5 7.5 15 12S17.5 14.5 19 12" fill="none" stroke-width="1.6" stroke-linecap="round" style="stroke:var(--ink)"/></svg>';

// The button that opens the chapter list over the page.
const MENU = '<svg aria-hidden="true" focusable="false" viewBox="0 0 16 16"><path d="M2 4h12M2 8h12M2 12h12" stroke="currentColor" stroke-width="1.6" stroke-linecap="round"/></svg>';

function layout({ title, body, current, prefix = '', panel = '', pager = '', license = TEXT_LICENSE }) {
  const links = chapters
    .map((c) => `${c.file.startsWith(FIRST_OPTIONAL) ? '<li class="part">Going further</li>\n        ' : ''}<li><a href="${prefix}${c.href}"${c.href === current ? ' aria-current="page"' : ''}>${escapeHtml(c.title)}</a></li>`)
    .join('\n      ');
  return `<!doctype html>
<html lang="en" data-site="glassworks">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>${escapeHtml(title)} | Glassworks</title>
<link rel="stylesheet" href="${prefix}../family.css">
<link rel="stylesheet" href="${prefix}../style.css">
<link rel="icon" href="${prefix}../favicon.ico" sizes="32x32">
<link rel="icon" href="${prefix}../favicon.svg" type="image/svg+xml">
<link rel="apple-touch-icon" href="${prefix}../apple-touch-icon.png">
<link rel="stylesheet" href="${prefix}../learn.css">${
    panel
      ? `
<script>try { const c = document.documentElement.classList; if (localStorage.getItem('code-panel') === 'off') c.add('code-off'); if (localStorage.getItem('code-known')) c.add('code-known'); } catch {}</script>
<script type="module" src="${prefix}../learn.js"></script>`
      : ''
  }
</head>
<body class="learn${panel ? ' with-code' : ''}">
<nav class="topnav">
  <div class="wrap">
    <div class="brand"><button class="chapters-open" type="button" popovertarget="chapters" aria-label="Chapters" title="Chapters">${MENU}</button><!--apex--><a class="home" href="${prefix}../index.html">${MARK}Glassworks</a></div>
    <div class="navright">
      <div class="links">
        <a href="${prefix}../index.html">Solver</a>
        <a href="${prefix}index.html" aria-current="${current === 'index.html' ? 'page' : 'true'}">Tutorial</a>
      </div>
      <button class="theme-toggle" hidden></button>
    </div>
  </div>
</nav>
<div class="learn-layout page">
  <nav class="chapters" id="chapters" popover aria-label="Tutorial chapters">
    <div class="chapters-head">Chapters<button type="button" popovertarget="chapters" popovertargetaction="hide" aria-label="Close">&times;</button></div>
    <ul>
      <li><a href="${prefix}index.html"${current === 'index.html' ? ' aria-current="page"' : ''}>Start here</a></li>
      ${links}
    </ul>
    <script>{ const nav = document.currentScript.parentElement; nav.addEventListener('toggle', (e) => e.newState === 'open' && nav.querySelector('[aria-current=page]')?.scrollIntoView({ block: 'center' })); }</script>
  </nav>
  <main>
${body}${pager}
    <footer class="learn-footer">${license}</footer>
  </main>${panel ? `\n  ${panel}` : ''}
</div>
</body>
</html>
`;
}

const write = (path, html) => check || writeFileSync(join(out, path), html);
if (!check) {
  rmSync(out, { recursive: true, force: true });
  mkdirSync(join(out, 'source'), { recursive: true });
}

// Renders a chapter, giving each block the code its following <!-- code: --> comments name.
function chapter(file, markdown) {
  const files = new Set([...markdown.matchAll(/\]\((\.\.\/[^)#]+)/g)].map((m) => posix.normalize(posix.join('docs', m[1]))));
  const used = [];
  const marked = renderer(file, files, used, markdown);
  const tokens = marked.lexer(markdown);
  let html = '';
  tokens.forEach((token, i) => {
    if (token.type === 'codeAnchor') return;
    let block = marked.parser([token]);
    const refs = [];
    for (let j = i + 1; tokens[j]?.type === 'codeAnchor' || (refs.length && tokens[j]?.type === 'space'); j++) if (tokens[j].ref) refs.push(tokens[j].ref);
    if (refs.length) {
      if (refs.some((r) => r.file !== refs[0].file)) throw new Error(`${file}: one passage can only point into one file`);
      block = block.replace(/^<(\w+)/, `<$1 data-code="${refAttr(refs)}"`);
    }
    html += block;
  });
  return { html, used };
}

// The pages in reading order, for the links at the foot of each one.
const order = [{ href: 'index.html', title: 'Start here' }, ...chapters];
function pager(href) {
  const at = order.findIndex((c) => c.href === href);
  const link = (c, rel, label) => c ? `<a href="${c.href}" rel="${rel}"><span>${label}</span>${escapeHtml(c.title)}</a>` : '<span></span>';
  return `\n    <nav class="pager" aria-label="Previous and next chapter">${link(order[at - 1], 'prev', 'Previous')}${link(order[at + 1], 'next', 'Next')}</nav>`;
}

// The index's chapter lists become a contents: a number, the title, and its gloss under it.
const contents = (html) =>
  html
    .replace(/<ol( start="\d+")?>\n(?=<li><a href="\d\d-)/g, '<ol class="contents"$1>\n')
    .replace(/<li><a href="((\d\d)-[^"]+\.html)">([\s\S]*?)<\/a>: ([\s\S]*?)<\/li>/g, (_, href, n, title, gloss) => `<li><span class="n" aria-hidden="true">${+n}</span><a href="${href}">${title}</a><span>${gloss}</span></li>`);

for (const file of ['README.md', ...chapters.map((c) => c.file)]) {
  const markdown = readFileSync(join(docs, file), 'utf8');
  const href = file === 'README.md' ? 'index.html' : file.replace(/\.md$/, '.html');
  const title = /^# (.+)$/m.exec(markdown)[1];
  let { html, used } = chapter(file, markdown);
  const panel = file === 'README.md' || !used.length ? '' : codePanel(used);
  const hint = panel ? '<p class="code-hint">This chapter\'s code is beside the text. Point at a name in <code>code</code> type, or a passage marked in the margin, to see the lines that do it.</p>\n' : '';
  // The closing "Next: ..." line is for GitHub; here the links at the foot do its job.
  html = file === 'README.md' ? contents(html) : html.replace(/<p>Next: [\s\S]*?<\/p>\n?$/, '');
  const body = `<article>${html.replace(/<\/h1>\n/, (h) => h + hint)}</article>`;
  write(href, layout({ title, body, current: href, panel, pager: pager(href) }));
}

for (const file of index) {
  const body = `<article class="source"><h1>${file.name}</h1><p>Solver source code, referenced by the tutorial.</p><pre><code>${codeLines(file, true, '')}</code></pre>${LIGHT_LINES}</article>`;
  write(join('source', `${file.name}.html`), layout({ title: file.name, body, prefix: '../', license: CODE_LICENSE }));
}

console.log(`${check ? 'Checked' : 'Wrote'} ${chapters.length + 1} tutorial pages and ${sources.length} source pages${check ? '' : ' to web/learn'}`);
