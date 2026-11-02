// Turns docs/*.md into static pages in web/learn, so docs/ stays the only copy of
// the tutorial. Equations become MathML here, at build time, so visitors download
// no math library or fonts. A broken equation or link fails the build.

import { mkdirSync, readFileSync, readdirSync, rmSync, writeFileSync } from 'node:fs';
import { basename, dirname, join, posix } from 'node:path';
import { fileURLToPath } from 'node:url';
import katex from 'katex';
import { Marked } from 'marked';

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
  ...['mesh', 'space', 'assemble', 'banded', 'eigen', 'problem', 'output', 'solve1d', 'sparse', 'sparse/symmetric', 'sparse/transposed', 'solve2d', 'solve2d/closed', 'solve2d/open', 'solve2d/driven', 'solve2d/leaky', 'solve2d/adjoint', 'solve2d/grating', 'solve2d/dipole', 'solve2d/validate'].map((n) => `crates/eigenmenno-core/src/${n}.rs`),
  ...['analytic', 'analytic2d', 'gratings', 'hz_design', 'cross_check'].map((n) => `crates/eigenmenno-core/tests/${n}.rs`),
];

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

function renderer(page) {
  const marked = new Marked({
    extensions: [
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
      link({ href, tokens }) {
        return `<a href="${escapeHtml(rewriteLink(href, page))}">${this.parser.parseInline(tokens)}</a>`;
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

function layout({ title, body, current, prefix = '' }) {
  const links = chapters
    .map((c) => `${c.file.startsWith(FIRST_OPTIONAL) ? '<li class="part">Going further</li>\n        ' : ''}<li><a href="${prefix}${c.href}"${c.href === current ? ' aria-current="page"' : ''}>${escapeHtml(c.title)}</a></li>`)
    .join('\n        ');
  return `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>${escapeHtml(title)} | Eigenmenno</title>
<link rel="stylesheet" href="${prefix}../style.css">
<link rel="stylesheet" href="${prefix}../learn.css">
</head>
<body class="learn">
<header class="learn-header">
  <a class="brand" href="${prefix}index.html">Eigenmenno tutorial</a>
  <a class="open-solver" href="${prefix}../index.html">Open the solver</a>
</header>
<div class="learn-layout">
  <nav aria-label="Tutorial chapters">
    <details open>
      <summary>Chapters</summary>
      <ul>
        <li><a href="${prefix}index.html"${current === 'index.html' ? ' aria-current="page"' : ''}>Start here</a></li>
        ${links}
      </ul>
    </details>
  </nav>
  <main>
${body}
  </main>
</div>
</body>
</html>
`;
}

rmSync(out, { recursive: true, force: true });
mkdirSync(join(out, 'source'), { recursive: true });

for (const file of ['README.md', ...chapters.map((c) => c.file)]) {
  const markdown = readFileSync(join(docs, file), 'utf8');
  const href = file === 'README.md' ? 'index.html' : file.replace(/\.md$/, '.html');
  const title = /^# (.+)$/m.exec(markdown)[1];
  const body = `<article>${renderer(file).parse(markdown)}</article>`;
  writeFileSync(join(out, href), layout({ title, body, current: href }));
}

for (const path of sources) {
  const name = basename(path);
  const lines = readFileSync(join(root, path), 'utf8').trimEnd().split('\n');
  const code = lines.map((l) => `<span class="line">${escapeHtml(l) || ' '}</span>`).join('');
  const body = `<article class="source"><h1>${name}</h1><p>Solver source code, referenced by the tutorial.</p><pre><code>${code}</code></pre></article>`;
  writeFileSync(join(out, 'source', `${name}.html`), layout({ title: name, body, prefix: '../' }));
}

console.log(`Wrote ${chapters.length + 1} tutorial pages and ${sources.length} source pages to web/learn`);
