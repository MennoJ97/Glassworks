// The code panel beside a tutorial chapter. Pointing at a name or a marked passage shows
// the lines it means, and so does reading: the passage crossing the upper third of the
// screen picks the code. A passage's `data-code` is "file.rs:12-30 file.rs:40-41", or
// just "file.rs" for the whole file. On a narrow screen the lines open under the text.

const root = document.documentElement;
const article = document.querySelector('.learn article');
const panel = document.querySelector('.code-panel');
const body = panel.querySelector('.code-body');
const full = panel.querySelector('.code-full');
const tabs = [...panel.querySelectorAll('[role=tab]')];
const files = new Map([...panel.querySelectorAll('pre')].map((pre) => [pre.dataset.file, { pre, lines: pre.querySelectorAll('.line') }]));
const wide = matchMedia('(min-width: 1100px)');
const still = matchMedia('(prefers-reduced-motion: reduce)');
const showing = () => wide.matches && !root.classList.contains('code-off');

function parse(spec) {
  const parts = spec.split(' ');
  return { file: parts[0].split(':')[0], ranges: parts.filter((p) => p.includes(':')).map((p) => p.split(':')[1].split('-').map(Number)) };
}
const fragment = ([a, b]) => (a === b ? `#L${a}` : `#L${a}-${b}`);

function selectFile(name) {
  for (const tab of tabs) tab.setAttribute('aria-selected', tab.dataset.file === name);
  for (const [file, { pre }] of files) pre.hidden = file !== name;
  full.href = `source/${name}.html`;
  return files.get(name);
}

let active = null;
let lit = [];

function show(el) {
  if (el === active) return;
  active?.classList.remove('code-active');
  el.classList.add('code-active');
  active = el;
  const { file, ranges } = parse(el.dataset.code);
  const { lines } = selectFile(file);
  for (const line of lit) line.classList.remove('hl');
  lit = ranges.flatMap(([a, b]) => [...lines].slice(a - 1, b));
  for (const line of lit) line.classList.add('hl');
  if (!lit.length) return body.scrollTo({ top: 0 });
  full.href += fragment(ranges[0]);
  // Leave the view alone if the lines are already in it; otherwise centre them, or show
  // the top of a range taller than the panel.
  const top = lit[0].offsetTop;
  const height = lit.at(-1).offsetTop + lit.at(-1).offsetHeight - top;
  const view = body.clientHeight;
  if (top >= body.scrollTop && top + height <= body.scrollTop + view) return;
  const target = height < view * 0.8 ? top - (view - height) / 2 : top - 16;
  body.scrollTo({ top: Math.max(0, target), behavior: still.matches ? 'auto' : 'smooth' });
}

// A name shows its code at once. A passage waits a moment, so a pointer crossing the
// text on its way to the panel doesn't drag the code along with it.
let pending = 0;
for (const el of article.querySelectorAll('[data-code]')) {
  const inline = el.matches('.code-ref');
  el.addEventListener('mouseenter', () => {
    clearTimeout(pending);
    if (showing()) pending = setTimeout(() => show(el), inline ? 0 : 150);
  });
  el.addEventListener('mouseleave', () => clearTimeout(pending));
  if (inline) {
    el.addEventListener('focus', () => showing() && show(el));
    el.addEventListener('click', (e) => {
      if (e.metaKey || e.ctrlKey || e.shiftKey || e.button) return;
      if (showing()) {
        e.preventDefault();
        show(el);
      } else if (!wide.matches && el.dataset.code.includes(':')) {
        e.preventDefault();
        peek(el.closest('li') ?? topLevel(el), el.dataset.code);
      }
    });
  }
}

// Narrow screens: no room beside the text, so the lines open under it. A marked passage
// gets a button for its own lines, and a name opens its lines under the paragraph or
// list item it sits in. A name for a whole file stays a link to that file's page.
const opened = new Map();
const buttons = new Map();
const topLevel = (el) => {
  while (el.parentElement !== article) el = el.parentElement;
  return el;
};

function snippet(spec) {
  const { file, ranges } = parse(spec);
  const { lines } = files.get(file);
  const box = document.createElement('div');
  box.className = 'code-snippet';
  const head = document.createElement('div');
  head.className = 'code-snippet-head';
  const name = document.createElement('span');
  name.textContent = `${file} · ${ranges.map(([a, b]) => (a === b ? a : `${a}–${b}`)).join(', ')}`;
  const link = document.createElement('a');
  link.href = `source/${file}.html${fragment(ranges[0])}`;
  link.textContent = 'Full page';
  head.append(name, link);
  const code = document.createElement('code');
  ranges.forEach(([a, b], i) => {
    if (i) code.append(Object.assign(document.createElement('span'), { className: 'line gap', textContent: '⋯' }));
    for (let n = a; n <= b; n++) {
      const line = lines[n - 1].cloneNode(true);
      line.classList.remove('hl');
      line.dataset.n = n;
      code.append(line);
    }
  });
  // The indent the lines share is only width the screen doesn't have.
  const shown = [...code.children].filter((line) => line.textContent.trim() && !line.matches('.gap'));
  const indent = Math.min(...shown.map((line) => line.textContent.search(/\S/)));
  for (const line of shown) {
    let left = indent;
    for (const walk = document.createTreeWalker(line, NodeFilter.SHOW_TEXT); left && walk.nextNode(); ) {
      const cut = Math.min(left, walk.currentNode.length);
      walk.currentNode.deleteData(0, cut);
      left -= cut;
    }
  }
  const pre = document.createElement('pre');
  pre.append(code);
  box.append(head, pre);
  return box;
}

// Opens `spec` under `host`, or closes it if that's what's open there.
function peek(host, spec) {
  const was = opened.get(host);
  was?.box.remove();
  opened.delete(host);
  if (was?.spec !== spec) {
    const box = snippet(spec);
    if (host.matches('li')) host.append(box);
    else (buttons.get(host) ?? host).after(box);
    opened.set(host, { box, spec });
    // A paragraph can run past the bottom of the screen, taking its code with it.
    const top = box.getBoundingClientRect().top;
    if (top > innerHeight - 96) scrollBy({ top: top - innerHeight / 2, behavior: still.matches ? 'auto' : 'smooth' });
  }
  const button = buttons.get(host);
  if (button) {
    const open = opened.get(host)?.spec === host.dataset.code;
    button.setAttribute('aria-expanded', open);
    button.textContent = open ? 'Hide the code' : 'Show the code';
  }
}

for (const el of article.querySelectorAll('[data-code]:not(a)')) {
  const button = document.createElement('button');
  button.type = 'button';
  button.className = 'code-peek';
  button.setAttribute('aria-expanded', false);
  button.textContent = 'Show the code';
  button.addEventListener('click', () => peek(el, el.dataset.code));
  el.after(button);
  buttons.set(el, button);
}

for (const tab of tabs) tab.addEventListener('click', () => selectFile(tab.dataset.file));

// Reading: each top-level passage with code, its own or that of the first name in it.
const passages = [...article.children].flatMap((el) => {
  const target = el.matches('[data-code]') ? el : el.querySelector('.code-ref');
  return target ? [{ el, target }] : [];
});
let read = null;
function follow() {
  if (!showing()) return;
  const line = innerHeight / 3;
  const current = passages.findLast((p) => p.el.getBoundingClientRect().top <= line);
  if (current && current !== read) show(current.target);
  read = current;
}
addEventListener('scroll', follow, { passive: true });

const toggle = document.createElement('button');
toggle.type = 'button';
toggle.className = 'code-toggle';
const label = () => (toggle.textContent = root.classList.contains('code-off') ? 'Show code' : 'Hide code');
toggle.addEventListener('click', () => {
  root.classList.toggle('code-off');
  try {
    localStorage.setItem('code-panel', root.classList.contains('code-off') ? 'off' : 'on');
  } catch {}
  label();
  read = null;
  follow();
});
label();
document.querySelector('.topnav .links').before(toggle);

follow();
