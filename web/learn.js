// The code panel beside a tutorial chapter. Pointing at a name or a marked passage shows
// the lines it means, and so does reading: the passage crossing the upper third of the
// screen picks the code. A passage's `data-code` is "file.rs:12-30 file.rs:40-41", or
// just "file.rs" for the whole file.

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
  const parts = el.dataset.code.split(' ');
  const { lines } = selectFile(parts[0].split(':')[0]);
  for (const line of lit) line.classList.remove('hl');
  lit = parts.filter((p) => p.includes(':')).flatMap((p) => {
    const [a, b] = p.split(':')[1].split('-').map(Number);
    return [...lines].slice(a - 1, b);
  });
  for (const line of lit) line.classList.add('hl');
  if (!lit.length) return body.scrollTo({ top: 0 });
  full.href += `#L${[...lines].indexOf(lit[0]) + 1}`;
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
      if (!showing() || e.metaKey || e.ctrlKey || e.shiftKey || e.button) return;
      e.preventDefault();
      show(el);
    });
  }
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

// The chapter list moves above the text when the panel needs its column.
const nav = document.querySelector('.learn nav details');
const narrow = matchMedia('(max-width: 1399px)');
const fold = () => (nav.open = !(showing() && narrow.matches));
wide.addEventListener('change', fold);
narrow.addEventListener('change', fold);

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
  fold();
  read = null;
  follow();
});
label();
document.querySelector('.learn-header .open-solver').before(toggle);

fold();
follow();
