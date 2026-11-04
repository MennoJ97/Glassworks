// Colours source code the way an editor does: keywords, control flow, types, functions,
// constants, numbers and variables, on top of the comments and strings that classify()
// in code-index.mjs already finds. It goes by the shape of each name, not by what the
// name resolves to, so `Some(x)` is a type (capitalised) and `f(x)` a function (called).

const KEYWORDS = {
  rs: 'as const crate dyn enum extern false fn impl in let mod move mut pub ref self Self static struct super trait true type unsafe use where',
  js: 'class const delete export extends false function import in instanceof let new null of static super this true typeof undefined var void',
};
const CONTROL = {
  rs: 'async await break continue else for if loop match return while yield',
  js: 'async await break case catch continue default do else finally for if return switch throw try while yield',
};
const PRIMITIVES = 'bool char str f32 f64 i8 i16 i32 i64 i128 isize u8 u16 u32 u64 u128 usize';
const words = (s) => new Set(s.split(' '));
const sets = {
  rs: { keyword: words(KEYWORDS.rs), control: words(CONTROL.rs), primitive: words(PRIMITIVES) },
  js: { keyword: words(KEYWORDS.js), control: words(CONTROL.js), primitive: new Set() },
};

const WORD = /[A-Za-z_$][\w$]*/y;
const NUMBER = /(?:0[xob][\da-fA-F_]+|\d[\d_]*(?:\.\d[\d_]*)?(?:[eE][+-]?\d+)?)(?:_?[iuf](?:8|16|32|64|128|size))?/y;
const LIFETIME = /'[A-Za-z_]\w*/y;

// The class of each character: '' for punctuation and plain text, otherwise kw (keyword),
// ctl (control flow), ty (type), fn (function, method or macro), co (constant), nu
// (number), va (variable or field), cm (comment) or st (string).
export function highlight(src, kind, lang) {
  const { keyword, control, primitive } = sets[lang];
  const cls = new Array(src.length).fill('');
  const mark = (from, to, c) => cls.fill(c, from, to);
  const nextCode = (j) => {
    while (j < src.length && !kind[j] && /\s/.test(src[j])) j++;
    return j;
  };
  let previous = '';
  for (let i = 0; i < src.length; ) {
    if (kind[i]) {
      mark(i, i + 1, kind[i] === 1 ? 'cm' : 'st');
      i++;
      continue;
    }
    const before = src[i - 1] ?? '';
    const match = (re) => ((re.lastIndex = i), re.exec(src)?.[0]);
    let token;
    if (/[A-Za-z_$]/.test(src[i]) && !/[\w$]/.test(before) && (token = match(WORD))) {
      const after = nextCode(i + token.length);
      const called = src[after] === '(' || (lang === 'rs' && src.startsWith('::<', after));
      const macro = lang === 'rs' && src[after] === '!' && src[after + 1] !== '=';
      const c = keyword.has(token) ? 'kw'
        : control.has(token) ? 'ctl'
        : primitive.has(token) ? 'ty'
        : /^[A-Z][A-Z\d_]+$/.test(token) ? 'co'
        : /^[A-Z]/.test(token) ? 'ty'
        : called || macro || previous === 'fn' || previous === 'function' ? 'fn'
        : 'va';
      mark(i, i + token.length + (macro ? 1 : 0), c);
      previous = token;
      i += token.length;
    } else if (/\d/.test(src[i]) && !/[\w$]/.test(before) && (token = match(NUMBER))) {
      mark(i, i + token.length, 'nu');
      i += token.length;
    } else if (lang === 'rs' && src[i] === "'" && (token = match(LIFETIME))) {
      mark(i, i + token.length, 'kw');
      i += token.length;
    } else {
      if (!/\s/.test(src[i])) previous = '';
      i++;
    }
  }
  return cls;
}
