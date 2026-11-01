// Canvas drawing: the field drawn over its structure, and the line charts.
// A 2D field view would replace FieldPlot with a WebGL triangle renderer that reads
// the same channels; nothing outside this file draws anything.

const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
export const seriesColors = () => ['--series-1', '--series-2', '--series-3'].map(css);

function fit(canvas) {
  const ratio = window.devicePixelRatio || 1;
  const { width, height } = canvas.getBoundingClientRect();
  if (canvas.width !== Math.round(width * ratio) || canvas.height !== Math.round(height * ratio)) {
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  ctx.clearRect(0, 0, width, height);
  ctx.font = '12px ui-sans-serif, system-ui, sans-serif';
  return { ctx, width, height };
}

function ticks(min, max, count) {
  const raw = (max - min) / count;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((m) => m >= raw);
  const out = [];
  for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-6; v += step) out.push(Math.abs(v) < step * 1e-6 ? 0 : v);
  return out;
}

const label = (v) => String(+v.toPrecision(3));

function drawYAxis(ctx, area, scale, values, format = label) {
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (const v of values) {
    const y = scale(v);
    ctx.strokeStyle = css('--grid');
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(area.left, y);
    ctx.lineTo(area.right, y);
    ctx.stroke();
    ctx.fillStyle = css('--muted');
    ctx.fillText(format(v), area.left - 8, y);
  }
}

function drawXAxis(ctx, area, scale, values, title, format = label) {
  ctx.fillStyle = css('--muted');
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (const v of values) ctx.fillText(format(v), scale(v), area.bottom + 7);
  ctx.textAlign = 'right';
  ctx.fillText(title, area.right, area.bottom + 22);
}

// The structure and its field share one x-axis, and the boundaries are draggable.
export class FieldPlot {
  constructor(canvas, handlers) {
    this.canvas = canvas;
    this.handlers = handlers;
    this.layers = [];
    this.field = null;
    this.visible = new Set();
    this.selected = -1;
    this.source = null;
    // A locked plot shows a design the optimizer owns: no handles, no editing.
    this.locked = false;
    this.phase = 0;
    this.hover = null;
    this.drag = null;
    canvas.addEventListener('pointerdown', (e) => this.onDown(e));
    canvas.addEventListener('pointermove', (e) => this.onMove(e));
    canvas.addEventListener('pointerup', () => { this.drag = null; });
    canvas.addEventListener('pointerleave', () => { this.hover = null; this.draw(); });
    canvas.addEventListener('dblclick', (e) => {
      if (this.locked) return;
      const hit = this.layerAt(this.toX(e));
      if (hit) handlers.onSplit(hit.index, hit.offset);
    });
  }

  area(width, height) {
    return { left: 44, right: width - 14, top: 8, bottom: height - 40 };
  }

  length() {
    return this.layers.reduce((s, l) => s + l.thickness, 0) || 1;
  }

  toX(e) {
    const rect = this.canvas.getBoundingClientRect();
    const a = this.area(rect.width, rect.height);
    return ((e.clientX - rect.left - a.left) / (a.right - a.left)) * this.length();
  }

  layerAt(x) {
    let start = 0;
    for (const [index, l] of this.layers.entries()) {
      if (x >= start && x <= start + l.thickness) return { index, offset: x - start };
      start += l.thickness;
    }
    return null;
  }

  // What a pointer at position `x` would grab: an interface or the source marker.
  grab(x) {
    if (this.locked) return null;
    const rect = this.canvas.getBoundingClientRect();
    const a = this.area(rect.width, rect.height);
    const reach = (10 / (a.right - a.left)) * this.length();
    if (this.source !== null && Math.abs(x - this.source) < reach) return { type: 'source' };
    let edge = 0;
    for (let i = 0; i < this.layers.length - 1; i++) {
      edge += this.layers[i].thickness;
      if (Math.abs(x - edge) < reach) return { type: 'interface', index: i };
    }
    return null;
  }

  onDown(e) {
    const x = this.toX(e);
    this.drag = this.grab(x);
    if (this.drag) {
      this.canvas.setPointerCapture(e.pointerId);
    } else {
      const hit = this.layerAt(x);
      if (hit && !this.locked) this.handlers.onSelect(hit.index);
    }
  }

  onMove(e) {
    const x = this.toX(e);
    if (this.drag?.type === 'interface') this.handlers.onInterface(this.drag.index, x);
    else if (this.drag?.type === 'source') this.handlers.onSource(x);
    this.hover = this.drag ?? this.grab(x);
    this.canvas.style.cursor = this.hover ? 'ew-resize' : 'default';
    this.draw();
  }

  draw() {
    const { ctx, width, height } = fit(this.canvas);
    const a = this.area(width, height);
    const L = this.length();
    const sx = (x) => a.left + (x / L) * (a.right - a.left);

    const peak = this.field ? this.peak() * 1.12 || 1 : 1;
    const sy = (v) => a.top + ((peak - v) / (2 * peak)) * (a.bottom - a.top);
    this.drawLayers(ctx, a, sx);
    drawYAxis(ctx, a, sy, ticks(-peak, peak, 4));
    drawXAxis(ctx, a, sx, ticks(0, L, Math.max(3, Math.floor(width / 110))), 'Position (µm)');
    if (this.field) this.drawField(ctx, sx, sy);
    this.drawHandles(ctx, a, sx);
  }

  drawLayers(ctx, a, sx) {
    const glass = css('--glass');
    const ink = css('--ink');
    let start = 0;
    for (const [i, l] of this.layers.entries()) {
      const [x0, x1] = [sx(start), sx(start + l.thickness)];
      // Deeper tint means higher index, so the stack reads like a cross-section.
      ctx.globalAlpha = Math.min(0.55, 0.3 * Math.max(0, l.n - 1) / 1.2);
      ctx.fillStyle = glass;
      ctx.fillRect(x0, a.top, x1 - x0, a.bottom - a.top);
      ctx.globalAlpha = 1;
      if (l.k > 0) this.hatch(ctx, x0, x1, a, ink, Math.max(5, 14 - l.k * 2));
      if (i === this.selected) {
        ctx.strokeStyle = ink;
        ctx.lineWidth = 1.5;
        ctx.strokeRect(x0 + 0.75, a.top + 0.75, x1 - x0 - 1.5, a.bottom - a.top - 1.5);
      }
      if (x1 - x0 > 34) {
        ctx.fillStyle = css('--muted');
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.fillText(l.k > 0 ? `${label(l.n)}+${label(l.k)}i` : label(l.n), (x0 + x1) / 2, a.top + 6);
      }
      start += l.thickness;
    }
  }

  // Lossy layers get a hatch, so absorption does not rely on colour alone.
  hatch(ctx, x0, x1, a, color, gap) {
    ctx.save();
    ctx.beginPath();
    ctx.rect(x0, a.top, x1 - x0, a.bottom - a.top);
    ctx.clip();
    ctx.strokeStyle = color;
    ctx.globalAlpha = 0.28;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let x = x0 - (a.bottom - a.top); x < x1; x += gap) {
      ctx.moveTo(x, a.bottom);
      ctx.lineTo(x + (a.bottom - a.top), a.top);
    }
    ctx.stroke();
    ctx.restore();
  }

  peak() {
    let peak = 0;
    for (const name of this.visible) {
      const { re, im } = this.field.components[name];
      for (let i = 0; i < re.length; i++) peak = Math.max(peak, Math.hypot(re[i], im[i]));
    }
    return peak;
  }

  drawField(ctx, sx, sy) {
    const colors = seriesColors();
    const [c, s] = [Math.cos(this.phase), Math.sin(this.phase)];
    const { x, components, names } = this.field;
    names.forEach((name, n) => {
      if (!this.visible.has(name)) return;
      const { re, im } = components[name];
      ctx.strokeStyle = colors[n];
      // The envelope is the reach of the oscillation, drawn thin behind it.
      ctx.globalAlpha = 0.35;
      ctx.lineWidth = 1;
      for (const sign of [1, -1]) {
        ctx.beginPath();
        for (let i = 0; i < x.length; i++) ctx.lineTo(sx(x[i]), sy(sign * Math.hypot(re[i], im[i])));
        ctx.stroke();
      }
      ctx.globalAlpha = 1;
      ctx.lineWidth = 2;
      ctx.lineJoin = 'round';
      ctx.beginPath();
      // Re[(re + i im) e^{-i phase}], the field at this instant.
      for (let i = 0; i < x.length; i++) ctx.lineTo(sx(x[i]), sy(re[i] * c + im[i] * s));
      ctx.stroke();
    });
  }

  drawHandles(ctx, a, sx) {
    if (this.locked) return;
    const ink = css('--ink');
    let edge = 0;
    for (let i = 0; i < this.layers.length - 1; i++) {
      edge += this.layers[i].thickness;
      const active = this.hover?.type === 'interface' && this.hover.index === i;
      ctx.strokeStyle = ink;
      ctx.globalAlpha = active ? 1 : 0.45;
      ctx.lineWidth = active ? 2 : 1;
      ctx.beginPath();
      ctx.moveTo(sx(edge), a.top);
      ctx.lineTo(sx(edge), a.bottom);
      ctx.stroke();
      ctx.globalAlpha = 1;
    }
    if (this.source !== null) {
      ctx.strokeStyle = ink;
      ctx.lineWidth = 2;
      ctx.setLineDash([5, 4]);
      ctx.beginPath();
      ctx.moveTo(sx(this.source), a.top);
      ctx.lineTo(sx(this.source), a.bottom);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = ink;
      ctx.textAlign = 'center';
      ctx.fillText('Source', sx(this.source), a.bottom - 16);
    }
  }
}

// A scalar field over a triangle mesh, with shape outlines on top. The mesh is
// rasterized once per size, so animating only recombines two stored images.
export class FieldMap {
  constructor(canvas) {
    this.canvas = canvas;
    this.mesh = null;
    this.shapes = [];
    this.marks = {};
    this.phase = 0;
    this.raster = null;
    this.image = document.createElement('canvas');
  }

  // `marks` may hold `bloch`, the phase per period across x and y; `open` for edges that
  // absorb; `source`, a point; `incident`, a direction of travel in degrees; and `monitor`,
  // a line that power is measured across; and `density`, a rectangle of pixels with `rho`.
  setMesh(mesh, shapes, marks) {
    this.mesh = mesh;
    this.shapes = shapes;
    this.marks = marks;
    this.raster = null;
    this.field = null;
    this.values = null;
  }

  setField(re, im) {
    this.field = { re, im };
    this.values = null;
  }

  clear() {
    this.mesh = null;
    this.draw();
  }

  // The domain drawn at one scale in both directions, inside the axes.
  box(width, height) {
    const [w, h] = [this.mesh.width, this.mesh.height];
    const scale = Math.min((width - 58) / w, (height - 48) / h);
    return { left: 44, right: 44 + w * scale, top: 8, bottom: 8 + h * scale, scale };
  }

  // For every pixel, the triangle holding it and its weights on the three corners.
  rasterize(a) {
    const { x, y, triangles } = this.mesh;
    const [w, h] = [Math.ceil(a.right - a.left), Math.ceil(a.bottom - a.top)];
    const owner = new Int32Array(w * h).fill(-1);
    const weights = new Float32Array(2 * w * h);
    const px = Array.from(x, (v) => v * a.scale);
    const py = Array.from(y, (v) => (this.mesh.height - v) * a.scale);
    for (let t = 0; t < triangles.length; t += 3) {
      const [i, j, k] = [triangles[t], triangles[t + 1], triangles[t + 2]];
      const det = (py[j] - py[k]) * (px[i] - px[k]) + (px[k] - px[j]) * (py[i] - py[k]);
      if (det === 0) continue;
      const x0 = Math.max(0, Math.floor(Math.min(px[i], px[j], px[k])));
      const x1 = Math.min(w - 1, Math.ceil(Math.max(px[i], px[j], px[k])));
      const y0 = Math.max(0, Math.floor(Math.min(py[i], py[j], py[k])));
      const y1 = Math.min(h - 1, Math.ceil(Math.max(py[i], py[j], py[k])));
      for (let row = y0; row <= y1; row++) {
        for (let col = x0; col <= x1; col++) {
          const [cx, cy] = [col + 0.5 - px[k], row + 0.5 - py[k]];
          const a0 = ((py[j] - py[k]) * cx + (px[k] - px[j]) * cy) / det;
          const a1 = ((py[k] - py[i]) * cx + (px[i] - px[k]) * cy) / det;
          if (a0 >= -1e-6 && a1 >= -1e-6 && a0 + a1 <= 1 + 1e-6) {
            const at = row * w + col;
            owner[at] = t;
            weights[2 * at] = a0;
            weights[2 * at + 1] = a1;
          }
        }
      }
    }
    this.raster = { w, h, owner, weights, width: a.right - a.left };
    this.image.width = w;
    this.image.height = h;
    this.values = null;
  }

  sample() {
    const { w, h, owner, weights } = this.raster;
    const { triangles } = this.mesh;
    const [re, im] = [new Float32Array(w * h), new Float32Array(w * h)];
    const { re: fr, im: fi } = this.field;
    let peak = 0;
    for (let at = 0; at < w * h; at++) {
      const t = owner[at];
      if (t < 0) continue;
      const [a, b] = [weights[2 * at], weights[2 * at + 1]];
      const c = 1 - a - b;
      const [i, j, k] = [triangles[t], triangles[t + 1], triangles[t + 2]];
      re[at] = a * fr[i] + b * fr[j] + c * fr[k];
      im[at] = a * fi[i] + b * fi[j] + c * fi[k];
      peak = Math.max(peak, Math.hypot(re[at], im[at]));
    }
    // A line source's field diverges at the point, which would wash out everything else.
    if (this.marks.source) {
      const sizes = Float32Array.from(re, (v, at) => Math.hypot(v, im[at])).sort();
      peak = sizes[Math.floor(0.99 * (sizes.length - 1))];
    }
    this.values = { re, im, peak: peak || 1 };
    const rgb = (name) => [1, 3, 5].map((i) => parseInt(css(name).slice(i, i + 2), 16));
    const [low, mid, high] = [rgb('--series-2'), rgb('--surface'), rgb('--series-1')];
    this.colors = new Uint8ClampedArray(3 * 511);
    for (let n = -255; n <= 255; n++) {
      const [end, f] = n < 0 ? [low, -n / 255] : [high, n / 255];
      for (let c = 0; c < 3; c++) this.colors[3 * (n + 255) + c] = mid[c] + f * (end[c] - mid[c]);
    }
  }

  draw() {
    const { ctx, width, height } = fit(this.canvas);
    if (!this.mesh || width < 80) return;
    const a = this.box(width, height);
    const sx = (x) => a.left + x * a.scale;
    const sy = (y) => a.top + (this.mesh.height - y) * a.scale;
    drawYAxis(ctx, a, sy, ticks(0, this.mesh.height, 4));
    drawXAxis(ctx, a, sx, ticks(0, this.mesh.width, Math.max(3, Math.floor((a.right - a.left) / 110))), 'Position (µm)');
    if (this.raster?.width !== a.right - a.left) this.rasterize(a);
    if (this.field) {
      if (!this.values) this.sample();
      this.paint();
      ctx.drawImage(this.image, a.left, a.top, a.right - a.left, a.bottom - a.top);
    }
    ctx.strokeStyle = css('--ink');
    ctx.lineWidth = 1.5;
    // Shapes can run on past the plot, into the absorbing layer.
    ctx.save();
    ctx.beginPath();
    ctx.rect(a.left, a.top, a.right - a.left, a.bottom - a.top);
    ctx.clip();
    if (this.marks.density) {
      const { x0, y0, width: w, height: h, nx, ny, rho } = this.marks.density;
      ctx.fillStyle = css('--ink');
      rho.forEach((v, i) => {
        if (v < 0.02) return;
        ctx.globalAlpha = 0.35 * v;
        // Whole device pixels, because overlapping edges would show as a grid.
        const [px, py] = [i % nx, Math.floor(i / nx)];
        const [left, right] = [Math.round(sx(x0 + (px * w) / nx)), Math.round(sx(x0 + ((px + 1) * w) / nx))];
        const [top, bottom] = [Math.round(sy(y0 + ((py + 1) * h) / ny)), Math.round(sy(y0 + (py * h) / ny))];
        ctx.fillRect(left, top, right - left, bottom - top);
      });
      ctx.globalAlpha = 1;
    }
    for (const s of this.shapes) {
      ctx.beginPath();
      if (s.type === 'circle') ctx.arc(sx(s.x), sy(s.y), (s.w / 2) * a.scale, 0, 2 * Math.PI);
      else ctx.rect(sx(s.x - s.w / 2), sy(s.y + s.h / 2), s.w * a.scale, s.h * a.scale);
      ctx.stroke();
    }
    ctx.restore();
    const { bloch, open, source, incident, monitor } = this.marks;
    if (bloch) this.drawBloch(ctx, a);
    else if (open) this.drawOpen(ctx, a);
    else ctx.strokeRect(a.left, a.top, a.right - a.left, a.bottom - a.top);
    if (source) this.drawSource(ctx, sx(source[0]), sy(source[1]));
    if (monitor) {
      ctx.strokeStyle = css('--ink');
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      ctx.moveTo(sx(monitor.x0), sy(monitor.y0));
      ctx.lineTo(sx(monitor.x1), sy(monitor.y1));
      ctx.stroke();
      this.drawLabel(ctx, 'Monitor', sx((monitor.x0 + monitor.x1) / 2), sy((monitor.y0 + monitor.y1) / 2) - 15);
    }
    if (incident !== undefined) {
      const angle = (incident * Math.PI) / 180;
      const size = 0.16 * Math.min(a.right - a.left, a.bottom - a.top);
      // Beside the corner the wave enters from, clear of shapes in the middle.
      const [cx, cy] = [a.left + 0.5 * size + 12 + (Math.cos(angle) < 0 ? a.right - a.left - size - 24 : 0), a.top + 0.5 * size + 12];
      this.drawArrow(ctx, cx, cy, Math.cos(angle), -Math.sin(angle), size, 'k');
    }
  }

  // Waves leave through these edges, into an absorbing layer outside the plot.
  drawOpen(ctx, a) {
    ctx.strokeStyle = css('--muted');
    ctx.lineWidth = 1.5;
    ctx.setLineDash([2, 4]);
    ctx.strokeRect(a.left, a.top, a.right - a.left, a.bottom - a.top);
    ctx.setLineDash([]);
  }

  drawSource(ctx, x, y) {
    ctx.fillStyle = css('--ink');
    ctx.strokeStyle = css('--surface');
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(x, y, 5, 0, 2 * Math.PI);
    ctx.fill();
    ctx.stroke();
    this.drawLabel(ctx, 'Source', x, y - 15);
  }

  drawLabel(ctx, text, x, y) {
    ctx.font = '600 13px ui-sans-serif, system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.lineWidth = 4;
    ctx.strokeStyle = css('--surface');
    ctx.strokeText(text, x, y);
    ctx.fillStyle = css('--ink');
    ctx.fillText(text, x, y);
  }

  // An arrow centred on (cx, cy), with its direction in screen coordinates.
  drawArrow(ctx, cx, cy, ux, uy, length, name) {
    const [tipX, tipY] = [cx + (ux * length) / 2, cy + (uy * length) / 2];
    const shape = () => {
      ctx.beginPath();
      ctx.moveTo(cx - (ux * length) / 2, cy - (uy * length) / 2);
      ctx.lineTo(tipX, tipY);
      ctx.moveTo(tipX - 10 * ux - 6 * uy, tipY - 10 * uy + 6 * ux);
      ctx.lineTo(tipX, tipY);
      ctx.lineTo(tipX - 10 * ux + 6 * uy, tipY - 10 * uy - 6 * ux);
      ctx.stroke();
    };
    ctx.lineCap = 'round';
    ctx.lineWidth = 6;
    ctx.strokeStyle = css('--surface');
    shape();
    ctx.lineWidth = 2.5;
    ctx.strokeStyle = css('--ink');
    shape();
    ctx.lineCap = 'butt';
    this.drawLabel(ctx, name, tipX + ux * 12, tipY + uy * 12);
  }

  // Matched chevrons mark each pair of tied edges, and the arrow shows the wavevector.
  drawBloch(ctx, a) {
    const ink = css('--ink');
    ctx.strokeStyle = ink;
    ctx.fillStyle = ink;
    ctx.lineWidth = 2;
    ctx.setLineDash([7, 5]);
    ctx.strokeRect(a.left, a.top, a.right - a.left, a.bottom - a.top);
    ctx.setLineDash([]);
    const [cx, cy] = [(a.left + a.right) / 2, (a.top + a.bottom) / 2];
    const chevrons = (x, y, dx, dy, count) => {
      for (let n = 0; n < count; n++) {
        const [px, py] = [x + dx * 7 * n, y + dy * 7 * n];
        ctx.beginPath();
        ctx.moveTo(px - dx * 6 - dy * 6, py - dy * 6 - dx * 6);
        ctx.lineTo(px, py);
        ctx.lineTo(px - dx * 6 + dy * 6, py - dy * 6 + dx * 6);
        ctx.stroke();
      }
    };
    // A quarter of the way along, clear of the axis labels at the middle of each side.
    const [markX, markY] = [(3 * a.left + a.right) / 4, (3 * a.top + a.bottom) / 4];
    chevrons(a.left, markY, 1, 0, 1);
    chevrons(a.right, markY, 1, 0, 1);
    chevrons(markX, a.bottom, 0, -1, 2);
    chevrons(markX, a.top, 0, -1, 2);

    const { bloch } = this.marks;
    const [kx, ky] = [bloch[0] / this.mesh.width, bloch[1] / this.mesh.height];
    const size = Math.hypot(kx, ky);
    if (size < 1e-9) {
      this.drawLabel(ctx, 'k = 0', cx, cy);
      return;
    }
    // Screen y points down. The longest k on the zone edge gets the longest arrow.
    const longest = Math.hypot(Math.PI / this.mesh.width, Math.PI / this.mesh.height);
    const length = (0.14 + 0.2 * (size / longest)) * Math.min(a.right - a.left, a.bottom - a.top);
    this.drawArrow(ctx, cx, cy, kx / size, -ky / size, length, 'k');
  }

  paint() {
    const { w, h, owner } = this.raster;
    const { re, im, peak } = this.values;
    const [c, s] = [Math.cos(this.phase), Math.sin(this.phase)];
    const image = new ImageData(w, h);
    const scale = 255 / peak;
    for (let at = 0; at < w * h; at++) {
      if (owner[at] < 0) continue;
      const n = 3 * (Math.max(-255, Math.min(255, Math.round((re[at] * c + im[at] * s) * scale))) + 255);
      image.data[4 * at] = this.colors[n];
      image.data[4 * at + 1] = this.colors[n + 1];
      image.data[4 * at + 2] = this.colors[n + 2];
      image.data[4 * at + 3] = 255;
    }
    this.image.getContext('2d').putImageData(image, 0, 0);
  }
}

// A small picture of a design region's densities, with `margin` pixels around them
// showing the `shapes` that lead into it. Strong colour is the high-index material.
export function drawDensity(canvas, { x0, y0, width, height, nx, ny }, shapes, rho, margin) {
  const [w, h] = [nx + 2 * margin, ny + 2 * margin];
  Object.assign(canvas, { width: w, height: h });
  const rgb = (name) => [1, 3, 5].map((i) => parseInt(css(name).slice(i, i + 2), 16));
  const [from, to] = [rgb('--surface'), rgb('--glass')];
  const image = new ImageData(w, h);
  for (let py = 0; py < h; py++) {
    for (let px = 0; px < w; px++) {
      const [ix, iy] = [px - margin, py - margin];
      const [x, y] = [x0 + ((ix + 0.5) * width) / nx, y0 + ((iy + 0.5) * height) / ny];
      const inside = ix >= 0 && iy >= 0 && ix < nx && iy < ny;
      const solid = inside ? rho[ix + iy * nx] : +shapes.some((s) => Math.abs(x - s.x) < s.w / 2 && Math.abs(y - s.y) < s.h / 2);
      const mix = 0.1 + 0.9 * solid;
      // The plots' y runs upward, and an image's downward.
      const at = 4 * (px + (h - 1 - py) * w);
      for (let c = 0; c < 3; c++) image.data[at + c] = from[c] + mix * (to[c] - from[c]);
      image.data[at + 3] = 255;
    }
  }
  canvas.getContext('2d').putImageData(image, 0, 0);
}

// Every pixel's density at every iteration, with iterations on the same x positions
// as a LineChart above it. Columns are painted once, so long runs stay cheap to draw.
export class EvolutionMap {
  constructor(canvas, onPick) {
    this.canvas = canvas;
    this.onPick = onPick;
    this.marker = null;
    this.total = 0;
    this.drawn = 0;
    const pick = (e) => {
      if (!this.drawn) return;
      const rect = canvas.getBoundingClientRect();
      const a = this.area(rect.width, rect.height);
      const frame = Math.round(((e.clientX - rect.left - a.left) / (a.right - a.left)) * (this.total - 1));
      onPick(Math.min(this.drawn - 1, Math.max(0, frame)));
    };
    canvas.addEventListener('pointerdown', (e) => { canvas.setPointerCapture(e.pointerId); pick(e); });
    canvas.addEventListener('pointermove', (e) => { if (e.buttons) pick(e); });
  }

  area(width, height) {
    return { left: 44, right: width - 14, top: 16, bottom: height - 40 };
  }

  reset(pixels, total, length, stages) {
    Object.assign(this, { pixels, total, length, stages, drawn: 0, marker: null });
    this.image = Object.assign(document.createElement('canvas'), { width: total, height: pixels });
    const rgb = (name) => [1, 3, 5].map((i) => parseInt(css(name).slice(i, i + 2), 16));
    this.colors = [rgb('--surface'), rgb('--glass')];
  }

  add(rho) {
    const column = new ImageData(1, this.pixels);
    const [from, to] = this.colors;
    rho.forEach((r, k) => {
      const at = 4 * (this.pixels - 1 - k);
      const mix = 0.1 + 0.9 * r;
      for (let c = 0; c < 3; c++) column.data[at + c] = from[c] + mix * (to[c] - from[c]);
      column.data[at + 3] = 255;
    });
    this.image.getContext('2d').putImageData(column, this.drawn++, 0);
  }

  draw() {
    const { ctx, width, height } = fit(this.canvas);
    if (!this.total) return;
    const a = this.area(width, height);
    const step = (a.right - a.left) / (this.total - 1);
    const sx = (frame) => a.left + frame * step;
    ctx.save();
    ctx.beginPath();
    ctx.rect(a.left, a.top, a.right - a.left, a.bottom - a.top);
    ctx.clip();
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(this.image, sx(-0.5), a.top, step * this.total, a.bottom - a.top);
    ctx.restore();
    const sy = (x) => a.bottom - (x / this.length) * (a.bottom - a.top);
    drawYAxis(ctx, a, sy, ticks(0, this.length, 3));
    ctx.fillStyle = css('--muted');
    ctx.textAlign = 'right';
    ctx.textBaseline = 'top';
    ctx.fillText('Iteration, with position (µm) upward', a.right, a.bottom + 7);
    // The sharpness schedule explains the jumps in the chart above.
    ctx.textAlign = 'left';
    ctx.textBaseline = 'bottom';
    for (const { frame, name } of this.stages) {
      ctx.strokeStyle = css('--muted');
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(sx(frame), a.top - 12);
      ctx.lineTo(sx(frame), a.bottom);
      ctx.stroke();
      ctx.fillText(name, sx(frame) + 3, a.top - 1);
    }
    if (this.marker !== null) {
      ctx.strokeStyle = css('--ink');
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.moveTo(sx(this.marker), a.top);
      ctx.lineTo(sx(this.marker), a.bottom);
      ctx.stroke();
    }
  }
}

// A line chart with a crosshair: spectra and band diagrams. Clicking picks a point.
export class LineChart {
  constructor(canvas, onPick) {
    this.canvas = canvas;
    this.onPick = onPick;
    this.spec = null;
    this.hover = null;
    canvas.addEventListener('pointermove', (e) => {
      this.hover = this.nearest(e);
      this.draw();
    });
    canvas.addEventListener('pointerleave', () => { this.hover = null; this.draw(); });
    canvas.addEventListener('click', (e) => {
      const hit = this.nearest(e);
      if (hit) this.onPick(hit);
    });
  }

  area(width, height) {
    return { left: 44, right: width - (this.spec?.endLabels ? 96 : 14), top: 12, bottom: height - 40 };
  }

  nearest(e) {
    if (!this.spec) return null;
    const rect = this.canvas.getBoundingClientRect();
    const a = this.area(rect.width, rect.height);
    const { xMin, xMax, yMin, yMax, series } = this.spec;
    const x = xMin + ((e.clientX - rect.left - a.left) / (a.right - a.left)) * (xMax - xMin);
    const y = yMax - ((e.clientY - rect.top - a.top) / (a.bottom - a.top)) * (yMax - yMin);
    const xs = series[0].x;
    let index = 0;
    for (let i = 1; i < xs.length; i++) if (Math.abs(xs[i] - x) < Math.abs(xs[index] - x)) index = i;
    let seriesIndex = 0;
    series.forEach((s, n) => {
      if (Math.abs(s.y[index] - y) < Math.abs(series[seriesIndex].y[index] - y)) seriesIndex = n;
    });
    return { x: xs[index], seriesIndex, index };
  }

  draw() {
    const { ctx, width, height } = fit(this.canvas);
    if (!this.spec) return;
    const { xMin, xMax, yMin, yMax, series, xLabel, gaps, marker, format, endLabels, yTicks, yFormat, levels, xTicks, xFormat } = this.spec;
    const a = this.area(width, height);
    const sx = (x) => a.left + ((x - xMin) / (xMax - xMin)) * (a.right - a.left);
    const sy = (y) => a.bottom - ((y - yMin) / (yMax - yMin)) * (a.bottom - a.top);

    ctx.fillStyle = css('--glass');
    ctx.globalAlpha = 0.16;
    for (const [lo, hi] of gaps ?? []) ctx.fillRect(a.left, sy(hi), a.right - a.left, sy(lo) - sy(hi));
    ctx.globalAlpha = 1;
    drawYAxis(ctx, a, sy, yTicks ?? ticks(yMin, yMax, 4), yFormat);
    drawXAxis(ctx, a, sx, xTicks ?? ticks(xMin, xMax, Math.max(3, Math.floor(width / 110))), xLabel, xFormat);

    // Dashed reference values, labelled at the right edge.
    for (const { y, name } of levels ?? []) {
      ctx.strokeStyle = css('--muted');
      ctx.lineWidth = 1;
      ctx.setLineDash([5, 4]);
      ctx.beginPath();
      ctx.moveTo(a.left, sy(y));
      ctx.lineTo(a.right, sy(y));
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = css('--muted');
      ctx.textAlign = 'right';
      ctx.textBaseline = 'bottom';
      ctx.fillText(name, a.right - 4, sy(y) - 2);
    }

    ctx.save();
    ctx.beginPath();
    ctx.rect(a.left, a.top - 2, a.right - a.left, a.bottom - a.top + 4);
    ctx.clip();
    ctx.lineWidth = 2;
    ctx.lineJoin = 'round';
    for (const s of series) {
      ctx.strokeStyle = s.color;
      ctx.beginPath();
      s.x.forEach((x, i) => ctx.lineTo(sx(x), sy(s.y[i])));
      ctx.stroke();
    }
    ctx.restore();

    if (endLabels) {
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';
      const taken = [];
      for (const s of series) {
        let y = sy(s.y[s.y.length - 1]);
        while (taken.some((t) => Math.abs(t - y) < 15)) y += 15;
        taken.push(y);
        ctx.fillStyle = s.color;
        ctx.fillRect(a.right + 8, y - 1, 12, 2);
        ctx.fillStyle = css('--ink');
        ctx.fillText(s.name, a.right + 25, y);
      }
    }
    if (marker !== undefined) this.vertical(ctx, sx(marker), a, css('--ink'), 1.5);
    if (this.hover) this.drawHover(ctx, a, sx, sy, format, width);
  }

  vertical(ctx, x, a, color, lineWidth) {
    ctx.strokeStyle = color;
    ctx.lineWidth = lineWidth;
    ctx.beginPath();
    ctx.moveTo(x, a.top);
    ctx.lineTo(x, a.bottom);
    ctx.stroke();
  }

  drawHover(ctx, a, sx, sy, format, width) {
    const { series } = this.spec;
    const { index, seriesIndex } = this.hover;
    // The chart may have changed under a resting pointer.
    if (!series[seriesIndex] || index >= series[0].x.length) return;
    const x = sx(series[0].x[index]);
    ctx.globalAlpha = 0.5;
    this.vertical(ctx, x, a, css('--muted'), 1);
    ctx.globalAlpha = 1;
    const s = series[seriesIndex];
    ctx.fillStyle = s.color;
    ctx.strokeStyle = css('--surface');
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(x, sy(s.y[index]), 5, 0, 2 * Math.PI);
    ctx.fill();
    ctx.stroke();

    const lines = format(index, seriesIndex);
    const boxWidth = Math.max(...lines.map((l) => ctx.measureText(l).width)) + 16;
    const boxX = x + boxWidth + 20 > width ? x - boxWidth - 10 : x + 10;
    ctx.fillStyle = css('--tooltip');
    ctx.fillRect(boxX, a.top + 4, boxWidth, lines.length * 17 + 10);
    ctx.fillStyle = css('--tooltip-ink');
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    lines.forEach((l, i) => ctx.fillText(l, boxX + 8, a.top + 10 + i * 17));
  }
}
