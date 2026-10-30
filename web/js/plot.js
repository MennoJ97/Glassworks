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

function drawXAxis(ctx, area, scale, values, title) {
  ctx.fillStyle = css('--muted');
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (const v of values) ctx.fillText(label(v), scale(v), area.bottom + 7);
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
    const { xMin, xMax, yMin, yMax, series, xLabel, gaps, marker, format, endLabels, yTicks, yFormat, levels } = this.spec;
    const a = this.area(width, height);
    const sx = (x) => a.left + ((x - xMin) / (xMax - xMin)) * (a.right - a.left);
    const sy = (y) => a.bottom - ((y - yMin) / (yMax - yMin)) * (a.bottom - a.top);

    ctx.fillStyle = css('--glass');
    ctx.globalAlpha = 0.16;
    for (const [lo, hi] of gaps ?? []) ctx.fillRect(a.left, sy(hi), a.right - a.left, sy(lo) - sy(hi));
    ctx.globalAlpha = 1;
    drawYAxis(ctx, a, sy, yTicks ?? ticks(yMin, yMax, 4), yFormat);
    drawXAxis(ctx, a, sx, ticks(xMin, xMax, Math.max(3, Math.floor(width / 110))), xLabel);

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
