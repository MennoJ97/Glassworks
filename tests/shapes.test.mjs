// What dragging does to a shape, without a page.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { dragShape, grabShape } from '../web/js/model.js';

const circle = () => ({ type: 'circle', x: 1, y: 1, w: 0.4, h: 0.4, n: 2 });
const rect = () => ({ type: 'rect', x: 1, y: 1, w: 0.6, h: 0.2, n: 2 });
const reach = 0.02;

test('a pointer grabs the outline, the inside or nothing', () => {
  assert.equal(grabShape([circle()], 1.2, 1, reach).part, 'rim');
  assert.equal(grabShape([circle()], 1.05, 1, reach).part, 'body');
  assert.equal(grabShape([circle()], 1.3, 1, reach), null);
  assert.deepEqual(grabShape([rect()], 1.3, 1, reach), { index: 0, part: 'edge', sideX: 1, sideY: 0 });
  assert.deepEqual(grabShape([rect()], 0.7, 0.9, reach), { index: 0, part: 'edge', sideX: -1, sideY: -1 });
  assert.equal(grabShape([rect()], 1, 1, reach).part, 'body');
  assert.equal(grabShape([rect()], 1, 1.2, reach), null);
});

test('later shapes are grabbed first, and even a tiny one can be moved', () => {
  assert.equal(grabShape([circle(), rect()], 1, 1, reach).index, 1);
  const tiny = { ...circle(), w: 0.01 };
  assert.equal(grabShape([tiny], 1.001, 1, reach).part, 'body');
});

test('moving keeps the pointer on the same spot of the shape', () => {
  const s = circle();
  const grab = grabShape([s], 1.05, 1.02, reach);
  dragShape(s, grab, 1.55, 0.52);
  assert.deepEqual([s.x, s.y, s.w], [1.5, 0.5, 0.4]);
});

test('resizing a circle follows the pointer around its rim', () => {
  const s = circle();
  dragShape(s, grabShape([s], 1.2, 1, reach), 1, 1.5);
  assert.deepEqual([s.x, s.y, s.w], [1, 1, 1]);
});

test('dragging an edge leaves the opposite edge in place', () => {
  const s = rect();
  dragShape(s, grabShape([s], 1.3, 1, reach), 1.5, 1);
  assert.deepEqual([s.x - s.w / 2, s.x + s.w / 2, s.y, s.h].map((v) => +v.toFixed(9)), [0.7, 1.5, 1, 0.2]);
  dragShape(s, { part: 'edge', sideX: -1, sideY: 1 }, 0.5, 1.4);
  assert.deepEqual([s.x - s.w / 2, s.x + s.w / 2, s.y - s.h / 2, s.y + s.h / 2].map((v) => +v.toFixed(9)), [0.5, 1.5, 0.9, 1.4]);
});

test('an edge can\'t cross its opposite edge', () => {
  const s = rect();
  dragShape(s, { part: 'edge', sideX: 1, sideY: 0 }, 0, 1);
  assert.equal(s.w, 0.005);
  assert.ok(Math.abs(s.x - s.w / 2 - 0.7) < 1e-9);
  const c = circle();
  dragShape(c, { part: 'rim' }, 1, 1);
  assert.equal(c.w, 0.005);
});
