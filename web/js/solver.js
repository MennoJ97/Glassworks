// Client for the solver worker. Only the newest request matters: while dragging,
// older solves are dropped, and a slow solve is abandoned by restarting the worker.

const ABANDON_AFTER_MS = 400;

export class Solver {
  constructor() {
    this.nextId = 0;
    this.active = null;
    this.queued = null;
    this.spawn();
  }

  spawn() {
    this.worker = new Worker(new URL('../worker.js', import.meta.url), { type: 'module' });
    this.worker.onmessage = ({ data }) => this.finish(data);
    this.worker.onerror = () =>
      this.finish({ error: 'The solver stopped unexpectedly. Try a coarser mesh.', fatal: true });
  }

  restart() {
    this.worker.terminate();
    this.spawn();
  }

  // Resolves with the results, or with null if a newer request replaced this one.
  run(requests) {
    return new Promise((resolve) => {
      this.queued?.resolve(null);
      this.queued = { id: ++this.nextId, requests, resolve };
      if (!this.active) {
        this.dispatch();
      } else if (performance.now() - this.active.started > ABANDON_AFTER_MS) {
        this.active.resolve(null);
        this.active = null;
        this.restart();
        this.dispatch();
      }
    });
  }

  dispatch() {
    if (!this.queued) return;
    this.active = { ...this.queued, started: performance.now() };
    this.queued = null;
    this.worker.postMessage({ id: this.active.id, requests: this.active.requests });
  }

  finish(data) {
    const job = this.active;
    this.active = null;
    if (data.fatal) this.restart();
    job?.resolve(this.queued ? null : data);
    this.dispatch();
  }
}
