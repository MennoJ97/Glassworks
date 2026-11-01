// Client for the solver workers. Only the newest request matters: while dragging,
// older solves are dropped, and a slow solve is abandoned by restarting the workers.
// Requests sent together are independent, so they share out over a few workers. Each
// worker's memory is its own, and a large mesh can take a gigabyte, so there aren't many.

const ABANDON_AFTER_MS = 400;
const MAX_WORKERS = Math.min(3, Math.max(1, (navigator.hardwareConcurrency ?? 2) - 1));

export class Solver {
  constructor() {
    this.nextId = 0;
    this.active = null;
    this.queued = null;
    this.workers = [];
  }

  // Starts workers only when a request first needs them.
  worker(index) {
    if (!this.workers[index]) {
      const worker = new Worker(new URL('../worker.js', import.meta.url), { type: 'module' });
      worker.onmessage = ({ data }) => this.finish(data, index);
      worker.onerror = () =>
        this.finish({ id: this.active?.id, error: 'The solver stopped unexpectedly. Try a coarser mesh.', fatal: true }, index);
      this.workers[index] = worker;
    }
    return this.workers[index];
  }

  restart() {
    this.workers.forEach((worker) => worker.terminate());
    this.workers = [];
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
    const { requests } = this.queued;
    const parts = Math.min(MAX_WORKERS, Math.max(1, requests.length));
    this.active = { ...this.queued, started: performance.now(), parts, waiting: parts, results: [], ms: 0 };
    this.queued = null;
    for (let part = 0; part < parts; part++) {
      const share = requests.filter((_, i) => i % parts === part);
      this.worker(part).postMessage({ id: this.active.id, requests: share });
    }
  }

  finish(data, part) {
    const job = this.active;
    // A worker's answer may already be on its way when its job is dropped.
    if (data.id !== job?.id) return;
    data.results?.forEach((result, i) => { job.results[part + i * job.parts] = result; });
    job.ms = Math.max(job.ms, data.ms ?? 0);
    job.error ??= data.error;
    job.fatal ||= data.fatal;
    if (--job.waiting > 0 && !job.fatal) return;
    this.active = null;
    if (job.fatal) this.restart();
    const reply = job.error ? { id: job.id, error: job.error, fatal: job.fatal } : { id: job.id, results: job.results, ms: job.ms };
    job.resolve(this.queued ? null : reply);
    this.dispatch();
  }
}
