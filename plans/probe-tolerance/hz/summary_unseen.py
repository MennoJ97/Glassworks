# Paired plain and rc2 on unseen starts, and the crystal with all its draws.
import json, glob, os, statistics as st
here = os.path.dirname(os.path.abspath(__file__)); med = st.median
def read(f): return [json.loads(l) for l in open(f) if l.startswith("{")]
J = {}
for f in glob.glob(f"{here}/log/judge2/*.log"):
    t, s = os.path.basename(f)[:-4].rsplit("-", 1); r = read(f)
    if any(x.get("finished") for x in r): J.setdefault(t, {})[int(s)] = r
crystal = read(f"{here}/log/tol/crystal.log"); extra = read(f"{here}/log/xdraws.log"); cb = next(x for x in crystal if x.get("base"))
cm = {0: cb["F0"]}
for g in (2, 5):
    v = [x["Fx"] for x in crystal + extra if x.get("sigma") == g]; cm[g] = med(v); print(f"crystal {g} nm: n={len(v)} median {med(v):.0f}, range {min(v):.0f} to {max(v):.0f}")
def base(r): return next(x for x in r if x.get("base"))
def ate(r, g, k): return med(x[k] for x in r if x.get("sigma") == g and x.get(k) is not None)
common = sorted(set(J["p"]) & set(J["rc2"])); print("common", len(common))
for t in ("p", "rc2"):
    print(f"{t}: F {min(base(J[t][s])['F0'] for s in common):.0f}/{med(base(J[t][s])['F0'] for s in common):.0f}/{max(base(J[t][s])['F0'] for s in common):.0f}, beats crystal {sum(base(J[t][s])['F0'] > cm[0] for s in common)}; Q {med(base(J[t][s])['Q'] for s in common):.0f}, shift n/a")
    for g in (2, 5):
        v = [ate(J[t][s], g, "Fx") for s in common]; q = [ate(J[t][s], g, "Qx") for s in common]; k = [ate(J[t][s], g, "Fx") / base(J[t][s])["F0"] for s in common]; o = [ate(J[t][s], g, "Fx_orig") for s in common]
        print(f"   {g} nm: F {min(v):.0f}/{med(v):.0f}/{max(v):.0f}, beats crystal {cm[g]:.0f} on {sum(x > cm[g] for x in v)}; kept {med(k):.2f}; Q {med(q):.0f}; at old w {med(o):.0f}")
a = {s: base(J["rc2"][s])["F0"] / base(J["p"][s])["F0"] for s in common}; print(f"paired no disorder: better {sum(x > 1 for x in a.values())}/{len(a)}, ratio {med(a.values()):.2f}")
aq = {s: base(J["rc2"][s])["Q"] / base(J["p"][s])["Q"] for s in common}; print(f"paired Q: better {sum(x > 1 for x in aq.values())}/{len(aq)}, ratio {med(aq.values()):.2f}")
for g in (2, 5):
    a = {s: ate(J["rc2"][s], g, "Fx") / ate(J["p"][s], g, "Fx") for s in common}; k = {s: (ate(J["rc2"][s], g, "Fx") / base(J["rc2"][s])["F0"]) / (ate(J["p"][s], g, "Fx") / base(J["p"][s])["F0"]) for s in common}
    print(f"paired {g} nm: better {sum(x > 1 for x in a.values())}/{len(a)}, ratio {med(a.values()):.2f}; share kept better {sum(x > 1 for x in k.values())}/{len(k)}, ratio {med(k.values()):.2f}")
