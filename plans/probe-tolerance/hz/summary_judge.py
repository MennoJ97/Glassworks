# Compares tags in log/judge (and log/corr as plain) at correlated disorder, against the crystal in log/tol.
import json, glob, os, statistics as st, sys
here = os.path.dirname(os.path.abspath(__file__)); med = st.median
def read(f): return [json.loads(l) for l in open(f) if l.startswith("{")]
plain = {}
for f in glob.glob(f"{here}/log/corr/*.log"): plain[int(os.path.basename(f)[:-4])] = read(f)
J = {}
for f in glob.glob(f"{here}/log/judge/*.log"):
    n = os.path.basename(f)[:-4]; t, s = n.rsplit("-", 1); r = read(f)
    if any(x.get("finished") for x in r): J.setdefault(t, {})[int(s)] = r
J["plain"] = plain
crystal = read(f"{here}/log/tol/crystal.log"); cb = next(x for x in crystal if x.get("base"))
cm = {0: cb["F0"], **{s: med(x["Fx"] for x in crystal if x.get("sigma") == s) for s in (2, 5)}}
tags = sys.argv[1:] or sorted(J)
def seeds_of(t, only): return sorted(s for s in J[t] if s in only) if only else sorted(J[t])
def row(t, only=None):
    S = seeds_of(t, only); b = {s: next(x for x in J[t][s] if x.get("base")) for s in S}
    per = {g: {s: med(x["Fx"] for x in J[t][s] if x.get("sigma") == g) for s in S} for g in (2, 5)}
    q = {g: {s: med(x["Qx"] for x in J[t][s] if x.get("sigma") == g and x.get("Qx")) for s in S} for g in (2, 5)}
    o = {g: {s: med(x["Fx_orig"] for x in J[t][s] if x.get("sigma") == g) for s in S} for g in (2, 5)}
    print(f"{t} n={len(S)}: F {med(b[s]['F0'] for s in S):.0f} (min {min(b[s]['F0'] for s in S):.0f}, max {max(b[s]['F0'] for s in S):.0f}) Q {med(b[s]['Q'] for s in S):.0f}")
    for g in (2, 5):
        v = per[g]; print(f"   {g} nm own F median {med(v.values()):.0f}, best {max(v.values()):.0f}, worst {min(v.values()):.0f}; beats crystal {cm[g]:.0f} on {sum(x > cm[g] for x in v.values())}; kept {med(v[s] / b[s]['F0'] for s in S):.2f}; Q {med(q[g].values()):.0f}; at old w {med(o[g].values()):.0f}")
    return S, b, per
print("crystal", {k: round(v) for k, v in cm.items()})
for t in tags: row(t)
for t in tags:
    if t == "plain": continue
    common = sorted(set(J[t]) & set(plain))
    bt = {s: next(x for x in J[t][s] if x.get("base"))["F0"] for s in common}; bp = {s: next(x for x in plain[s] if x.get("base"))["F0"] for s in common}
    print(f"paired {t} over plain, n={len(common)}: no disorder better on {sum(bt[s] > bp[s] for s in common)}, ratio {med(bt[s] / bp[s] for s in common):.2f}", end="")
    for g in (2, 5):
        a = {s: med(x["Fx"] for x in J[t][s] if x.get("sigma") == g) for s in common}; c = {s: med(x["Fx"] for x in plain[s] if x.get("sigma") == g) for s in common}
        print(f"; {g} nm better on {sum(a[s] > c[s] for s in common)}, ratio {med(a[s] / c[s] for s in common):.2f}", end="")
    print()
