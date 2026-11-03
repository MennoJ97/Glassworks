# All tags in log/judge2 on the unseen starts, paired with the page's current robust goal, and against the crystal.
import json, glob, os, statistics as st, sys
here = os.path.dirname(os.path.abspath(__file__)); med = st.median
def read(f): return [json.loads(l) for l in open(f) if l.startswith("{")]
J = {}
for f in glob.glob(f"{here}/log/judge2/*.log"):
    t, s = os.path.basename(f)[:-4].rsplit("-", 1); r = read(f)
    if any(x.get("finished") for x in r): J.setdefault(t, {})[int(s)] = r
crystal = read(f"{here}/log/tol/crystal.log"); extra = read(f"{here}/log/xdraws.log"); cb = next(x for x in crystal if x.get("base"))
cm = {0: cb["F0"], **{g: med(x["Fx"] for x in crystal + extra if x.get("sigma") == g) for g in (2, 5)}}
print("crystal", {k: round(v) for k, v in cm.items()})
def base(r): return next(x for x in r if x.get("base"))
def ate(r, g, k): return med(x[k] for x in r if x.get("sigma") == g and x.get(k) is not None)
for t in sys.argv[1:]:
    S = sorted(J[t]); F = [base(J[t][s])["F0"] for s in S]; Q = [base(J[t][s])["Q"] for s in S]
    print(f"{t} n={len(S)}: F {min(F):.0f}/{med(F):.0f}/{max(F):.0f} beats {sum(x > cm[0] for x in F)}; Q {med(Q):.0f}")
    for g in (2, 5):
        v = [ate(J[t][s], g, "Fx") for s in S]; k = [ate(J[t][s], g, "Fx") / base(J[t][s])["F0"] for s in S]; q = [ate(J[t][s], g, "Qx") for s in S]; o = [ate(J[t][s], g, "Fx_orig") for s in S]
        print(f"   {g} nm: F {min(v):.0f}/{med(v):.0f}/{max(v):.0f} beats {cm[g]:.0f} on {sum(x > cm[g] for x in v)}; kept {med(k):.2f}; Q {med(q):.0f}; at old w {med(o):.0f}")
    if t != "rc2":
        C = sorted(set(S) & set(J["rc2"])); a = [base(J[t][s])["F0"] / base(J["rc2"][s])["F0"] for s in C]
        line = f"   against rc2 n={len(C)}: none better {sum(x > 1 for x in a)}, ratio {med(a):.2f}"
        for g in (2, 5):
            a = [ate(J[t][s], g, "Fx") / ate(J["rc2"][s], g, "Fx") for s in C]; line += f"; {g} nm better {sum(x > 1 for x in a)}, ratio {med(a):.2f}"
        print(line)
for t in sys.argv[1:]:
    r = []
    for f in glob.glob(f"{here}/log/chk-{t}-*.s"):
        try: r.append(json.loads(open(f).read()))
        except Exception: pass
    if r: print(t, "shift", round(med(x["shift"] for x in r), 2), "range", round(min(x["shift"] for x in r), 2), round(max(x["shift"] for x in r), 2), "source in material", sum(x["source"] == "material" for x in r), "of", len(r))
