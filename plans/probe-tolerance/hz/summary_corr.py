# Summarizes log/corr: correlated disorder on the designs, against log/tol's crystal and uncorrelated results.
import json, glob, os, statistics as st
here = os.path.dirname(os.path.abspath(__file__)); med = st.median
def rng(v): v = sorted(v); return f"{v[0]:.4g}/{med(v):.4g}/{v[-1]:.4g}" if v else "-"
def rows(sub):
    R = {}
    for f in glob.glob(f"{here}/log/{sub}/*.log"):
        R[os.path.basename(f)[:-4]] = [json.loads(l) for l in open(f) if l.startswith("{")]
    return R
C, T = rows("corr"), rows("tol")
ids = sorted(C, key=int)
base = {i: next(x for x in C[i] if x.get("base")) for i in ids}
crystal = T["crystal"]; cb = next(x for x in crystal if x.get("base"))
print("crystal base", cb["F0"], "w", cb["w"], "Q", cb["Q"])
for s in [2, 5]:
    d = [(x, base[i]) for i in ids for x in C[i] if x.get("sigma") == s]
    cr = [x for x in crystal if x.get("sigma") == s]
    u = [(i, x) for i in ids for x in T[i] if x.get("sigma") == s]
    def moved(r): return (sum((((x["wx"] + x["wy"]) / 2 - b["w"]) / b["w"] * 100) ** 2 for x, b in r) / len(r)) ** .5
    print(f"sigma {s} correlated designs n={len(d)}: own F {rng([x['Fx'] for x, b in d])}, share {rng([x['Fx'] / b['F0'] for x, b in d])}, F at old w {rng([x['Fx_orig'] for x, b in d])}, apart {rng([abs(x['wx'] - x['wy']) * 100 for x, b in d])}, moved rms {moved(d):.3f}, Q ratio {rng([x['Qx'] / b['Q'] for x, b in d if x.get('Qx')])}")
    print(f"   crystal n={len(cr)}: own F {rng([x['Fx'] for x in cr])}, F at old w {rng([x['Fx_orig'] for x in cr])}, Q {rng([x['Qx'] for x in cr if x.get('Qx')])}")
    per = {i: med([x['Fx'] for x in C[i] if x.get('sigma') == s]) for i in ids}
    cm = med([x['Fx'] for x in cr])
    print(f"   designs whose median own F beats the crystal's {cm:.0f}: {sum(v > cm for v in per.values())} of 18; best {max(per.values()):.0f}, median {med(per.values()):.0f}")
    pu = {i: med([x['Fx'] for x in T[i] if x.get('sigma') == s]) for i in ids}
    print(f"   same designs, uncorrelated: beats {sum(v > cm for v in pu.values())} of 18; best {max(pu.values()):.0f}, median {med(pu.values()):.0f}")
s1 = [x for x in crystal if x.get("sigma") == 1]; cm = med([x['Fx'] for x in s1]); pu = {i: med([x['Fx'] for x in T[i] if x.get('sigma') == 1]) for i in ids}
print(f"sigma 1 uncorrelated: crystal {cm:.0f}, beats {sum(v > cm for v in pu.values())} of 18; best {max(pu.values()):.0f}, median {med(pu.values()):.0f}")
b = {i: next(x for x in T[i] if x.get('base'))['F0'] for i in ids}; print(f"no disorder: crystal {cb['F0']:.0f}, beats {sum(v > cb['F0'] for v in b.values())} of 18; best {max(b.values()):.0f}, median {med(b.values()):.0f}")
