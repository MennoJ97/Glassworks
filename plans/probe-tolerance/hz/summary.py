# Summarizes log/tol: baselines, uniform moves and disorder, for the designs and for the crystal.
import json, glob, os, statistics as st, sys
here = os.path.dirname(os.path.abspath(__file__))
med = st.median
def rng(v): v = sorted(v); return f"{v[0]:.4g}/{med(v):.4g}/{v[-1]:.4g}" if v else "-"
R = {}
for f in glob.glob(here + "/log/tol/*.log"):
    id = os.path.basename(f)[:-4]; R[id] = [json.loads(l) for l in open(f) if l.startswith("{")]
def base(r): return next(x for x in r if x.get("base"))
def nm(r, v): return next((x for x in r if x.get("nm") == v and "w" in x), None)
def describe(name, ids):
    print(f"\n== {name}, n={len(ids)}")
    B = {i: base(R[i]) for i in ids}
    print(" base  w", rng([B[i]['w'] for i in ids]), "Q", rng([B[i]['Q'] for i in ids]), "F along x", rng([B[i]['F0'] for i in ids]), "F along y", rng([B[i]['F90'] for i in ids]))
    print(" below Q 1000:", sum(B[i]['Q'] < 1000 for i in ids))
    for v in [0.5, 1, 2, 5]:
        up = [(nm(R[i], v), nm(R[i], -v), B[i]) for i in ids]
        ok = [(u, d, b) for u, d, b in up if u and d]
        print(f" +-{v} nm  n={len(ok)} shift %/nm", rng([(u['w'] - d['w']) / (2 * v) / b['w'] * 100 for u, d, b in ok]), "F own/base +", rng([u['Fown'] / b['Fmode'] for u, d, b in ok]), "-", rng([d['Fown'] / b['Fmode'] for u, d, b in ok]), "F at the old wavelength +", rng([u['Forig'] for u, d, b in ok]), "-", rng([d['Forig'] for u, d, b in ok]), "Q ratio +", rng([u['Q'] / b['Q'] for u, d, b in ok]), "-", rng([d['Q'] / b['Q'] for u, d, b in ok]))
    for s in [1, 2, 5]:
        rows = [(x, B[i]) for i in ids for x in R[i] if x.get('sigma') == s]
        fx = [x['Fx'] / b['F0'] for x, b in rows]; fy = [x['Fy'] / b['F90'] for x, b in rows]
        split = [abs(x['wx'] - x['wy']) / b['w'] * 100 for x, b in rows]
        moved = [((x['wx'] + x['wy']) / 2 - b['w']) / b['w'] * 100 for x, b in rows]
        q = [min(x.get('Qx', 1e12), x.get('Qy', 1e12)) for x, b in rows if 'Qx' in x or 'Qy' in x]
        qr = [x['Qx'] / b['Q'] for x, b in rows if x.get('Qx')]
        print(f" sigma {s} draws={len(rows)} F own, x, abs", rng([x['Fx'] for x, b in rows]), "share of base", rng(fx), "y share", rng(fy), "| best of x, y over base", rng([max(x['Fx'] / b['F0'], x['Fy'] / b['F90']) for x, b in rows]), "| peaks apart %", rng(split), "| mean moved rms %", f"{(sum(m * m for m in moved) / len(moved)) ** .5:.3f}", "| F at old w, abs", rng([x['Fx_orig'] for x, b in rows]), "| Q x", rng([x['Qx'] for x, b in rows if x.get('Qx')]), "ratio", rng(qr))
designs = sorted((i for i in R if i != 'crystal'), key=int)
describe("Designs", designs)
describe("Crystal", ['crystal'])
print("\nby design at sigma 2: base F, Q, shift, median own F")
for i in designs:
    b = base(R[i]); u, d = nm(R[i], 1), nm(R[i], -1); r = [x for x in R[i] if x.get('sigma') == 2]
    print(i, f"F {b['F0']:.0f} Q {b['Q']:.0f} shift {(u['w'] - d['w']) / 2 / b['w'] * 100:.2f}", "own", rng([x['Fx'] for x in r]), "Q", rng([x['Qx'] for x in r if x.get('Qx')]))
