import json,statistics as st,sys
R=[json.loads(l) for f in sys.argv[1:] for l in open(f) if l.startswith('{')]
tags=[]
for r in R:
    if r['tag'] not in tags: tags.append(r['tag'])
print('fails',[(r['tag'],r['seed']) for r in R if isinstance(r.get('fail'),str)], 'draw fails',sum(r.get('fail')==1 for r in R))
for t in tags:
    v=[r for r in R if r['tag']==t]; b=[r for r in v if r.get('base')]; seeds=sorted(r['seed'] for r in b)
    f=lambda k:[r[k] for r in b]
    print(f"\n{t}: n={len(b)} LDOS {min(f('F')):.0f}/{st.median(f('F')):.0f}/{max(f('F')):.0f} shift {min(f('slope')):.3f}/{st.median(f('slope')):.3f}/{max(f('slope')):.3f} Q {min(f('Q')):.0f}/{st.median(f('Q')):.0f}/{max(f('Q')):.0f} below1e4 {sum(x<1e4 for x in f('Q'))} fill {st.median(f('fill')):.2f} time {st.median(f('secs')):.0f}s")
    fu=[st.median(r['F1'] for r in v if r.get('nm') in (1,-1) and r['seed']==s) for s in seeds]
    print(f"   uniform 1 nm LDOS at original wavelength: {st.median(fu):.1f}, 2 nm: {st.median([st.median(r['F1'] for r in v if r.get('nm') in (2,-2) and r['seed']==s) for s in seeds]):.1f}")
    for g in (1,2,5):
        q=[];fo=[];w=[];fn=[]
        for s in seeds:
            d=[r for r in v if r.get('sigma')==g and r['seed']==s and 'Q' in r]
            if not d: continue
            q.append(st.median(r['Q'] for r in d));fo.append(st.median(r['Fown'] for r in d));w.append(st.pstdev(r['w'] for r in d)*100);fn.append(st.median(r['F1'] for r in v if r.get('sigma')==g and r['seed']==s))
        print(f"   sigma {g}: Q {min(q):.0f}/{st.median(q):.0f}/{max(q):.0f} | retuned {min(fo):.0f}/{st.median(fo):.0f} | spread {st.median(w):.3f}% | at original wavelength {min(fn):.0f}/{st.median(fn):.0f}/{max(fn):.0f}")
