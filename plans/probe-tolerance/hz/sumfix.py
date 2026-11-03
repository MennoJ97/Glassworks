import json,glob,os,statistics as st,sys
here=os.path.dirname(os.path.abspath(__file__))
med=st.median
for t in sys.argv[1:]:
    r=[json.loads(l) for l in open(f'{here}/log/fix/{t}.log') if l.startswith('{"seed"')]
    F=sorted(x['F'] for x in r); Q=sorted(x['Q'] for x in r if x.get('Q'))
    off=[];err=[];swing=[]
    for f in glob.glob(f'{here}/log/fix/trace-{t}-*.txt'):
        rows=[json.loads(l) for l in open(f) if l.startswith('{')]
        v=[x for x in rows if 'off' in x]; p={x['iteration']:x['predicted'] for x in rows if x.get('pred')}
        off+=[abs(x['off'])>1 for x in v if x['iteration']>100]
        late=[x['F'] for x in v if x['iteration']>120]
        if late: swing.append(max(late)/max(min(late),1e-9))
        err+=[abs(x['pole']-p[x['iteration']-1])/x['pole']*x['Q'] for x in v if x['iteration']-1 in p and x['iteration']>100]
    print(f"{t}: n={len(r)} F {F[0]:.0f}/{med(F):.0f}/{F[-1]:.0f} Q {Q[0]:.0f}/{med(Q):.0f}/{Q[-1]:.0f} over-371 {sum(x>371 for x in F)} sec {med(x['sec'] for x in r):.0f} | late off>1 {100*sum(off)/max(1,len(off)):.0f}% | prediction error, bandwidths {med(err) if err else float('nan'):.2f} | swing x{med(swing) if swing else 0:.0f}")
