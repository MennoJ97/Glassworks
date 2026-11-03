import json,glob,os,statistics as st,sys
# Summary of log/mix runs: LDOS, Q, how well each step was predicted, and the drop at each stage's first iteration.
here=os.path.dirname(os.path.abspath(__file__))
med=st.median
for t in sys.argv[1:]:
    r=[json.loads(l) for l in open(f'{here}/log/mix/{t}.log') if l.startswith('{"seed"')]
    F=sorted(x['F'] for x in r); Q=sorted(x['Q'] for x in r if x.get('Q'))
    err=[];hold=[];drops=[];swing=[]
    for f in glob.glob(f'{here}/log/mix/trace-{t}-*.txt'):
        v={x['iteration']:x for x in (json.loads(l) for l in open(f) if l.startswith('{'))}
        for i,x in v.items():
            bw=x['pole']/x['Q']
            if i>100 and i-1 in v and v[i-1]['Q']==x['Q']: err.append(abs(x['pole']-v[i-1]['predicted'])/bw)
            if i>100: hold.append(abs(x['held']-x['want'])/(2*3.141592653589793/x['pole']/x['Q']))
        starts=[i for i in v if i-1 in v and v[i-1]['Q']!=v[i]['Q'] and i>40]
        if starts: drops.append(min(v[i]['F']/max(v[i-1]['F'],1e-9) for i in starts))
        late=[v[i]['F'] for i in v if i>120 and i not in starts]
        if late: swing.append(max(late)/max(min(late),1e-9))
    n=len(r)
    print(f"{t}: n={n} F {F[0]:.0f}/{med(F):.0f}/{F[-1]:.0f} Q {Q[0]:.0f}/{med(Q):.0f}/{Q[-1]:.0f} over-371 {sum(x>371 for x in F)} sec {med(x['sec'] for x in r):.0f} | prediction error {med(err):.2f} bw | hold residual {med(hold):.3f} bw | worst stage-start ratio {med(drops):.3f} | late swing x{med(swing):.1f}")
