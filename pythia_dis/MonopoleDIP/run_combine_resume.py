"""Resume combine.py: computes only the (spin,m,T) keys missing from res_comb.json."""
import json,os,combine
f='res_comb.json'; out=json.load(open(f)) if os.path.exists(f) else {}
P={1e3:'1e3',1e4:'1e4',1e5:'1e5'}
for spin in (0,0.5):
    for m in (1e3,1e4,1e5):
        for tag in ('11','14'):
            key=f"s{spin}_m{m:g}_T{tag}"
            if key in out: continue
            out[key]=combine.run(tag,m,spin,php_file=f"php/php_T{tag}_m{P[m]}.npz")
            json.dump(out,open(f,'w')); print('done',key,flush=True)
