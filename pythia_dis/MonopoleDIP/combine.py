"""DIS (Q2>1 GeV^2, Pythia DIS slices + analytic mass/spin flux weight) + photoproduction (Q2<1 GeV^2, Pythia lepton2gamma flux with the real heavy mass).
Returns weighted E dsigma/dE spectra, sigma and neutrino-tail statistics, with the DIS and photoproduction parts separately."""
import numpy as np, json, sys, glob
import nu
THR=2e8; Z2=3600.0
L11=[1,4,1e2,1e4,1e6,1e8,1e10]; L14=L11+[1e12]
def flux_w(spin,y,Q2,m):
    Q2min=m*m*y*y/np.maximum(1-y,1e-300)
    core=(1-y)*(1-Q2min/Q2)+(y*y/2 if spin==0.5 else 0.0)
    return np.where(Q2>=Q2min,np.maximum(core,0),0)/((1-y)+y*y/2)
class Acc:
    def __init__(s,nb): s.h=dict(pi=np.zeros(nb),mu=np.zeros(nb),e=np.zeros(nb)); s.sig=0.;s.sumw2=0.;s.Wn1=0.;s.Wn=0.;s.Wmu=0.;s.We=0.;s.Wnpi=0.;s.e1=0.;s.en=0.
def add(acc,E,ev,w,n,edges,R,rng):
    E=np.where(E>1e30,0.5,E); we=w[ev]
    acc.sig+=w.sum(); acc.sumw2+=(w**2).sum(); acc.Wnpi+=we.sum()
    acc.h['pi']+=np.histogram(np.log10(E),bins=edges,weights=we)[0]
    ind=np.zeros(n); cnt=np.zeros(n)
    for k in range(R):
        f1,f2,f3=nu.decay(E,rng); e1,e2,e3=E*f1,E*f2,E*f3
        for arr,key in ((e1,'mu'),(e2,'mu'),(e3,'e')):
            msk=arr>0; acc.h[key]+=np.histogram(np.log10(arr[msk]),bins=edges,weights=(we/R)[msk])[0]
        c=(e1>THR).astype(int)+(e2>THR)+(e3>THR); per=np.bincount(ev,weights=c,minlength=n)
        cm=np.bincount(ev,weights=(e1>THR).astype(int)+(e2>THR),minlength=n); ce=np.bincount(ev,weights=(e3>THR).astype(int),minlength=n)
        acc.Wn+=(w*per).sum()/R; acc.Wn1+=(w*(per>=1)).sum()/R; acc.Wmu+=(w*cm).sum()/R; acc.We+=(w*ce).sum()/R
        ind+=(per>=1)/R; cnt+=per/R
    acc.e1+=((w*ind)**2).sum(); acc.en+=((w*cnt)**2).sum()
def run(tag,m,spin,R=6,seed=21,php_file=None,with_dis=True,with_php=True):
    rng=np.random.default_rng(seed); T=1e11 if tag=='11' else 1e14
    edges=np.arange(0,np.log10(T)+0.001,0.5); nb=len(edges)-1
    L=L11 if tag=='11' else L14
    A_dis=Acc(nb); A_php=Acc(nb)
    if with_dis:
        for k,lo in enumerate(L):
            d=np.load(f"events/h{tag}_{k}.npz"); hi=L[k+1] if k+1<len(L) else np.inf; n=int(d['n'])
            keep=(d['Q2']>=lo)&(d['Q2']<hi)
            w=np.where(keep,float(d['sigma_mb'])*1e3/n*flux_w(spin,d['y'],d['Q2'],m),0.)*Z2
            add(A_dis,d['E'],d['ev'],w,n,edges,R,rng)
    if with_php:
        d=np.load(php_file); ok=d['ok']; nok=int(ok.sum()); n=int(d['n'])
        w=np.where(ok,float(d['sigma_mb'])*1e3/nok,0.)*Z2
        add(A_php,d['E'],d['ev'],w,n,edges,R,rng)
    tot=Acc(nb)
    for A in (A_dis,A_php):
        for k in tot.h: tot.h[k]=tot.h[k]+A.h[k]
        for a in ('sig','sumw2','Wn1','Wn','Wmu','We','Wnpi','e1','en'): setattr(tot,a,getattr(tot,a)+getattr(A,a))
    dl=np.diff(edges)*np.log(10); S=tot.sig
    out=dict(T=T,m=m,spin=spin,edges=edges.tolist(),sigma_ub=S,sigma_dis_ub=A_dis.sig,sigma_php_ub=A_php.sig,npi=tot.Wnpi/S,
        P_ge1=tot.Wn1/S,P_err=float(np.sqrt(tot.e1)/S),mean_N=tot.Wn/S,N_err=float(np.sqrt(tot.en)/S),sigP_ub=tot.Wn1,sigP_dis_ub=A_dis.Wn1,sigP_php_ub=A_php.Wn1,
        sigP_err_ub=float(np.sqrt(tot.e1)),mean_mu=tot.Wmu/S,mean_e=tot.We/S,
        pi=(tot.h['pi']/dl).tolist(),nu_mu=(tot.h['mu']/dl).tolist(),nu_e=(tot.h['e']/dl).tolist(),
        pi_php=(A_php.h['pi']/dl).tolist(),pi_dis=(A_dis.h['pi']/dl).tolist())
    return out
if __name__=="__main__":
    out={}
    for spin in (0,0.5):
        for m in (1e3,1e4,1e5):
            for tag in ('11','14'):
                key=f"s{spin}_m{m:g}_T{tag}"
                out[key]=run(tag,m,spin,php_file=f"php/php_T{tag}_m{ {1e3:'1e3',1e4:'1e4',1e5:'1e5'}[m] }.npz")
                v=out[key]; print(key,{k:(float('%.4g'%x) if isinstance(x,float) else x) for k,x in v.items() if not isinstance(x,list)},flush=True)
                json.dump(out,open('res_comb.json','w'))
