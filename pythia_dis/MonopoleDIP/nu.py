import numpy as np
mpi=0.13957039; mmu=0.1056583755; r=(mmu/mpi)**2
def rej(kind,P,rng):
    N=len(P); x=np.empty(N); c=np.empty(N); todo=np.arange(N)
    while len(todo):
        k=len(todo); xs=rng.random(k); cs=rng.uniform(-1,1,k); u=rng.random(k)*2.0
        p=P[todo]
        if kind=='mu':   # partner of mu (nubar_mu from mu+): x^2[(3-2x)+P c(2x-1)]
            d=xs**2*((3-2*xs)+p*cs*(2*xs-1))
        else:            # nu_e from mu+: 6x^2(1-x)(1-Pc)
            d=6*xs**2*(1-xs)*(1-p*cs)
        ok=u<d
        x[todo[ok]]=xs[ok]; c[todo[ok]]=cs[ok]; todo=todo[~ok]
    return x,c
def decay(E,rng):
    """E: pion energies. returns fractions of E_pi for (nu_mu direct, nu_mu-type from mu, nu_e)"""
    N=len(E)
    xm=r+(1-r)*rng.random(N)
    P=(2*r/xm-(1+r))/(1-r)
    P=np.clip(P,-1,1)
    f1=1-xm
    x2,c2=rej('mu',P,rng); f2=xm*x2*(1+c2)/2
    x3,c3=rej('e',P,rng);  f3=xm*x3*(1+c3)/2
    return f1,f2,f3
