import pythia8mc, math, sys, time
import numpy as np
M=0.9382720813
def gen(T, lo, n, seed, out, m=1e5):
    E=T+m                       # total energy of the heavy projectile
    S=2*M*E                     # s-M^2 with the lepton-mass term removed (massless-equivalent run)
    s=M**2+S
    Lam=math.sqrt((2*E+M)/M)    # CM(massless) -> target rest frame, e^eta
    py=pythia8mc.Pythia("",False)
    for st in ["Beams:frameType = 1","Beams:idA = 13","Beams:idB = 2212",
      f"Beams:eCM = {math.sqrt(s)!r}",
      "WeakBosonExchange:ff2ff(t:gmZ) = on","WeakZ0:gmZmode = 1",
      f"PhaseSpace:Q2Min = {lo}","PDF:lepton = off",
      "HadronLevel:Decay = on","ParticleDecays:limitTau0 = on","ParticleDecays:tau0Max = 1.0",
      "SpaceShower:dipoleRecoil = on","SpaceShower:QEDshowerByL = off","TimeShower:QEDshowerByL = off","Next:numberCount = 0","Print:quiet=on","Random:setSeed = on",f"Random:seed = {seed}"]:
        py.readString(st)
    assert py.init()
    Elab=[];evid=[];Q2=np.zeros(n);X=np.zeros(n);Y=np.zeros(n); t0=time.time()
    for i in range(n):
        if not py.next(): continue
        info=py.infoPython()
        pr=py.process
        l_in=pr[3]; q_in=pr[4]
        l_out=pr[5] if abs(pr[5].id())==13 else pr[6]
        pp=l_in.e()+l_in.pz()
        if l_out.pz()>0: pm=l_out.pT()**2/(l_out.e()+l_out.pz())
        else: pm=l_out.e()-l_out.pz()
        Q2[i]=pp*pm
        X[i]=(q_in.e()-q_in.pz())/(pr[2].e()-pr[2].pz())
        Y[i]=Q2[i]/(X[i]*S)
        ev=py.event
        for j in range(ev.size()):
            pt=ev[j]
            if pt.isFinal() and abs(pt.id())==211:
                kp=pt.e()+pt.pz(); km=(pt.m()**2+pt.pT()**2)/kp if kp>0 else 1e300
                Elab.append(0.5*(Lam*kp+km/Lam)); evid.append(i)
    info=py.infoPython()
    np.savez(out,E=np.array(Elab),ev=np.array(evid),Q2=Q2,x=X,y=Y,n=n,sigma_mb=info.sigmaGen(),lo=lo,T=T,m=m,S=S)
    print(f"T={T:g} lo={lo:g} sigma_gen={info.sigmaGen()*1e3:.4g} ub  ymax={Y.max():.3g}  time={time.time()-t0:.0f}s",flush=True)
if __name__=="__main__":
    T=float(sys.argv[1]); n=int(sys.argv[2]); tag=sys.argv[3]
    los=[float(v) for v in sys.argv[4:]]
    for k,lo in enumerate(los): gen(T,lo,n,100+k,f"h{tag}_{k}.npz")
