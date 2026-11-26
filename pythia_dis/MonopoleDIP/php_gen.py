import pythia8mc, math, sys, time
import numpy as np
M=0.9382720813
def gen(T, m, n, seed, out, Wmin=5.0, q2max=1.0, verbose=False):
    E=T+m; p=math.sqrt((E-m)*(E+m)); s=M**2+m**2+2*M*E
    Lam=math.sqrt((E+M+p)/(M+m**2/(E+p)))          # CM -> target rest frame
    py=pythia8mc.Pythia("",False)
    for st in ["Beams:frameType = 1","Beams:idA = 15","Beams:idB = 2212",f"Beams:eCM = {math.sqrt(s)!r}",
      f"15:m0 = {m}","15:mWidth = 0.0","15:mayDecay = off",
      "PDF:lepton = off","PDF:lepton2gamma = on",f"Photon:Q2max = {q2max}",f"Photon:Wmin = {Wmin}",
      "SoftQCD:all = on",
      "HadronLevel:Decay = on","ParticleDecays:limitTau0 = on","ParticleDecays:tau0Max = 1.0",
      "Next:numberCount = 0","Print:quiet = on","Random:setSeed = on",f"Random:seed = {seed}"]:
        py.readString(st)
    if not py.init(): raise RuntimeError("init failed")
    Elab=[];evid=[];Y=np.zeros(n);Q2=np.zeros(n);W2=np.zeros(n); ok=np.zeros(n,bool); t0=time.time()
    for i in range(n):
        if not py.next(): continue
        ev=py.event; ok[i]=True
        lep=None; kp_had=0.0; kp_lep_in=None; kpB=0.0
        for j in range(ev.size()):
            pt=ev[j]
            if j==1 or (kp_lep_in is None and pt.id()==15 and pt.status()==-12): kp_lep_in=pt.e()+pt.pz()
            if j==2: kpB=pt.e()+pt.pz()
            if not pt.isFinal(): continue
            if pt.id()==15 and lep is None: lep=pt; continue
            kp=pt.e()+pt.pz(); kp_had+=kp
            if abs(pt.id())==211:
                km=(pt.m()**2+pt.pT()**2)/kp if kp>0 else 1e300
                Elab.append(0.5*(Lam*kp+km/Lam)); evid.append(i)
        if lep is not None and kp_lep_in:
            xi=(lep.e()+lep.pz())/kp_lep_in; y=1-xi
            Y[i]=max(kp_had-kpB,1e-30)/kp_lep_in      # y from hadronic k+ sum (accurate for tiny y)
            Q2[i]=(m*m*Y[i]**2+lep.pT()**2)/max(1-Y[i],1e-300)
    info=py.infoPython()
    np.savez(out,E=np.array(Elab),ev=np.array(evid),y=Y,Q2=Q2,ok=ok,n=n,sigma_mb=info.sigmaGen(),T=T,m=m,Wmin=Wmin,q2max=q2max)
    print(f"T={T:g} m={m:g} Wmin={Wmin:g} sigma_gen={info.sigmaGen()*1e3:.5g} ub  nok={ok.sum()}/{n}  <Npi>={len(Elab)/max(ok.sum(),1):.1f}  y[min,med,max]=({Y[ok].min():.2g},{np.median(Y[ok]):.2g},{Y[ok].max():.2g}) Q2max={Q2[ok].max():.3g}  time={time.time()-t0:.1f}s",flush=True)
if __name__=="__main__":
    gen(float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5])
