# Regenerate Daniele's event set (run_generation.sh) in parallel.
import sys, time, multiprocessing as mp
L11=[1,4,1e2,1e4,1e6,1e8,1e10]; L14=L11+[1e12]
def job(a):
    kind=a[0]; t0=time.time()
    if kind=='dis':
        import gen2; _,T,k,lo,n,tag=a
        gen2.gen(T,lo,n,100+k,f"events/h{tag}_{k}.npz")
    else:
        import php_gen; _,T,m,n,seed,out=a
        php_gen.gen(T,m,n,seed,out)
    return a, time.time()-t0
jobs=[]
for tag,T,L in (('11',1e11,L11),('14',1e14,L14)):
    for k,lo in enumerate(L):
        jobs.append(('dis',T,k,lo,30000 if k==0 else 25000,tag))
for m,ms in ((1e3,'1e3'),(1e4,'1e4'),(1e5,'1e5')):
    jobs.append(('php',1e14,m,100000,int(m/1e2)+14,f"php/php_T14_m{ms}.npz"))
    jobs.append(('php',1e11,m,100000,int(m/1e2)+11,f"php/php_T11_m{ms}.npz"))
if __name__=='__main__':
    t0=time.time()
    with mp.Pool(10) as pool:
        for a,dt in pool.imap_unordered(job,jobs):
            print(f"done {a[:3]} in {dt:.0f}s (elapsed {time.time()-t0:.0f}s)",flush=True)
    print("ALL DONE",flush=True)
