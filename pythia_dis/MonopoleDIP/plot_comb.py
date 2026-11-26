import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=json.load(open('results/res_comb.json'))
for spin,name in ((0,'scalar'),(0.5,'spin-1/2')):
    fig,ax=plt.subplots(2,2,figsize=(11,8))
    for j,T in enumerate(('11','14')):
        for m,c in ((1000,'C0'),(10000,'C1'),(100000,'C2')):
            r=d.get(f's{spin}_m{m}_T{T}')
            if not r: continue
            x=0.5*(np.array(r['edges'][1:])+np.array(r['edges'][:-1])); lab=f'm={m/1e3:g} TeV'
            ax[0,j].semilogy(x,np.array(r['pi'])*1e-3,c,label=lab+' total'); ax[0,j].semilogy(x,np.array(r['pi_php'])*1e-3,c+'--',lw=1)
            ax[1,j].semilogy(x,(np.array(r['nu_mu'])+np.array(r['nu_e']))*1e-3,c,label=lab)
        ax[0,j].set_title(f'{name}: pi+- spectrum, T=1e{T} GeV (dashed: photoproduction)'); ax[1,j].set_title(f'{name}: neutrinos (all flavours), T=1e{T} GeV')
        ax[1,j].axvline(np.log10(2e8),color='k',ls=':')
        for i in (0,1): ax[i,j].set_xlabel('log10(E/GeV)');ax[i,j].set_ylabel('E dsigma/dE [mb]');ax[i,j].legend()
        ax[0,j].set_ylim(bottom=1e-2)
    plt.tight_layout();plt.savefig('plots/spectra_combined_%s.png'%('scalar' if spin==0 else 'spin_half'),dpi=120);plt.close()
