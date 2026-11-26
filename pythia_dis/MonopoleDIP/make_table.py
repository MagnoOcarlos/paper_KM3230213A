import json
d=json.load(open('results/res_comb.json'))
L=['spin | m [TeV] | T [GeV] | sigma_tot [mb] | sigma_DIS | sigma_php | <Npi> | <N_nu>(>200PeV)/collision | P(>=1 nu) | sigma_nu [mb] | sigma(>=1nu) [mb] | N per DIS collision','-|-|-|-|-|-|-|-|-|-|-|-']
for k,r in d.items():
    mb=1e-3; sn=r['mean_N']*r['sigma_ub']*mb
    L.append(f"{r['spin']} | {r['m']/1e3:g} | {r['T']:.0e} | {r['sigma_ub']*mb:.4g} | {r['sigma_dis_ub']*mb:.4g} | {r['sigma_php_ub']*mb:.4g} | {r['npi']:.3g} | {r['mean_N']:.3g} +- {r['N_err']:.2g} | {r['P_ge1']:.3g} +- {r['P_err']:.2g} | {sn:.4g} | {r['sigP_ub']*mb:.4g} | {sn/(r['sigma_dis_ub']*mb):.3g}")
open('results/summary_table.md','w').write('\n'.join(L)+'\n');print('\n'.join(L))
