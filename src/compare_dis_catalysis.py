"""Compare Daniele's Pythia result (electromagnetic DIS + photoproduction of a heavy charged
projectile, z = 60, on a proton at rest; pythia_dis/MonopoleDIP/results/res_comb.json) with our
Callan-Rubakov catalysis kinematics (mm_kinematics.py) at the same T = E_M = 1e11 GeV and
the same masses m_M = 1, 10, 100 TeV.

Everything is expressed per event: per collision (sigma_tot-normalised) for DIS,
per catalysis for our channels.  Bins: 0.5 dex in log10(E/GeV), as in combine.py.

usage:  python src/compare_dis_catalysis.py [path/to/res_comb.json] [path/to/MonopoleDIP]
"""
import sys, json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import mm_kinematics as mk

RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'pythia_dis', 'MonopoleDIP', 'results', 'res_comb.json')
DIP = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'pythia_dis', 'MonopoleDIP')
sys.path.insert(0, DIP)
import nu as dnu                                  # Daniele's collinear pi -> mu -> e decay chain

T = 1e11; THR = 2e8; N = 200_000
MASSES = [1e3, 1e4, 1e5]
CHANS = ['p -> e+ pi+ pi-', 'p -> mu+ K0', 'p -> nubar_e pi+']
d = json.load(open(RES))
edges = np.array(d['s0.5_m1000_T11']['edges']); x = 0.5 * (edges[1:] + edges[:-1]); dl = np.diff(edges) * np.log(10)

def dis_hist(key):
    r = d[key]; S = r['sigma_ub']
    nu = (np.array(r['nu_mu']) + np.array(r['nu_e'])) * dl / S     # neutrinos per collision per bin
    pi = np.array(r['pi']) * dl / S
    return nu, pi, r

def hist(E):
    h, _ = np.histogram(np.log10(np.clip(E, 1e-3, None)), bins=edges)
    return h

def stats(E, n_ev):
    E = np.asarray(E)
    return dict(mult=E.size / n_ev, median=np.median(E) if E.size else np.nan,
                mean=E.mean() if E.size else np.nan, frac_gt=(E > THR).sum() / n_ev)

rows = []
fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True)
cat = {}
for j, m in enumerate(MASSES):
    key = f's0.5_m{m:g}_T11'
    nu_d, pi_d, r = dis_hist(key)
    # DIS stats from the binned spectrum
    c = np.cumsum(nu_d) / nu_d.sum(); med = 10 ** np.interp(0.5, c, x)
    rows.append(dict(process='EM DIS+photoprod (Pythia, z=60)', m_M=m, mult=nu_d.sum(), median=med,
                     mean=(nu_d * 10 ** x).sum() / nu_d.sum(), frac_gt=r['mean_N'], sigma_mb=r['sigma_ub'] / 1e3))
    ax = axes[0, j]; ax.step(x, nu_d, where='mid', color='k', lw=2, label='EM DIS + photoproduction (Daniele, Pythia)')
    axp = axes[1, j]; axp.step(x, pi_d, where='mid', color='k', lw=2, label='pi+- EM DIS + photoprod.')
    for ch, col in zip(CHANS, ['C0', 'C1', 'C2']):
        res = mk.catalysis(ch, T, m, n=N, seed=7)
        cat[(ch, m)] = res
        hn = hist(res['nu']['E']) / N
        ax.step(x, hn, where='mid', color=col, lw=1.6, label=f'catalysis {ch}')
        pl = [p['E'] for p in res['vertex'] if p['name'] in ('pi+', 'pi-')]
        pions = np.concatenate(pl) if pl else np.array([])
        if pions.size:
            axp.step(x, hist(pions) / N, where='mid', color=col, lw=1.6, label=f'pi+- catalysis {ch}')
        s = stats(res['nu']['E'], N); s.update(process=f'catalysis {ch}', m_M=m, sigma_mb=np.nan); rows.append(s)
    for a in (ax, axp):
        a.axvline(np.log10(THR), color='grey', ls=':'); a.set_yscale('log'); a.set_ylim(1e-6, 30); a.grid(alpha=.3)
    ax.set_title(f'$m_M$ = {m/1e3:g} TeV,  T = $E_M$ = $10^{{11}}$ GeV'); axp.set_xlabel('log10(E / GeV)')
    ax.text(np.log10(THR) + 0.1, 10, '220 PeV', rotation=90, va='top', fontsize=8)
axes[0, 0].set_ylabel('neutrinos per event per 0.5 dex'); axes[1, 0].set_ylabel('charged pions per event per 0.5 dex')
axes[0, 0].legend(fontsize=8, loc='lower left'); axes[1, 0].legend(fontsize=8, loc='lower left')
fig.suptitle('Neutrinos (top) and charged pions (bottom) per event: electromagnetic DIS of a heavy charged projectile (Pythia) vs Callan-Rubakov catalysis')
fig.tight_layout(); fig.savefig(os.path.join(ROOT, 'figures', 'compare_DIS_vs_catalysis_E1e11.png'), dpi=130)

# ---- energy-transfer fraction y = (E_M - E_M') / E_M  (catalysis: E_W / E_M) ------------------------
fig2, ax2 = plt.subplots(figsize=(8, 5))
yedges = np.arange(-13, 0.01, 0.5)
res = cat[('p -> e+ pi+ pi-', 1e3)]
ycat = res['EW_lab'] / T
ax2.hist(np.log10(ycat), bins=yedges, weights=np.full(N, 1 / N), histtype='step', lw=2, color='C0',
         label='catalysis p -> e+ pi+ pi- (any m_M <= 1e5 GeV): y = E_W / E_M')
events_dir = os.path.join(DIP, 'events')
if os.path.isdir(events_dir) and os.path.exists(os.path.join(events_dir, 'h11_0.npz')):
    sys.path.insert(0, DIP); import combine as cmb
    for m, col in zip(MASSES, ['k', 'grey', 'silver']):
        ys, ws = [], []
        for k, lo in enumerate(cmb.L11):
            e = np.load(os.path.join(events_dir, f'h11_{k}.npz')); hi = cmb.L11[k + 1] if k + 1 < len(cmb.L11) else np.inf
            keep = (e['Q2'] >= lo) & (e['Q2'] < hi); n = int(e['n'])
            w = np.where(keep, float(e['sigma_mb']) * 1e3 / n * cmb.flux_w(0.5, e['y'], e['Q2'], m), 0.) * cmb.Z2
            ys.append(e['y']); ws.append(w)
        ms = {1e3: '1e3', 1e4: '1e4', 1e5: '1e5'}[m]
        p = np.load(os.path.join(DIP, 'php', f'php_T11_m{ms}.npz')); ok = p['ok']
        ys.append(p['y'][ok]); ws.append(np.full(ok.sum(), float(p['sigma_mb']) * 1e3 / ok.sum() * cmb.Z2))
        y = np.concatenate(ys); w = np.concatenate(ws); S = w.sum()
        ax2.hist(np.log10(np.clip(y, 1e-13, 1)), bins=yedges, weights=w / S, histtype='step', lw=2, color=col,
                 label=f'EM DIS + photoprod., m = {m/1e3:g} TeV (Pythia, per collision)')
        rows.append(dict(process='EM DIS+photoprod y-median', m_M=m, mult=np.nan, median=np.nan, mean=np.nan,
                         frac_gt=np.nan, sigma_mb=S / 1e3, y_median=float(np.interp(0.5, np.cumsum(w[np.argsort(y)]) / S, np.sort(y)))))
ax2.set_xlabel('log10 of the energy fraction y transferred to the light system'); ax2.set_ylabel('fraction of events per 0.5 dex')
ax2.set_yscale('log'); ax2.set_ylim(1e-5, 1.5); ax2.grid(alpha=.3); ax2.legend(fontsize=8)
ax2.set_title('How much of $E_M$ = $10^{11}$ GeV the proton side receives')
fig2.tight_layout(); fig2.savefig(os.path.join(ROOT, 'figures', 'energy_transfer_DIS_vs_catalysis.png'), dpi=130)

# ---- decay-chain cross check: our pions through Daniele's nu.py ---------------------------------
rng = np.random.default_rng(3)
res = cat[('p -> e+ pi+ pi-', 1e3)]
pions = np.concatenate([p['E'] for p in res['vertex'] if p['name'] in ('pi+', 'pi-')])
f1, f2, f3 = dnu.decay(pions, rng)
E_dan = np.concatenate([pions * f1, pions * f2, pions * f3])
E_our = res['nu']['E']
chk = dict(n_pions=pions.size, nu_per_pion_daniele=E_dan.size / pions.size, nu_per_pion_ours=E_our.size / pions.size,
           median_daniele=np.median(E_dan), median_ours=np.median(E_our), mean_daniele=E_dan.mean(), mean_ours=E_our.mean(),
           frac_gt_daniele=(E_dan > THR).mean(), frac_gt_ours=(E_our > THR).mean())
fig3, ax3 = plt.subplots(figsize=(7, 4.5))
ax3.step(x, hist(E_our) / N, where='mid', lw=2, label='our 3-D decay chain (mm_kinematics)')
ax3.step(x, hist(E_dan) / N, where='mid', lw=2, ls='--', label="Daniele's collinear chain with mu polarisation (nu.py)")
ax3.set_yscale('log'); ax3.set_ylim(1e-6, 10); ax3.set_xlabel('log10(E_nu / GeV)'); ax3.set_ylabel('neutrinos per catalysis per 0.5 dex')
ax3.set_title('Same pions (p -> e+ pi+ pi-, E_M = 1e11, m_M = 1 TeV), two decay codes'); ax3.grid(alpha=.3); ax3.legend(fontsize=8)
fig3.tight_layout(); fig3.savefig(os.path.join(ROOT, 'figures', 'decay_chain_crosscheck.png'), dpi=130)

import csv
keys = ['process', 'm_M', 'mult', 'median', 'mean', 'frac_gt', 'sigma_mb', 'y_median']
with open(os.path.join(ROOT, 'results', 'compare_DIS_vs_catalysis_E1e11.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
    for r in rows: w.writerow({k: r.get(k, '') for k in keys})
print(f"{'process':40s} {'m_M':>7s} {'nu/event':>9s} {'median E_nu':>12s} {'mean E_nu':>10s} {'N(>220PeV)/event':>17s}")
for r in rows:
    if 'y_median' in r: print(f"{r['process']:40s} {r['m_M']:7.0f}  y_median = {r['y_median']:.3g}"); continue
    print(f"{r['process']:40s} {r['m_M']:7.0f} {r['mult']:9.2f} {r['median']:12.3g} {r['mean']:10.3g} {r['frac_gt']:17.3g}")
print('\ndecay-chain cross check:'); [print(f'  {k:22s} {v:.4g}') for k, v in chk.items()]
