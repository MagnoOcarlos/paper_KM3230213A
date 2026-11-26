"""Callan-Rubakov catalysis with our vertex and Pythia 8 for hadronization + decays.

Pythia has no catalysis process.  We therefore
  1. generate the vertex  M + t -> M + W  ourselves (same kinematics as mm_kinematics.py:
     exact s, gamma_cm, isotropic theta_cm for q_J = 1/2, m_W = m_l + m_had,min + Q, Q ~ U(0, Q_max));
  2. in the W rest frame split W -> lepton + X, where X is the quark-level hadronic system
     (d dbar for p -> e+ X,  d sbar for p -> mu+ X,  u dbar for p -> nubar_e X,
      (uu)bar antidiquark + dbar string for e- -> X), with m_X ~ U(m_X,min, m_W - m_l);
  3. hand lepton + quarks to Pythia in the W rest frame (ProcessLevel off = "standalone hadronization"):
     Lund string fragmentation, resonance decays, and pi/K/mu decays in flight;
  4. boost every final particle to the lab with the exact W 4-velocity (cancellation-free, as in
     mm_kinematics.lab_energy), and collect the neutrinos.

Pythia works at GeV energies in the W rest frame; the 1e11 GeV boost is done by us in a precision-safe way.

usage: python pythia_catalysis.py            (needs the venv with pythia8mc 8.318)
"""
import sys, os, time, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
import mm_kinematics as mk
import pythia8mc

THR = 2e8
NU_NAME = {12: 'nu_e', -12: 'nubar_e', 14: 'nu_mu', -14: 'nubar_mu', 16: 'nu_tau', -16: 'nubar_tau'}

# inclusive quark-level channels.  m_had_min: lightest hadronic final state (fixes the threshold as in the
# exclusive channels of mm_kinematics); quarks: Pythia ids of the string ends; junction: antibaryon system.
CHANNELS = {
    'p -> e+ X(d dbar)':        dict(target='p',  lepton=-11, quarks=[1, -1],      m_had_min=mk.MASS['pi0'], junction=False,
                                     ours=['p -> e+ pi0', 'p -> e+ pi+ pi-']),
    'p -> mu+ X(d sbar)':       dict(target='p',  lepton=-13, quarks=[1, -3],      m_had_min=mk.MASS['K0'],  junction=False,
                                     ours=['p -> mu+ K0']),
    'p -> nubar_e X(u dbar)':   dict(target='p',  lepton=-12, quarks=[2, -1],      m_had_min=mk.MASS['pi+'], junction=False,
                                     ours=['p -> nubar_e pi+']),
    # antibaryon system as an antidiquark-antiquark string (Pythia's own representation of baryon strings);
    # the three-antiquark junction topology does not hadronize in standalone mode (tested, silently ignored)
    'e- -> X(ubar ubar dbar)':  dict(target='e-', lepton=None, quarks=[-2203, -1],   m_had_min=mk.MASS['p'],   junction=False,
                                     ours=['e- -> pbar pi+ pi-', 'e- -> pbar K0']),
}

def make_pythia(m_M, seed):
    py = pythia8mc.Pythia("", False)
    for s in ["ProcessLevel:all = off", "Print:quiet = on", "Next:numberCount = 0",
              "HadronLevel:Decay = on", "ParticleDecays:limitTau0 = off",          # everything decays (ISM case)
              "211:mayDecay = on", "321:mayDecay = on", "130:mayDecay = on", "310:mayDecay = on", "13:mayDecay = on",
              "2112:mayDecay = off", "3122:mayDecay = on",
              f"4110000:new = MM MMbar 2 0 0 {m_M}", "4110000:mayDecay = off",
              "Random:setSeed = on", f"Random:seed = {seed}"]:
        py.readString(s)
    assert py.init()
    return py

def unit_vectors(n, rng):
    c = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n); s = np.sqrt(1 - c ** 2)
    return np.stack([s * np.cos(ph), s * np.sin(ph), c], axis=1)

def boost_vec(p4, bvec):
    """ordinary boost of 4-vectors p4 (N,4: E,px,py,pz) by velocity bvec (N,3); fine at GeV energies"""
    b2 = np.sum(bvec ** 2, axis=1); g = 1 / np.sqrt(1 - b2)
    bp = np.sum(p4[:, 1:] * bvec, axis=1)
    g2 = np.where(b2 > 0, (g - 1) / np.maximum(b2, 1e-300), 0.0)
    E = g * (p4[:, 0] + bp)
    p = p4[:, 1:] + (g2 * bp + g * p4[:, 0])[:, None] * bvec
    return np.concatenate([E[:, None], p], axis=1)

def run(channel, E_kin, m_M, n=50_000, Q_max=1.0, seed=11):
    ch = CHANNELS[channel]; rng = np.random.default_rng(seed)
    py = make_pythia(m_M, seed); ev = py.event; pd = py.particleData
    m_t = mk.MASS[ch['target']]
    m_l = 0.0 if ch['lepton'] is None else pd.m0(abs(ch['lepton']))
    mq = np.array([pd.m0(abs(q)) for q in ch['quarks']])                 # Pythia constituent masses
    m_X_min = max(ch['m_had_min'], mq.sum() * 1.001)
    fr = mk.cm_frame(E_kin, m_M, m_t)
    mW_max = fr['s_minus_mM2'] / (fr['sqrt_s'] + m_M)
    if mW_max <= m_l + m_X_min: return None
    # threshold of the inclusive channel = lepton + lightest string (Pythia constituent masses), slightly above
    # the physical hadron threshold m_had_min; Q is the excess above it, as in mm_kinematics
    Q_eff = min(Q_max, 0.999 * (mW_max - m_l - m_X_min))
    Q = rng.uniform(0, Q_eff, n); mW = m_l + m_X_min + Q
    # ---- vertex in the CM frame, boosted to the lab (as in mm_kinematics.catalysis) -----------------
    pcm, EW, EM = mk.vertex_cm(fr, m_M, mW)
    cth = rng.uniform(-1, 1, n); pzW = -pcm * cth; pT = pcm * np.sqrt(1 - cth ** 2)
    EW_lab = mk.lab_energy(EW, pzW, mW ** 2, pT ** 2, fr['beta_cm'], fr['gamma_cm'], fr['one_minus_beta_cm'])
    pzW_lab = fr['gamma_cm'] * (pzW + fr['beta_cm'] * EW)
    pW_lab = np.sqrt(np.maximum(EW_lab ** 2 - mW ** 2, pT ** 2))          # |p_W| in the lab
    nW = np.stack([pT, np.zeros(n), pzW_lab], axis=1) / pW_lab[:, None]  # W direction in the lab (phi = 0)
    bW, gW, ombW = mk.parent_boost(EW_lab, mW)                           # exact beta, gamma, 1-beta of W
    # ---- W rest frame: W -> lepton + X ---------------------------------------------------------------
    if ch['lepton'] is None:                                              # no lepton: X is the whole W
        mX = mW.copy(); pl = np.zeros(n)
    else:
        mX = rng.uniform(m_X_min, np.maximum(mW - m_l, m_X_min * 1.0001))
        mX = np.minimum(mX, mW - m_l - 1e-6)
        pl = mk.pstar(mW, m_l, mX)
    u = unit_vectors(n, rng)
    El = np.sqrt(pl ** 2 + m_l ** 2); EX = mW - El
    pX = -pl[:, None] * u
    bX = pX / EX[:, None]                                                 # velocity of X in the W frame
    # ---- event loop -----------------------------------------------------------------------------------
    out_E, out_pz, out_pT2, out_m2, out_id, out_evt = [], [], [], [], [], []
    n_prim_pi, n_prim_had, n_fail = np.zeros(n), np.zeros(n), 0
    t0 = time.time()
    for i in range(n):
        ev.reset()
        if ch['lepton'] is not None:
            ev.append(ch['lepton'], 1, 0, 0, pl[i] * u[i, 0], pl[i] * u[i, 1], pl[i] * u[i, 2], El[i], m_l)
        if not ch['junction']:
            # q qbar back to back in the X rest frame, then boosted to the W frame
            pq = mk.pstar(mX[i], mq[0], mq[1]); v = unit_vectors(1, rng)[0]
            q4 = np.array([[np.hypot(pq, mq[0]), *(pq * v)], [np.hypot(pq, mq[1]), *(-pq * v)]])
            q4 = boost_vec(q4, np.repeat(bX[i][None], 2, axis=0))
            ev.append(ch['quarks'][0], 23, 101, 0, q4[0, 1], q4[0, 2], q4[0, 3], q4[0, 0], mq[0])
            ev.append(ch['quarks'][1], 23, 0, 101, q4[1, 1], q4[1, 2], q4[1, 3], q4[1, 0], mq[1])
        else:
            vecs = mk.three_body_rest(mX[i], mq[0], mq[1], mq[2], np.random.default_rng(rng.integers(1 << 31)))
            q4 = np.array([[np.sqrt(mq[k] ** 2 + np.sum(vecs[k][0] ** 2)), *vecs[k][0]] for k in range(3)])
            q4 = boost_vec(q4, np.repeat(bX[i][None], 3, axis=0))
            for k, tag in enumerate((101, 102, 103)):
                ev.append(ch['quarks'][k], 23, 0, tag, q4[k, 1], q4[k, 2], q4[k, 3], q4[k, 0], mq[k])
            ev.appendJunction(2, 101, 102, 103)
        if not py.next():
            n_fail += 1; continue
        for j in range(ev.size()):
            pt = ev[j]; st = pt.statusAbs()
            if 81 <= st <= 89:                                            # primary hadrons from the string
                n_prim_had[i] += 1; n_prim_pi[i] += (abs(pt.id()) == 211)
            if pt.isFinal():
                px, pyy, pz, e = pt.px(), pt.py(), pt.pz(), pt.e()
                pn = px * nW[i, 0] + pyy * nW[i, 1] + pz * nW[i, 2]      # momentum along the W direction
                out_E.append(e); out_pz.append(pn); out_pT2.append(max(px * px + pyy * pyy + pz * pz - pn * pn, 0.0))
                out_m2.append(pt.m() ** 2); out_id.append(pt.id()); out_evt.append(i)
    dt = time.time() - t0
    E = np.array(out_E); pz = np.array(out_pz); pT2 = np.array(out_pT2); m2 = np.array(out_m2)
    ids = np.array(out_id); evt = np.array(out_evt)
    E_lab = mk.lab_energy(E, pz, m2, pT2, bW[evt], gW[evt], ombW[evt])
    isnu = np.isin(np.abs(ids), (12, 14, 16))
    n_ok = n - n_fail
    res = dict(channel=channel, E_kin=E_kin, m_M=m_M, n=n, n_ok=n_ok, n_fail=n_fail, Q_max=Q_max, time_s=dt,
               sqrt_s=fr['sqrt_s'], gamma_cm=fr['gamma_cm'],
               EW_lab_mean=float(EW_lab.mean()), mX_mean=float(mX.mean()),
               prim_had=float(n_prim_had.sum() / n_ok), prim_pi=float(n_prim_pi.sum() / n_ok),
               nu=dict(id=ids[isnu], E=E_lab[isnu], evt=evt[isnu]),
               final=dict(id=ids, E=E_lab, evt=evt))
    Enu = E_lab[isnu]
    res.update(nu_mult=Enu.size / n_ok, nu_median=float(np.median(Enu)) if Enu.size else np.nan,
               nu_mean=float(Enu.mean()) if Enu.size else np.nan,
               nu_p10=float(np.percentile(Enu, 10)) if Enu.size else np.nan,
               nu_p90=float(np.percentile(Enu, 90)) if Enu.size else np.nan,
               nu_gt=float((Enu > THR).sum() / n_ok))
    # energy check: sum of final lab energies vs E_W_lab (per event)
    sumE = np.bincount(evt, weights=E_lab, minlength=n)
    ok = np.ones(n, bool); ok[:] = True
    if n_fail: ok = np.bincount(evt, minlength=n) > 0
    res['energy_check'] = float(np.max(np.abs(sumE[ok] / EW_lab[ok] - 1)))
    return res

def flavour_table(res):
    out = {}
    ids = res['nu']['id']; E = res['nu']['E']
    for pid, name in NU_NAME.items():
        m = ids == pid
        if m.any(): out[name] = dict(mult=m.sum() / res['n_ok'], median=float(np.median(E[m])), gt=float((E[m] > THR).sum() / res['n_ok']))
    return out

if __name__ == '__main__':
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    E_kin = 1e11; N = 50_000
    edges = np.arange(3, 11.26, 0.25); x = 0.5 * (edges[1:] + edges[:-1])
    def h(E, nev): return np.histogram(np.log10(np.clip(E, 1e-3, None)), bins=edges)[0] / nev
    runs = []
    plan = [(c, 1e3, q) for c in list(CHANNELS)[:3] for q in (1.0, 3.0, 10.0)] + \
           [('e- -> X(ubar ubar dbar)', m, q) for m in (1e3, 1e4, 1e5) for q in (1.0, 3.0)]
    results = {}
    for c, m, q in plan:
        r = run(c, E_kin, m, n=N, Q_max=q)
        if r is None: print(c, m, q, 'closed'); continue
        results[(c, m, q)] = r
        print(f"{c:26s} m_M={m:7.0f} Q_max={q:4.1f}  fail={r['n_fail']:5d}  prim.had={r['prim_had']:.2f} prim.pi+-={r['prim_pi']:.2f}  "
              f"nu/ev={r['nu_mult']:.2f}  median={r['nu_median']:.3g}  N(>220PeV)/ev={r['nu_gt']:.2f}  "
              f"<m_X>={r['mX_mean']:.2f}  dE/E={r['energy_check']:.1e}  {r['time_s']:.0f}s", flush=True)
        runs.append({k: v for k, v in r.items() if k not in ('nu', 'final')} | dict(flavours=flavour_table(r)))
    json.dump(runs, open(os.path.join(ROOT, 'results', 'pythia_catalysis_E1e11.json'), 'w'), indent=1, default=float)
    # ---- figure: Pythia inclusive channels vs our exclusive phase-space channels ----------------------
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)
    for ax, (c, m) in zip(axes.flat, [(c, 1e3) for c in list(CHANNELS)[:3]] + [('e- -> X(ubar ubar dbar)', 1e3)]):
        for q, ls in ((1.0, '-'), (3.0, '--'), (10.0, ':')):
            r = results.get((c, m, q))
            if r: ax.step(x, h(r['nu']['E'], r['n_ok']), where='mid', color='k', ls=ls, lw=1.8, label=f'Pythia string, Q_max = {q:g} GeV')
        for oc, col in zip(CHANNELS[c]['ours'], ('C0', 'C1')):
            o = mk.catalysis(oc, E_kin, m, n=N, seed=3)
            if o is not None and o['nu']['E'].size:
                ax.step(x, h(o['nu']['E'], N), where='mid', color=col, lw=1.6, label=f'phase space {oc} (Q_max = 1 GeV)')
            elif o is not None: ax.plot([], [], color=col, label=f'phase space {oc}: no neutrinos')
        ax.axvline(np.log10(THR), color='grey', ls=':'); ax.set_yscale('log'); ax.set_ylim(1e-5, 10); ax.grid(alpha=.3)
        ax.set_title(f'{c},  m_M = {m/1e3:g} TeV'); ax.legend(fontsize=7.5, loc='upper left')
    for ax in axes[1]: ax.set_xlabel('log10(E_nu / GeV)')
    for ax in axes[:, 0]: ax.set_ylabel('neutrinos per catalysis per 0.25 dex')
    fig.suptitle('Our catalysis vertex + Pythia hadronization/decays (black) vs our exclusive phase-space model (colour), E_M = 1e11 GeV, all flavours')
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, 'figures', 'pythia_catalysis_vs_phasespace.png'), dpi=130)
    # mass dependence of the electron channel
    fig2, ax2 = plt.subplots(figsize=(7.5, 5))
    for m, col in ((1e3, 'C0'), (1e4, 'C2'), (1e5, 'C3')):
        r = results.get(('e- -> X(ubar ubar dbar)', m, 1.0))
        if r: ax2.step(x, h(r['nu']['E'], r['n_ok']), where='mid', color=col, lw=1.8, label=f'Pythia, m_M = {m/1e3:g} TeV')
        o = mk.catalysis('e- -> pbar pi+ pi-', E_kin, m, n=N, seed=3)
        if o is not None: ax2.step(x, h(o['nu']['E'], N), where='mid', color=col, ls='--', lw=1.2, label=f'phase space pbar pi+ pi-, m_M = {m/1e3:g} TeV')
    ax2.axvline(np.log10(THR), color='grey', ls=':'); ax2.set_yscale('log'); ax2.set_ylim(1e-5, 10); ax2.grid(alpha=.3)
    ax2.set_xlabel('log10(E_nu / GeV)'); ax2.set_ylabel('neutrinos per catalysis per 0.25 dex'); ax2.legend(fontsize=8)
    ax2.set_title('e- -> M + (ubar ubar dbar): Pythia (uu)bar-dbar string vs phase space, E_M = 1e11 GeV')
    fig2.tight_layout(); fig2.savefig(os.path.join(ROOT, 'figures', 'pythia_catalysis_electron_mass.png'), dpi=130)
    print('done')
