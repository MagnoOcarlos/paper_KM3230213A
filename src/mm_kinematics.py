"""
Kinematics of neutrino production by a relativistic magnetic monopole (MM)
through Callan-Rubakov catalysis.

Two stages, both done exactly with Lorentz invariants (no "products at rest in
the MM frame" assumption):

  Stage 1  (vertex)   M + t  ->  M + W ,  t = p or e- at rest in the lab,
                       W = the light final-state system (e+ pi0, pbar pi+ pi-, ...)
           - s = m_M^2 + m_t^2 + 2 E_M m_t          (E_M = total lab energy)
           - MM scattered isotropically in the CM frame (lowest partial wave,
             q_J = 1/2: dsigma/dOmega independent of theta_cm, Candela et al.
             arXiv:2504.14918 Eq. 2.3 / Khoze arXiv:2405.18689 Eq. 6.5, 6.12)
           - W decays isotropically (flat phase space) into its constituents
  Stage 2  (chains)   pi, K, mu decay in flight (ISM case: nothing is absorbed)
                       until only stable particles remain; neutrinos are collected
                       with their flavour, lab energy and "generation"
                       (0 = produced directly at the catalysis vertex).

Only lab-frame ENERGIES are tracked after each decay.  Because every decay is
isotropic in the parent rest frame, the daughter lab energy depends only on the
parent lab energy and the daughter's rest-frame momentum component along the
parent direction, so directions need not be propagated.

All energies in GeV.
"""
import numpy as np
from scipy.spatial.transform import Rotation

# ----------------------------------------------------------------------------
# particle data (PDG 2024)
# ----------------------------------------------------------------------------
MASS = {
    'p': 0.938272, 'pbar': 0.938272, 'n': 0.939565,
    'e-': 0.000511, 'e+': 0.000511,
    'mu-': 0.105658, 'mu+': 0.105658,
    'pi+': 0.139570, 'pi-': 0.139570, 'pi0': 0.134977,
    'K+': 0.493677, 'K-': 0.493677, 'K0': 0.497611, 'K_S': 0.497611, 'K_L': 0.497611,
    'gamma': 0.0,
    'nu_e': 0.0, 'nubar_e': 0.0, 'nu_mu': 0.0, 'nubar_mu': 0.0,
}
NEUTRINOS = ('nu_e', 'nubar_e', 'nu_mu', 'nubar_mu')
STABLE = {'p', 'pbar', 'n', 'e-', 'e+', 'gamma'} | set(NEUTRINOS)

# decay modes: species -> list of (branching ratio, daughters)
# 'michel' modes are handled with the V-A matrix element, the rest with flat phase space
DECAYS = {
    'pi+': [(1.0, ['mu+', 'nu_mu'])],
    'pi-': [(1.0, ['mu-', 'nubar_mu'])],
    'pi0': [(1.0, ['gamma', 'gamma'])],
    'mu+': [(1.0, ['e+', 'nu_e', 'nubar_mu'])],      # Michel
    'mu-': [(1.0, ['e-', 'nubar_e', 'nu_mu'])],      # Michel
    'K+':  [(0.6356, ['mu+', 'nu_mu']), (0.2067, ['pi+', 'pi0']),
            (0.0558, ['pi+', 'pi+', 'pi-']), (0.0507, ['pi0', 'e+', 'nu_e']),
            (0.0335, ['pi0', 'mu+', 'nu_mu']), (0.0176, ['pi+', 'pi0', 'pi0'])],
    'K-':  [(0.6356, ['mu-', 'nubar_mu']), (0.2067, ['pi-', 'pi0']),
            (0.0558, ['pi-', 'pi-', 'pi+']), (0.0507, ['pi0', 'e-', 'nubar_e']),
            (0.0335, ['pi0', 'mu-', 'nubar_mu']), (0.0176, ['pi-', 'pi0', 'pi0'])],
    'K0':  [(0.5, ['K_S']), (0.5, ['K_L'])],         # strangeness eigenstate -> mass eigenstates
    'K_S': [(0.6920, ['pi+', 'pi-']), (0.3069, ['pi0', 'pi0'])],
    'K_L': [(0.20275, ['pi+', 'e-', 'nubar_e']), (0.20275, ['pi-', 'e+', 'nu_e']),
            (0.1352, ['pi+', 'mu-', 'nubar_mu']), (0.1352, ['pi-', 'mu+', 'nu_mu']),
            (0.1952, ['pi0', 'pi0', 'pi0']), (0.1254, ['pi+', 'pi-', 'pi0'])],
}
MICHEL = {'mu+': 'nu_e', 'mu-': 'nubar_e'}   # the neutrino whose energy carries the (1-x) factor

# the catalysis channels requested by Daniele (Candela Eq. 2.1 / 2.2; Super-K Eq. 1;
# Hu et al. 2201.02386 Eq. 3.2).  'products' = constituents of the light system W.
CHANNELS = {
    'p -> e+ pi0':          dict(target='p',  products=['e+', 'pi0']),
    'p -> e+ pi+ pi-':      dict(target='p',  products=['e+', 'pi+', 'pi-']),   # via rho0/eta/omega -> pions
    'p -> mu+ K0':          dict(target='p',  products=['mu+', 'K0']),
    'p -> nubar_e pi+':     dict(target='p',  products=['nubar_e', 'pi+']),     # prompt neutrino, BR ~ 1e-4
    'e- -> pbar pi+ pi-':   dict(target='e-', products=['pbar', 'pi+', 'pi-']),
    'e- -> pbar K0':        dict(target='e-', products=['pbar', 'K0']),
    'e- -> pbar':           dict(target='e-', products=['pbar']),               # pure 2->2 (Candela check)
}


# ----------------------------------------------------------------------------
# relativistic helpers
# ----------------------------------------------------------------------------
def kallen(x, y, z):
    return (x - y - z) ** 2 - 4.0 * y * z


def pstar(M, m1, m2):
    """two-body momentum of m1, m2 in the rest frame of M (arrays allowed)."""
    lam = kallen(M ** 2, m1 ** 2, m2 ** 2)
    return np.sqrt(np.clip(lam, 0.0, None)) / (2.0 * M)


def cm_frame(E_kin, m_M, m_t):
    """Lab: target at rest, MM with kinetic energy E_kin.
    Returns dict with s, sqrt_s, gamma_cm, beta_cm, one_minus_beta_cm, E_M (total), p_M,
    and s_minus_mM2 = s - m_M^2 computed without cancellation (needed for GUT-scale masses)."""
    E = E_kin + m_M
    p = np.sqrt(E_kin * (E_kin + 2.0 * m_M))             # exact sqrt(E^2 - m^2)
    s_minus_mM2 = 2.0 * E * m_t + m_t ** 2                # = s - m_M^2, exact
    s = m_M ** 2 + s_minus_mM2
    rs = np.sqrt(s)
    Etot = E + m_t
    gamma = Etot / rs
    beta = p / Etot
    one_minus_beta = s / (Etot * (Etot + p))              # exact, no cancellation
    return dict(s=s, sqrt_s=rs, gamma_cm=gamma, beta_cm=beta, s_minus_mM2=s_minus_mM2,
                one_minus_beta_cm=one_minus_beta, E_M=E, p_M=p, E_tot=Etot)


def vertex_cm(fr, m_M, mW):
    """CM momentum and energies of W and M in  M + t -> M + W, cancellation-free."""
    rs, d = fr['sqrt_s'], fr['s_minus_mM2']
    lam = (d - mW ** 2) ** 2 - 4.0 * m_M ** 2 * mW ** 2     # Kallen(s, mM^2, mW^2)
    p = np.sqrt(np.clip(lam, 0.0, None)) / (2.0 * rs)
    EW = (d + mW ** 2) / (2.0 * rs)
    EM = rs - EW
    return p, EW, EM


def lab_energy(Estar, pz, m2, pT2, beta, gamma, one_minus_beta):
    """E_lab = gamma*(E* + beta*pz) computed without catastrophic cancellation
    for backward emission (pz<0) at gamma ~ 1e5-1e10."""
    Estar, pz, m2, pT2 = map(np.asarray, (Estar, pz, m2, pT2))
    beta, gamma, one_minus_beta = (np.broadcast_to(np.asarray(a, float), Estar.shape)
                                   for a in (beta, gamma, one_minus_beta))
    apz = np.abs(pz)
    diff = (m2 + pT2) / (Estar + apz)             # = E* - |pz|  (exact)
    back = gamma * (diff + one_minus_beta * apz)  # pz < 0
    fwd = gamma * (Estar + beta * pz)             # pz >= 0
    return np.where(pz < 0, back, fwd)


def parent_boost(E, m):
    """beta, gamma, 1-beta of a particle of lab energy E and mass m (arrays)."""
    E = np.asarray(E, float)
    p = np.sqrt(np.clip(E ** 2 - m ** 2, 0, None))
    gamma = E / m
    beta = p / E
    one_minus_beta = m ** 2 / (E * (E + p))
    return beta, gamma, one_minus_beta


# ----------------------------------------------------------------------------
# rest-frame decay generators  (return list of (E*, pz*, pT2) per daughter, random orientation)
# ----------------------------------------------------------------------------
def _random_dirs(n, rng):
    cos = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    sin = np.sqrt(1 - cos ** 2)
    return np.column_stack([sin * np.cos(phi), sin * np.sin(phi), cos])


def two_body_rest(M, m1, m2, rng):
    M = np.asarray(M, float)
    n = M.size
    p = pstar(M, m1, m2)
    d = _random_dirs(n, rng)
    v1 = p[:, None] * d
    v2 = -v1
    return [v1, v2]


def three_body_rest(M, m1, m2, m3, rng):
    """flat phase space, massive daughters; returns three (n,3) rest-frame momenta."""
    M = np.asarray(M, float)
    n = M.size
    lo = np.full(n, m1 + m2)
    hi = M - m3
    # envelope for rejection in m12
    grid = lo[:, None] + (hi - lo)[:, None] * np.linspace(0, 1, 64)[None, :]
    wmax = 1.05 * np.max(pstar(M[:, None], grid, m3) * pstar(grid, m1, m2), axis=1)
    m12 = np.empty(n)
    todo = np.ones(n, bool)
    while todo.any():
        k = todo.sum()
        trial = lo[todo] + (hi - lo)[todo] * rng.uniform(size=k)
        w = pstar(M[todo], trial, m3) * pstar(trial, m1, m2)
        acc = rng.uniform(size=k) * wmax[todo] < w
        idx = np.flatnonzero(todo)[acc]
        m12[idx] = trial[acc]
        todo[idx] = False
    p3 = pstar(M, m12, m3)
    d3 = _random_dirs(n, rng)
    v3 = p3[:, None] * d3
    # (12) system recoils with -p3; decay it in its own rest frame then boost (mild boost)
    E12 = np.sqrt(m12 ** 2 + p3 ** 2)
    b = -p3 / E12
    g = E12 / m12
    p12 = pstar(m12, m1, m2)
    d12 = _random_dirs(n, rng)
    q1 = p12[:, None] * d12
    q2 = -q1
    out = [v3]
    for q, m in ((q1, m1), (q2, m2)):
        Eq = np.sqrt(m ** 2 + p12 ** 2)
        qpar = np.sum(q * d3, axis=1)            # component along d3
        qperp = q - qpar[:, None] * d3
        qpar_new = g * (qpar + b * Eq)
        out.append(qperp + qpar_new[:, None] * d3)
    return [out[1], out[2], out[0]]              # order m1, m2, m3


def michel_rest(m_mu, which_nu_idx, rng, n):
    """mu -> e nu nu with the V-A matrix element.  Daughters ordered (e, nu_a, nu_b) where
    nu_a is the one with the x^2(1-x) spectrum (nu_e for mu+, nubar_e for mu-).
    All daughters treated massless.  Returns three (n,3) rest-frame momenta."""
    x = np.empty((n, 3))
    todo = np.ones(n, bool)
    while todo.any():
        k = todo.sum()
        x1 = rng.uniform(size=k); x2 = rng.uniform(size=k)
        ok = (x1 + x2) >= 1.0
        x3 = 2.0 - x1 - x2
        # x1 = electron, x2 = nu_a (x^2(1-x) type), x3 = nu_b ; weight for nu_a:
        w = x2 * (1.0 - x2) / 0.25
        acc = ok & (rng.uniform(size=k) < w)
        idx = np.flatnonzero(todo)[acc]
        x[idx, 0] = x1[acc]; x[idx, 1] = x2[acc]; x[idx, 2] = x3[acc]
        todo[idx] = False
    p = 0.5 * m_mu * x                           # magnitudes
    # planar configuration: p1 along x, p2 at angle th12, p3 = -(p1+p2)
    c12 = np.clip((p[:, 2] ** 2 - p[:, 0] ** 2 - p[:, 1] ** 2) / (2 * p[:, 0] * p[:, 1] + 1e-300), -1, 1)
    s12 = np.sqrt(1 - c12 ** 2)
    v1 = np.column_stack([p[:, 0], np.zeros(n), np.zeros(n)])
    v2 = np.column_stack([p[:, 1] * c12, p[:, 1] * s12, np.zeros(n)])
    v3 = -(v1 + v2)
    R = Rotation.random(n, random_state=rng.integers(2**31 - 1)).as_matrix()
    return [np.einsum('nij,nj->ni', R, v) for v in (v1, v2, v3)]


# ----------------------------------------------------------------------------
# decay chains in the lab (energies only)
# ----------------------------------------------------------------------------
def decay_all(particles, rng, max_gen=20, at_rest=False):
    """particles: list of dicts {name, E (array), evt (array), gen (int)}.
    Decays everything unstable; returns list of stable particles (same format)."""
    stable, queue = [], list(particles)
    while queue:
        part = queue.pop()
        name, E, evt, gen = part['name'], part['E'], part['evt'], part['gen']
        if name in STABLE or gen > max_gen:
            stable.append(part)
            continue
        modes = DECAYS[name]
        brs = np.array([b for b, _ in modes]); brs = brs / brs.sum()
        choice = rng.choice(len(modes), size=E.size, p=brs)
        m = MASS[name]
        if at_rest:                                        # dense medium: parent stops before decaying
            beta, gamma, omb = np.zeros_like(E), np.ones_like(E), np.ones_like(E)
        else:
            beta, gamma, omb = parent_boost(E, m)
        for k, (_, daughters) in enumerate(modes):
            sel = choice == k
            if not sel.any():
                continue
            Ek, evk, n = E[sel], evt[sel], int(sel.sum())
            masses = [MASS[d] for d in daughters]
            if len(daughters) == 1:                       # K0 -> K_S / K_L
                queue.append(dict(name=daughters[0], E=Ek, evt=evk, gen=gen))
                continue
            if name in MICHEL:
                vecs = michel_rest(m, 1, rng, n)
                # order daughters so that index 1 is the x^2(1-x) neutrino
                order = [daughters[0], MICHEL[name],
                         [d for d in daughters[1:] if d != MICHEL[name]][0]]
                daughters, masses = order, [MASS[d] for d in order]
            elif len(daughters) == 2:
                vecs = two_body_rest(np.full(n, m), masses[0], masses[1], rng)
            else:
                vecs = three_body_rest(np.full(n, m), *masses, rng)
            for d, md, v in zip(daughters, masses, vecs):
                Estar = np.sqrt(md ** 2 + np.sum(v ** 2, axis=1))
                pz = v[:, 2]
                pT2 = v[:, 0] ** 2 + v[:, 1] ** 2
                Elab = lab_energy(Estar, pz, md ** 2, pT2, beta[sel], gamma[sel], omb[sel])
                queue.append(dict(name=d, E=Elab, evt=evk, gen=gen + 1))
    return stable


# ----------------------------------------------------------------------------
# Stage 1: the catalysis vertex
# ----------------------------------------------------------------------------
def catalysis(channel, E_kin, m_M, n=100_000, Q_max=1.0, seed=1, cos_theta_cm=None, decay_at_rest=False):
    """Simulate n catalysis events  M + t -> M + W,  W -> products.

    E_kin   : MM kinetic energy in the lab (GeV)
    m_M     : MM mass (GeV)
    Q_max   : the invariant mass of W is  sum(m_i) + Q,  Q ~ U(0, Q_max)  (strong-scale assumption);
              Q_max = 0 gives W exactly at threshold (2->2 for a single product).
    cos_theta_cm : None -> isotropic (q_J = 1/2); or a number to force a fixed CM angle.
    decay_at_rest: True -> every unstable secondary is stopped before it decays (dense medium,
                   e.g. the solar core); False -> all decays in flight (ISM).

    Returns dict with the lab energies of all final particles after the decay chains,
    the vertex-level quantities and the neutrino table.
    """
    rng = np.random.default_rng(seed)
    ch = CHANNELS[channel]
    m_t = MASS[ch['target']]
    prods = ch['products']
    masses = np.array([MASS[q] for q in prods])
    fr = cm_frame(E_kin, m_M, m_t)
    rs = fr['sqrt_s']

    mW_max = fr['s_minus_mM2'] / (fr['sqrt_s'] + m_M)     # = sqrt(s) - m_M, exact
    if mW_max <= masses.sum():
        return None                                         # channel kinematically closed
    Q_eff = min(Q_max, 0.999 * (mW_max - masses.sum())) if len(prods) > 1 else 0.0
    Q = rng.uniform(0, Q_eff, n) if Q_eff > 0 else np.zeros(n)
    mW = masses.sum() + Q
    pcm, EW, EMs = vertex_cm(fr, m_M, mW)
    cth = rng.uniform(-1, 1, n) if cos_theta_cm is None else np.full(n, float(cos_theta_cm))
    # MM at angle theta_cm w.r.t. the boost axis; W opposite
    pzM, pzW = pcm * cth, -pcm * cth
    pT2 = pcm ** 2 * (1 - cth ** 2)
    EW_lab = lab_energy(EW, pzW, mW ** 2, pT2, fr['beta_cm'], fr['gamma_cm'], fr['one_minus_beta_cm'])
    EM_lab = lab_energy(EMs, pzM, m_M ** 2, pT2, fr['beta_cm'], fr['gamma_cm'], fr['one_minus_beta_cm'])

    evt = np.arange(n)
    parts = []
    if len(prods) == 1:
        parts.append(dict(name=prods[0], E=EW_lab, evt=evt, gen=0))
    else:
        bW, gW, ombW = parent_boost(EW_lab, mW)
        if len(prods) == 2:
            vecs = two_body_rest(mW, masses[0], masses[1], rng)
        else:
            vecs = three_body_rest(mW, *masses, rng)
        for q, mq, v in zip(prods, masses, vecs):
            Estar = np.sqrt(mq ** 2 + np.sum(v ** 2, axis=1))
            Elab = lab_energy(Estar, v[:, 2], mq ** 2, v[:, 0] ** 2 + v[:, 1] ** 2, bW, gW, ombW)
            parts.append(dict(name=q, E=Elab, evt=evt, gen=0))
    vertex = [dict(p) for p in parts]
    final = decay_all(parts, rng, at_rest=decay_at_rest)

    # neutrino table
    rows = [(p['name'], p['E'], p['evt'], np.full(p['E'].size, p['gen'])) for p in final if p['name'] in NEUTRINOS]
    if rows:
        flav = np.concatenate([np.full(r[1].size, r[0]) for r in rows])
        E = np.concatenate([r[1] for r in rows]); ev = np.concatenate([r[2] for r in rows])
        gen = np.concatenate([r[3] for r in rows])
    else:
        flav, E, ev, gen = np.array([], str), np.array([]), np.array([], int), np.array([], int)
    return dict(channel=channel, E_kin=E_kin, m_M=m_M, n=n, frame=fr, mW=mW, cos_theta_cm=cth,
                EW_lab=EW_lab, EM_lab=EM_lab, vertex=vertex, final=final,
                nu=dict(flavour=flav, E=E, evt=ev, gen=gen))


# ----------------------------------------------------------------------------
# summaries
# ----------------------------------------------------------------------------
def summarize(res, prompt_only=False):
    """per-flavour: multiplicity per catalysis, mean, median, 10-90 % range of E_nu (GeV)."""
    out = {}
    if res is None:
        return out
    nu = res['nu']
    sel_all = np.ones(nu['E'].size, bool) if not prompt_only else (nu['gen'] == 0)
    for f in NEUTRINOS + ('all',):
        sel = sel_all if f == 'all' else (sel_all & (nu['flavour'] == f))
        E = nu['E'][sel]
        if E.size == 0:
            out[f] = dict(mult=0.0, mean=np.nan, median=np.nan, p10=np.nan, p90=np.nan)
            continue
        out[f] = dict(mult=E.size / res['n'], mean=E.mean(), median=np.median(E),
                      p10=np.percentile(E, 10), p90=np.percentile(E, 90))
    return out


def draft_estimate(E_kin, m_M, E0=(0.030, 0.050)):
    """the draft's Eq. 5 assumption: E_nu = gamma_M * E_nu^(0), products at rest in the MM frame."""
    g = 1 + E_kin / m_M
    return g * E0[0], g * E0[1]


def analytic_limits(E_kin, m_M, m_t, mW):
    """E_W,lab for theta_cm = 0 (MM straight on), pi (MM bounced back) and the isotropic mean."""
    fr = cm_frame(E_kin, m_M, m_t)
    rs = fr['sqrt_s']
    p, EW, _ = vertex_cm(fr, m_M, mW)
    g, b, omb = fr['gamma_cm'], fr['beta_cm'], fr['one_minus_beta_cm']
    fwd = g * (mW ** 2 / (EW + p) + omb * p)        # W backward in CM
    bwd = g * (EW + b * p)                           # W forward in CM
    mean = g * EW
    return dict(forward=fwd, backward=bwd, mean=mean, sqrt_s=rs, gamma_cm=g)
