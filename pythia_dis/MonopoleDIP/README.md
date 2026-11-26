# Heavy charged projectile on a proton at rest: pi+- and neutrinos above 200 PeV (DIS + photoproduction)

Projectile: scalar or spin-1/2, mass m = 1, 10, 100 TeV, charge q = 60 x muon charge (z = 60; cross section scaled by z^2 = 3600, one-photon exchange).
Kinetic energies T = 1e11 and 1e14 GeV (proton at rest). Primary-vertex pi+- (Pythia), then pi -> mu nu -> e nu nu decays in vacuum (all neutrinos counted).

## Folder content
| path | content |
|---|---|
| gen2.py | DIS generator (Q2 > lo slice), Pythia 8.318 |
| php_gen.py | photoproduction generator (Q2 < 1 GeV^2), Pythia 8.318 |
| nu.py | pion -> mu nu -> e nu nu chain in the collinear limit, with muon polarization |
| combine.py | combines DIS slices (analytic massive-projectile flux weight) + photoproduction, decays, spectra, tail statistics |
| run_combine_resume.py | same, but computes only missing keys of results/res_comb.json |
| plot_comb.py, make_table.py | plots/ and results/summary_table.md |
| run_generation.sh | commands used to produce the events (see comments) |
| events/ | DIS events h{11,14}_{k}.npz (k = Q2 slice; slice boundaries below) |
| php/ | photoproduction events php_T{11,14}_m{1e3,1e4,1e5}.npz |
| results/ | res_comb.json (all numbers and spectra), summary_table.md |
| plots/ | spectra_combined_scalar.png, spectra_combined_spin_half.png |
| logs/ | generation and combination logs |

Run from this folder: `python3 combine.py` (needs numpy, pythia8mc 8.318). Pythia is not needed to recombine, only to regenerate events.

## Settings
### DIS (gen2.py), Q2 > 1 GeV^2
- Pythia: `Beams:frameType=1`, idA = 13 (muon used as massless carrier), idB = 2212, eCM = sqrt(M^2 + 2 M E) with E = T + m ("massless-equivalent" run, S = 2ME).
- `WeakBosonExchange:ff2ff(t:gmZ)=on`, `WeakZ0:gmZmode=1` (photon exchange only), `PhaseSpace:Q2Min = lo`, `PDF:lepton=off`, built-in NNPDF2.3 (frozen below x ~ 1e-9).
- `SpaceShower:dipoleRecoil=on` (default recoil distorts the lepton kinematics in ~20% of events and violates energy conservation), QED showers by leptons off.
- `HadronLevel:Decay=on`, `ParticleDecays:limitTau0=on`, `tau0Max=1 mm`: prompt resonances decay, pi+-, K+-, mu stay stable.
- Events are generated in the CM frame and boosted by hand to the target rest frame with light-cone variables (precision), E_lab = 0.5(Lam k+ + k-/Lam), Lam = sqrt((2E+M)/M).
- Q2-stratified: slices with Q2 > [1, 4, 1e2, 1e4, 1e6, 1e8, 1e10] (and 1e12 at T = 1e14); each slice weight sigma_gen/n; an event is assigned to the slice of its hard-process Q2. 25000 events/slice (30000 for Q2 > 1).
- Heavy mass and spin enter through an analytic reweighting of the photon flux (combine.py, flux_w), with y = Q2/(xS), Q2min = m^2 y^2/(1-y), weight 0 for Q2 < Q2min:
  - spin 1/2: [(1-y)(1-Q2min/Q2) + y^2/2] / [(1-y) + y^2/2]
  - scalar:   [(1-y)(1-Q2min/Q2)] / [(1-y) + y^2/2]
- Charge: multiplicative factor Z2 = 3600.

### Photoproduction (php_gen.py), Q2 < 1 GeV^2
- Beam id 15 (tau) carrying the real mass: `15:m0 = m`, `15:mWidth=0`, `15:mayDecay=off`; true kinematics s = M^2 + m^2 + 2ME.
- `PDF:lepton2gamma=on` (equivalent-photon flux with the correct m-dependent Q2min), `Photon:Q2max=1.0`, `Photon:Wmin=5.0` (minimum allowed by Pythia), `SoftQCD:all=on` (VMD + MPI), same decay settings. 100000 events per (m, T); 4-5% of events fail at T=1e14 and are normalised away.
- Flux is the fermion one; for y <~ 1/m the scalar flux is identical (photoproduction events have y <~ 1e-3, 1e-4, 1e-5 for m = 1, 10, 100 TeV at 1e14).
- Charge factor Z2 = 3600 as for DIS.

### Decays (nu.py)
Collinear limit. Pion lab energy fraction x_mu uniform in [r,1], r = (m_mu/m_pi)^2 = 0.573; muon polarization P_L = (2r/x - (1+r))/(1-r); muon-decay neutrino densities with polarization; 3 neutrinos per pion; R = 6 decay replicas per event.

## Results (results/summary_table.md, all with z^2 = 3600)
- mean_N = expected number of neutrinos with E > 200 PeV (2e8 GeV) per collision (collision = any DIS Q2>1 or photoproduction event, i.e. normalised to sigma_tot); P_ge1 = probability of at least one.
- sigma_nu = mean_N x sigma_tot: absolute cross section for neutrinos > 200 PeV (all flavours).
- Scalar and spin-1/2 agree to better than 1% (the spin difference only matters at large y).

| m | T [GeV] | sigma_tot [mb] | N(>200 PeV) per collision | P(>=1) | N per DIS collision |
|---|---|---|---|---|---|
| 1 TeV | 1e11 | 424 | 1.2e-4 | 7.3e-5 | 4.7e-4 |
| 10 TeV | 1e11 | 253 | 6.2e-7 | 4.0e-7 | 4.0e-6 |
| 100 TeV | 1e11 | 150 | 1.05e-9 | 6.9e-10 | 1.3e-8 |
| 1 TeV | 1e14 | 1106 | 0.95 | 0.155 | 3.1 |
| 10 TeV | 1e14 | 838 | 0.24 | 0.070 | 0.76 |
| 100 TeV | 1e14 | 601 | 0.016 | 0.0080 | 0.050 |

res_comb.json keys: s{0|0.5}_m{GeV}_T{11|14}; `edges` = log10(E/GeV) bin edges (0.5 dex); `pi`, `nu_mu`, `nu_e`, `pi_php`, `pi_dis` = E dsigma/dE [microbarn per bin]; sigma_*_ub in microbarn; sigP_* = sigma(>=1 nu above 200 PeV).

## Known limitations / systematics
- Matching at Q2 = 1 GeV^2 is discontinuous: per unit charge^2, dsigma/dlnQ2 is ~13 ub just below (Pythia VMD/photoproduction) vs ~31-66 ub just above (PDF-based DIS) at T=1e14, and ~6.7 vs ~10-23 ub at T=1e11. Expect a factor 2-3 systematic on sigma_tot, hence on per-collision numbers, and on the T=1e14 tail cross section.
- Q2 < 1 and W < 5 GeV (and elastic/diffractive-elastic) are not included.
- PDFs frozen at x < ~1e-9; one-photon exchange only (z alpha = 0.44, multi-photon/ Coulomb corrections not included); no energy loss/EW (W, Z) corrections; no kaon/charm decay neutrinos (K+-, only pi+- decays are followed, pi0 etc. ignored).
- Earlier DIS-only runs (Q2 > 4) gave larger per-collision numbers because sigma_tot was smaller; the absolute sigma_nu at T=1e11 is unchanged by the photoproduction.
