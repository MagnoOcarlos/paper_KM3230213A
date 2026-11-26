# Catalysis neutrinos with Pythia: what we want, why Pythia cannot do it alone, and what we did

Carlos Costa, 21 November 2026. Code: `src/pythia_catalysis.py` (uses `src/mm_kinematics.py` and Pythia 8.318 through the `pythia8mc` package). Results: `results/pythia_catalysis_E1e11.json`, `figures/pythia_catalysis_vs_phasespace.png`, `figures/pythia_catalysis_electron_mass.png`. Companion notes: `kinematics_note.md` (our phase-space calculation), `pythia_replication_note.md` (Daniele's electromagnetic run).

## 1. The physics we want to simulate

A relativistic magnetic monopole of mass m_M = 10²–10⁵ GeV and kinetic energy E_M ≈ 10¹¹ GeV crosses the interstellar medium and hits protons and electrons at rest. Through Callan–Rubakov catalysis the monopole core turns the target into a lepton plus a hadronic system without being changed itself:

| target | quark level | hadron level (Super-K 1203.0940 Eq. 1, Hu et al. 2201.02386 Eq. 3.2, Candela et al. 2504.14918 Eq. 2.2) |
|---|---|---|
| proton | u + u + M → d̄ + e⁺ + M, the spectator d and the d̄ form a neutral meson | p → e⁺π⁰, e⁺(ρ⁰, η, ω) → e⁺π⁺π⁻… |
| proton | u + u + M → s̄ + μ⁺ + M | p → μ⁺K⁰ |
| proton | u + d + M → d̄ + ν̄ₑ + M | p → ν̄ₑπ⁺ (prompt neutrino) |
| electron | e⁻ + M → ū + ū + d̄ + M, the three antiquarks form an antibaryon plus mesons | e⁻ → p̄π⁺π⁻, p̄K⁰ |

The quantity the paper needs is the lab-frame energy spectrum of the neutrinos from these final states after the pions, kaons and muons decay in flight, as a function of E_M and m_M, in particular the number of neutrinos above the KM3-230213A energy of 220 PeV per catalysis. Two pieces of physics enter:

1. **The vertex.** A contact interaction at the monopole core, lowest partial wave only. For the minimal charge combination q_J = ½ the Callan–Rubakov cross section is isotropic in the centre-of-mass frame (Candela et al. Eq. 2.3). With √s ≈ 433 TeV on a proton and γ_cm ≈ 2.3×10⁵, this isotropy is what decides the energy of the light system in the lab: on average half of E_M, with a flat distribution from GeV (forward) to E_M (backward).
2. **Hadronization and decays.** How the hadronic system turns into pions and kaons, with what multiplicity, and how those decay into neutrinos. This is where a generator like Pythia is the tool of choice.

Our previous calculation (`mm_kinematics.py`) did piece 1 exactly and piece 2 with a minimal model: the exclusive final states of the table, flat phase space, an invariant mass of the light system of threshold plus Q with Q uniform in [0, Q_max = 1 GeV], and our own decay chains with PDG branching ratios.

## 2. Why Pythia cannot do the catalysis natively

- **No process in the library.** Pythia's hard processes are perturbative matrix elements: photon, W, Z, gluon exchange and so on, listed in its process database. Callan–Rubakov catalysis is non-perturbative. It is not a Feynman diagram but a boundary condition that the monopole's core imposes on the fermion fields; the cross section is of geometric size, 1/GeV², with no coupling constant to expand in. There is no matrix element to plug into Pythia's phase-space machinery, and no switch of the form `Catalysis:all = on`.
- **The monopole is not a known particle.** PDG code 4110000 is reserved for monopoles, but Pythia 8.318 has no entry for it (checked: `particleData.isParticle(4110000)` is false). It has to be added by hand and Pythia cannot use it as a beam particle, because beams need parton densities or a lepton structure that Pythia knows.
- **Baryon-number violation at the string level.** The final state of the electron channel is three antiquarks that must be colour-connected through a junction. In the "standalone hadronization" mode Pythia accepts the junction without error but silently does not hadronize it (tested with kinds 1 and 2 and both parton status codes). The standard remedy, used by Pythia itself for proton remnants, is to represent two of the antiquarks as an antidiquark, which gives an ordinary antidiquark–antiquark string.
- **Numerical range.** Pythia works in double precision in the frame it is given. A 10¹¹ GeV system with GeV masses loses the mass information entirely (E² − p² cannot be resolved), and Lorentz factors of 10⁵–10¹⁰ give catastrophic cancellations in backward boosts. Pythia must therefore work in a frame where energies are of order GeV, and the boost to the lab must be done outside with cancellation-free formulas.

The electromagnetic run in `pythia_dis/MonopoleDIP` does not have these problems because photon exchange is a standard Pythia process; that is why it could be done with Pythia alone, and also why it is a different physical process from catalysis (see `pythia_replication_note.md`).

## 3. What we did: our vertex, Pythia's hadronization

The pipeline in `pythia_catalysis.py` is the standard way external processes are fed to Pythia, except that we fill the event record directly instead of writing a Les Houches file:

1. **Vertex, as before.** Exact s, γ_cm, β_cm; isotropic cos θ_cm; invariant mass of the light system m_W = m_ℓ + m_X,min + Q with Q uniform in [0, Q_max]. The light system's lab energy and direction come from the same cancellation-free boost as in `mm_kinematics.py`.
2. **Quark-level final state in the W rest frame.** W → lepton + X, back to back, isotropic, with the mass of the hadronic system m_X uniform between its threshold and m_W − m_ℓ. X is a colour-singlet string: d d̄ for the e⁺ channel, d s̄ for the μ⁺ channel, u d̄ for the ν̄ₑ channel, (uu)̄ antidiquark + d̄ for the electron channel. The string ends receive Pythia's constituent masses (0.33 GeV for u, d; 0.50 GeV for s; 0.77 GeV for the uu antidiquark), which set the thresholds m_X,min = 0.66, 0.83, 0.66 and 1.10 GeV.
3. **Pythia in the W rest frame.** `ProcessLevel:all = off`; the lepton (status 1) and the two string ends (status 23, colour tags 101) are appended to the event record and `pythia.next()` runs the Lund string fragmentation and all decays. Every unstable particle is allowed to decay, including π±, K±, K_L, K_S and μ, which is the interstellar-medium case. Pythia's own decay tables and branching ratios are used, including the muon's V−A matrix element (without the spin correlation from the pion decay, a 2 % effect as checked against Daniele's `nu.py`).
4. **Boost to the lab.** Each final particle's momentum component along the W direction is boosted with the exact β_W, γ_W and 1 − β_W of W, using `mm_kinematics.lab_energy`. Neutrinos are collected by flavour.

Checks: no event failed in Pythia; the sum of the lab energies of all final particles equals the lab energy of W to 10⁻¹⁴ in every event; with Q_max → 0 the string has two primary hadrons and the result tends to the exclusive channels.

## 4. Results at E_M = 10¹¹ GeV

Per catalysis, all flavours, decays in flight, m_M = 1 TeV unless stated. "Phase space" is our previous exclusive model with Q_max = 1 GeV.

| channel | model | ν per catalysis | median E_ν [GeV] | N_ν above 220 PeV | primary hadrons |
|---|---|---|---|---|---|
| p → e⁺ + X | Pythia string, Q_max = 1 GeV | 4.3 | 2.2×10⁹ | 3.9 | 2.0 |
| | Pythia string, Q_max = 3 GeV | 5.4 | 1.4×10⁹ | 4.6 | 2.2 |
| | Pythia string, Q_max = 10 GeV | 8.9 | 0.63×10⁹ | 6.5 | 3.2 |
| | phase space p → e⁺π⁺π⁻ | 6.0 | 2.6×10⁹ | 5.5 | 2 |
| | phase space p → e⁺π⁰ | 0 | – | 0 | 1 |
| p → μ⁺ + X | Pythia string, Q_max = 1 GeV | 7.8 | 1.7×10⁹ | 6.8 | 2.0 |
| | Pythia string, Q_max = 10 GeV | 12.2 | 0.71×10⁹ | 9.1 | 3.2 |
| | phase space p → μ⁺K⁰ | 6.1 | 2.3×10⁹ | 5.5 | 1 |
| p → ν̄ₑ + X | Pythia string, Q_max = 1 GeV | 5.0 | 2.3×10⁹ | 4.5 | 2.0 |
| | Pythia string, Q_max = 10 GeV | 9.8 | 0.76×10⁹ | 7.4 | 3.2 |
| | phase space p → ν̄ₑπ⁺ | 4.0 | 5.2×10⁹ | 3.8 | 1 |
| e⁻ → X, m_M = 1 TeV | Pythia string, Q_max = 1 GeV | 3.3 | 1.2×10⁹ | 2.8 | 2.0 |
| | Pythia string, Q_max = 3 GeV | 5.9 | 0.92×10⁹ | 4.7 | 2.7 |
| | phase space e⁻ → p̄π⁺π⁻ | 6.0 | 1.2×10⁹ | 5.1 | 3 |
| | phase space e⁻ → p̄K⁰ | 4.1 | 1.2×10⁹ | 3.4 | 2 |
| e⁻ → X, m_M = 10 TeV | Pythia string, Q_max = 1 GeV | 3.3 | 0.61×10⁹ | 2.4 | 2.0 |
| e⁻ → X, m_M = 100 TeV | Pythia string, Q_max = 1 GeV | 3.3 | 1.2×10⁷ | 0.01 | 2.0 |

Flavour composition with Pythia (Q_max = 1 GeV, 1 TeV): the e⁺ channel gives ν_e : ν̄_e : ν_μ : ν̄_μ = 0.7 : 0.7 : 1.4 : 1.4 per catalysis; the μ⁺ channel 2.0 : 1.0 : 1.9 : 2.9, the extra ν_e and ν̄_μ coming from the μ⁺ itself; the ν̄ₑ channel 1.2 : 1.2 : 1.3 : 1.3, with the prompt ν̄ₑ at a median of 3.5×10⁹ GeV; the electron channel 0.5 : 0.6 : 1.1 : 1.1.

## 5. Reading the figures

`figures/pythia_catalysis_vs_phasespace.png`: one panel per inclusive channel, neutrinos per catalysis per 0.25 dex, black for Pythia with Q_max = 1, 3, 10 GeV (solid, dashed, dotted), colour for our exclusive phase-space channels, the dotted vertical line at 220 PeV.

- In every panel the Pythia curve with Q_max = 1 GeV lies on top of the phase-space curve, with the same peak at 10⁹–10¹⁰ GeV, the same cutoff at E_M and the same low-energy slope. Hadronization does not move the spectrum; the vertex kinematics does.
- Pythia's inclusive d d̄ string sits between our two exclusive e⁺ channels: it contains π⁰ (no neutrinos), ρ⁰, η and ω with Pythia's own rates, so it gives 4.3 neutrinos per catalysis against 0 and 6.
- Raising Q_max increases the multiplicity and lowers the median in proportion, about a factor 3.5 in the median between 1 and 10 GeV, while the number of neutrinos above 220 PeV changes by less than a factor 2. This replaces the 30 % sensitivity to Q_max quoted in the previous note with a proper statement: the median scales as 1 over the multiplicity, the count above threshold is robust.
- Bottom right, the electron channel: Pythia's antidiquark string mostly gives antiproton plus one meson (two primary hadrons), hence fewer neutrinos than the p̄π⁺π⁻ phase-space channel, but with the same median.

`figures/pythia_catalysis_electron_mass.png`: the electron channel for three masses, Pythia solid and phase space dashed. The mass dependence is identical in the two models because it comes from the vertex: √s on an electron is 10 TeV, so a 100 TeV monopole has its centre of mass at rest with itself and the products get 10⁷ GeV, below the event.

## 6. What this changes and what it does not

- The median neutrino energy, about 2 % of E_M for the proton channels, and the required E_M ≈ 10¹⁰ GeV for 220 PeV are confirmed with Pythia's hadronization and decay tables in place of our minimal model. The result of `kinematics_note.md` stands.
- The multiplicity per catalysis is now better founded: 4 to 8 neutrinos for the proton channels and 3 for the electron channel at the strong-scale invariant mass, rising with Q_max.
- What Pythia cannot settle, because it is our input: the isotropy of the vertex in the centre of mass (if the final state is forward peaked, every spectrum moves down), the invariant mass of the light system (Q_max), and whether the target is the whole proton or a single parton (Khoze's picture for β > 0.2, which would reduce the collision energy by the parton momentum fraction and is not implemented here).
- Small technical caveats: the inclusive thresholds are set by Pythia's constituent masses, slightly above the physical π⁰ and p̄ thresholds; muon decays are unpolarised; the electron-channel antibaryon is built as an antidiquark string rather than a junction.

## 7. How to run

```
./pyenv/bin/python src/pythia_catalysis.py        # ~2 min for 15 runs of 50 000 events; needs pythia8mc 8.318.3
```

`run(channel, E_kin, m_M, n, Q_max, seed)` returns the lab energies and ids of every final particle together with the neutrino table and the vertex quantities, so new channels or energies are one call away. A new channel needs only its lepton, its string ends and its physical threshold in the `CHANNELS` dictionary.
