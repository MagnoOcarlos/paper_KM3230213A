# Catalysis neutrinos for the KM3-230213A monopole paper — full explanation (10 Nov 2026)

Working notes for Carlos. Everything below refers to the draft *MM 200 PeV nu* (file `MM_200_PeV_nu (1).pdf`, sections and equation numbers as in that version) and to Daniele's request. Companion files in this folder: `mm_kinematics.py` (code), `neutrino_kinematics.ipynb` (executed notebook), `kinematics_note.md` (note for the coauthors), `calculation_walkthrough.md` (step-by-step with references), `figures/`, `results/`.

Papers cited: **[C]** Candela, Khoze, Turner, arXiv:2504.14918; **[K]** Khoze, arXiv:2405.18689; **[H]** Hu et al., arXiv:2201.02386; **[SK]** Super-Kamiokande, arXiv:1203.0940; **[HLS]** Houston, Li, Sun, arXiv:1803.02835; **[D]** the draft.

---

## 1. What Daniele asked

Two messages define the task.

* **Daniele's direct email to Carlos** (late October): *what is the kinetic energy of the neutrinos produced by catalysis, as a function of the monopole kinetic energy, channel by channel?* Benchmark E_M = 10^11 GeV. Use [C] Eqs. 2.1 and 2.2 for the list of processes, [SK] Eqs. 1–3 for the proton-decay channels, [H] Sects. 3–4 for how the neutrinos emerge. For the antiproton side take only e⁻ → p̄ π⁺π⁻ and e⁻ → p̄ K⁰.
* **Daniele to the group, 4 Nov 2026** (thread p. 71): "we need to know from Carlos the details about the kinematic and update the draft." Pedro agreed the same day and asked to include also the deep-inelastic case. Meeting Friday 13 Nov, 9:00.

So the deliverable is a **calculation**: numbers and formulas that replace Sect. 2.1 / Eq. 5 of [D], with the reasoning that justifies them.

Why he needs it: everything in [D] rests on Eq. 5,

> γ ∼ E_ν / E_ν⁽⁰⁾ ∼ 10⁹–10¹⁰ with E_ν⁽⁰⁾ ∼ 30–50 MeV,

obtained by assuming "catalysis products are dominated by pion decay at rest in the MM frame" ([D] Sect. 2.1). That number fixes the mass window 10²–10⁵ GeV ([D] Sect. 2.1, 4.1, 6, Fig. 1, abstract, conclusions). In October the group found (thread, 4–19 Oct) that (i) neutrinos are almost never produced at the vertex, (ii) the catalysis cross section changes for fast monopoles ([K] Sect. 6.2), (iii) at these energies pions, kaons and muons do not decay inside the atmosphere or the detector (Pedro, 2 Nov). All three make the kinematics the open point.

---

## 2. What catalysis is (background for the whole problem)

**Mechanism.** A GUT monopole is a knot of the Higgs field. Outside it is a magnetic charge; inside a core of radius ~1/M_GUT the unified symmetry is restored and the X, Y bosons that turn quarks into leptons are massless. Ordinary proton decay goes through a virtual X and is suppressed by 1/M_GUT⁴; a fermion that reaches the core gets the transition for free. Rubakov (1981) and Callan (1982); the elementary process is u + u + M → d̄ + e⁺ + M ([K] Eq. 1.1), ΔB = ΔL = −1, B−L conserved.

**Why the cross section is large.** The fermion–monopole system carries an extra angular momentum eg/4π in the field. In the lowest partial wave it cancels the centrifugal barrier, so the fermion is pulled into the core whatever the core size. The cross section is then set by the strong-interaction scale: σ ∼ q_J²/GeV² ([C] p. 4, "not suppressed by the monopole mass or by any Grand Unification scale"), i.e. 10⁻²⁷–10⁻²⁴ cm² as [D] Sect. 2 uses. Only the lowest partial wave contributes ([C] p. 4, [K] Sects. 4–6).

**The channels.** Hadronising d̄ + e⁺ with the spectator quark gives p → e⁺π⁰, p → e⁺ + (ρ⁰, η, ω) → pions, p → μ⁺K⁰ with the second family ([SK] Eq. 1, [H] Eq. 3.2), and the rare prompt p → ν̄ₑπ⁺ with branching ~(m_p/m_W)² ∼ 10⁻⁴ ([H] Eq. 3.2, [HLS]). For a neutron, n → e⁺π⁻ and n → ν̄ₑπ⁰ ([D] Eq. 3). The monopole is unchanged: a catalyst.

**The reverse process.** The same core turns an electron into antiquarks: e⁻ → p̄ + mesons ([C] Eq. 2.2). It needs enough centre-of-mass energy to make a proton, which for a monopole on an electron at rest means γ_M ≳ 1800 ([C] Eq. 2.8), and it yields an antiproton of at least 861 GeV in the lab ([C] Sect. 2.3, Fig. 2).

**Velocity dependence.** For β < 0.2 the proton falls on the core as a whole and σ ∝ 1/β² ([C] Eq. 2.5, [K] Eq. 6.13). Above that the valence quarks miss the core and only sea quarks with small impact parameter reach it: geometric factor ξ = (m_p/m_M)² ([K] Sect. 6.2, p. 21–22). The electron has no substructure, so e⁻ channels keep the full cross section ([K] Sect. 6.1). This is Daniele's 15 Oct message and the reason the antiproton channels are in the paper.

**Where neutrinos come from.** Not from the vertex (except the 10⁻⁴ channel): from π⁺ → μ⁺ν_μ, μ⁺ → e⁺ν_eν̄_μ, kaon decays, and their conjugates. e⁺π⁰ gives photons only. In dense matter the mesons stop first and give the 29.8 MeV line and the Michel spectrum ([SK] Fig. 1, [H] Sect. 3); in the ISM they decay in flight; in the atmosphere/rock at PeV–EeV energies they are absorbed before decaying (Pedro, 2 Nov).

**Angles.** Lowest partial wave only ⇒ for the minimal monopole (q_J = 1/2) the amplitude has no angular dependence in the centre-of-mass frame: [C] Eq. 2.3 and 2.5, [K] Eq. 6.5 and 6.12, exponent 4|q_J|−2 = 0 ([C] p. 8; [K] p. 19). This isotropy is the single physics input that drives the result.

---

## 3. The concept that replaces Eq. 5 of the draft

[D] Sect. 2.1 assumes the pions and muons are at rest in the monopole frame, E_ν⁽⁰⁾ = 30–50 MeV, and multiplies by γ_M. In the monopole frame, however, the target proton arrives with energy γ_M m_p ≈ 10⁹–10¹⁰ GeV, far above the monopole's own mass of 10²–10⁵ GeV. Energy conservation alone forbids products at rest there. The monopole frame is **not** the centre-of-mass frame; the monopole is the *light* partner of the collision.

### Key equations (lab frame: target at rest, monopole with kinetic energy E_M)

Invariant mass:

    s = m_M² + m_t² + 2 (E_M + m_M) m_t                          [C] Eq. 2.6

Centre-of-mass boost:

    γ_cm = (E_M + m_M + m_t) / √s ,    β_cm = p_M / (E_M + m_M + m_t)    [C] Eq. A.9

Two-body M + t → M + W in the CM frame (W = light system of mass m_W):

    E_W* = (s + m_W² − m_M²) / 2√s ,   p* = λ^{1/2}(s, m_M², m_W²) / 2√s    [C] Eqs. A.2, A.3, A.11

Back to the lab, monopole scattered by θ_cm, W opposite:

    E_W^lab = γ_cm ( E_W* − β_cm p* cos θ_cm )

Isotropic θ_cm (lowest partial wave, q_J = 1/2):

    ⟨E_W^lab⟩ = γ_cm E_W*   ≈ E_M / 2        light monopole, m_M ≪ √(2E_M m_t)   (flat between ~m_W²/2m_t and ~E_M)
                            ≈ γ_M² m_t       heavy monopole, m_M ≫ √(2E_M m_t)   ("relativistic mirror")
                            → m_t shared     β_M → 0                              (this is [D] Eq. 5)

Then W decays with flat phase space into its constituents, and π, K, μ decay in flight (PDG branching ratios, V−A Michel spectrum for the muon). Each neutrino ends with ~1–3 % of E_M.

### What γ_cm is

The Lorentz factor of the frame where monopole and target have equal and opposite momenta, as seen from the detector. It is **not** γ_M = 1 + E_M/m_M. At E_M = 10¹¹ GeV on a proton: √s = 4.3×10⁵ GeV, γ_cm = 2.3×10⁵, while γ_M = 10⁸ for m_M = 10³ GeV. In the CM frame the proton and the monopole each carry ~216 TeV; after an isotropic scattering the light products carry up to that, and boosting back with γ_cm spreads their lab energy uniformly between ~0 and ~E_M. That is the origin of "half the monopole energy" and of the mass independence.

### Variables

| Symbol | Meaning | Value here |
|---|---|---|
| E_M | monopole kinetic energy (lab) | benchmark 10¹¹ GeV; scan 10⁸–10¹⁴ |
| m_M | monopole mass | 10²–10⁵ GeV ([D] window) |
| γ_M, β_M | monopole Lorentz factor and velocity, γ_M = 1 + E_M/m_M | 10⁶–10⁹ |
| t, m_t | target at rest (p or e⁻) and its mass | 0.938 / 0.000511 GeV |
| s, √s | invariant mass squared / invariant mass of the collision | √s ≈ 433 TeV (proton, benchmark) |
| γ_cm, β_cm | boost of the CM frame relative to the lab | 2.3×10⁵ |
| W, m_W | light final-state system and its invariant mass = Σ masses + Q | ~0.3–2 GeV |
| Q, Q_max | internal energy of W, uniform in [0, Q_max] | Q_max = 1 GeV default |
| E_W*, p* | energy and momentum of W in the CM frame | ≈ √s/2 (light MM) |
| θ_cm | monopole scattering angle in the CM frame | cos θ_cm uniform |
| E_W^lab | lab energy of W | flat 0 … E_M, mean E_M/2 |
| λ | Källén function (x−y−z)² − 4yz | |
| q_J | half the magnetic charge in Dirac units; 1/2 = minimal monopole | 1/2 |
| ξ | Khoze's suppression of proton channels for β > 0.2, (m_p/m_M)² | 10⁻⁴ at 10² GeV, 10⁻¹⁰ at 10⁵ GeV |
| σ_cat | catalysis cross section (not computed here) | 10⁻²⁷–10⁻²⁴ cm² ([D] Sect. 2) |
| E_ν | lab energy of one neutrino | medians 1–3×10⁹ GeV at benchmark |
| E_ν⁽⁰⁾ | [D]'s rest-frame neutrino energy, 30–50 MeV | valid only for β_M → 0 |
| E_KM3 | event energy | 220 PeV = 2.2×10⁸ GeV |
| gen | 0 = neutrino made at the vertex; ≥1 = from a decay chain | only gen 0 counts in the atmosphere |
| κ_ν | neutrinos per catalysis ([D] Sect. 2) | 4–6 in the ISM; ~10⁻⁴ ξ in the atmosphere |
| x | Michel variable 2E/m_μ | 0–1 |

---

## 4. Channels

E_M = 10¹¹ GeV, m_M = 10³ GeV, all decays in flight (ISM). From `results/benchmark_E1e11.csv`.

| Channel | Target | W | Neutrino sources | ν per catalysis | Median E_ν | 10–90 % | Weight (slow MM, literature) | Relativistic factor | Counts where |
|---|---|---|---|---|---|---|---|---|---|
| p → e⁺π⁰ | p | e⁺π⁰ | none (π⁰ → γγ) | 0 | – | – | ≈ 0.5 of e⁺+meson ([H] Eq. 3.2) | ξ | gamma rays only |
| p → e⁺π⁺π⁻ (via ρ⁰,η,ω) | p | e⁺π⁺π⁻ | π⁺, π⁻ chains | 6 | 2.6×10⁹ GeV | 2.4×10⁸ – 1.2×10¹⁰ | ≈ 0.5 | ξ | ISM |
| p → μ⁺K⁰ | p | μ⁺K⁰ | μ⁺; K_S → ππ; K_L → πℓν | 6.1 | 2.3×10⁹ | 2.0×10⁸ – 1.2×10¹⁰ | 1/400–1/2 ([H]) | ξ | ISM |
| p → ν̄ₑπ⁺ | p | ν̄ₑπ⁺ | **prompt ν̄ₑ** + π⁺ chain | 1 + 3 | prompt: 1.6×10¹⁰ | 1.6×10⁹ – 5.3×10¹⁰ | ~10⁻⁴ ([H], [HLS]) | ξ | prompt ν̄ₑ is the only neutrino in the atmosphere |
| e⁻ → p̄π⁺π⁻ | e⁻ | p̄π⁺π⁻ | π⁺, π⁻ chains (p̄ stable) | 6 | 1.2×10⁹ | 1.1×10⁸ – 5.9×10⁹ | [C] Eq. 2.3 | none | ISM |
| e⁻ → p̄K⁰ | e⁻ | p̄K⁰ | K⁰ chain (p̄ stable) | 4.1 | 1.2×10⁹ | 1.1×10⁸ – 6.5×10⁹ | [C] Eq. 2.3 | none | ISM |
| e⁻ → p̄ | e⁻ | p̄ | none; check against [C] 861 GeV | 0 | – | – | | | check only |

Flavour split of the pion chains 1 : 1 : 2 : 2 (ν_e : ν̄_e : ν_μ : ν̄_μ); kaon channels add ν_e from K_L. Proton-channel numbers are identical for m_M = 10²–10⁴ GeV and 5 % lower at 10⁵. Electron channels hold up to m_M = 10³ GeV and fall for heavier monopoles (heavy regime on an electron: E_W ≈ γ_M² m_e; median 6×10⁸ GeV at 10⁴, 1.2×10⁷ GeV at 10⁵).

Monopole energy at which the **median** neutrino has 220 PeV (`results/required_EM_for_220PeV.csv`):

| Channel | m_M = 10² | 10³ | 10⁴ | 10⁵ GeV |
|---|---|---|---|---|
| p → e⁺π⁺π⁻ | 8.6×10⁹ (10–90 %: 1.9×10⁹–9.3×10¹⁰) | 8.6×10⁹ | 8.6×10⁹ | 1.2×10¹⁰ |
| p → μ⁺K⁰ | 9.7×10⁹ | 9.7×10⁹ | 9.8×10⁹ | 1.4×10¹⁰ |
| prompt ν̄ₑ (atmosphere) | 1.4×10⁹ | 1.4×10⁹ | 1.4×10⁹ | 1.5×10⁹ |
| e⁻ → p̄π⁺π⁻ | 1.8×10¹⁰ | 1.9×10¹⁰ | 5.3×10¹⁰ | 4.3×10¹¹ |
| e⁻ → p̄K⁰ | 1.8×10¹⁰ | 1.9×10¹⁰ | 5.3×10¹⁰ | 4.3×10¹¹ |
| [D] Eq. 5 (γ_M = 220 PeV / 40 MeV) | 5.5×10¹¹ | 5.5×10¹² | 5.5×10¹³ | 5.5×10¹⁴ |

---

## 5. How to read the plots

**`figures/spectra_E1e11.png`** — six panels, one per channel, E_M = 10¹¹ GeV. x: E_ν (log). y: neutrinos per catalysis per log bin (sum over bins = multiplicity). Line styles = the four masses. Black vertical line = 220 PeV; grey dotted = E_M (no neutrino can exceed it). Curves rise with slope 1 in log–log (flat dN/dE, the signature of isotropic CM scattering), peak at a few 10⁹ GeV, cut off below E_M. 220 PeV sits on the rising flank. Proton panels: all masses coincide. Electron panels: 10⁴ shifts slightly left, 10⁵ by two decades (heavy regime). Top-left panel empty (π⁰).

**`figures/Enu_vs_EM.png`** — left p → e⁺π⁺π⁻, right e⁻ → p̄π⁺π⁻. x: E_M, y: E_ν, both log. Solid line = median; pale band = 10–90 %; hatched narrow bands = [D] Eq. 5 for each mass (γ_M × 30–50 MeV), far below and spaced by a decade per decade of mass; black horizontal = 220 PeV. Where the solid line crosses 220 PeV you read the required E_M (≈ 10¹⁰ GeV, left); where the band crosses, the range. Slope 1 in the light regime; slope 2 at low E_M for the heaviest masses on the right (E_ν ∝ γ_M² m_e).

---

## 6. Relation to the draft, section by section

* **Abstract, Sect. 7.** "γ ∼ 10⁹–10¹⁰, achievable for masses 10²–10⁵ GeV; Galactic acceleration alone cannot" → the event needs E_M ≈ 10¹⁰ GeV; the mass enters only through acceleration, the flux bounds and ξ; the Galactic gain of Eq. 6 (~10¹¹ GeV) is sufficient.
* **Sect. 2, Eqs. 2–4, κ_ν.** κ_ν ≈ 1–2 → 4–6 per catalysis in the ISM (0 for e⁺π⁰); in the atmosphere κ_ν = BR(prompt) × ξ ≈ 10⁻⁴ ξ. State the two separately.
* **Sect. 2.1, Eq. 5.** The paragraph to rewrite: "pion decay at rest in the MM frame" holds only for β_M → 0. Replace with the invariant-mass treatment above. Daniele's "[DP: check also direct production]" → prompt ν̄ₑ carries ~0.15 E_M. The sentence on which masses reach the boost disappears; GUT-scale monopoles remain excluded, but because at E_M = 10¹¹ GeV they are non-relativistic (back in the 30–50 MeV regime), not for lack of γ.
* **Sect. 2.2, Sect. 3.1, Appendix.** Losses computed at γ = 10⁹. With E_M ≈ 10¹⁰ GeV, γ ≤ 10⁸ and Eq. 9 gives ~20× smaller losses; the Appendix table must be recomputed. The conclusion that the monopole cannot cross 147 km of rock stands.
* **Sect. 3.2 (atmosphere).** "114 neutrinos per crossing with κ_ν = 2" → only the prompt ν̄ₑ, BR ~10⁻⁴ × ξ: many orders of magnitude fewer; the argument for a smaller σ_cat weakens. Its energy 0.15 E_M means 220 PeV needs E_M ≈ 1.4×10⁹ GeV in this scenario.
* **Sect. 3.3.** Unchanged, confirmed by Pedro's decay lengths.
* **Sect. 4 (directionality).** The Galactic-only regime was discarded for lack of γ. With E_M ≈ 10¹⁰ GeV it is open; the field-line-following case must be discussed again.
* **Sect. 5.** Formula unchanged; use the new multiplicities and the flat spectrum from ~10⁸ GeV to E_M. Gamma-ray comparison: e⁺π⁰ puts comparable energy into photons at similar energies.
* **Sect. 6, Fig. 1.** The mass range was justified by the 30–50 MeV uncertainty; that goes. Mass axis from acceleration to 10¹⁰ GeV, bounds, and ξ. The required-flux curves do not move (E₀ cancels in Eq. 17).
* **Pedro's closing remark (PDL).** Sharper now: the ISM scenario gives several neutrinos per catalysis with a flat spectrum reaching far above 220 PeV, which IceCube should also see.

---

## 7. How the calculations were made

The whole thing is a Monte Carlo of the kinematics, written in `mm_kinematics.py` and driven by `neutrino_kinematics.ipynb`. No cross sections are evaluated: every catalysis event is given weight one, and the physics input is only (a) the frames, (b) the isotropic CM angle, (c) flat phase space for the decays, (d) PDG branching ratios. For each configuration (channel, E_M, m_M) I generate N events and read the neutrino energies off the result.

**Step by step for one event** (function `catalysis`, line ≈283):

1. **Frames.** From E_M, m_M and the target mass compute s, √s, γ_cm, β_cm (`cm_frame`, line ≈93). Two numerical precautions: s − m_M² is formed as 2(E_M + m_M)m_t + m_t² instead of subtracting two numbers of order 10³² (needed for the GUT-mass checks), and 1 − β_cm is formed as s / [E_tot (E_tot + p)], which is exact, instead of 1 − p/E_tot, which would be 0 to double precision when γ_cm ~ 10⁵–10⁷.
2. **Mass of the light system.** m_W = Σ(constituent masses) + Q, Q uniform in [0, Q_max], Q_max = 1 GeV, clipped so that m_W ≤ √s − m_M (lines ≈305–310). For the single-particle check e⁻ → p̄, m_W = m_p exactly.
3. **CM two-body kinematics.** E_W*, E_M*, p* from the Källén function (`vertex_cm`, line ≈110, again written with s − m_M²).
4. **Angle.** cos θ_cm drawn uniformly in [−1, 1] (line ≈312); the monopole goes at θ_cm, W opposite, so p_z,W* = −p* cos θ_cm and p_T² = p*²(1 − cos²θ_cm).
5. **Boost of W to the lab.** E_W^lab = γ_cm (E_W* + β_cm p_z,W*) (`lab_energy`, line ≈120). For backward emission the bracket is a difference of nearly equal numbers; it is rewritten exactly as (E* − |p_z|) + (1 − β)|p_z| with E* − |p_z| = (m² + p_T²)/(E* + |p_z|). The same routine is used for every later decay. Also E_M^lab is computed so that energy conservation can be checked event by event.
6. **Decay of W.** In the W rest frame, isotropic two-body (`two_body_rest`) or flat three-body phase space (`three_body_rest`: m₁₂ sampled by rejection with weight p*(W → (12)+3) × p*((12) → 1+2), then the standard sequential construction). Each constituent is boosted to the lab with β_W, γ_W obtained from E_W^lab and m_W (`parent_boost`, line ≈133, with 1 − β = m²/[E(E + p)]). Only the momentum component along the parent direction is needed, because every decay is isotropic in its parent's rest frame; directions are therefore never stored.
7. **Decay chains** (`decay_all`, line ≈234). A queue of particles; each unstable one picks a mode according to the PDG table `DECAYS` (line ≈47), decays in its rest frame, and the daughters are boosted with the parent's lab energy. Muons use the V−A matrix element: on flat massless three-body phase space the weight x_νe (1 − x_νe) (for μ⁺) reproduces the Michel spectra (`michel_rest`, line ≈203). π⁰ → γγ stops the chain. p̄, e±, γ, ν are stable. Each neutrino is stored with flavour, lab energy, event number and **generation** (0 if it came out of W directly, otherwise the number of decays above it). Option `decay_at_rest=True` sets every parent velocity to zero before it decays, which is the dense-medium (Sun) case used only for checking.
8. **Summary** (`summarize`, line ≈352). Per flavour: multiplicity = neutrinos / N, mean, median, 10th and 90th percentile. `prompt_only=True` keeps generation 0 only (the atmosphere case).

**Runs made in the notebook.**

* Sect. 1: analytic table (`analytic_limits`) for the forward, backward and isotropic-mean E_W, and [D] Eq. 5 for comparison.
* Sect. 3: benchmark E_M = 10¹¹ GeV, six channels × four masses, N = 2×10⁵ events each → `results/benchmark_E1e11.csv`, spectra figure.
* Sect. 4: scan of 13 energies from 10⁸ to 10¹⁴ GeV × four masses × six channels, N = 4×10⁴ each → `results/scan_Enu_vs_EM.csv`, figure `Enu_vs_EM.png`; then the E_M at which the median (and the 10 %, 90 % quantiles) reach 220 PeV by log–log interpolation → `results/required_EM_for_220PeV.csv`.
* Sect. 5: prompt ν̄ₑ only, from the p → ν̄ₑπ⁺ run → `results/prompt_nubar_e_E1e11.csv`.
* Sect. 6: sensitivity, Q_max = 0.3/1/3 GeV and cos θ_cm forced to 0, 0.9, 0.99.
* Sect. 7: the five checks (energy conservation; [C] 861 GeV; γ_M² m_p; √s vs m_M; stopped-decay lines).

**Running time.** About 0.2 s per 10⁵ events; the full notebook executes in a few minutes (`jupyter nbconvert --to notebook --execute --inplace neutrino_kinematics.ipynb`).

**What a hand check looks like.** Take E_M = 10¹¹ GeV, proton target, m_M = 10³ GeV: s ≈ 2×10¹¹×0.938 = 1.88×10¹¹ GeV², √s = 4.33×10⁵ GeV, γ_cm = 10¹¹/4.33×10⁵ = 2.3×10⁵; E_W* ≈ p* ≈ 2.17×10⁵ GeV; isotropic mean E_W^lab = γ_cm E_W* = 5×10¹⁰ GeV = E_M/2. Three constituents share it, roughly a third each; a π⁺ of 1.7×10¹⁰ GeV gives a ν_μ with on average 0.21 of that (≈ 3.5×10⁹ GeV) and a muon whose three daughters take about a third each of the rest. The medians in the tables (2–3×10⁹ GeV) are consistent with this estimate once the flat distribution in cos θ_cm is folded in.

## 8. Is the calculation correct? What is solid and what is a choice

Solid (verified in notebook Sect. 7): the relativistic kinematics (energy conservation to 10⁻¹⁶; [C]'s 861 GeV forward antiproton for m_M = 10³ and 10⁹ GeV; heavy-MM limit ⟨E_W⟩ = γ_M² m_p; √s independent of m_M); the decay chains (PDG BRs, V−A Michel spectrum; with stopped secondaries the 29.79 MeV line and the Michel mean 31.7 MeV / endpoint 52.8 MeV); and the qualitative conclusion that [D] Eq. 5 fails for light relativistic monopoles (energy conservation alone).

Choices to agree with the group: (1) **isotropy in the CM frame** extended from the two-body formulas to the whole light system — the assumption that matters; cos θ_cm fixed at 0.9 (0.99) lowers all energies by 10 (100); (2) invariant mass of W, Σm + U(0, 1 GeV) — ±30 % on the medians for Q_max = 0.3–3 GeV; (3) branching ratios taken from the slow-monopole literature; ξ listed but not applied; channels not combined.

Not included: cross sections and rates, fluxes, oscillations, detector response, Daniele's DIS case (same Steps 1–2 and 5 apply; only the final-state model differs).

---

## 9. The answer to send (sent/ready for Daniele, Pedro, Michele in copy)

See `kinematics_note.md` for the full note. Core of the message: first pass at the kinematics, channel by channel; Eq. 5 is the β → 0 picture and does not apply; with invariants and isotropic CM scattering the light system takes E_M/2 and each neutrino ~1–3 % of E_M; numbers at the benchmark as in Sect. 4 above; a 220 PeV neutrino needs E_M ≈ 10¹⁰ GeV independent of the mass (electron channels: 2×10¹⁰ for m_M ≲ 10³ GeV); checks done; the isotropy assumption is the one to agree on; DIS case can go through the same chain; Friday 13 Nov at 9 is fine.

---

## 10. Conclusion to Daniele's question

Daniele asked for the neutrino energy as a function of the monopole kinetic energy, channel by channel, at E_M = 10¹¹ GeV. The answer is: for every mass in the draft's window the catalysis products share the centre-of-mass energy, not the monopole-frame rest energy, so the light system carries on average E_M/2 and each neutrino 1–3 % of E_M — median 1–3×10⁹ GeV at the benchmark, 4–6 neutrinos per catalysis in the ISM for the channels with charged pions or kaons, none for e⁺π⁰, and in the atmosphere only the prompt ν̄ₑ of p → ν̄ₑπ⁺ at ~0.15 E_M. The neutrino energy is set by E_M and not by γ_M, so KM3-230213A needs E_M ≈ 10¹⁰ GeV whatever the mass, below the Galactic acceleration estimate of the draft. Eq. 5 and the mass-window argument built on it must be replaced; the result rests on the isotropy of the lowest-partial-wave scattering in the centre-of-mass frame, which is the point to settle at Friday's meeting.
