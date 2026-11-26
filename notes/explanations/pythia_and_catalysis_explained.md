# Catalysis neutrinos, the vertex, and Pythia: everything needed to understand the problem and its limits

Carlos Costa, 26 November 2026. Written for myself, to be able to explain the work to the coauthors and to answer questions about it. Technical companions: `notes/kinematics_note.md` (the kinematics), `notes/pythia_replication_note.md` (Daniele's electromagnetic run), `notes/pythia_catalysis_note.md` (the Pythia-with-our-vertex run). Code: `src/mm_kinematics.py`, `src/pythia_catalysis.py`.

---

## 0. The story in twelve lines

1. KM3NeT saw a neutrino of about 220 PeV = 2.2×10⁸ GeV (event KM3-230213A). The paper asks whether a magnetic monopole crossing the Galaxy could have made it.
2. A monopole makes neutrinos by **catalysis**: when it touches a proton or an electron, its core turns that particle into a lepton plus hadrons (Callan–Rubakov effect). The pions and kaons then decay into neutrinos.
3. The question for me was: **what is the energy of those neutrinos** as a function of the monopole kinetic energy E_M and its mass m_M? Benchmark: E_M = 10¹¹ GeV.
4. The draft assumed the products are made at rest in the monopole frame and then boosted by γ_M. That is wrong for m_M ≤ 10⁵ GeV, because the monopole frame is not the centre-of-mass frame.
5. The right way: compute the invariants (√s, γ_cm), decide the angular distribution at the vertex, and boost. For the minimal monopole the vertex is isotropic in the centre of mass, and on average **half of E_M** goes to the light products.
6. Result: the neutrinos carry about 2 % of E_M each, independent of m_M for proton targets. So E_M ≈ 10¹⁰ GeV reproduces the event.
7. Daniele ran Pythia for monopole–proton collisions. That run is **not** catalysis: it is ordinary photon exchange (electromagnetic scattering) of a heavy charged projectile. It transfers 10⁻⁷ of E_M to the proton, so its neutrinos are at tens of GeV, not hundreds of PeV.
8. I replicated Daniele's run exactly, to be sure what it contains, and then asked: can Pythia do the catalysis itself?
9. Answer: **no, not natively.** Pythia has no catalysis process, does not know the monopole, cannot hadronize the electron-channel final state in the mode we need, and cannot handle 10¹¹ GeV numerically.
10. What Pythia *can* do is the second half of the problem: turn quarks into hadrons and decay them. So I generated the vertex myself and gave Pythia only the quark-level final state, in a frame where energies are of order GeV.
11. Pythia's hadronization gives the same neutrino spectrum as my simple phase-space model: same median, same shape, 4–8 neutrinos per catalysis. The result of step 6 stands.
12. What no generator can decide, because it is physics input: whether the vertex is really isotropic, how heavy the hadronic system is, and whether the monopole sees the whole proton or one quark inside it.

---

## 1. The problem in plain words

### 1.1 What we need to compute

One number family: **N_ν(E_ν > 220 PeV) per catalysis**, and behind it the full neutrino spectrum dN/dE_ν for each channel, as a function of E_M and m_M. From that the paper reads off which (E_M, m_M) can produce the event, and multiplies by the catalysis cross section σ_cat and the column of matter crossed to get a rate.

### 1.2 The actors

| symbol | meaning | typical value here |
|---|---|---|
| m_M | monopole mass | 10²–10⁵ GeV (paper's window); 10¹⁶ GeV for a GUT monopole |
| E_M | monopole **kinetic** energy in the lab (Galaxy frame) | 10¹¹ GeV after one crossing of the Galactic magnetic field; independent of m_M |
| γ_M = 1 + E_M/m_M | Lorentz factor of the monopole | 10⁶–10⁹ |
| target | proton or electron of the interstellar medium, at rest | m_p = 0.938 GeV, m_e = 0.000511 GeV |
| √s | total energy in the centre-of-mass (CM) frame of monopole + target | 433 TeV on a proton, 10 TeV on an electron, at E_M = 10¹¹ GeV |
| γ_cm | Lorentz factor of the CM frame seen from the lab | 2.3×10⁵ on a proton |
| W | the "light system": everything produced except the monopole | lepton + hadrons, mass m_W of order 1 GeV |
| y | fraction of E_M handed to W | catalysis: flat, median 0.5; photon exchange: median 4×10⁻⁷ |

Why E_M is independent of the mass: a monopole with magnetic charge g gains energy g B L in a field B over a path L, and that does not involve m_M. Lighter monopoles are faster (bigger γ_M) but carry the same energy.

### 1.3 Catalysis and the vertex

The core of a GUT monopole contains the heavy X bosons of grand unification. A quark that enters the core can change into a lepton or an antiquark without any suppression by the X mass, because the quark is inside the region where the X field lives. This is the Callan–Rubakov effect. In practice the monopole acts like a catalyst: it is unchanged, and the target is converted.

The **vertex** is the point of the diagram where this happens: the monopole line goes straight through, the target line comes in, and the lepton plus quark lines go out. The vertex is characterised by:

- which particles come out (the channel, see 1.4);
- the cross section σ_cat (order of the geometric size 1/GeV², i.e. 10⁻²⁷ cm² = 1 mb, with a large uncertainty, 10⁻³¹ to 10⁻²⁴ cm²);
- the angular distribution of the products in the CM frame. For the lowest partial wave and the minimal charge combination q_J = ½ it is **isotropic**: the products come out in any direction with equal probability.

The word is vertex, not vortex. A vortex is a whirlpool.

### 1.4 The channels and their particles

| target | what happens at the quark level | observed final state | gives neutrinos how |
|---|---|---|---|
| proton (uud) | u + u → d̄ + e⁺; the spectator d and the d̄ form a neutral meson | p → e⁺π⁰, e⁺ρ⁰, e⁺η, e⁺ω → e⁺ + pions | π⁰ gives none; ρ⁰, η, ω → π⁺π⁻… → ν |
| proton | u + u → s̄ + μ⁺; the d and s̄ form a K⁰ | p → μ⁺K⁰ | μ⁺ decay (3 ν); K⁰ decays |
| proton | u + d → d̄ + ν̄ₑ; the u and d̄ form a π⁺ | p → ν̄ₑπ⁺ | one **prompt** ν̄ₑ from the vertex itself, plus 3 ν from π⁺ → μ⁺ → e⁺ |
| electron | e⁻ → ū ū d̄; three antiquarks form an antibaryon and mesons | e⁻ → p̄π⁺π⁻, p̄K⁰ | the mesons; the p̄ is stable |

Particle glossary: e⁺ positron; μ⁺ antimuon (decays in 2.2 µs at rest to e⁺ ν_e ν̄_μ); π⁺ charged pion (decays to μ⁺ν_μ); π⁰ neutral pion (decays to two photons, no neutrinos); ρ⁰, η, ω neutral mesons that decay to pions; K⁰ neutral kaon (a d s̄ pair, which appears as K_S → π⁺π⁻ or K_L → π l ν); p̄ antiproton, the antiparticle of the proton, stable; ν̄ₑ electron antineutrino.

Suppression: Khoze argued that for β > 0.2 the proton no longer falls onto the core as a whole; the monopole hits individual quarks and the proton channels are suppressed by ξ ≈ (m_p/m_M)². The electron channels are not suppressed, which is why they matter for heavier monopoles.

---

## 2. Why the kinematics decides everything

### 2.1 The draft's picture and why it fails

The draft took each product to be born at rest in the monopole frame with E⁰ ≈ 30–50 MeV (a pion decay at rest), then boosted it by γ_M: E_ν = γ_M E⁰. That gives E_ν ∝ 1/m_M at fixed E_M and produced the mass window of the draft's Fig. 1.

This is correct only if the monopole frame is the CM frame, i.e. m_M ≫ γ_M m_p. For a GUT monopole that holds. For m_M = 10³ GeV and E_M = 10¹¹ GeV, γ_M m_p ≈ 10⁸ GeV ≫ m_M. The proton, seen from the monopole, carries 10⁸ GeV; nothing is "at rest".

### 2.2 The invariants

With the target at rest and the monopole of total energy E = E_M + m_M:

- s = m_M² + m_t² + 2 E m_t. At E_M = 10¹¹ GeV on a proton: s ≈ 1.88×10¹¹ GeV², √s ≈ 4.33×10⁵ GeV = 433 TeV, for any m_M ≪ 4×10⁵ GeV.
- γ_cm = (E + m_t)/√s ≈ 2.3×10⁵. β_cm ≈ 1 − 10⁻¹¹.
- In the CM frame the light system W and the monopole go back to back with momentum p* ≈ (s − m_M²)/(2√s) ≈ 2.2×10⁵ GeV.

### 2.3 Why half of E_M goes to the products

Boost W back to the lab: E_W = γ_cm (E*_W + β_cm p* cos θ*), where θ* is the CM angle of W.

- θ* = 180° (W goes backward in the CM, monopole goes straight on): E_W ≈ m_W. The products are essentially at rest. This is the only configuration where the draft's picture is right.
- θ* = 0° (W goes forward, monopole deflected backward in the CM): E_W ≈ 2 γ_cm p* ≈ E_M. All the energy goes to the products.
- Isotropic θ*: E_W is uniform between those extremes, mean E_M/2.

Intuition: deflecting a very heavy, very fast particle requires the light particle to carry away a lot of momentum. Any non-forward configuration is expensive, and the isotropic vertex spends that price in most events. The photon exchange of Daniele's run, by contrast, is peripheral: it strongly prefers tiny deflections, hence y ~ 10⁻⁷.

### 2.4 The two halves of the calculation

1. **The vertex**: s, γ_cm, θ*, the mass m_W of the light system. This fixes how much energy and in which direction W moves in the lab. It is pure two-body kinematics plus one physics assumption (isotropy).
2. **Hadronization and decays**: how W turns into a lepton plus pions and kaons, and how those decay to neutrinos in flight. This is where a generator like Pythia has real knowledge: the Lund string model, measured decay tables, branching ratios.

My first calculation did half 1 exactly and half 2 with a minimal model (the exclusive final states of the table, flat phase space, our own decay chains). The Pythia exercise replaces half 2 with Pythia and keeps half 1.

---

## 3. What Pythia is and what it does

Pythia 8 is a Monte Carlo event generator for high-energy collisions. One event goes through:

1. **Hard process**: a perturbative matrix element (quark–quark scattering, photon exchange, Z production…) chosen from a library and sampled with the proper phase space and parton densities of the beams.
2. **Parton showers**: extra gluon and photon radiation.
3. **Hadronization**: colour-connected quarks and gluons are joined by "strings" (the Lund model) that break into hadrons. A colour singlet string between a quark and an antiquark gives mesons; a string ending on a diquark gives a baryon.
4. **Decays**: every unstable hadron and lepton decays according to measured tables (branching ratios, matrix elements, lifetimes).

Everything is written in the **event record**: one line per particle with an id (the PDG code: 211 is π⁺, −11 is e⁺, 12 is ν_e, 1 is d, −3 is s̄, 2203 is the uu diquark), a status (1 = final, 23 = outgoing from the hard process, 81–89 = primary hadrons from a string), **colour tags** (two integers saying which other parton a quark's colour line connects to), and the four-momentum and mass.

**Standalone hadronization**: Pythia can be told to skip steps 1 and 2 (`ProcessLevel:all = off`). You then fill the event record yourself with partons and leptons and call `pythia.next()`, which runs steps 3 and 4. This is the standard way to feed Pythia a process it does not have, either directly, as we did, or through a Les Houches Event (LHE) file.

**What Daniele's run uses.** `MonopoleDIP` uses a normal hard process: one-photon exchange between a charged lepton beam and a proton (deep inelastic scattering for Q² > 1 GeV², photoproduction for Q² < 1 GeV²). The heavy charge of the monopole is mimicked by an electric charge z = 60 and an analytic reweighting for its mass and spin. All of that exists in Pythia, so Pythia could do it alone. It describes the electromagnetic energy loss of a fast monopole and its soft-pion background. It does not contain baryon-number violation, so it is not catalysis.

---

## 4. Why Pythia cannot do catalysis by itself

Each point below is something I tested or checked, not a guess.

### 4.1 There is no catalysis process in the library

Pythia's hard processes are perturbative: a Feynman diagram with a coupling constant, squared, integrated over phase space. Callan–Rubakov catalysis is not a diagram. It is a boundary condition that the monopole core imposes on the fermion fields; the cross section is geometric, of order 1/GeV², and does not come from any small coupling. There is no matrix element to give to Pythia and no setting like `Catalysis:all = on`. Pythia cannot know σ_cat, its β dependence, or the Khoze suppression. Those have to come from theory papers and be multiplied in afterwards.

### 4.2 Pythia does not know the monopole

PDG reserves the code 4110000 for a magnetic monopole. Pythia 8.318 has no entry for it: `particleData.isParticle(4110000)` returns false. One can add it by hand (`4110000:new = MM MMbar 2 0 0 1000.0`, i.e. a stable spin-½ neutral heavy particle of 1000 GeV), and that is what the code does, but a hand-added particle has no parton densities and no interactions. Pythia cannot use it as a **beam**. Beams must be particles Pythia can describe: protons, nuclei, leptons, photons. So even if the process existed, "monopole beam on proton target" is not a configuration Pythia accepts.

### 4.3 The electron channel needs a junction that standalone mode does not hadronize

In the electron channel three antiquarks (ū ū d̄) must end up in an antibaryon. In QCD three colour lines meet at a **junction**. Pythia has a junction object (`event.appendJunction`) and hadronizes junction systems when they come from its own processes (proton remnants, baryon-number-violating supersymmetry). In standalone hadronization mode, however, I found that it accepts the junction without any error message and then does nothing: with junction kind 1 `next()` returns false, with kind 2 it returns true but the partons are left unchanged, and `forceHadronLevel` also fails silently. There is no documented way to make it work in this mode.

The remedy is Pythia's own way of describing baryons in strings: two of the antiquarks are combined into an **antidiquark** (ūū, PDG −2203), and the string is drawn between the antidiquark and the remaining d̄. That is an ordinary two-ended string and it hadronizes normally, giving an antibaryon (p̄, n̄, Δ̄…) plus mesons. The price is that the diquark is treated as a point-like string end, which slightly raises the threshold (constituent masses 0.77 + 0.33 GeV instead of m_p = 0.94 GeV).

### 4.4 Numerical precision at 10¹¹ GeV

Pythia uses double precision, about 16 significant digits. Consider a pion with E = 10¹¹ GeV:

- E² ≈ 10²² GeV². Its mass squared is m² = E² − p² ≈ 0.02 GeV². To obtain 0.02 from a difference of two numbers of order 10²² you need 24 digits. In double precision m² comes out as garbage (zero or negative). All the internal mass checks and decay kinematics break.
- The boost from the CM frame to the lab has γ ≈ 10⁵ to 10¹⁰. A backward-moving product has E_lab = γ(E* − β p*), a difference of two nearly equal numbers. With β = 1 − 10⁻¹¹ or smaller, 1 − β is below what double precision can represent next to 1 (β rounds to exactly 1.0 for γ > 10⁸).

The solution is to let Pythia work only in the W rest frame, where every energy is of order 1 GeV, and to do the boost to the lab ourselves with formulas that never subtract nearly equal numbers (they use 1 − β and m²/(E + p) instead of E − p). This is what `mm_kinematics.lab_energy` does, and it was checked to conserve energy to 10⁻¹⁴ in every event.

### 4.5 What Pythia cannot know even in principle

Pythia is a tool for the second half of the problem. The first half is input. Nothing in Pythia tells us:

- whether the catalysis vertex is isotropic in the CM;
- how heavy the hadronic system is (how much energy the core releases into the strong sector, our parameter Q_max);
- whether the monopole interacts with the whole proton or with one parton carrying a fraction x of its momentum (Khoze's picture for β > 0.2, which reduces √s by √x);
- the value and β dependence of σ_cat.

Those decide the result far more than hadronization details do.

---

## 5. How the vertex was built and fed to Pythia

The code is `src/pythia_catalysis.py`. One event, with numbers for E_M = 10¹¹ GeV, m_M = 1000 GeV, proton target, channel p → e⁺ + X(dd̄):

**Step 1, the vertex (our code).**
- s = 1.88×10¹¹ GeV², √s = 433 TeV, γ_cm = 2.31×10⁵.
- Draw the mass of the light system: m_W = m_e + m_X,min + Q with Q uniform in [0, Q_max]. m_X,min = 0.66 GeV is Pythia's threshold for a d d̄ string (two constituent masses of 0.33 GeV). With Q_max = 1 GeV, m_W ≈ 1.2 GeV on average.
- Two-body CM kinematics M + W: p* ≈ 2.165×10⁵ GeV.
- Draw cos θ* uniformly in [−1, 1] (isotropy). Boost W to the lab: E_W runs from about m_W (θ* = 180°) to about 10¹¹ GeV (θ* = 0°), uniformly. Record the exact β_W, γ_W and 1 − β_W of W without cancellations. For a typical event with E_W = 5×10¹⁰ GeV, γ_W = 4×10¹⁰ and 1 − β_W = 3×10⁻²², a number double precision could never obtain from β itself.

**Step 2, the quark-level final state in the W rest frame (our code).**
- W → e⁺ + X back to back, isotropic direction, with m_X drawn uniformly between m_X,min and m_W − m_e.
- X → d + d̄ back to back in the X frame, each with its constituent mass, then boosted into the W frame. All energies are now of order 1 GeV.

**Step 3, Pythia (its job).**
- `event.reset()`, then `event.append(-11, 1, 0, 0, px, py, pz, E, m_e)` for the positron (status 1, no colour), `event.append(1, 23, 101, 0, …)` for the d with colour tag 101, and `event.append(-1, 23, 0, 101, …)` for the d̄ with anticolour tag 101. The matching tags tell Pythia the two quarks form one colour-singlet string.
- `pythia.next()`: the string fragments (typically into π⁺π⁻, or a ρ⁰ or ω that then decays, or a single π⁰ when m_X is at threshold), and every unstable particle decays, including π±, K±, K_L, K_S and μ±, because in the interstellar medium everything decays before interacting. This is set with `ParticleDecays:limitTau0 = off` and `mayDecay = on` for those ids.
- Loop over the final particles (status 1) and read id, energy, and the momentum component along the W direction.

**Step 4, back to the lab (our code).**
- Each final particle is boosted with the stored β_W, γ_W, 1 − β_W of its event, using the cancellation-free formula. Neutrinos (|id| = 12, 14, 16) are collected with their flavour.
- Check: the sum of lab energies of all final particles equals E_W to 10⁻¹⁴, in every one of the 50 000 events per run, and Pythia rejected no event.

For the electron channel there is no outgoing lepton: W is the antidiquark–antiquark string itself. For the μ⁺ channel the string is d s̄, for the ν̄ₑ channel u d̄ (and the ν̄ₑ is appended as a final particle like the positron, so it appears directly in the lab spectrum as the prompt neutrino).

---

## 6. What came out

At E_M = 10¹¹ GeV, m_M = 1 TeV, Q_max = 1 GeV, per catalysis, all neutrino flavours:

| channel | ν per catalysis: Pythia / phase space | median E_ν [GeV]: Pythia / phase space | ν above 220 PeV: Pythia / phase space |
|---|---|---|---|
| p → e⁺ + X | 4.3 / 6.0 (e⁺π⁺π⁻) | 2.2×10⁹ / 2.6×10⁹ | 3.9 / 5.5 |
| p → μ⁺ + X | 7.8 / 6.1 (μ⁺K⁰) | 1.7×10⁹ / 2.3×10⁹ | 6.8 / 5.5 |
| p → ν̄ₑ + X | 5.0 / 4.0 (ν̄ₑπ⁺) | 2.3×10⁹ / 5.2×10⁹ | 4.5 / 3.8 |
| e⁻ → X | 3.3 / 6.0 (p̄π⁺π⁻) | 1.2×10⁹ / 1.2×10⁹ | 2.8 / 5.1 |

How to read it:

- **The spectra coincide.** In `figures/pythia_catalysis_vs_phasespace.png` the black Pythia curve lies on the coloured phase-space curve in every panel: same peak at 10⁹–10¹⁰ GeV, same cutoff at E_M, same low-energy slope. The spectrum is set by the vertex, not by how the string breaks.
- **Multiplicities differ a little** because Pythia's inclusive string contains all the mesons at their proper rates (π⁰ gives no neutrinos, ρ⁰/ω give four), while the phase-space model used one exclusive final state per channel. Pythia's numbers are the better founded ones: 4 to 8 neutrinos per proton catalysis, 3 per electron catalysis.
- **Q_max** (how much energy the vertex puts into the strong sector) is the one knob that moves things: from 1 to 10 GeV the multiplicity doubles and the median drops by 3.5, but the number of neutrinos above 220 PeV changes by less than a factor 2. Energy is conserved, so more neutrinos means softer neutrinos, and the count above threshold is robust.
- **Mass dependence** (`figures/pythia_catalysis_electron_mass.png`): proton channels do not depend on m_M up to 100 TeV because √s = 433 TeV ≫ m_M. The electron channel does: √s is only 10 TeV, so for m_M = 100 TeV the CM frame is the monopole frame and the products get only 10⁷ GeV. Pythia and phase space agree exactly here too, because this is vertex kinematics.
- **Decay chains**: our own pion decay code and Daniele's `nu.py`, applied to the same pions, agree to 2 % (`figures/decay_chain_crosscheck.png`). Pythia's decay tables agree as well. The decay stage is not a source of uncertainty.

---

## 7. Daniele's run next to ours

| | electromagnetic run (Pythia alone) | catalysis (our vertex + Pythia) |
|---|---|---|
| process | one-photon exchange, charge z = 60 | Callan–Rubakov, baryon number violated |
| cross section at 10¹¹ GeV, 1 TeV | 424 mb | σ_cat ≈ 1 mb reference, times ξ = (m_p/m_M)² for protons if Khoze applies |
| energy to the light system, median y | 4×10⁻⁷ | 0.5 |
| neutrinos per event | 41 | 4–8 |
| median E_ν | 44 GeV | 2×10⁹ GeV |
| neutrinos above 220 PeV per event | 1.2×10⁻⁴ | 4–7 |

They are two different physical processes, not two estimates of the same thing. The electromagnetic one is the energy loss and the soft background; catalysis is the source of hard neutrinos. Both should appear in the draft, named as such.

---

## 8. What is solid and what is an assumption

Solid (checked numerically, independent of models):

- The invariants s, γ_cm, and the fact that the monopole frame is not the CM frame for m_M ≤ 10⁵ GeV.
- Energy and momentum conservation event by event in both codes.
- The decay chains (two independent codes plus Pythia agree).
- The fact that photon exchange cannot deliver 220 PeV at E_M = 10¹¹ GeV (Daniele's own run, replicated).

Assumptions (physics input, to be stated in the paper):

- Isotropic vertex in the CM (q_J = ½, lowest partial wave). If the final state is forward peaked, every catalysis spectrum moves down toward the electromagnetic one.
- The mass of the light system: threshold plus Q ≤ 1 GeV ("strong-scale" release). Robust for the count above 220 PeV, not for the median.
- Whole-proton target. If the monopole interacts with a single parton of momentum fraction x (Khoze, β > 0.2), √s falls by √x and the spectrum softens; not implemented.
- σ_cat and its β dependence, and whether the Khoze suppression applies.

Technical caveats of the Pythia pipeline, all small: thresholds set by constituent masses instead of physical hadron masses; unpolarised muon decay (2 % effect); antidiquark string instead of a junction for the antibaryon.

---

## 9. Questions I expect, and the answers

**Did you simulate the monopole in Pythia?** No. Pythia never sees the monopole move or interact. We computed the collision kinematics, and Pythia received only the lepton and the quarks produced, in their own rest frame.

**Then what did Pythia add?** The hadronization (which hadrons, how many, with what momenta) and all decays with measured tables. That replaces our hand-made exclusive final states and confirms them.

**Why not write a Les Houches file?** An LHE file is just another way to fill the event record. Filling it directly in Python is simpler and avoids writing 10¹¹ GeV momenta into a text file, where the precision problem of 4.4 would appear immediately.

**Why does Daniele's code work with Pythia alone?** Because photon exchange between a charged lepton and a proton is a standard process. He replaced the monopole by a charged lepton with an analytic reweighting for mass, spin and charge. That is legitimate for electromagnetic scattering and impossible for catalysis.

**Is 4110000 a real PDG code?** Yes, reserved for a monopole, but Pythia has no particle attached to it. We add one by hand only so that the event record can mention it; it plays no dynamical role.

**Why do the proton channels not depend on the monopole mass?** Because √s = 433 TeV is far above m_M, so the monopole is a negligible part of the CM energy balance. Only when m_M approaches √s (as in the electron channel at 100 TeV) does the mass matter.

**What is y?** The fraction of E_M given to the light system. It is what distinguishes the two processes: 10⁻⁷ for photon exchange, 0.5 for catalysis.

**Could Geant4 do this?** Geant4 has a monopole physics list for ionisation and tracking, not for catalysis. Same situation: the vertex would have to be written by hand.

**What should the paper say about Pythia?** That Pythia was used for hadronization and decays of a user-supplied catalysis vertex, and separately for the electromagnetic energy-loss channel. Not that Pythia simulated catalysis.

---

## 10. Glossary

| term | meaning |
|---|---|
| catalysis, Callan–Rubakov effect | conversion of a proton or electron by a monopole core, with baryon number violation and no suppression by the GUT mass |
| vertex | the point of the diagram where the conversion happens; characterised by channel, cross section and CM angular distribution |
| channel | one specific final state, e.g. p → e⁺π⁺π⁻ |
| inclusive / exclusive | inclusive: a class of final states (e⁺ plus any hadrons from a d d̄ string); exclusive: one specific final state |
| CM frame | frame where the total momentum of monopole + target is zero |
| √s | total energy in the CM frame |
| γ_cm, β_cm | Lorentz factor and velocity of the CM frame seen from the lab |
| W, light system | all products except the monopole |
| y | fraction of E_M transferred to W |
| Q_max | maximum energy above threshold put into the hadronic system; our knob for its mass |
| hadronization | formation of hadrons from quarks; in Pythia the Lund string model |
| string | colour field between a quark and an antiquark (or diquark); breaks into hadrons |
| junction | point where three colour lines meet; needed for a baryon made of three separate quarks |
| diquark | two quarks treated as one string end; Pythia's way to make baryons in strings |
| constituent mass | effective quark mass inside a hadron used by Pythia: 0.33 GeV (u, d), 0.50 GeV (s), 0.77 GeV (uu diquark) |
| event record | Pythia's list of particles in an event with id, status, colour tags, momentum |
| status code | 1 final particle; 23 outgoing parton from the hard process; 81–89 primary hadron from a string |
| colour tag | integer label; a quark and an antiquark with the same tag are connected by a string |
| standalone hadronization | Pythia mode with `ProcessLevel:all = off`: the user fills the event record and Pythia only hadronizes and decays |
| LHE | Les Houches Event file, text format for passing parton-level events between generators |
| PDG code | integer identifying a particle type (211 π⁺, 13 μ⁻, 12 ν_e, 2212 p, 4110000 monopole) |
| DIS, photoproduction | deep inelastic scattering (Q² > 1 GeV²) and quasi-real photon exchange (Q² < 1 GeV²): the two parts of Daniele's electromagnetic run |
| Khoze suppression | factor ξ ≈ (m_p/m_M)² for proton catalysis at β > 0.2, when the monopole sees individual quarks |
| 220 PeV | 2.2×10⁸ GeV, the energy of KM3-230213A |
