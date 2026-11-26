# From zero to the result: a beginner's guide to the monopole-neutrino problem

*Written for Carlos, who knows CRPropa and cosmic-ray propagation but has not worked on monopoles, catalysis or collision kinematics before. Read it top to bottom once; the later files (`explanation_2026-11-10.md`, `calculation_walkthrough.md`, `kinematics_note.md`) then make sense.*

---

## 0. The story in ten lines

1. KM3NeT saw one muon of about 70 PeV, implying a neutrino of about 220 PeV. IceCube, which is bigger, saw nothing like it. That tension makes people look for exotic explanations.
2. Our paper (Doro, Perri, de la Torre Luque, you) asks: could a **magnetic monopole** have made that neutrino?
3. A monopole does something ordinary particles cannot: when it touches a proton it makes the proton decay into a positron and pions. This is called **catalysis** (Callan–Rubakov effect). The pions decay into neutrinos.
4. If the monopole moves very fast, the neutrinos it makes are very energetic. The draft said: to get 220 PeV you need a monopole Lorentz factor of 10⁹ to 10¹⁰, which only light monopoles (100 GeV to 100 TeV) can reach. That number is Eq. 5 of the draft and it decides the whole paper.
5. In October the group realised Eq. 5 was derived for a *slow* monopole and asked whether it survives for a fast one. Daniele asked you to redo the kinematics: **what neutrino energy comes out, as a function of the monopole energy, channel by channel?**
6. You did it with a small Monte Carlo (like a CRPropa run, but for one collision instead of a propagation).
7. Result: for a fast monopole the neutrino energy is a few percent of the monopole **kinetic energy**, whatever the monopole mass. The event needs a monopole energy of about 10¹⁰ GeV, not a Lorentz factor of 10⁹.
8. Eq. 5 of the draft is off by orders of magnitude and the mass window built on it has to be rethought.
9. The result rests on one physics assumption (the collision is isotropic in the centre-of-mass frame) which the group must accept or refine.
10. You send the note and the plots, and discuss on Friday.

---

## 1. The pieces of physics you need

### 1.1 The event
KM3-230213A: a near-horizontal muon track in KM3NeT/ARCA, 13 Feb 2023, reconstructed muon energy ≈ 70 PeV. A muon cannot cross 147 km of rock and sea from the atmosphere at that energy, so it must come from a neutrino interacting near the detector; with the usual energy sharing the neutrino had ≈ 220 PeV ([D] Sect. 1). Nothing comparable in IceCube: 2–3.5σ tension. [D] Eq. 1 is the flux KM3NeT assigns to the event.

### 1.2 Magnetic monopoles
A hypothetical particle with a single magnetic pole, predicted by Grand Unified Theories (GUT). Properties that matter here:
* **Mass m_M**: unknown. GUT monopoles ~10¹⁶ GeV; the draft considers lighter ones, 10²–10⁵ GeV, because only those can be fast enough (see below). Collider limits exclude below ~75–100 GeV.
* **Acceleration**: a magnetic field pushes a magnetic charge along the field line, the way an electric field pushes an electron. (In CRPropa a charged cosmic ray is *bent* by B; a monopole is *accelerated* by B.) The energy gain over the Galaxy is g B √(R λ) ≈ 10¹¹ GeV ([D] Eq. 6), independent of mass; intergalactic fields can give up to 10¹⁴ GeV ([D] Eq. 7, Perri et al. 2025). So a monopole arrives with kinetic energy E_M ~ 10¹¹–10¹⁴ GeV whatever its mass, and its Lorentz factor is γ_M = 1 + E_M/m_M: huge for a light monopole, ~1 for a GUT one.
* **Interaction with matter**: it ionises like a particle of charge 68.5e, radiates pairs and photonuclear showers at high γ ([D] Eqs. 8–9), and, if it is a GUT-type monopole, **catalyses nucleon decay**.

### 1.3 Catalysis in plain words
Inside the monopole there is a tiny core where the GUT symmetry is unbroken. There, quarks can turn into leptons freely. Normally a proton would need to tunnel through a 10⁻²⁹ cm barrier to do that (which is why protons live > 10³⁴ yr). But a charged particle approaching a monopole in its lowest angular-momentum state has **no centrifugal barrier**: it is sucked straight into the core. So the proton decays with a cross section of nuclear size, ~10⁻²⁷ cm² (not 10⁻⁵⁸ cm²), and the monopole comes out unchanged, ready to do it again. Hence "catalysis". Rubakov 1981, Callan 1982.

At the quark level: u + u + M → d̄ + e⁺ + M. Baryon number and lepton number each drop by one, their difference B − L is conserved, charge is conserved. Dress it with the spectator quark and you get the **channels**:

| channel | what happens | gives neutrinos? |
|---|---|---|
| p → e⁺ π⁰ | the d̄ pairs with the spectator d into a π⁰ | no (π⁰ → γγ) |
| p → e⁺ (ρ⁰, η, ω) → e⁺ π⁺ π⁻ | same, with a heavier meson that decays to charged pions | yes |
| p → μ⁺ K⁰ | second family: s̄ instead of d̄, μ⁺ instead of e⁺ | yes |
| p → ν̄ₑ π⁺ | the lepton is a neutrino; rare, ~10⁻⁴ | yes, directly |
| e⁻ → p̄ π⁺ π⁻ , e⁻ → p̄ K⁰ | the **reverse** process: an electron becomes an antiproton | yes (from the pions/kaon) |

The last line matters because of velocity: a fast monopole does not "see" a whole proton any more (its quarks fly past the core); Khoze (arXiv:2405.18689, Sect. 6.2) estimates a suppression (m_p/m_M)² for the proton channels above β ≈ 0.2. The electron has no inner structure, so the e⁻ channels are not suppressed. This is why Daniele asked for the antiproton channels too.

### 1.4 Where the neutrinos come from — you already know this from CRPropa
In CRPropa a UHECR proton hits a CMB photon, makes a Δ⁺ → π⁺ n, then π⁺ → μ⁺ ν_μ and μ⁺ → e⁺ ν_e ν̄_μ. Each neutrino takes roughly 5 % of the proton energy. **Catalysis neutrinos are produced by exactly the same decay chains**; only the first step (what makes the pion) is different. So:

    π⁺ → μ⁺ ν_μ            (ν_μ takes on average 21 % of the pion energy when the pion decays in flight)
    μ⁺ → e⁺ ν_e ν̄_μ        (each takes about a third of the muon energy)
    π⁻ → μ⁻ ν̄_μ ,  μ⁻ → e⁻ ν̄_e ν_μ
    K⁰ → K_S (→ ππ) or K_L (→ π e ν, π μ ν, 3π)
    π⁰ → γ γ                (no neutrinos)

One difference from the solar monopole searches (Super-K, Hu et al.): in the Sun the pions **stop** before decaying, so the neutrinos have fixed low energies (29.8 MeV line, 53 MeV Michel endpoint). In the interstellar medium the pions decay **in flight** and carry a fraction of their energy to the neutrinos, exactly as in CRPropa. In the atmosphere or rock at these energies they are **absorbed before decaying** (Pedro checked, 2 Nov), so there only the neutrino made directly at the vertex (the "prompt" ν̄ₑ of p → ν̄ₑπ⁺) exists.

---

## 2. What the draft assumed, and why it is wrong for our monopoles

[D] Sect. 2.1 says: the pions and muons decay at rest *in the monopole frame*, producing 30–50 MeV neutrinos there; seen from Earth these are boosted by γ_M, so E_ν = γ_M × 40 MeV, and to reach 220 PeV you need γ_M ≈ 10⁹–10¹⁰ ([D] Eq. 5).

Picture it from the monopole's point of view. The monopole sits still; the proton of the interstellar gas comes at it with Lorentz factor γ_M, i.e. with energy γ_M m_p ≈ 10⁹–10¹⁰ GeV. Now compare with the monopole's own mass, 10²–10⁵ GeV. **The incoming proton carries ten thousand to a hundred million times more energy than the monopole weighs.** The proton cannot be absorbed into pions "at rest": there is nowhere for 10⁹ GeV to go except into the products. The picture of Eq. 5 is the picture of a *slow* monopole, where the proton arrives with just its rest energy 0.94 GeV. It is right for monopoles trapped in the Sun and wrong for ours.

The correct way is to go to the frame where the total momentum is zero — the **centre-of-mass (CM) frame** — do the collision there, and transform back.

---

## 3. The relativity you need (five formulas)

If you know how CRPropa computes the pion energy in photopion production, you know all of this; here it is written out.

**(a) Energy, momentum, mass.** For any particle E² = p² + m², γ = E/m, β = p/E. Kinetic energy is E − m. For the monopole, E_M in all our files is the kinetic energy and the total is E_M + m_M.

**(b) Invariant mass of a collision.** For two particles the quantity s = (E₁+E₂)² − (p⃗₁+p⃗₂)² is the same in every frame. With the target at rest:

    s = m_M² + m_t² + 2 (E_M + m_M) m_t

√s is the total energy available in the CM frame. At E_M = 10¹¹ GeV on a proton: s ≈ 2 × 10¹¹ × 0.94 ≈ 1.9 × 10¹¹ GeV², √s ≈ 4.3 × 10⁵ GeV = 430 TeV. Notice m_M² = 10⁶ GeV² for m_M = 10³ GeV does not even register: **√s does not depend on the monopole mass** as long as m_M ≪ √(2E_M m_t) ≈ 4 × 10⁵ GeV.

**(c) The CM frame moves.** Its Lorentz factor relative to the lab is

    γ_cm = (total lab energy) / √s = (E_M + m_M + m_t) / √s ≈ 2.3 × 10⁵

Compare with γ_M = 10⁸ for a 1 TeV monopole: the CM frame is 400 times *slower* than the monopole. The monopole is the light partner; think of a ping-pong ball (monopole) hitting a bowling ball (proton) that, in the CM frame, both move with 216 TeV each.

**(d) Two-body kinematics in the CM frame.** After the collision there is the monopole and "W", the bundle of light products (e⁺π⁰, p̄π⁺π⁻, …) with some invariant mass m_W of order 1 GeV. In the CM frame W has

    E_W* = (s + m_W² − m_M²) / 2√s ≈ √s / 2 ≈ 216 TeV,   momentum p* ≈ E_W*

and flies off opposite to the monopole, at some angle θ_cm.

**(e) Back to the lab.** The lab energy of W is

    E_W^lab = γ_cm ( E_W* − β_cm p* cos θ_cm )

If the monopole continues straight on (θ_cm = 0), W goes backwards in the CM frame and gets almost nothing in the lab (this is the situation of the draft, and gives m_W²/2m_t ≈ 1 GeV). If the monopole is bounced back (θ_cm = 180°), W gets essentially the whole monopole energy. **The crucial fact**: the Callan–Rubakov process happens only in the lowest partial wave, and for the minimal monopole that means the angle θ_cm is **isotropic** (Candela et al. 2504.14918 Eq. 2.3 and 2.5: the angular factor is exactly constant for q_J = 1/2). Isotropic in cos θ_cm plus the formula above being linear in cos θ_cm means the lab energy of W is **uniformly distributed between ~0 and ~E_M, with mean E_M/2**.

That single sentence is the whole result. Everything else is sharing E_W among the pions and then among the neutrinos, exactly like in CRPropa.

**Analogy.** A light ball thrown at a heavy wall moving towards you comes back with 2γ² times its energy (relativistic mirror). Here the roles are: the *monopole* is the light ball, the *proton* is the moving wall (in the monopole frame). The products of the collision fly back towards Earth with a large fraction of the monopole's energy.

**Three regimes to remember**

| if | then ⟨E_W^lab⟩ ≈ | when |
|---|---|---|
| m_M ≪ √(2E_M m_t) ("light") | E_M / 2 | our window, proton target, all E_M of interest |
| m_M ≫ √(2E_M m_t) ("heavy") | γ_M² m_t | GUT monopoles; also m_M ≥ 10⁴ GeV on an *electron* |
| β_M → 0 | m_t shared among products | monopoles in the Sun; **this is [D] Eq. 5** |

---

## 4. The method, told as if it were a CRPropa run

In CRPropa you create a candidate, push it through modules (propagation, energy loss, interactions, decays) and observers record what arrives. Here a "candidate" is **one catalysis event**, and the modules are:

1. **Source module – the collision.** Fix channel, E_M, m_M. Compute s, γ_cm, β_cm. Draw the invariant mass of W as Σ(masses) + Q with Q uniform in 0–1 GeV (a "strong scale" guess; Candela et al. ignore the extra mesons altogether). Draw cos θ_cm uniform in [−1, 1]. Compute E_W^lab with formula (e). Also compute the monopole's outgoing energy so that energy conservation can be checked.
2. **Decay module – W → its constituents.** In W's rest frame, split W into 2 or 3 particles with flat phase space and random orientation (like a decay in CRPropa). Boost each to the lab with W's velocity.
3. **Decay module – chains.** Every unstable particle (π±, μ±, K⁰, K_S, K_L, K±) picks a decay mode by its PDG branching ratio, decays in its rest frame (muons with the proper V−A spectrum), daughters are boosted with the parent's velocity. Repeat until only stable particles remain. Nothing is absorbed — this is the ISM assumption. There is a switch to force "decay at rest" to reproduce the solar case for checking.
4. **Observer – the neutrino list.** Each neutrino is stored with flavour, lab energy, event number, and a *generation* tag: 0 if it came straight out of the vertex, ≥1 if from a decay chain. The atmosphere case = generation 0 only.
5. **Analysis.** Per channel and per flavour: number of neutrinos per catalysis, mean, median, 10th–90th percentile of E_ν. Repeat over a grid of E_M (10⁸–10¹⁴ GeV) and masses (10²–10⁵ GeV). Find the E_M at which the median neutrino has 220 PeV.

Two technical points that cost a few hours: with γ ~ 10⁵–10⁷ the usual formulas lose all digits in double precision (1 − β is 10⁻¹¹; s − m_M² for a GUT mass is 10¹⁶ out of 10³²). The code computes those differences exactly instead of by subtraction. Everything is vectorised with numpy; 10⁵ events take 0.2 s.

Files: `mm_kinematics.py` (the modules), `neutrino_kinematics.ipynb` (the run and the plots), `results/*.csv` (the numbers), `figures/*.png`.

**Checks that prove the machinery works** (notebook Sect. 7): energy conserved to 10⁻¹⁶; the forward-scattering antiproton of Candela et al., 861.4 GeV, reproduced for two monopole masses; the heavy-monopole formula γ_M² m_p reproduced; with "decay at rest" the 29.79 MeV pion line and the Michel spectrum (mean 31.7 MeV, endpoint 52.8 MeV) come out exactly.

---

## 5. The results in plain words

At E_M = 10¹¹ GeV (Daniele's benchmark), monopole mass 10³ GeV, everything decaying in flight:

| channel | neutrinos per catalysis | typical (median) E_ν | 80 % of neutrinos between |
|---|---|---|---|
| p → e⁺π⁰ | 0 | – | – |
| p → e⁺π⁺π⁻ | 6 | 2.6 × 10⁹ GeV | 2 × 10⁸ – 1 × 10¹⁰ |
| p → μ⁺K⁰ | 6 | 2.3 × 10⁹ | 2 × 10⁸ – 1 × 10¹⁰ |
| p → ν̄ₑπ⁺, the prompt ν̄ₑ alone | 1 | 1.6 × 10¹⁰ | 2 × 10⁹ – 5 × 10¹⁰ |
| e⁻ → p̄π⁺π⁻ | 6 | 1.2 × 10⁹ | 1 × 10⁸ – 6 × 10⁹ |
| e⁻ → p̄K⁰ | 4 | 1.2 × 10⁹ | 1 × 10⁸ – 6 × 10⁹ |

Read: a neutrino carries ~1–3 % of the monopole kinetic energy (the prompt one ~15 %, because it does not share with a chain). The proton-channel numbers are the same for 10², 10³, 10⁴ GeV and 5 % lower at 10⁵ GeV. The electron channels keep their numbers up to 10³ GeV and drop for heavier monopoles (on an electron target, 10⁴–10⁵ GeV is already the "heavy" regime).

Turn it around — **which monopole energy gives a 220 PeV neutrino?** About 9 × 10⁹ GeV for the proton channels, independent of the mass; about 2 × 10¹⁰ GeV for the electron channels if m_M ≲ 10³ GeV; about 1.4 × 10⁹ GeV for the prompt ν̄ₑ. The draft's Eq. 5 would say 5.5 × 10¹¹ GeV × (m_M / 100 GeV), i.e. 10¹²–10¹⁵ GeV for the masses in the window.

**Consequences for the paper**
* Eq. 5 must be replaced: E_ν is a fraction of E_M, not γ_M × 40 MeV.
* The requirement becomes "E_M ≈ 10¹⁰ GeV" and no longer selects a mass. Galactic acceleration (10¹¹ GeV, [D] Eq. 6) is enough — the draft said it was not.
* The mass window must come from other things: what fields can accelerate, the flux bounds in Fig. 1, and the (m_p/m_M)² suppression that kills the proton channels for heavy monopoles.
* In the atmosphere only the prompt ν̄ₑ survives (BR 10⁻⁴, times the suppression): the draft's "114 neutrinos per crossing" collapses.
* In the ISM there are 4–6 neutrinos per catalysis with a broad spectrum up to E_M; relevant to why IceCube saw nothing (Pedro's remark).
* The energy-loss numbers of Sect. 2.2/3.1 were computed at γ = 10⁹ and need to be redone at the new γ.

---

## 6. How to read the two plots

**`figures/spectra_E1e11.png`.** One panel per channel. Horizontal: neutrino energy (log). Vertical: neutrinos per catalysis in each log bin. Four line styles = four masses. Black line = 220 PeV. The curves rise with slope 1 (that *is* the uniform distribution of Sect. 3e), peak at a few 10⁹ GeV, stop at E_M. 220 PeV is on the rising side, i.e. a 220 PeV neutrino from a 10¹¹ GeV monopole is a low-energy outlier. All masses overlap in the proton panels; in the electron panels the 10⁵ GeV curve sits two decades lower (heavy regime). Top-left panel empty (π⁰).

**`figures/Enu_vs_EM.png`.** Horizontal: monopole kinetic energy. Vertical: neutrino energy. Solid line = median, pale band = 10–90 % range, hatched narrow bands = the draft's Eq. 5 for each mass (far below), black horizontal = 220 PeV. Where the solid line crosses 220 PeV, read the required E_M. Straight lines of slope 1 = "a fixed fraction of E_M". On the right the heavy masses bend to slope 2 at low E_M (γ_M² m_e regime).

---

## 7. What is certain and what is an assumption (to say honestly on Friday)

Certain: the relativity (Sect. 3), the decay chains, the fact that products cannot be at rest in the monopole frame, the checks.

Assumptions, in decreasing importance:
1. **Isotropy in the CM frame** of the whole final state. It is what the lowest-partial-wave formulas say for the two-body process; extending it to the inclusive final state is a step. If the real distribution were forward-peaked, all energies drop: forcing cos θ_cm = 0.9 lowers them tenfold. This is the point to settle with Daniele and Pedro.
2. Invariant mass of the light system (Q up to 1 GeV): ±30 % on the medians.
3. Branching ratios between channels are the slow-monopole literature values; the suppression (m_p/m_M)² is listed but not applied; channels are not combined.
Not done at all: cross sections, rates, fluxes, oscillations, detector response, the deep-inelastic case Daniele simulated with Pythia (the same frames apply to it).

---

## 8. Glossary (the symbols you will see)

E_M monopole kinetic energy · m_M monopole mass · γ_M = 1 + E_M/m_M · t, m_t the target (p or e⁻) at rest and its mass · s, √s invariant mass squared / invariant mass of the collision · γ_cm, β_cm boost of the centre-of-mass frame · W the bundle of light products, m_W its mass, Q its extra internal energy · E_W*, p* its CM energy and momentum · θ_cm the monopole scattering angle in the CM frame · E_W^lab its lab energy · ξ = (m_p/m_M)² Khoze's suppression · q_J half the magnetic charge (1/2 = minimal monopole, isotropic) · κ_ν neutrinos per catalysis · gen 0 = prompt, ≥1 = from a chain · E_ν⁽⁰⁾ the draft's 30–50 MeV · E_KM3 = 220 PeV.

Particles: p proton, e⁻ electron (targets); p̄ antiproton, e⁺ positron (stable products); π⁺ π⁻ π⁰ pions; μ± muons; K⁰ (= half K_S + half K_L), K± kaons; ρ⁰ η ω light mesons that decay to pions; ν_e ν̄_e ν_μ ν̄_μ the four neutrino species; γ photon; M the monopole.

---

## 9. Short FAQ

**Why does the answer not depend on the monopole mass?** Because √s ≈ √(2E_M m_t) as long as m_M is small compared with that, and all lab energies follow from √s and γ_cm = E_M/√s. The mass only re-enters when it approaches √(2E_M m_t) (4 × 10⁵ GeV on a proton, 10⁴ GeV on an electron at 10¹¹ GeV).

**Why half the monopole energy?** The mean of a uniform distribution between 0 and E_M. Uniform because the amplitude is isotropic in the CM frame and the lab energy is linear in cos θ_cm.

**Why does the draft's Eq. 5 fail?** It assumes the products are at rest in the monopole frame, which requires the incoming proton to have negligible energy there — true for β_M → 0, false when the proton carries 10⁹ GeV and the monopole weighs 10³.

**Why no neutrinos from p → e⁺π⁰?** π⁰ decays to two photons and the positron is stable. It gives gamma rays instead (relevant to the KASCADE-Grande comparison in the draft).

**What is "prompt"?** A neutrino produced directly at the catalysis vertex (p → ν̄ₑπ⁺), not from a later decay. It is the only one that survives when the medium absorbs the pions and muons before they decay (atmosphere, rock). Branching ratio ~10⁻⁴.

**Why are the electron channels different for heavy monopoles?** On an electron √s ≈ √(2E_M m_e) ≈ 10 TeV at 10¹¹ GeV; a 10⁴–10⁵ GeV monopole is then heavier than √s, the CM frame rides with the monopole, and the products get γ_M² m_e instead of E_M/2.

**What should I answer Daniele?** The kinematics channel by channel are in `kinematics_note.md`; the one-paragraph version is Sect. 10 of `explanation_2026-11-10.md`; the email text is in the chat of 10 Nov. The message: Eq. 5 is the slow-monopole formula; with proper frames each neutrino carries 1–3 % of E_M; a 220 PeV neutrino needs E_M ≈ 10¹⁰ GeV for any mass; the result rests on CM isotropy, to be agreed; happy to run Daniele's DIS output through the same chain.
