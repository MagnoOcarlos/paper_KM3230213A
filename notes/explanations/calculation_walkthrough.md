# Walkthrough of the calculation, step by step, with the source of every ingredient

*For Carlos, to check the chain of reasoning in `mm_kinematics.py` / `neutrino_kinematics.ipynb`.*
Line numbers refer to `mm_kinematics.py`. Papers:
**[C]** Candela, Khoze, Turner, arXiv:2504.14918 (*Monopoles at Future Neutrino Detectors*);
**[K]** Khoze, arXiv:2405.18689 (*Monopoles and Fermions in the Standard Model*);
**[H]** Hu, Cheng, Guo, Wang, arXiv:2201.02386 (solar monopole neutrinos);
**[SK]** Super-Kamiokande, arXiv:1203.0940; **[D]** our draft *MM 200 PeV nu*; **[PDG]** PDG 2024.
Everything marked **standard** is textbook special relativity / particle kinematics (e.g. PDG kinematics review, Byckling & Kajantie) and not a model choice.

---

## Step 0. The question and the set-up

Lab frame: target at rest (proton or electron of the medium), monopole with kinetic energy $E_M$ and mass $m_M$.
Reaction: $M + t \to M + W$, then $W \to$ the listed products, then the products decay.

* Which reactions: [C] Eq. 2.1 $M+p\to M+e^+ +X$ and Eq. 2.2 $M+e\to M+\bar p+X$, with $X$ = pions ($\pi^0$, $\pi^+\pi^-$) or kaons/muons when other families are included ([C] text after Eq. 2.2). The proton channels with explicit mesons are [SK] Eq. 1 $p\to(\rho^0,\omega,\eta,K^+,\dots)+e^+(\mu^+)$ and [H] Eq. 3.2 ($p\to e^+\tilde\pi$, $p\to\mu^+K^0$, $p\to\bar\nu_e\pi^+$). Daniele's email: for the antiproton side use $e^-\to\bar p\pi^+\pi^-$ and $e^-\to\bar pK^0$.
* Code: the channel list `CHANNELS`, line 69.

---

## Step 1. Invariant mass of the collision  — **standard**

$$ s = (p_M + p_t)^2 = m_M^2 + m_t^2 + 2E_M^{\rm tot}m_t ,\qquad E_M^{\rm tot}=E_M+m_M $$

* This is [C] Eq. 2.6 written for a general target (they write it for the electron: $s = 2m_eE_M^{\rm lab}+m_e^2+m_M^2$).
* Code: `cm_frame`, line 93. I compute $s-m_M^2 = 2E^{\rm tot}m_t+m_t^2$ separately (`s_minus_mM2`) because for $m_M=10^{16}$ GeV the subtraction $s-m_M^2$ would lose all digits in double precision.
* Hand check at $E_M=10^{11}$ GeV, proton target, $m_M=10^3$: $s \simeq 2\times10^{11}\times0.938 = 1.88\times10^{11}$ GeV², $\sqrt s = 4.33\times10^5$ GeV. The $m_M^2=10^6$ term is negligible. This is why the result does not depend on $m_M$ below $\sim10^5$ GeV.

## Step 2. Boost from the lab to the centre-of-mass (CM) frame — **standard**

$$ \beta_{cm} = \frac{p_M^{\rm lab}}{E_M^{\rm tot}+m_t},\qquad \gamma_{cm}=\frac{E_M^{\rm tot}+m_t}{\sqrt s} $$

* [C] Eq. A.9 for $\beta_{cm}$. $\gamma_{cm}$ = total lab energy / invariant mass.
* Code: `cm_frame`, lines 93–107. `one_minus_beta = s/(E_tot (E_tot+p))` is the exact identity $1-\beta = (E_{\rm tot}-p)/E_{\rm tot}$ with $E_{\rm tot}-p = s/(E_{\rm tot}+p)$, used to avoid computing $1-\beta$ as a difference of two numbers equal to 1 within $10^{-11}$.
* Hand check: $\gamma_{cm} = 10^{11}/4.33\times10^5 = 2.3\times10^5$.

## Step 3. Two-body kinematics $M + t \to M + W$ in the CM frame — **standard**

With $W$ of invariant mass $m_W$:
$$ E_W^* = \frac{s + m_W^2 - m_M^2}{2\sqrt s},\qquad p^* = \frac{\lambda^{1/2}(s,m_M^2,m_W^2)}{2\sqrt s},\qquad \lambda(x,y,z)=(x-y-z)^2-4yz $$

* [C] Eq. A.2 (momenta via the Källén function), Eq. A.3 (definition of $\lambda$), Eq. A.11 (energy of the light particle in the CM, written there for $\bar p$).
* Code: `kallen` line 83, `pstar` line 87, `vertex_cm` line 110 (again written with $s-m_M^2$ to keep precision).
* Hand check (light MM): $E_W^* \simeq p^* \simeq \sqrt s/2 = 2.2\times10^5$ GeV.

## Step 4. The angular distribution in the CM frame — **the key physics input**

Callan–Rubakov scattering happens only in the lowest partial wave $j_0=|q_J|-1/2$ ([C] p. 4, "the Callan–Rubakov processes occur only in the lowest available partial wave"; [K] Sect. 4–6). The CM differential cross sections are

* $Me^-\to M\bar p$: [C] Eq. 2.3 / [K] Eq. 6.5, $\;d\sigma/d\Omega \propto q_J^2/|p_e^{cm}|^2\,[\sin(\theta_{cm}/2)]^{4|q_J|-2}$;
* $Mp\to Me^+X$: [C] Eq. 2.5 / [K] Eq. 6.12–6.13, same angular factor.

For the minimal monopole $q_J=1/2$ the exponent $4|q_J|-2=0$, so **$d\sigma/d\Omega$ is independent of $\theta_{cm}$** ([C] p. 8: "In this instance, the cross-section of Eq. (2.3) is independent of the scattering angle"; [K] p. 19: "the dependence on $\sin(\theta/2)$ disappears"). Hence I draw $\cos\theta_{cm}$ uniformly in $[-1,1]$.

* Code: line 312 `cth = rng.uniform(-1, 1, n)`; option `cos_theta_cm=` to force a fixed angle (used for the checks and the sensitivity toy).
* **What to scrutinise:** the formulas are for the two-body final states ($\bar p$ alone, $e^+$ alone, "inclusive" in the sense of [K] Sect. 6). I apply the same isotropy to the whole light system $W$. This is the assumption that drives the result (Sect. 6 of the notebook: forcing $\cos\theta_{cm}=0.9$ lowers all energies by 10). It should be discussed with Daniele and Pedro.
* Also note [C] Sect. 2.3, p. 8: "when $\cos\theta_{cm}\neq1$ ... altering the monopole's trajectory requires the antiproton to carry a substantial amount of momentum" — that sentence is the physics of everything that follows.

## Step 5. Energy of $W$ back in the lab — **standard**

The MM is scattered at $\theta_{cm}$ with respect to the boost axis, $W$ goes opposite, so $p_{z,W}^* = -p^*\cos\theta_{cm}$ and

$$ E_W^{\rm lab} = \gamma_{cm}\left(E_W^* + \beta_{cm}\,p^*_{z,W}\right) = \gamma_{cm}\left(E_W^* - \beta_{cm}p^*\cos\theta_{cm}\right). $$

* Code: lines 314–317 and `lab_energy` line 120. For $p_z<0$ the bracket is a difference of two nearly equal numbers; line 127–128 rewrite it exactly as $(E^*-|p_z|) + (1-\beta)|p_z|$ with $E^*-|p_z| = (m^2+p_T^2)/(E^*+|p_z|)$ — same numbers, no cancellation.
* Consequences you can verify by hand (notebook Sect. 1 and `analytic_limits` line 376):
  * isotropic average: $\langle E_W^{\rm lab}\rangle=\gamma_{cm}E_W^*$. Light MM: $=\gamma_{cm}\sqrt s/2 = (E_M^{\rm tot}+m_t)/2\simeq E_M/2$. The distribution in $E_W^{\rm lab}$ is **flat** because it is linear in $\cos\theta_{cm}$ and $\cos\theta_{cm}$ is uniform.
  * $\theta_{cm}=0$ (MM goes straight on): $E_W^{\rm lab}\to\gamma_{cm}\,m_W^2/(E_W^*+p^*)\simeq m_W^2/(2m_t)$. For $m_W=m_p$, $m_t=m_e$: $0.880/0.001022 = 861$ GeV = the plateau of [C] Fig. 2 (right) and p. 9 ("plateaus at around 861.4 GeV"). Reproduced numerically in notebook Sect. 7, check 2.
  * $\theta_{cm}=\pi$: $E_W^{\rm lab}\simeq E_M$ (the MM is bounced back and loses everything).
  * heavy MM ($m_M\gg\sqrt{2E_Mm_t}$): $\gamma_{cm}\to\gamma_M$, $E_W^*\to\gamma_Mm_t$, so $\langle E_W^{\rm lab}\rangle=\gamma_M^2m_t$ (check 3 in the notebook, $\gamma_M=11$: 113.5 GeV).
  * $\beta_M\to0$: $E_W^*\to m_t$, products share the target rest energy: this, and only this, is the regime of [D] Eq. 5 and of [SK]/[H].

## Step 6. Invariant mass of the light system $W$ — **model choice**

$m_W=\sum_i m_i + Q$, $Q$ uniform in $[0,Q_{\max}]$, default $Q_{\max}=1$ GeV, clipped to the kinematic limit $\sqrt s-m_M$.

* Motivation: Callan–Rubakov is "governed by the strong interactions scale ∼1 GeV" ([C] p. 4) and [C] p. 5 argue that the cross section "receive[s] dominant contributions at energies not much above the kinematic threshold" so they "ignore contributions of $X$". There is no calculation of the $X$ spectrum in the literature; $Q\lesssim1$ GeV is my reading of "strong scale".
* Code: lines 305–310.
* Effect: notebook Sect. 6, $Q_{\max}=0.3/1/3$ GeV changes the medians by ≈ ±30 %. It only changes how $E_W$ is shared among the constituents, not $E_W$ itself (Step 5 is independent of $m_W$ except in the $\theta_{cm}\to0$ corner).

## Step 7. Decay of $W$ into its constituents — **standard**

Flat (Lorentz-invariant) phase space in the $W$ rest frame, isotropic orientation, then boost to the lab with $\beta_W,\gamma_W$ from $E_W^{\rm lab}$ and $m_W$.

* Two-body: momenta from Step 3's formula with $M\to m_W$ (`two_body_rest`, line 153).
* Three-body: $d\Phi_3\propto p^*_{(12)3}\,p^*_{12}\,dm_{12}$ sampled by rejection, then the standard sequential construction $W\to(12)+3$, $(12)\to1+2$ (`three_body_rest`, line 163). Reference: PDG kinematics review, Sect. "three-body decays" (Dalitz plot uniform in $m_{12}^2,m_{23}^2$ ⇔ this weight in $m_{12}$).
* Boost: `parent_boost` line 133 gives $\beta=p/E$, $\gamma=E/m$, $1-\beta=m^2/(E(E+p))$ (exact); lab energy via `lab_energy` using only the component of the daughter momentum along the parent direction. Because each decay is isotropic in its parent's rest frame, no directions need to be kept (lines 324–332).

## Step 8. Decay chains to neutrinos — **standard, PDG numbers**

Table `DECAYS`, line 47, with [PDG] branching ratios:
$\pi^\pm\to\mu\nu$ (100 %); $\pi^0\to\gamma\gamma$ (no neutrinos); $K^\pm$: $\mu\nu$ 63.56 %, $\pi\pi^0$ 20.67 %, $3\pi$ 5.58 %+1.76 %, $\pi^0e\nu$ 5.07 %, $\pi^0\mu\nu$ 3.35 %; $K^0\to K_S,K_L$ 50/50; $K_S$: $\pi^+\pi^-$ 69.2 %, $\pi^0\pi^0$ 30.7 %; $K_L$: $\pi e\nu$ 40.55 %, $\pi\mu\nu$ 27.04 %, $3\pi^0$ 19.52 %, $\pi^+\pi^-\pi^0$ 12.54 %.

* Muon decay uses the V–A matrix element: $|{\cal M}|^2\propto(p_\mu\!\cdot\!p_{\nu_e})(p_e\!\cdot\!p_{\bar\nu_\mu})$ for $\mu^+$, which in the muon rest frame is the weight $x_{\nu_e}(1-x_{\nu_e})$ on flat massless phase space (`michel_rest`, line 203, weight at line 215). Resulting spectra: $\nu_e\propto x^2(1-x)$ (mean $0.30\,m_\mu=31.7$ MeV), $\bar\nu_\mu,e^+\propto x^2(3-2x)$ (mean $0.35\,m_\mu$), endpoint $m_\mu/2=52.8$ MeV — the standard Michel result (PDG muon decay parameters review). Verified in notebook Sect. 7, check 5.
* Pion at rest: $E_\nu=(m_\pi^2-m_\mu^2)/2m_\pi=29.79$ MeV ([D] Sect. 2.1, [SK] Fig. 1 caption). Verified, check 5.
* All decays **in flight** (ISM): [H] Sect. 3 and [SK] footnote 1 explain that the lines of the solar case exist only because pions stop in dense matter before decaying; in the ISM the interaction length is ~Mpc and every pion/kaon/muon decays first. Option `decay_at_rest=True` (line 248) reproduces the stopped case for checking.
* Atmosphere / rock: Pedro (thread, 2 Nov 2026) computed that $\pi,K,\mu$ at these energies do not decay before being absorbed, so only neutrinos produced **at the vertex** count. Those are tagged `gen = 0`; `summarize(..., prompt_only=True)` (line 352) selects them. The only such channel is $p\to\bar\nu_e\pi^+$ with branching $\sim10^{-4}$ ([H] Eq. 3.2, from Houston–Li–Sun 1803.02835; Daniele's message of 19 Oct).
* $\bar p$ is stable; it would annihilate only after ~Mpc in the ISM, so it gives no neutrinos here (it is the AMS-02 antiproton handle Pedro mentioned on 14–16 Oct).

## Step 9. Weights between channels — **taken from the literature, not recomputed**

* $B(p\to e^+\tilde\pi\to\pi^+)\approx0.5$, $B(p\to\mu^+K^0)\approx(m_d/m_s)^2\sim1/400$–$1/2$, $B(p\to\bar\nu_e\pi^+)\approx(m_p/m_W)^2\sim10^{-4}$: [H] Eq. 3.2 and text (slow monopoles).
* Relativistic suppression of the proton channels: [K] Sect. 6.2, Eq. 6.11 ($\beta<0.2$ needed for the proton to fall on the core as a whole) and p. 21–22: for faster protons the process goes through sea quarks with "the geometric suppression factor $\xi=(R_M/R_p)^2\sim(m_p/M_M)^2$". The electron channel has no such factor ([K] Sect. 6.1; Daniele's message of 15 Oct).
* I report per-channel results and list $\xi$ in the channel table; I did **not** combine channels into one spectrum, because the relativistic branching ratios are not known.

## Step 10. What the numbers mean — reading the outputs

* `results/benchmark_E1e11.csv`: at $E_M=10^{11}$ GeV each neutrino carries ~1–3 % of $E_M$ (median $1$–$3\times10^9$ GeV), 4–6 neutrinos per catalysis for any channel with a charged pion or a kaon, 0 for $e^+\pi^0$.
* `results/required_EM_for_220PeV.csv`: the $E_M$ at which the median neutrino has 220 PeV: $\approx 9\times10^9$ GeV for the proton channels, independent of $m_M$; $\approx2\times10^{10}$ GeV for the electron channels when $m_M\lesssim10^3$ GeV, rising to $4\times10^{11}$ at $10^5$ GeV (heavy regime on an electron, Step 5).
* Compare [D] Eq. 5: $\gamma_M=E_\nu/E_\nu^{(0)}\sim10^9$–$10^{10}$, i.e. $E_M=\gamma_Mm_M\sim 5\times10^{11}(m_M/100\,{\rm GeV})$ GeV. That formula is Step 5's $\beta_M\to0$ corner applied outside its validity.

## Step 11. Checks you can rerun (notebook Sect. 7)

1. $|E_M'+E_W-E_{\rm tot}|/E_{\rm tot}<10^{-15}$ event by event.
2. [C] Fig. 2: 861.4 GeV forward antiproton for $m_M=10^3$ and $10^9$ GeV.
3. $\langle E_W\rangle=\gamma_M^2m_p$ for a heavy monopole.
4. $\sqrt s$ equal for $m_M=10^2$–$10^4$ GeV.
5. Stopped secondaries: 29.79 MeV line and Michel spectrum.

## What is *not* in the calculation

Cross sections and rates (Candela Eq. 2.3/2.5 give them; Daniele/Pedro are handling $\sigma_{\rm cat}$ and fluxes), neutrino oscillations, detector effects, the deep-inelastic channel Daniele simulated with Pythia (same Steps 1–2 and 5 apply; only the final-state model differs), and any statement about whether catalysis itself survives at $\sqrt s\sim10^2$–$10^3$ TeV beyond [K]'s geometric argument.
