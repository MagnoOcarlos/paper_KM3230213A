# Neutrino energies from Callan–Rubakov catalysis by a relativistic monopole

*Carlos Magno Ribeiro da Costa — note for Daniele, Pedro, Michele — 5 Nov 2026*
*Code: `notebooks/neutrino_kinematics.ipynb` + `src/mm_kinematics.py` (repository `paper_KM3230213A`, figures in `figures/`, tables in `results/`).*

## 1. The point

Eq. 5 of the draft takes the catalysis products (π⁺, μ⁺) **at rest in the monopole frame**, E_ν⁽⁰⁾ ≈ 30–50 MeV, and boosts them by γ_M. That is the β_M → 0 picture (monopole stopped in the Sun, Super-K / Hu et al.). For our monopoles it does not apply: in the monopole frame the target proton arrives with γ_M m_p ≈ 10⁹–10¹⁰ GeV, far above m_M = 10²–10⁵ GeV, so the monopole frame is not the centre-of-mass frame and nothing is at rest in it.

Doing the kinematics with invariants (target at rest in the lab, E_M = monopole kinetic energy):

    s = m_M² + m_t² + 2 E_M^tot m_t ,     γ_cm = (E_M^tot + m_t)/√s ,
    M + t → M + W ,  W = light system (e⁺π⁰, μ⁺K⁰, p̄π⁺π⁻, …)
    E_W^lab(θ_cm) = γ_cm ( E_W* − β_cm p* cos θ_cm )

The Callan–Rubakov amplitude is pure lowest partial wave; for q_J = 1/2 it is isotropic in the CM frame (Candela–Khoze–Turner 2504.14918 Eq. 2.3/2.5, Khoze 2405.18689 Eq. 6.5/6.12). Averaging over cos θ_cm:

| regime | ⟨E_W^lab⟩ | distribution |
|---|---|---|
| light MM, m_M ≪ √(2E_M m_t) | **E_M / 2**, independent of m_M | flat between m_W²/(2m_t) (MM goes straight on) and ≈ E_M (MM bounced back) |
| heavy MM, m_M ≫ √(2E_M m_t) | γ_M² m_t ("relativistic mirror") | |
| β_M → 0 | m_t shared among the products | the draft's Eq. 5 |

At E_M = 10¹¹ GeV on a proton √s = 433 TeV and γ_cm = 2.3×10⁵ for every mass in our window. On an electron √s ≈ 10 TeV, so m_M = 10⁴–10⁵ GeV is already in the heavy regime for the e⁻ channels.
Check: for M e⁻ → M p̄ at θ_cm = 0 this gives the 861.4 GeV plateau of Candela et al. (their Fig. 2).

## 2. Channels and decay treatment

| channel | target | neutrinos from | weight for slow MM (literature) | relativistic suppression |
|---|---|---|---|---|
| p → e⁺π⁰ | p | none | ≈ 0.5 of e⁺π̃ | ξ = (m_p/m_M)² |
| p → e⁺(ρ⁰,η,ω) → e⁺π⁺π⁻ | p | π⁺, π⁻ chains | ≈ 0.5 of e⁺π̃ | ξ |
| p → μ⁺K⁰ | p | μ⁺, K_S → ππ, K_L → πℓν | 1/400–1/2 of e⁺π̃ | ξ |
| p → ν̄ₑπ⁺ | p | **prompt ν̄ₑ** + π⁺ chain | ≈ 10⁻⁴ | ξ |
| e⁻ → p̄π⁺π⁻ | e⁻ | π⁺, π⁻ chains (p̄ stable) | Candela Eq. 2.3 | none |
| e⁻ → p̄K⁰ | e⁻ | K⁰ chain (p̄ stable) | Candela Eq. 2.3 | none |

ξ = (m_p/m_M)² is Khoze's geometric factor for β > 0.2 (sea quarks): 10⁻⁴ at 100 GeV, 10⁻¹⁰ at 10⁵ GeV.
W has invariant mass Σmᵢ + Q with Q uniform in [0, 1 GeV] (strong-scale assumption; Candela ignore X altogether) and decays with flat phase space. Then, **in the ISM**, π, K, μ decay in flight with PDG branching ratios (Michel spectrum for μ). **In the atmosphere / rock** nothing decays (Pedro, 2 Nov), so only the prompt ν̄ₑ counts.

## 3. Results at the benchmark E_M = 10¹¹ GeV (m_M = 10³ GeV, all flavours, ISM case)

| channel | ν per catalysis | median E_ν | 10–90 % | mean |
|---|---|---|---|---|
| p → e⁺π⁰ | 0 | – | – | – |
| p → e⁺π⁺π⁻ | 6 | 2.6×10⁹ GeV | 2.4×10⁸ – 1.2×10¹⁰ | 4.6×10⁹ |
| p → μ⁺K⁰ | 6.1 | 2.3×10⁹ | 2.0×10⁸ – 1.2×10¹⁰ | 4.6×10⁹ |
| p → ν̄ₑπ⁺ (all ν) | 4 | 5.2×10⁹ | 4.3×10⁸ – 2.8×10¹⁰ | 1.1×10¹⁰ |
| p → ν̄ₑπ⁺, **prompt ν̄ₑ only** | 1 | 1.6×10¹⁰ | 1.6×10⁹ – 5.3×10¹⁰ | 2.2×10¹⁰ |
| e⁻ → p̄π⁺π⁻ | 6 | 1.2×10⁹ | 1.1×10⁸ – 5.9×10⁹ | 2.3×10⁹ |
| e⁻ → p̄K⁰ | 4.1 | 1.2×10⁹ | 1.1×10⁸ – 6.5×10⁹ | 2.5×10⁹ |

Each neutrino carries ~1–3 % of E_M (the prompt ν̄ₑ ~15 %, because it does not share W's energy with a chain). The proton-channel numbers are the same for m_M = 10²–10⁴ GeV and 5 % lower at 10⁵ GeV. The electron channels drop for m_M ≥ 10⁴ GeV (heavy regime on an electron: E_W ≈ γ_M² m_e): median 6×10⁸ GeV at 10⁴ GeV and 1.2×10⁷ GeV at 10⁵ GeV. Flavour ratios are the usual pion-chain 1 : 1 : 2 : 2 (νₑ : ν̄ₑ : ν_μ : ν̄_μ); the kaon channels add a νₑ excess from K_L.

The draft's Eq. 5 at the same point would give E_ν = γ_M × 30–50 MeV = 3–5×10⁶ GeV for m_M = 10³ GeV: three orders of magnitude too low, and wrongly mass-dependent.

## 4. Monopole energy needed for a 220 PeV neutrino

E_M at which the median (10 %, 90 %) neutrino of the channel has 220 PeV:

| channel | m_M = 10² GeV | 10³ | 10⁴ | 10⁵ |
|---|---|---|---|---|
| p → e⁺π⁺π⁻ | 8.6×10⁹ (1.9×10⁹ – 9.3×10¹⁰) | 8.6×10⁹ | 8.6×10⁹ | 1.2×10¹⁰ |
| p → μ⁺K⁰ | 9.7×10⁹ (1.8×10⁹ – 1.1×10¹¹) | 9.7×10⁹ | 9.8×10⁹ | 1.4×10¹⁰ |
| p → ν̄ₑπ⁺, prompt ν̄ₑ (atmosphere) | 1.4×10⁹ (1.4×10⁸ – 1.4×10¹⁰) | 1.4×10⁹ | 1.4×10⁹ | 1.5×10⁹ |
| e⁻ → p̄π⁺π⁻ | 1.8×10¹⁰ (3.7×10⁹ – 2.0×10¹¹) | 1.9×10¹⁰ | 5.3×10¹⁰ | 4.3×10¹¹ |
| e⁻ → p̄K⁰ | 1.8×10¹⁰ (3.4×10⁹ – 2.0×10¹¹) | 1.9×10¹⁰ | 5.3×10¹⁰ | 4.3×10¹¹ |
| draft Eq. 5 (γ_M = 220 PeV / 40 MeV) | 5.5×10¹¹ | 5.5×10¹² | 5.5×10¹³ | 5.5×10¹⁴ |

So the event needs E_M ≈ 10⁹–10¹⁰ GeV, essentially **independent of the monopole mass** for the proton channels, and ≈ 2×10¹⁰ GeV for the electron channels as long as m_M ≲ 10³ GeV. The condition is on the energy, not on γ_M. (Figure `figures/Enu_vs_EM.png`: solid = this calculation, hatched = draft Eq. 5.)

## 5. What I think changes in the draft (to discuss)

1. Sect. 2.1 / Eq. 5: replace "products at rest in the MM frame, E_ν = γ_M E_ν⁽⁰⁾" by the CM treatment above. E_ν ≈ (0.01–0.03) E_M for light monopoles; the heavy-MM limit E_ν ∝ γ_M² m_t interpolates to the old picture only when β_M → 0.
2. The required monopole energy (≈ 10¹⁰ GeV) is below the Galactic gain of Eq. 6 (~10¹¹ GeV, mass-independent). The mass window is then no longer fixed by the neutrino energy; it is fixed by acceleration (Perri et al. 2025), by the flux bounds of Fig. 1, and by ξ = (m_p/m_M)², which kills the proton channels above ~10³–10⁴ GeV and leaves the electron channels (which need m_M ≲ 10³ GeV to be efficient).
3. The statement "Galactic acceleration alone cannot reach the required boost" (Sect. 2.1, 4, 7) should be re-examined with the numbers in Sect. 4 of this note: at E_M = 10¹¹ GeV the median neutrino is ~2 EeV and 220 PeV sits at the ~10 % quantile.
4. Atmosphere scenario (Sect. 3.2): only the prompt ν̄ₑ of p → ν̄ₑπ⁺ survives, BR ~10⁻⁴ and suppressed by ξ; its energy is ~0.15 E_M.
5. ISM scenario: 4–6 neutrinos per catalysis, each ~1–3 % of E_M; κ_ν of Sect. 2 becomes ≈ 6 for the channels with charged pions (and 0 for e⁺π⁰).

## 6. Assumptions I would like us to agree on

* Isotropy of the inclusive final state in the CM frame at √s ~ 10²–10³ TeV (it is what the lowest-partial-wave formulas say). This is the assumption that matters: forcing cos θ_cm = 0.9 (0.99) instead of isotropic lowers every energy by a factor 10 (100), because the light system only gets energy when the monopole is deflected.
* Invariant mass of the light system W (I use Σm + U(0, 1 GeV)). Taking Q_max = 0.3 or 3 GeV instead moves the medians by about ±30 % (e⁻ → p̄π⁺π⁻: 1.0 / 1.2 / 1.4×10⁹ GeV; p → e⁺π⁺π⁻: 3.2 / 2.5 / 2.1×10⁹ GeV), so this choice is secondary.
* Whether to keep the proton channels at all for m_M ≥ 10⁴ GeV given ξ.
* Daniele's DIS case can go through exactly the same frame bookkeeping (same s, γ_cm); only the final-state model differs.

## 7. Checks done

Energy conservation at the vertex to 10⁻¹⁶; Candela's 861.4 GeV forward p̄ for m_M = 10³ and 10⁹ GeV; heavy-MM limit ⟨E_W⟩ = γ_M² m_p; √s independent of m_M in the light regime; with stopped secondaries the code gives the 29.79 MeV ν_μ line and the Michel spectrum (mean 31.7 MeV, endpoint 52.8 MeV), i.e. the Super-K / Hu et al. regime.
