# Daniele's Pythia run (MonopoleDIP) replicated and compared with the catalysis kinematics

Carlos Costa, 21 November 2026. Files: `pythia_dis/MonopoleDIP/` (Daniele's original), `pythia_dis/MonopoleDIP_rerun/` (my rerun), `src/compare_dis_catalysis.py`, `figures/compare_DIS_vs_catalysis_E1e11.png`, `figures/energy_transfer_DIS_vs_catalysis.png`, `figures/decay_chain_crosscheck.png`, `results/compare_DIS_vs_catalysis_E1e11.csv`.

## 1. What the MonopoleDIP code actually computes

It is **not** Callan–Rubakov catalysis. It is the ordinary electromagnetic scattering of a heavy charged projectile on a proton at rest, with Pythia 8.318:

| piece | what it does | script |
|---|---|---|
| DIS, Q² > 1 GeV² | one-photon exchange ff → ff, muon beam used as a massless carrier, heavy mass and spin put back by an analytic reweighting of the photon flux, seven Q²-stratified slices | `gen2.py`, weights in `combine.py` |
| photoproduction, Q² < 1 GeV² | equivalent-photon flux of a "tau" given the real mass m, SoftQCD events with W > 5 GeV | `php_gen.py` |
| charge | z = 60 times the electron charge, cross section scaled by z² = 3600 | `combine.py` |
| decays | every charged pion → μν → eνν in the collinear limit with muon polarisation, three neutrinos per pion | `nu.py` |
| output | E dσ/dE spectra of π± and ν, σ_tot, number of neutrinos above 200 PeV per collision | `results/res_comb.json` |

So it answers: "if the monopole behaves like a very heavy, highly charged particle interacting through photon exchange, what neutrinos come out of its collisions with interstellar protons?" The magnetic charge g = 68.5 e is approximated by an electric charge z = 60 (for a relativistic monopole the field seen by the proton is that of a charge g β). It does not contain the baryon-number-violating vertex, the (m_p/m_M)² Khoze factor, or the electron channels.

## 2. Replication

Setup: Python 3.11 virtual environment, `pip install pythia8mc==8.318.3` (same Pythia version Daniele used, binary wheel, no compilation), numpy, scipy, matplotlib. Daniele's folder contains the scripts and the final numbers but not the events (about 800 MB), so I regenerated them with his commands from `run_generation.sh` at the original statistics (25 000–30 000 events per DIS slice, 100 000 photoproduction events per mass and energy), with the same seeds, using the driver `MonopoleDIP_rerun/drive_gen.py` on ten cores. Generation took 11 minutes, `combine.py` another 15.

Result (spin ½, Daniele's `summary_table.md` against mine):

| m | T [GeV] | σ_tot [mb] D / mine | ⟨N_π⟩ D / mine | N_ν(>200 PeV) per collision D / mine |
|---|---|---|---|---|
| 1 TeV | 10¹¹ | 423.5 / 423.1 | 17.2 / 17.2 | 1.22e-4 ± 0.05 / 1.21e-4 ± 0.05 |
| 10 TeV | 10¹¹ | 252.8 / 252.5 | 13.7 / 13.6 | 6.2e-7 ± 1.5 / 4.4e-7 ± 0.4 |
| 100 TeV | 10¹¹ | 150.2 / 150.2 | 10.6 / 10.6 | 1.05e-9 ± 0.10 / 0.91e-9 ± 0.10 |
| 1 TeV | 10¹⁴ | 1106 / 1107 | 35.3 / 35.3 | 0.952 ± 0.009 / 0.964 ± 0.009 |
| 10 TeV | 10¹⁴ | 837.5 / 835.8 | 27.5 / 27.5 | 0.236 ± 0.003 / 0.242 ± 0.003 |
| 100 TeV | 10¹⁴ | 601 / 598 | 22.1 / 22.2 | 0.0155 ± 0.0004 / 0.0155 ± 0.0004 |

The photoproduction cross sections in the logs agree digit for digit (87.158, 59.374, 38.354 µb at 10¹¹ GeV; 213.99, 160.08, 115.45 µb at 10¹⁴ GeV), the DIS cross sections to better than 1 %, and the tail numbers within the quoted statistical errors. The only visible difference, the 10 TeV tail at 10¹¹ GeV, is driven by a handful of events in the highest-Q² slice and is covered by Daniele's own error bar. The code is reproducible as distributed.

## 3. Comparison with the catalysis kinematics at T = E_M = 10¹¹ GeV

Both results are expressed per event: per collision (normalised to σ_tot) for the Pythia run, per catalysis for ours. `figures/compare_DIS_vs_catalysis_E1e11.png`, top row neutrinos, bottom row charged pions, one column per mass.

| process | m_M | ν per event | median E_ν [GeV] | N_ν(>220 PeV) per event |
|---|---|---|---|---|
| EM DIS + photoproduction (Pythia, z = 60) | 1 TeV | 41 | 44 | 1.2e-4 |
| | 10 TeV | 30 | 17 | 4.4e-7 |
| | 100 TeV | 21 | 7.5 | 0.9e-9 |
| catalysis p → e⁺π⁺π⁻ | 1–100 TeV | 6.0 | 2.5e9 | 5.5 |
| catalysis p → μ⁺K⁰ | 1–100 TeV | 6.1 | 2.2e9 | 5.5 |
| catalysis p → ν̄ₑπ⁺ (1 prompt + 3 from π⁺) | 1–100 TeV | 4.0 | 5.2e9 | 3.8 |
| catalysis e⁻ → p̄π⁺π⁻ | 1 / 10 / 100 TeV | 6.0 | 1.2e9 / 6.2e8 / 1.2e7 | 5.1 / 4.5 / 0.02 |
| catalysis e⁻ → p̄K⁰ | 1 / 10 / 100 TeV | 4.1 | 1.2e9 / 6.2e8 / 1.2e7 | 3.4 / 3.0 / 0.03 |

The two processes give opposite pictures. The electromagnetic collision produces many pions, 10 to 17 per collision, but they are soft: the median neutrino has tens of GeV and only one collision in 10⁴ (1 TeV) to 10⁹ (100 TeV) yields a neutrino above 220 PeV. Catalysis produces 4 to 6 neutrinos per event and almost all of them are above 220 PeV, independently of the mass for the proton channels.

## 4. Why: the energy transferred to the proton side

`figures/energy_transfer_DIS_vs_catalysis.png` shows the distribution of y, the fraction of the monopole energy given to the light system, event by event from the regenerated Pythia events (weighted as in `combine.py`) and from our catalysis.

| | median y | fraction of events with y > 10⁻² |
|---|---|---|
| EM, m = 1 TeV | 4e-7 | ~1e-3 |
| EM, m = 10 TeV | 4e-8 | ~1e-4 |
| EM, m = 100 TeV | 8e-9 | ~1e-5 |
| catalysis, proton target, any m ≤ 100 TeV | 0.5 | 0.98 |

Photon exchange is peripheral: the flux goes as dy/y and dQ²/Q², so most collisions transfer a tiny fraction of the energy, and the kinematic floor Q²_min = m² y²/(1 − y) suppresses large y more strongly the heavier the projectile, which is why the EM tail falls with the mass. Catalysis is a contact interaction in the lowest partial wave: for q_J = ½ the centre-of-mass angle is isotropic, the monopole is deflected at large angle, and the proton side receives on average half the energy. In Candela et al. (2504.14918, Sect. 2.3) this is the statement that only the exactly forward configuration gives a soft product.

Daniele's run is therefore the right description of the electromagnetic energy loss and of the soft-pion background of a fast monopole, and our run is the right description of the catalysis neutrinos. They are not two estimates of the same thing.

## 5. Cross-check of the decay chain

Our pions from p → e⁺π⁺π⁻ at E_M = 10¹¹ GeV, m_M = 1 TeV, pushed through Daniele's `nu.py` (collinear, with muon polarisation) and through our full three-dimensional decay chain (`figures/decay_chain_crosscheck.png`):

| | ν per pion | median E_ν [GeV] | mean E_ν [GeV] | fraction above 220 PeV |
|---|---|---|---|---|
| Daniele's `nu.py` | 3 | 2.60e9 | 4.64e9 | 0.913 |
| our `mm_kinematics.py` | 3 | 2.55e9 | 4.57e9 | 0.912 |

The two decay codes agree to 2 %. The only difference is below 10⁴ GeV, where the collinear approximation drops the small transverse tail. This validates the decay stage of both calculations independently.

## 6. Neutrinos above 220 PeV per unit column: putting the cross sections in

Per target proton the effective cross section for a neutrino above 220 PeV is σ × N_ν(>220 PeV):

| process at E_M = 10¹¹ GeV | 1 TeV | 10 TeV | 100 TeV |
|---|---|---|---|
| EM DIS + photoproduction (Pythia) | 0.05 mb | 1.6e-4 mb | 1.6e-7 mb |
| catalysis on protons, unsuppressed, σ_cat = 10⁻²⁷ cm² = 1 mb (draft's reference) | 5.5 mb | 5.5 mb | 5.5 mb |
| same with the Khoze factor (m_p/m_M)² | 5e-6 mb | 5e-8 mb | 5e-10 mb |
| catalysis on electrons, σ_cat = 1 mb, no suppression | ~5 mb | ~4 mb | ~0.02 mb |

So at Galactic energies: if the Khoze suppression applies, the electromagnetic channel beats catalysis on protons by four orders of magnitude for a 1 TeV monopole, yet both are far below one neutrino per crossing; the electron channels, which are not suppressed, remain the dominant source of 220 PeV neutrinos for m_M ≤ 10 TeV. At E_M = 10¹⁴ GeV the Pythia run gives about one neutrino above 220 PeV per collision for 1 TeV, so the electromagnetic channel becomes efficient there, but that energy needs a thousand Galactic crossings.

## 7. What this means for the draft

- Sect. 2 and 3 should present the two mechanisms as different things: electromagnetic scattering (Daniele's Pythia run, hundreds of mb, soft products, negligible above 220 PeV at 10¹¹ GeV) and catalysis (mb scale, hard products, every event above 220 PeV). The sentence "the monopole produces neutrinos by scattering on protons" needs to say which one.
- Equation 5 stays replaced by the catalysis kinematics of the previous note, E_ν ≈ 0.02 E_M. The Pythia run supports this picture indirectly: the decay chain agrees, and its own low y shows that a photon-mediated process cannot reach the event energy at 10¹¹ GeV.
- With the Khoze suppression the proton catalysis channels are below the electromagnetic channel, so the discussion of Sect. 3.2 and 4 should be made with the electron channels or with the unsuppressed cross section, stating which.
- Open point, as before: the isotropic centre-of-mass angle. If the catalysis final state is forward peaked, y drops and the catalysis spectrum moves towards the electromagnetic one.

## 8. How to rerun

```
python3 -m venv pyenv && ./pyenv/bin/pip install pythia8mc==8.318.3 numpy scipy matplotlib
cd pythia_dis/MonopoleDIP && mkdir -p events php && python ../MonopoleDIP_rerun/drive_gen.py   # 11 min on 10 cores
python combine.py && mv res_comb.json results/ && python make_table.py && python plot_comb.py
cd ../.. && python src/compare_dis_catalysis.py pythia_dis/MonopoleDIP/results/res_comb.json pythia_dis/MonopoleDIP
```

The regenerated events (800 MB) are not kept in the shared folder; `pythia_dis/MonopoleDIP_rerun/results/res_comb.json` holds everything derived from them.
