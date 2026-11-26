# paper_KM3230213A

Neutrino kinematics for the paper *Can Magnetic Monopoles Explain the KM3-230213A Event?*
(M. Doro, D. Perri, P. de la Torre Luque, C. Costa). The draft itself lives on Overleaf
("MM 200 PeV nu"); this repository holds the calculations, notes and figures behind it.

## Layout

```
src/                     Python code
  mm_kinematics.py         Callan–Rubakov catalysis kinematics: exact two-body vertex (lab → CM → lab),
                           light-system decay, pion/kaon/muon decay chains, neutrino spectra vs E_M and m_M
  pythia_catalysis.py      catalysis vertex fed to Pythia 8 standalone hadronization and decays
  compare_dis_catalysis.py comparison of Daniele's electromagnetic Pythia run with the catalysis kinematics
notebooks/
  neutrino_kinematics.ipynb  scans and figures behind notes/kinematics_note.md
pythia_dis/              electromagnetic scattering (not catalysis) of a heavy charged projectile on protons
  MonopoleDIP/             Daniele Perri's Pythia 8 code and results, as distributed
  MonopoleDIP_rerun/       driver, results and plots of the independent rerun
notes/
  kinematics_note.md           neutrino energy per channel at E_M = 1e11 GeV; proposed replacement of Eq. 5
  pythia_replication_note.md   replication of the Pythia run and comparison with catalysis
  pythia_catalysis_note.md     what Pythia can and cannot do for catalysis; results with its hadronization
  explanations/                step-by-step explanations of the calculation
results/                 CSV / JSON tables (benchmark at 1e11 GeV, scans, required E_M for 220 PeV, comparisons)
figures/                 all figures referenced in the notes
private/                 draft PDF and email thread; ignored by git
```

## Setup and use

```
python3 -m venv pyenv
./pyenv/bin/pip install -r requirements.txt
```

| task | command | time |
|---|---|---|
| catalysis tables and figures | `cd notebooks && jupyter notebook neutrino_kinematics.ipynb` | seconds |
| catalysis through Pythia | `./pyenv/bin/python src/pythia_catalysis.py` | ~2 min |
| regenerate the electromagnetic run | see `notes/pythia_replication_note.md`, Sect. 8 | ~25 min on 10 cores |
| compare the two | `./pyenv/bin/python src/compare_dis_catalysis.py` | ~1 min |

Scripts write to `figures/` and `results/` relative to the repository root, whatever the current directory.
The Pythia event files (about 800 MB) are not stored here.

## Main results (E_M = 1e11 GeV)

- Catalysis on protons gives neutrinos with median energy ≈ 2 % of E_M, independent of m_M for m_M ≤ 100 TeV;
  E_M ≈ 1e10 GeV reproduces the 220 PeV event. Electron channels need E_M ≈ 2e10 GeV and m_M ≤ 1 TeV.
- The electromagnetic (photon-exchange) collisions transfer a fraction ~1e-7 of E_M to the proton:
  many soft pions, essentially nothing above 220 PeV at 1e11 GeV.
- Pythia's hadronization and decay tables reproduce the phase-space catalysis spectra.
