#!/bin/sh
# 1) DIS slices (Q2 > lo), one file per slice: events/h{11|14}_{k}.npz
#    slices k>=1 are gen2.py runs (25000 events each, seeds 100+k'), slice k=0 is lo=1 GeV^2 with 30000 events:
#    python3 -c "import gen2; gen2.gen(1e11,1.0,30000,SEED,'events/h11_0.npz')"   (same for 1e14)
python3 gen2.py 1e11 25000 11 4 1e2 1e4 1e6 1e8 1e10     # -> h11_0..5  (then renamed to events/h11_1..6)
python3 gen2.py 1e14 25000 14 4 1e2 1e4 1e6 1e8 1e10 1e12 # -> h14_0..6 (then renamed to events/h14_1..7)
# 2) photoproduction (Q2<1, W>5 GeV), 100000 events per mass and energy
for m in 1e3 1e4 1e5; do
  python3 -c "import php_gen; php_gen.gen(1e14,$m,100000,int($m/1e2)+14,'php/php_T14_m${m}.npz')"
  python3 -c "import php_gen; php_gen.gen(1e11,$m,100000,int($m/1e2)+11,'php/php_T11_m${m}.npz')"
done
# 3) combination, decays, spectra, P(E_nu>200 PeV) -> res_comb.json
python3 combine.py        # or run_combine_resume.py (skips keys already in res_comb.json)
python3 plot_comb.py; python3 make_table.py
