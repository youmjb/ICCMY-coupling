# Reproducing the paper

All commands are run from the repository root. Simulation scripts skip runs whose output already exists, so a step can be repeated after
an interruption. A full rerun takes many hours on a desktop computer; the analyses and figures alone take minutes once `out/` exists.
Set `PYTHONIOENCODING=utf-8` on Windows.

## 1. Build

```
build.bat            → bin\gisim.exe   (MSVC 2022 x64, OpenMP; Eigen 3 in C:\eigen)
```

## 2. Simulations

Calibrated column and robustness analyses (Fig. 2A–C, Supplement S3–S4, Figs. S3–S8):

```
powershell -File run_v4win.ps1                 → out/v4win
powershell -File run_v4ijp.ps1                 → out/v4ijp
powershell -File paper_R1/run_paper_cal.ps1    → out/paper_cal
powershell -File paper_R1/run_paper_sens.ps1   → out/paper_sens
powershell -File paper_R1/run_paper_anat.ps1   → out/paper_anat   (and run_paper_anat2.ps1)
powershell -File paper_R1/run_paper_rev1.ps1   → out/paper_rev1   (and run_paper_rev1b.ps1, run_paper_rev1c.ps1, run_paper_rev2.ps1)
python paper_R1/run_paper_rev3.py              → out/paper_rev3
```

Models D and C, the original-conductance model and the tests against isolated cells:

```
python paper_R1/icc_hyperpol.py run            → out/paper_icc_hyp   isolated ICC-MY, upstroke and clock (Fig. 4, Tables S6–S7)
FM_MODEL=A python paper_R1/final_model.py run  → out/paper_final     original-conductance model (Table 2, Fig. 3)
FM_MODEL=C python paper_R1/final_model.py run                         model C
FM_MODEL=D python paper_R1/final_model.py run                         model D
python paper_R1/figA_runs.py run               → out/paper_figA      coupling and muscle variants (Fig. 5A, Table S8)
python paper_R1/icc_lowg.py run                → out/paper_lowg      isolated cell for each factor f (Table 3, Table S9)
python paper_R1/lowg_coupled.py run            (and lowg_rest.py, lowg_active.py)   column for each f (Fig. 5A, Table S9)
python paper_R1/gown_coupled.py run            own conductance in the column, ±5 pA (Results, transfers above 1)
python paper_R1/goto_clamp.py run              → out/paper_goto      voltage clamp as in Goto et al. 2004 (Table 3, Fig. 6A, Table S15)
python paper_R1/soce_hold.py run               → out/paper_soce      spontaneous currents at −80 mV (Table S16)
python paper_R1/cluster_sim.py run             → out/paper_cluster   isolated cells and clusters
python paper_R1/cluster_pulse.py run                                  time constant and input resistance (Fig. 6B, Table S10)
python paper_R1/k1_only.py run                 → out/paper_k1        inward rectifier reduced alone
python paper_R1/f_only.py run                  → out/paper_fonly     f with the published upstroke and clock
python paper_R1/alt_paths.py run               → out/paper_alt       direct ICC-MY–FLC coupling
```

(`FM_MODEL=…` is a shell environment variable; in PowerShell use `$env:FM_MODEL="D"; python …`.)

## 3. Analyses

Replace `run` by `analyze` in each command of step 2 (with the same `FM_MODEL`) to write the `paper_R1/*.out.txt` files. The 95%
intervals of the recorded transfers (Table 2) come from `python paper_R1/transfer_intervals.py` (no simulation needed).
`sh paper_R1/run_analyses.sh` reruns the analyses of the calibrated column.

## 4. Figures and supplementary tables

```
FM_MODEL=D python paper_R1/figA_fig12.py   → paper_R1/fig/D_fig1, D_fig2
FM_MODEL=D python paper_R1/figA.py         → D_fig3, D_fig4
FM_MODEL=D python paper_R1/figC_fig5.py    → D_fig5
python paper_R1/figD_fig6.py               → D_fig6
python paper_R1/figs_assemble.py           → figS1, figS2, figS4 (Figs. S1, S2, S8) and fig3–fig6 (Figs. S7, S5, S6, S4)
python paper_R1/figA_figS3.py              → A_figS3
```

In `paper_R1/fig/`, the supplementary figures are named as in the paper: A_figS1 = figS1, A_figS2 = figS2, A_figS3, A_figS4 = fig6,
A_figS5 = fig4, A_figS6 = fig5, A_figS7 = fig3, A_figS8 = figS4.

```
FM_MODEL=D python paper_R1/planA/make_supp_tables.py   → paper_R1/planA/supp_tables_D.md (Tables S5–S10)
python paper_R1/planA/make_supp_S10.py                 → paper_R1/planA/supp_S10_D.md (Tables S14–S16)
```

## 5. Where the numbers come from

| Paper item | File |
|---|---|
| Table 2, Fig. 3 (drug tests) | `final_model.out.txt` (original-conductance), `final_model_C.out.txt`, `final_model_D.out.txt`; intervals `transfer_intervals.out.txt` |
| Fig. 4 (isolated ICC-MY) | `icc_hyperpol.out.txt`, `figA.out.txt` |
| Fig. 5, coupling variants | `figA.out.txt`, `figD_fig5.out.txt`, `alt_paths.out.txt` |
| Own conductance and transfers above 1 | `gown_coupled.out.txt`, `f_only.out.txt` |
| Table 3, Fig. 6 | `goto_clamp.out.txt`, `figD_fig6.out.txt`, `cluster_pulse.out.txt` |
| Input resistance in situ (Fig. 5B) | `figD_fig5.out.txt`, `network_divider.out.txt` |
| Supplement S10 | `goto_clamp.out.txt`, `soce_hold.out.txt`, `k1_only.out.txt` |
| Calibration (Fig. 2, Supplement S3–S4) | `criteria.out.txt`, `input_conductance.out.txt`, `analyze_paper_anat.out.txt`, `figD_fig12.out.txt` |
| Literature data | `kito_data.md` (Supplement S9) |

`digitize_kito2005.py` documents how the pinacidil amplitudes were read from Fig. 9 of Kito et al. (2005); the figure itself is not
distributed (set `KITO2005_RENDER_DIR` to a folder with a rendering of it to rerun the script).
