# Pair functions: status

Spec: `src/work/hrg_pair_density_spec.md`. Reports: `reports/pair_functions_M0_M2/`,
`reports/pair_functions_M3/`, `reports/pair_functions_M6/`. Run directory: `build/t2/` (holds the
Integrator cache `fmultiplicity_E33B55BA.txt`, `--samples 1000000`); set `THERMINATOR_SHARE`.

| Milestone | Status | Report |
|---|---|---|
| M0–M2 | done, all tests pass | `reports/pair_functions_M0_M2/report.pdf` |
| M3 | done (T3.1 pass) | `reports/pair_functions_M3/report.pdf` |
| M4 | not started | |
| M6 | done, T6.1–T6.4 pass, both weak settings | `reports/pair_functions_M6/report.pdf` |
| M7 | done, T7.1–T7.4 pass, both weak settings | `reports/pair_functions_M7/report.pdf` |
| M8 | not started | |

## Programs

- `ex_pair_functions_table.exe` — M1 table and T1.x tests; `--level L1` (M7) uses the parent
  currents (`ToL1`, `CKCu3Current`) and runs T7.1–T7.3 instead; `--norm quad|mc` for the currents.
- `ex_pair_functions_mc.exe` — event MC, pair histograms at L0/L1x/L2x, T0.1. `--spectra` (M6) adds
  `hSpeciesN`, `hSpeciesN2` (per-event count and count² of every species, all particles, decayed or
  not; bins = ParticleDB index) and `hSpec_<pdg>` (dN/dp, 0–3 GeV, 10 MeV) for
  211, 321, 2212, 113, 313, 223, 3122, 2224.
- `ex_pair_functions_compare.exe` — M2/M3 comparisons; `--level L1` (M7) compares the L1 table with
  the cumulative MC level L1 (T7.4rho = T2.1, T7.4kap = T2.2 criteria; M3/T2.3 sections skipped).
- `ex_pair_functions_currents.exe` — M6 single-particle currents (`BuildCurrents` in
  `PairFunctions.h`) and T6.1–T6.4. Options: `--weak on|off`, `--stat qs|boltzmann`,
  `--norm quad|mc`, `--kcut` (grid cut / T, default 100), `--cutcheck` (second cut for T6.2,
  default 40, 0 = skip), `--threads`, `--mc <file>` (T6.4; use with `--norm mc`).

## M6 (26 October 2026)

```bash
cd build/t2
../work/ex_pair_functions_currents.exe --weak on  --output m6/currents_weakon_quad.root
../work/ex_pair_functions_currents.exe --weak off --output m6/currents_weakoff_quad.root
for s in 1 2 3 4; do ../work/ex_pair_functions_mc.exe --spectra --weak on --events 5000 \
  --seed 60$s --output m6/jobs/spectra_weakon_s$s.root & done   # same for --weak off
hadd m6/spectra_weakon.root m6/jobs/spectra_weakon_s*.root
../work/ex_pair_functions_currents.exe --weak on --norm mc --cutcheck 0 \
  --mc m6/spectra_weakon.root --output m6/currents_weakon_mc.root
```

Results (weak on; weak off the same or similar):
- T6.1: A vs F_1 6.2e-13, A+B vs F_1+F_u 1.5e-11 (nodes with E- > E_max - 50T excluded) — PASS
- T6.2: n_grid vs counting 1.5e-8 (photon); cut 40T vs 100T 2.8e-10 — PASS
- T6.3: no non-finite / non-positive / spacelike hadron nodes; k = 1e-6, 1e-8 finite — PASS
- T6.4: 2e4 events, species pulls rms 1.01 (weak off 0.97), max 2.8 (3.7); dN/dp χ² p ≥ 0.035 — PASS

Deviations from the spec:
- The two-body rule is integrated in the parent energy E (linear in w), with panels on the thermal
  scale starting at E-. A fixed w rule fails at large k. For narrow windows the y = E − E_c form is
  used on the exactly symmetric window to avoid the 1/k³ cancellation.
- Output TTree `currents` stores the total A and A+B (thermal + feed-down) on the GL nodes.
- Weak-on and weak-off MC samples share seeds 601–604, so they are not statistically independent.

## M7 (26 October 2026)

```bash
cd build/t2
../work/ex_pair_functions_table.exe --level L1 --weak on --output m7/pair_table_L1_weakon.root
../work/ex_pair_functions_compare.exe --level L1 --mc mc_weakon.root \
  --output m7/compare_L1_weakon.root          # same for --weak off / mc_weakoff.root
```

Results (weak on / off):
- T7.1: thermal-only currents vs M1 table: continua 5.9e-14, poles 4.4e-16 — PASS
- T7.2: rho(M) vs Dalitz with n_c 6.2e-14; poles carry B n_c 3.8e-15; int dM rho = c^L1 5.6e-14 — PASS
- T7.3: no non-finite / negative / spacelike cells, C_u > 0; small kappa finite and smooth — PASS
- T7.4: rho(M) all ordered pairs, min p 0.030 / 0.217; dP/dkappa at 9 masses, min p 0.109 / 0.028
  — PASS
- c^L1/c^L0: pi+pi- 1.448 (weak off 1.341), K+K- 1.009, like-sign 1.000.

Notes:
- L1 continuum C_u = L0 expression with f -> A_c plus the B_c term, so thermal-only currents
  reproduce the M1 arithmetic. rho_ch(M) at L1 is the Dalitz projection with n_c, so the M-binned
  predictions reuse `ContDensityM` with `na = n_c`.
- `--level L0` output of both programs is unchanged (compare log identical except its header; the
  table's T1.5 varies run to run because `--seed 0` gives a random ROOT seed).
