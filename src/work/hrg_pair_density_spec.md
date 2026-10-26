# Pair functions in a hadron resonance gas at rest: implementation spec

Implementation target: `federicacapellino/therminator2-hboltz` (THERMINATOR 2, hboltz fork).
Physics reference: `local_correlations_short.tex` (LCS) and `invariant_mass.tex` (IM) in this repo.
Equation numbers below are those of the compiled notes.

## 0. Scope and conventions

Goals:

1. Tabulate the pair functions $\mathcal C_K^{ss'}(M^2,\epsilon_K)$, $\mathcal C_u^{ss'}(M^2,\epsilon_K)$ for a gas at rest from the
   analytic formulas (LCS Eqs. 13–20), using the THERMINATOR hadron list.
2. Compare the pair density from the table with the connected pair density of THERMINATOR HRG events, at three levels
   of ancestry, to measure what the table leaves out (fed-down parents, cascades, $Q$ dependence).
3. Check the expectations of IM: the pair mass distribution (IM Eqs. 5–9) and the normalized pair spectrum (IM Eq. 12)
   against the thermal spectrum of a particle of mass $M$ (IM Eq. 13).
4. Check the charge sum rule.
5. (Last) Separate $\mathcal C_K$ from $\mathcal C_u$ with a source whose $d\Sigma$ is not comoving (BlastWave).

Conventions (as in the project instructions):

- Metric mostly plus. $K=p_1+p_2$, $Q=p_1-p_2$, $M^2=-K^2$ (invariant mass of the measured pair 1,2),
  $\epsilon_i=-u\cdot p_i$, $\epsilon_K=-u\cdot K$, $\epsilon_Q=-u\cdot Q$, $\kappa=\sqrt{\epsilon_K^2-M^2}$.
- Parent mass $m_a$, daughter masses $m_1,m_2$ (measured), $m_3$ (unobserved third daughter).
- Constant matrix elements, zero widths, one-step decays in the table.
- Gas at rest: $u=(1,\mathbf 0)$, $d\Sigma_\mu u^\mu=V$, so $\epsilon_K=E_K$ and $\kappa=|\mathbf K|$ in the lab.
- Units: GeV throughout; $V$ in GeV$^{-3}$ internally ($\hbar c=0.19733$ GeV fm).
- No integration over a freeze-out surface except in milestone M5.

## 1. What THERMINATOR does (constraints on the implementation)

Verified in the source of the fork:

- **HRG model** (`src/src/Model_HRG.cxx`): $T=0.155$ GeV, cubic box of side 30 fm, $\mu=0$. **Do not change these
  values.** They are private members; add read-only getters (`GetTemperature()`, `GetVolume()` returning
  $\text{size}^3$) so the new code reads them from the model instead of hard-coding them again.
- **Statistics**: primordials use Bose–Einstein (integer spin) or Fermi–Dirac (half-integer spin),
  $f_a(\epsilon)=1/(e^{\epsilon/T}\mp1)$, degeneracy $g_a=2J+1$; charge states are separate entries. The table uses the
  same $f_a$ and $g_a$.
- **Multiplicities**: Poisson, means from Monte Carlo integration (`Integrator`, cached in `fmultiplicity_<hash>.txt`).
  The table computes $n_a$ by quadrature and cross-checks against `GetMultiplicity()/V`.
- **Decays** (`src/src/ParticleDecayer.cxx`): pole masses (the `BreitWigner` function is never called; widths only set
  decay times). Two-body decays are isotropic in the parent rest frame. Three-body decays sample $(E_1,E_2)$ uniformly on
  the Dalitz plot (constant matrix element) followed by a random rotation. This matches the assumptions of the table.
- **Branching ratios** (`src/src/Parser.cxx`, `src/src/DecayTable.cxx`):
  - channels closed at the pole mass are dropped when the table is read;
  - for two-body channels flagged in `decays.data`, BR is multiplied by the squared isospin Clebsch–Gordan coefficient,
    and by 2 for daughters with equal mass, spin, strangeness and charm but different $I_3$ (so $\rho^0\to\pi^+\pi^-$ has
    BR 1);
  - BRs are **not** renormalized (`_PARTICLE_DECAYER_RESCALE_CHANNELS_` is not defined): if they sum to less than 1 the
    particle survives with the remaining probability (e.g. about half of the $K^0$).

  **Requirement:** the table code takes all channels and BRs from THERMINATOR's `ParticleDB` after parsing, never by
  re-parsing `decays.data`.
- **Ancestry** (`src/src/Particle.cxx`): primordials have `fathereid = -1` and `rooteid = eid`; decay products have
  `fathereid` = parent `eid` and `rooteid` = primordial ancestor `eid`. Decayed particles stay in the list with
  `decayed = 1`.
- **Weak decays** are on by default. Weakly decaying species with channels in the table: $K^0,\bar K^0$ (BR sum 0.5),
  $\Lambda$, $\Sigma^\pm$, $\Xi^0$, $\Xi^-$, $\Omega^-$ and antiparticles ($\Gamma\sim10^{-15}$ GeV). $\Sigma^0\to\Lambda\gamma$
  and $\eta$ decays are electromagnetic / isospin violating and always stay on. $\pi^\pm,K^\pm,p,n$ have no channels.
- Existing example `src/work/ex_mass_distribution_connected.cxx` groups final particles by `rooteid`; reuse its event
  loop. (Its comments mention `ex_balance_function.cxx` and `ex_connected_correlations_2d.cxx`, which are not in the
  synced repo.)

### THERMINATOR API notes

- Programs (in `src/work/`, registered in `src/work/CMakeLists.txt` with `add_executable` and
  `target_link_libraries(<name> therminator ${ROOT_LIBRARIES})`; macros via `configure_file`):
  `hrg_pair_table.exe` (tables, M1 tests), `hrg_pair_mc.exe` (events, histograms, sum-rule counters),
  `hrg_pair_compare.C` (plots and pass/fail numbers). Options via `CLI11.hpp`: `--share`, `--weak on|off`,
  `--output`; MC also `--events`, `--samples`, `--seed`.
- Setup: `ParticleDB(shareDir)`, `Integrator(nSamples, model)`, `integrator.SetMultiplicities(db)`,
  `Event(db, integrator)`.
- `ParticleType`: `GetName()`, `GetPDGCode()`, `GetMass()`, `GetSpin()` ($g=2J+1$, integer spin = boson),
  `GetCharge()`, `GetBarionN()`, `GetStrangeN()`, `GetGamma()`, `GetMultiplicity()` = mean number of primordials per
  event (density = `GetMultiplicity()/V`), `GetTable()` (null if no channels).
- `DecayTable::GetChannelCount()` returns the number of channels **minus one**: loop `c = 0 .. GetChannelCount()`
  inclusive. `DecayChannel`: `GetBranchingRatio()` (the effective $B$), `GetParticle1/2/3()` = indices into
  `ParticleDB`, `GetParticle3() == -1` for two-body channels.
- Photons (`gam000zer`) are never generated as primordials: exclude them from sums over primordial species.
- Event loop: `Reset()` (also restarts the `eid` counter, so `eid` is unique per event), `GeneratePrimordials()`,
  `DecayParticles()`, `GetParticleList()`. Per particle: `e, px, py, pz, mass, pid` (PDG), `eid`, `fathereid`,
  `rooteid`, `GetDecayed()`.
  - final: `GetDecayed() == 0` (includes particles that survived because their BRs sum to less than 1);
  - primordial: `fathereid == -1` (then `rooteid == eid`);
  - parent: particle with `eid == fathereid` (build an `eid -> particle` map per event; decayed particles stay in
    the list); primordial ancestor: particle with `eid == rooteid`.
  - L2: `p1.rooteid == p2.rooteid`; L1: `p1.fathereid == p2.fathereid != -1`; L0: L1 and
    `parent(p1).fathereid == -1`. Group by `rooteid` first, as in `ex_mass_distribution_connected.cxx`.
- Measured species by `pid`: ±211, ±321, ±2212.
- Weak switch: `Event::DecayParticles()` decays types with `GetGamma() >= 0` and a non-empty table. Implement
  `--weak off` as a set of PDG codes skipped in that loop (an `Event` option, default empty): ±311, ±3122, ±3222,
  ±3112, ±3322, ±3312, ±3334 (check against the `MC` column of `particles.data`). The table program treats the same
  set as stable.
- Pair kinematics: $\cos\theta^*$ from `TLorentzVector::Boost(-K.BoostVector())` applied to $p_1$, then
  $\hat{\mathbf p}_1^*\cdot\hat{\mathbf K}$.
- After each milestone, update `STATUS_pair_functions.md` (commands, file and histogram names, test results,
  deviations from this spec).

### Switches (same setting for table and MC in every comparison)

- `--weak on|off`. `off` makes the 8 species above and their antiparticles stable (no decay; they become final
  particles). Default `on`.
- `--stat qs|boltzmann` for the table parents. Default `qs` (as THERMINATOR). `boltzmann` is only for debugging; the MC
  always uses THERMINATOR's statistics.

## 2. Species and pairs

Measured species: $\pi^+,\pi^-,K^+,K^-,p,\bar p$. Same-species pairs only (no $\pi K$, $\pi p$, $Kp$ yet), all four
ordered charge combinations: $(\pi^+\pi^+)$, $(\pi^+\pi^-)$, $(\pi^-\pi^+)$, $(\pi^-\pi^-)$, and likewise for $K$ and $p$:
12 ordered pairs.

Ordered pairs of distinct particles: for $s=s'$ both orderings $(i,j)$ and $(j,i)$ are filled, consistent with the
counting $n_s(n_{s'}-\delta_{ss'})$ of LCS.

## 3. The table (Goal 1)

### 3.1 Channel loop

For each parent $a$ in `ParticleDB` and each channel with BR $B$ (as stored after parsing), for each ordered pair of
distinct daughter slots $(i,j)$ whose species are an ordered measured pair $(s,s')$:

- **two-body channel** ($a\to ij$): add a pole (LCS Eqs. 13–14)
  $$\mathcal C_K^{a\to ss'}=\frac{A}{m_a}\,f_a(\epsilon_K)\,\delta(M-m_a),\qquad
  A=\frac{B}{\Phi_2(m_a;m_i,m_j)}\frac{g_a}{(2\pi)^3},\qquad \Phi_2(M;m_1,m_2)=\frac{4\pi p^*}{M},\qquad
  \mathcal C_u^{a\to ss'}=0 .$$
- **three-body channel** ($a\to ijk$, $m_1=m_i$, $m_2=m_j$, $m_3=m_k$): add a continuum (LCS Eqs. 15–18),
  $$D=\frac{B}{\Phi_3(m_a)}\frac{g_a}{(2\pi)^3},\qquad
  E_P^*=\frac{m_a^2+M^2-m_3^2}{2M},\qquad p_3^*=\frac{\lambda^{1/2}(m_a^2,M^2,m_3^2)}{2M},\qquad
  \epsilon_P(t)=\frac{\epsilon_KE_P^*+\kappa p_3^*\,t}{M},\ t\in[-1,1],$$
  $$\boxed{\;\mathcal C_K=\frac{2\pi D\,p_3^*}{\kappa M^2}\int_{-1}^{1}dt\,\big(\kappa E_P^*+\epsilon_Kp_3^*\,t\big)\,f_a(\epsilon_P(t)),\qquad
  \mathcal C_u=-\frac{2\pi D\,p_3^{*2}}{\kappa\,\epsilon_K}\int_{-1}^{1}dt\;t\,f_a(\epsilon_P(t)),\;}$$
  for $m_1+m_2\le M\le m_a-m_3$, zero outside. This is LCS Eq. (18) with $\epsilon_P$ traded for $t$; the substitution
  keeps the integrand regular as $\kappa\to0$, where both functions are finite. Use Gauss–Legendre with $\ge 64$ points.

Phase space (LCS Sec. 2): $\Phi_3(m_a)=\int_{m_1+m_2}^{m_a-m_3}M\,dM\,\Phi_2(m_a;M,m_3)\,\Phi_2(M;m_1,m_2)$, and also
$\Phi_3=2\pi^2A_{\rm D}/m_a^2$ with $A_{\rm D}$ the Dalitz area in $(s_{12},s_{23})$,
$$A_{\rm D}=\int_{(m_1+m_2)^2}^{(m_a-m_3)^2}ds_{12}\;4\,p_2^\dagger p_3^\dagger,\qquad
E_2^\dagger=\frac{s_{12}-m_1^2+m_2^2}{2\sqrt{s_{12}}},\quad E_3^\dagger=\frac{m_a^2-s_{12}-m_3^2}{2\sqrt{s_{12}}},\quad
p_i^\dagger=\sqrt{E_i^{\dagger2}-m_i^2},$$
where $\dagger$ denotes the rest frame of particles 1 and 2, and $4p_2^\dagger p_3^\dagger=(s_{23})_{\max}-(s_{23})_{\min}$ (PDG kinematics
review, Sec. 49.4.3.1). Compute both (test T1.1).

The ordered-slot loop gives each channel $n_s(n_{s'}-\delta_{ss'})$ terms automatically (e.g. $\eta\to\pi^+\pi^-\pi^0$
contributes once to $\pi^+\pi^-$ and once to $\pi^-\pi^+$; a channel with two $\pi^+$ contributes twice to $\pi^+\pi^+$).

Channels with photons or neutral hadrons as the third daughter are included ($m_3=0$ for $\gamma$).

### 3.2 Grid and output

- $M$: from $m_1+m_2$ to 2 GeV, step 2–5 MeV, refined near threshold.
- $\kappa$: from 0 to 3 GeV, refined at small $\kappa$.
- Output ROOT file, per ordered pair $(s,s')$:
  - `CK_cont_<pair>`, `Cu_cont_<pair>`: `TH2D` in $(M,\kappa)$, summed over three-body channels; also per channel
    (`CK_<parent>_<channel>`) for the mixture weights;
  - a `TTree` of poles: parent name, $m_a$, $A/m_a$, BR, and $f_a$ on the $\kappa$ grid;
  - metadata: $T$, $V$, weak switch, statistics, `ParticleDB` hash, list of dropped channels.

### 3.3 Derived quantities (gas at rest)

From IM Eqs. (5), (12):

- Pair mass distribution per unit volume:
  $$\rho_{ss'}(M)=4\pi p^*(M)\int d^3K\,(\mathcal C_K+\mathcal C_u)=16\pi^2p^*\int_0^\infty\kappa^2d\kappa\,(\mathcal C_K+\mathcal C_u)(M^2,\epsilon_K).$$
  A pole contributes $B\,n_a\,\delta(M-m_a)$ (IM Eq. 7), with $n_a=g_a\int d^3k\,f_a/(2\pi)^3$.
- Normalized pair spectrum at fixed $M$, integrated over the direction of $\mathbf K$ and over $\Omega^*$:
  $$\frac{dP}{d\kappa}\Big|_M=\frac{16\pi^2p^*\,\kappa^2\,(\mathcal C_K+\mathcal C_u)(M^2,\epsilon_K)}{\rho_{ss'}(M)},\qquad \int d\kappa\,\frac{dP}{d\kappa}=1 .$$
- Boltzmann reference of a particle of mass $M$:
  $$\frac{dP_{\rm th}}{d\kappa}\Big|_M=\frac{\kappa^2e^{-\epsilon_K/T}}{M^2T\,K_2(M/T)} .$$
- All predictions are integrated over the same $M$ bins as the MC histograms. A pole falls entirely inside one bin; in
  that bin the prediction is the mixture of the pole and any continuum underneath, with weights
  $w_{\rm ch}=\rho_{\rm ch}/\rho_{ss'}$ (IM Eq. 14).

## 4. The MC (Goal 2)

New program, HRG model, standard event loop (`GeneratePrimordials`, `DecayParticles`). Final particles: `decayed == 0`.
Build a map `eid -> particle` per event to find fathers.

### 4.1 Ancestry levels

For each ordered pair of distinct final particles $(1,2)$ with species in the measured list:

| Level | Condition | Content |
|---|---|---|
| L0 | same `fathereid` $\ne-1$, and the father has `fathereid == -1` (father primordial) | exactly what the table contains |
| L1 | same `fathereid` $\ne-1$ (any father) | + siblings from fed-down parents (non-thermal parent spectrum, still isotropic) |
| L2 | same `rooteid`, and not the same particle | + cascades (pairs from different stages of a chain): $\mathcal C_Q\neq0$ |

L0 $\subset$ L1 $\subset$ L2. Store the three levels separately (or L0, L1$-$L0, L2$-$L1).

### 4.2 Observables per level and ordered pair

Gas at rest, so the lab is the fluid rest frame.

- $M$; $\kappa=|\mathbf K|$; $\cos\theta^*$ = cosine of the angle between $\mathbf p_1^*$ (particle 1 boosted to the pair
  rest frame) and $\mathbf K$. Skip $\kappa<10^{-6}$ GeV.
- `TH3D` in $(M,\kappa,\cos\theta^*)$, filled with weight $1/(N_{\rm ev}V)$: the pair density per unit volume. Same $M$
  and $\kappa$ binning as the table (coarser $\kappa$ is fine).
- `TProfile2D` in $(M,\kappa)$ of $P_1(\cos\theta^*)$ and $P_2(\cos\theta^*)$, plus entries, giving the Legendre moments
  $a_\ell=\langle P_\ell\rangle$ with errors $\sigma(a_1)\simeq1/\sqrt{3N}$, $\sigma(a_2)\simeq1/\sqrt{5N}$ for a nearly
  flat distribution (use the profile's own error in practice).
- Single-particle counts per species (for the sum rule and per-trigger normalizations).
- Per event: net charge $Q$, net baryon number $B$, net strangeness $S$ of all final particles (all species, not only
  the measured ones), and the counts $N_+$, $N_-$ of all positive and negative final particles.
- For the sum rule (Sec. 6): integrated L0, L1, L2 pair counts for every ordered pair of **all** charged final species
  (cross-species included, counts only, no histograms), from which the $(++)$, $(--)$, $(+-)$ sums are formed.

Run with both `--weak on` and `--weak off`. Number of events: enough for $\lesssim1\%$ statistical error on $\rho(M)$ at
the $\rho^0$ pole bin and on $a_{1,2}$ in the populated $(M,\kappa)$ bins (start with $10^4$ events and scale).

## 5. Comparisons (Goals 2 and 3)

ROOT macros in `src/work/`, in the style of the existing `ex_plot_*.C`. For each ordered pair, each level, and a few
masses: the $\rho^0$, $\eta$, $\omega$ (two-body $\omega\to\pi^+\pi^-$ if present), and continuum masses, e.g. 0.3, 0.45,
0.6 GeV for $\pi\pi$.

1. **Mass distribution** $\rho_{ss'}(M)$: MC L0, L1, L2 overlaid on the table prediction, with ratio panels. Also the
   per-channel decomposition from the table.
2. **Normalized spectrum** $dP/d\kappa|_M$: MC L0, L1, L2 and the table prediction, all divided by the Boltzmann
   reference $dP_{\rm th}/d\kappa|_M$; the line 1 shown. Expected: the table deviates from 1 by quantum statistics at the
   poles (order $e^{-m_a/T}$: ~0.7% for the $\rho$, a few % for the $\eta$, largest at small $\kappa$), by the recoil
   against $X$ on the continua, and by the mixture where channels overlap. L0 must follow the table.
3. **Angular structure**: $a_1(\kappa)$ and $a_2(\kappa)$ at fixed $M$, for $(+-)$ and $(++)$ separately. For a weight
   $w(x)=\epsilon_K(\mathcal C_K+\mathcal C_u)+\epsilon_Q\mathcal C_Q$ with $\epsilon_Q=2\kappa p^*x/M$ (equal masses) and
   $\mathcal C$'s independent of $\epsilon_Q$,
   $$a_1=\frac{2\kappa p^*}{3M\epsilon_K}\,\frac{\mathcal C_Q}{\mathcal C_K+\mathcal C_u},$$
   so $a_1$ at L2 measures $\mathcal C_Q/(\mathcal C_K+\mathcal C_u)$ versus $(M,\epsilon_K)$; $a_2\ne0$ signals an $\epsilon_Q^2$
   dependence. For like-sign pairs filled in both orderings, $a_1=0$ by construction and only $a_2$ is informative.
   Also show a few raw $\cos\theta^*$ histograms.

Expected pattern:

| Projection | L0 | L1 | L2 |
|---|---|---|---|
| $\rho(M)$ | = table | same pole positions and continuum shapes, larger weights | extra continua from chains |
| $dP/d\kappa\,\vert_M$ | = table | different (non-thermal parent spectrum) | different |
| $a_1$, $a_2$ | 0 | 0 | $\ne0$ in general |

Like-sign pairs: essentially empty at L0 and L1 (one-step decays with two same-sign hadrons are rare in the list); they
are mainly an L2 observable. This is expected.

## 6. Sum rule (Goal 4)

1. **Chemistry of the table** (one-step, no MC):
   $$\int dM\,\rho_{ss'}(M)=c^{\rm 1step}_{ss'}=\sum_a n_a\sum_{\rm ch}B_{\rm ch}\,n^{\rm ch}_s\big(n^{\rm ch}_{s'}-\delta_{ss'}\big),$$
   IM Eq. (6). The right-hand side needs no integration.
2. **Exact charge sum rule in the MC.** With Poisson primordials and charge-conserving decays,
   $$\frac{\mathrm{Var}(Q_{\rm final})}{V}=\sum_aq_a^2\,n^{\rm in}_a=\frac{\langle N_+\rangle+\langle N_-\rangle}{V}+c_{++}+c_{--}-c_{+-},$$
   with $n^{\rm in}_a$ the primordial densities (`GetMultiplicity()/V`), $\pm$ = sums over all positive / negative final
   species, $c_{+-}$ including both orderings, and the $c$'s from **L2** pair counts per unit volume (IM Eq. 16,
   `susceptibility_sumrule.tex`). Same for baryon number and strangeness with the corresponding charges. Compare with the
   Poisson value $\sum q^2n$, not with the thermodynamic $T\chi$ (THERMINATOR does not include the $f(1\pm f)$ factor in
   its fluctuations).
3. **What the table misses.** Evaluate $n_++n_-+c_{++}+c_{--}-c_{+-}$ with the $c$'s from L0, from L1, and from the
   one-step counting $c^{\rm 1step}$ of check 1. The sum runs over all charged final species, including cross-species
   pairs (e.g. $K^{*0}\to K^+\pi^-$; with `--weak off` also $\Sigma^\pm$, $\Xi^-$, $\Omega^-$). Use the counting, which needs no
   table, for these; the table is only built for the measured same-species pairs. Expected: $c^{\rm 1step}$ = L0
   within statistics; the shortfall of L0 and L1 with respect to check 2 is the share of charge balancing carried by
   fed-down parents (L1$-$L0) and by cascades (L2$-$L1), resolved in $M$ by Sec. 5.1 for the measured pairs.

## 7. Milestones and pass/fail tests

Each milestone ends with its tests; stop and report before starting the next.

**M0. Setup.**
- Getters for $T$ and $V$ in `Model_HRG` (no change of values); `--weak` switch; new targets in
  `src/work/CMakeLists.txt`; build succeeds.
- T0.1: with `--weak off`, the 8 weakly decaying species and their antiparticles appear as final particles and none of
  their decay products do.

**M1. Table and analytic checks (no MC).**
- T1.1: $\Phi_3$ from the $M$ integral and from the Dalitz area agree to $10^{-8}$ for every three-body channel.
- T1.2: for every three-body channel, $\rho^{a\to ss'X}(M)$ from the table equals the Dalitz projection
  $B\,n_a\,M\,\Phi_2(M;m_1,m_2)\,\Phi_2(m_a;M,m_3)/\Phi_3(m_a)$ (IM Eq. 9) to $10^{-6}$. The Python prototype of
  Appendix A reproduces it to $10^{-10}$ for $\omega\to\pi^+\pi^-\pi^0$; use it as the reference for the C++ functions.
- T1.3: every pole carries $\int dM\,\rho=B\,n_a$ (IM Eq. 7) to $10^{-6}$.
- T1.4: $\int dM\,\rho_{ss'}(M)=c^{\rm 1step}_{ss'}$ (Sec. 6.1) to $10^{-6}$.
- T1.5: $n_a$ from quadrature agrees with `GetMultiplicity()/V` within the Monte Carlo integration error.
- T1.6: $\mathcal C^\mu=\mathcal C_KK^\mu+\mathcal C_u\epsilon_Ku^\mu$ is future timelike everywhere on the grid:
  $\epsilon_K(\mathcal C_K+\mathcal C_u)>0$ (positive pair density) and $\epsilon_K(\mathcal C_K+\mathcal C_u)>\kappa|\mathcal C_K|$; also
  $\mathcal C_u>0$ (since $f_a$ decreases with energy). Both functions finite as $\kappa\to0$. **$\mathcal C_K$ alone is not
  sign-definite**: it is negative near threshold for channels where the parent is fast in the pair frame
  (Boltzmann, $\kappa\to0$: $\mathcal C_K\propto E_P^*-p_3^{*2}/3T$); e.g. $\omega\to\pi^+\pi^-\pi^0$ near $M=2m_\pi$. Do not flag
  this as an error.
- Report: list of dropped channels (closed at the pole mass) and the BR sums below 1.

**M2. MC at L0 versus the table.**
- T2.1: $\rho(M)$, L0 versus table, $\chi^2/\mathrm{ndf}\simeq1$ for each ordered pair, including the pole bins.
- T2.2: $dP/d\kappa|_M$, L0 versus table, $\chi^2/\mathrm{ndf}\simeq1$ at the chosen masses.
- T2.3: $a_1$, $a_2$ at L0 consistent with 0 (pull distribution of mean 0 and width 1 over populated bins). A failure
  for three-body channels points to the rotation in the three-body decayer.

**M3. L1, L2 and the thermal comparison.**
- The plots of Sec. 5 for both weak settings.
- Table of $a_1$, $a_2$ at L2 and the extracted $\mathcal C_Q/(\mathcal C_K+\mathcal C_u)$ in $(M,\epsilon_K)$ bins.
- T3.1: L1 $a_{1,2}$ consistent with 0.

**M4. Sum rule.**
- T4.1: Sec. 6.2 holds within statistical errors for $Q$, $B$, $S$, both weak settings.
- T4.2: $c^{\rm 1step}$ from the counting equals the L0 pair counts per unit volume within statistics, for every pair of
  charged final species.
- Report: the decomposition of Sec. 6.3 (table, L0, L1, L2).

**M5. BlastWave (later; outline only).**
- Same comparisons with the `BlastWave` model, where $d\Sigma$ is not parallel to $u$ and $\mathcal C_K$, $\mathcal C_u$
  enter separately: table integrated over the surface versus MC L0, in the pair variables. Only after M1–M4 pass.

## 8. Not in scope (now)

- Cross-species pairs ($\pi K$, $\pi p$, $Kp$).
- Finite widths, non-constant matrix elements, spin alignment.
- Saddle-point checks.
- Acceptance cuts and balance functions (they use these tables; see the acceptance notes).
- $\mathcal C_Q$ and the $\epsilon_Q^2$ dependence of chain pairs (Part II keeps the $\Omega^*$ average).

## 9. Part II: fed-down parents (L1) and decay chains (L2)

Physics reference: `pair_functions_chains.tex` (PFC) in the notes repo. Status before Part II: M0–M3 done (report of
24 October 2026), M4 not started. M4 can be done before, after or in parallel with Part II.

### 9.1 Idea

Every correlated pair has a unique last common ancestor (LCA) $c$. The pair density is
$\sum_c$ (current of $c$) $\otimes$ (pair kernel of the vertex $c$). Levels: L0 = thermal current, one-step vertex
(the existing table); L1 = full current, one-step vertex; L2 = full current, any number of further decays below the
vertex. $\mathcal C_Q=0$ and no $\epsilon_Q$ dependence are kept: at L1 this is exact; at L2 the table is the $\Omega^*$
average at fixed $(M,K)$ (exact for $\rho(M)$, $dP/d\kappa$ and the sum rule; it misses the $a_2$ found in M3).

The tables remain functions of $(M^2,\epsilon_K)$ at fixed $T$, $\mu$; flow (M5) uses them cell by cell and does not
require redoing any chain.

### 9.2 Single-particle currents (FastReso recursion, ported to THERMINATOR)

For every species $c$ in `ParticleDB`: $S_c^\mu=\mathcal A_c(\epsilon)p^\mu+\mathcal B_c(\epsilon)\epsilon u^\mu$; in the gas rest
frame the density per unit volume and $d^3p$ is $\mathcal A_c+\mathcal B_c$.

- Initialization: primordial species $\mathcal A_c=g_cf_c/(2\pi)^3$ (Bose–Einstein/Fermi–Dirac as `Model_HRG`),
  $\mathcal B_c=0$. Photons are not primordial.
- Order: decreasing mass. When species $a$ is reached its current is complete; each channel adds to every daughter
  slot $b$ (all species, not only the measured ones). Weak-stable species (switch off) are not decayed.
- Two-body rule $a\to bc$ (PFC Sec. 2; FastReso Eq. 12), regular form with $w\in[-1,1]$:
  $$\mathcal A_b(\epsilon)\mathrel{+}=B\,\frac{m_a^2}{m_b^2}\,\frac12\int_{-1}^1dw\,\mathcal A_a(E(w))\,\frac{Q(w)}{k},\qquad
  (\mathcal A_b+\mathcal B_b)(\epsilon)\mathrel{+}=B\,\frac{m_a^2}{m_b^2}\,\frac12\int_{-1}^1dw\,(\mathcal A_a+\mathcal B_a)(E(w))\,\frac{E(w)}{\epsilon},$$
  $E(w)=\frac{m_a}{m_b^2}(\epsilon E_b^*-wkp^*)$, $Q(w)=\frac{m_a}{m_b^2}(E_b^*k-wp^*\epsilon)$, $k=\sqrt{\epsilon^2-m_b^2}$, with $E_b^*$,
  $p^*$ the energy and momentum of $b$ in the rest frame of $a$. For a thermal parent this equals the two-body $F_1$,
  $F_1+F_u$ already used for the single-particle densities.
- Three-body rule $a\to bcd$: average the two-body rule over the effective partner mass
  $m_{\tilde c}\in[m_c+m_d,m_a-m_b]$ with weight $\propto p^a_{b|\tilde c}\,p^{\tilde c}_{c|d}$ (FastReso Eq. 4), normalized to 1.
- Grid: per species in $k$ (or kinetic energy) from 0 to a cut where the current is negligible (e.g. $\epsilon-m\le40T$);
  values of the parent beyond its grid are set to 0. Note that $E(w)$ reaches large values for light daughters
  (pions from heavy parents), so check that the cut does not matter (T6.2). Small $k$: the $1/k$ in $Q/k$ cancels in
  the integral; treat as for the three-body $\kappa\to0$ limit and test finiteness.
- Output: per species, $\mathcal A_c$, $\mathcal B_c$ on the grid, and $n_c=\int d^3p\,(\mathcal A_c+\mathcal B_c)$.

### 9.3 L1 table

Same channel loop, grid and output format as M1, with the thermal current replaced by $S_c$ (PFC Sec. 3):

- pole $c\to ss'$: $\mathcal C_K=\frac{B}{\Phi_2}\frac{\delta(M-m_c)}{m_c}\mathcal A_c(\epsilon_K)$,
  $\mathcal C_u=\frac{B}{\Phi_2}\frac{\delta(M-m_c)}{m_c}\mathcal B_c(\epsilon_K)$ (store $\mathcal A_c$, $\mathcal B_c$ on the $\kappa$ grid with the pole);
- continuum $c\to ss'X$, $\mathcal D=B/\Phi_3$, $\epsilon_P(t)$ as before:
  $$\mathcal C_K=\frac{2\pi\mathcal Dp_3^*}{\kappa M^2}\int_{-1}^1dt\,(\kappa E_P^*+\epsilon_Kp_3^*t)\,\mathcal A_c(\epsilon_P(t)),\qquad
  \mathcal C_K+\mathcal C_u=\frac{2\pi\mathcal Dp_3^*}{M\epsilon_K}\int_{-1}^1dt\,\epsilon_P(t)\,(\mathcal A_c+\mathcal B_c)(\epsilon_P(t)).$$

The table program gets `--level L0|L1|L2` (L0 = current behaviour).

### 9.4 Chain kernels and the L2 table

- Kernels by Monte Carlo decays at rest (an analytic two-particle recursion is a possible later replacement). For each
  species $c$ with channels (not weak-stable): create $N_c$ particles of type $c$ at rest in a particle list, decay them
  with `ParticleDecayer` until the end of the list is reached (as in `Event::DecayParticles`). Particles that do not
  decay (BR sum $<1$) still count in $N_c$.
- For each decay tree, the first-generation branch of a final particle is the daughter of $c$ it descends from
  (walk `fathereid` up to the daughter whose `fathereid` is the `eid` of $c$). Pairs: ordered pairs of distinct final
  particles of the measured species in **different** branches, of which **at least one is not a direct daughter** of
  $c$ (one-step pairs are in the L1 table). Pairs in the same branch belong to a deeper LCA and are not counted here.
- Fill $h_c^{\rm chain}(M,k_c)$, $k_c=|\mathbf K|$ in the rest frame of $c$ (the lab here), normalized per particle $c$:
  counts$/(N_c\,\Delta M\,\Delta k_c)$. $M$ binning = the table's.
- Convolution (PFC Sec. 4, master formula), per parent, with $E_P^*=\frac{m_c}{M}\sqrt{M^2+k_c^2}$,
  $p_P^*=\frac{m_c}{M}k_c$, $\epsilon_P(t)=(\epsilon_KE_P^*+\kappa p_P^*t)/M$:
  $$\bar{\mathcal C}_K\mathrel{+}=\frac{1}{4\pi p^*}\frac{m_c^2}{2M^2}\int dk_c\,h_c\int_{-1}^1dt\;\mathcal A_c(\epsilon_P)\frac{\kappa E_P^*+\epsilon_Kp_P^*t}{M\kappa},\qquad
  \bar{\mathcal C}_K+\bar{\mathcal C}_u\mathrel{+}=\frac{1}{4\pi p^*}\frac{m_c^2}{2M^2}\int dk_c\,h_c\int_{-1}^1dt\;(\mathcal A_c+\mathcal B_c)(\epsilon_P)\frac{\epsilon_P}{\epsilon_K}.$$
  The $k_c$ integral runs over the histogram bins (bin centre, or a few points per bin). $p^*$ is that of the measured
  pair at mass $M$.
- Process one parent at a time: build $h_c$, convolve, add to the table, discard $h_c$. Storing the kernels (sparse) is
  a later option (it would make a change of $T$ a re-convolution only).
- $L2$ table $=$ L1 table $+$ $\bar{\mathcal C}^{\rm chain}$. Statistics $N_c$: adaptive, so that the kernel's statistical error
  on the table is small compared with the event-MC error in the comparisons (record $N_c$ and the propagated
  uncertainty).
- Symmetry (may be used to save time, but check on a subset): at $\mu=0$ the kernel of $\bar c$ is the charge conjugate of
  that of $c$, and for the $\Omega^*$-averaged table the mirror pair $(\pi^-\pi^+)$ equals $(\pi^+\pi^-)$ and
  $(\pi^-\pi^-)$ equals $(\pi^+\pi^+)$.

### 9.5 Counting (exact references, PFC Sec. 5)

With $n_c$ the total and $n_c^{\rm th}$ the thermal density of $c$, and $\bar n_{d\to s}$ the mean number of final $s$
per particle $d$ ($=1$ if $d=s$ is final; recursion over the channels otherwise):
$c^{\rm L1}_{ss'}=\sum_cn_c\sum_{\rm ch}B\,n_s^{\rm ch}(n_{s'}^{\rm ch}-\delta_{ss'})$ and
$c^{\rm L2}_{ss'}=\sum_cn_c\sum_{\rm ch}B\sum_{i\ne j}\bar n_{d_i\to s}\bar n_{d_j\to s'}$ (ordered pairs of distinct daughter slots).
Also $n_c=n_c^{\rm th}+\sum_an_a\sum_{\rm ch}B\,n_c^{\rm ch}$.

### 9.6 Milestones and tests

**M6. Single-particle currents.**
- T6.1: one step from thermal parents reproduces the closed-form two-body $F_1$, $F_1+F_u$ (relative $10^{-10}$).
- T6.2: $n_c=\int d^3p(\mathcal A_c+\mathcal B_c)$ equals the counting of Sec. 9.5 for every species ($10^{-6}$); also checks the
  grid cut.
- T6.3: positivity: $\mathcal A_c+\mathcal B_c>0$ and $S_c$ future timelike, $\epsilon(\mathcal A_c+\mathcal B_c)>k|\mathcal A_c|$, finite at small $k$.
  ($\mathcal A_c$ alone need not be positive.)
- T6.4 (MC, new run of the event program with a `--spectra` option, e.g. $2\times10^4$ events): per-event mean count of
  every species $c$, decayed or not, divided by $V$, versus $n_c$ (pulls); momentum spectra $dN/dp$ of
  $\pi^+,K^+,p,\rho^0,K^{*0},\omega,\Lambda,\Delta^{++}$ (all, decayed or not) versus $4\pi p^2V(\mathcal A_c+\mathcal B_c)$,
  $\chi^2$ with $p>10^{-3}$. Both weak settings.

**M7. L1 table.**
- T7.1: with thermal-only currents, `--level L1` reproduces the M1 table (relative $10^{-12}$).
- T7.2: $\int dM\,\rho^{\rm L1}_{ss'}=c^{\rm L1}_{ss'}$; each pole carries $B\,n_c$ ($10^{-6}$).
- T7.3: $\mathcal C^\mu$ future timelike, $\mathcal C_K+\mathcal C_u>0$, finite as $\kappa\to0$ (as T1.6).
- T7.4: versus the MC at L1 (existing M2 files, cumulative L1): $\rho(M)$ and $dP/d\kappa|_M$ at the test masses, same
  criteria as T2.1, T2.2 (the like-sign variance $2\mu$ applies). Both weak settings.

**M8. Chain kernels and L2 table.**
- T8.1: for each parent, the kernel of all LCA pairs (one-step included, for this test only) integrates to
  $\sum_{\rm ch}B\sum_{i\ne j}\bar n_{d_i\to s}\bar n_{d_j\to s'}$ within its statistical error (pulls over parents).
- T8.2: $\int dM\,\rho^{\rm L2}_{ss'}=c^{\rm L2}_{ss'}$ within the kernel statistics.
- T8.3: versus the MC at L2 (existing M2 files): $\rho(M)$ and $dP/d\kappa|_M$ at the test masses, $\chi^2$ with the
  kernel uncertainty added in quadrature, $p>10^{-3}$. Both weak settings.
- Report: the chain contribution versus $M$ (fraction of L2), and the MC $a_2$ at L2 next to the table's $a_2=0$ as the
  size of the neglected structure.

After M8, M4 (sum rule) can use $c^{\rm L2}$: the L2 densities of all charged final species (counting, Sec. 9.5,
cross-species included) must close $\mathrm{Var}(Q)/V=\sum_aq_a^2n_a^{\rm th}$.

## Appendix A. Python reference prototype (three-body functions)

Run with `python3 prototype_3body.py` (needs numpy, scipy). Expected output: T1.1 relative difference about 1e-13, T1.2 worst
relative deviation about 2e-10 (PASS), finite positive values at small kappa. The C++ functions should reproduce
`CK_Cu` and both `Phi3` functions to machine precision for the same inputs.

```python
"""Reference prototype for the three-body pair functions (spec Sec. 3.1, tests T1.1-T1.2).

Gas at rest, one channel a -> 1 2 3 with constant matrix element, Boltzmann parent,
branching ratio B = 1 and degeneracy g_a = 1. Units: GeV.
"""
import numpy as np
from scipy import integrate

T = 0.155


def lam(x, y, z):
    return x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * z * x


def pstar(M, m1, m2):
    return np.sqrt(max(lam(M * M, m1 * m1, m2 * m2), 0.0)) / (2 * M)


def Phi2(M, m1, m2):
    return 4 * np.pi * pstar(M, m1, m2) / M


def Phi3_mass(ma, m1, m2, m3):
    f = lambda M: M * Phi2(ma, M, m3) * Phi2(M, m1, m2)
    return integrate.quad(f, m1 + m2, ma - m3, limit=200)[0]


def Phi3_dalitz(ma, m1, m2, m3):
    def width(s12):  # (s23)_max - (s23)_min = 4 p2 p3 in the 12 rest frame
        r = np.sqrt(s12)
        E2 = (s12 - m1**2 + m2**2) / (2 * r)
        E3 = (ma**2 - s12 - m3**2) / (2 * r)
        return 4 * np.sqrt(max(E2**2 - m2**2, 0)) * np.sqrt(max(E3**2 - m3**2, 0))
    AD = integrate.quad(width, (m1 + m2) ** 2, (ma - m3) ** 2, limit=200)[0]
    return 2 * np.pi**2 * AD / ma**2


f_parent = lambda e: np.exp(-e / T)            # Boltzmann parent distribution
tq, wq = np.polynomial.legendre.leggauss(64)   # Gauss-Legendre nodes on [-1, 1]


def CK_Cu(M, kap, ma, m3, D):
    """Spec Sec. 3.1 (LCS Eq. 18 with eps_P = (eK EP* + kap p3* t)/M)."""
    eK = np.sqrt(M * M + kap * kap)
    EP = (ma**2 + M**2 - m3**2) / (2 * M)
    p3 = np.sqrt(lam(ma**2, M**2, m3**2)) / (2 * M)
    fP = f_parent((eK * EP + kap * p3 * tq) / M)
    I0, I1 = np.sum(wq * fP), np.sum(wq * tq * fP)
    CK = 2 * np.pi * D * p3 / (kap * M**2) * (kap * EP * I0 + eK * p3 * I1)
    Cu = -2 * np.pi * D * p3**2 / (kap * eK) * I1
    return CK, Cu


if __name__ == "__main__":
    ma, m1, m2, m3 = 0.78266, 0.13957, 0.13957, 0.13498   # omega -> pi+ pi- pi0
    P3a, P3b = Phi3_mass(ma, m1, m2, m3), Phi3_dalitz(ma, m1, m2, m3)
    print(f"T1.1  Phi3: mass integral {P3a:.12g}, Dalitz area {P3b:.12g}, rel. diff {abs(P3a / P3b - 1):.1e}")

    D = 1.0 / ((2 * np.pi) ** 3 * P3a)                    # B = g_a = 1
    n_a = integrate.quad(lambda k: 4 * np.pi * k * k * f_parent(np.hypot(ma, k)), 0, 10)[0] / (2 * np.pi) ** 3
    worst = 0.0
    for M in [0.30, 0.40, 0.50, 0.60]:
        rho = 16 * np.pi**2 * pstar(M, m1, m2) * integrate.quad(
            lambda k: k * k * sum(CK_Cu(M, k, ma, m3, D)), 1e-6, 8, limit=200)[0]
        dalitz = n_a * M * Phi2(M, m1, m2) * Phi2(ma, M, m3) / P3a   # IM Eq. 9
        worst = max(worst, abs(rho / dalitz - 1))
        print(f"T1.2  M = {M:.2f} GeV: rho(M) from C = {rho:.10e}, Dalitz projection = {dalitz:.10e}")
    print(f"T1.2  worst relative deviation {worst:.1e}  ->  {'PASS' if worst < 1e-6 else 'FAIL'}")
    print("T1.6  small kappa, (C_K, C_u) finite; C_u > 0 (C_K may be negative near threshold):", [tuple(f"{c:.4e}" for c in CK_Cu(0.4, k, ma, m3, D)) for k in (1e-4, 1e-2, 0.5)])
```
