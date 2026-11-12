cd /home/alice/therminator2-hboltz/build/work
export THERMINATOR_SHARE=/home/alice/therminator2-hboltz/share

Run everything from build/work/ (that's where the .exe and the copied .C macro live), with the SHARE dir exported:

mkdir myrun && cd myrun
cp /home/alice/therminator2-hboltz/ini/blastwave.ini .
cp -r /path/to/hboltz/ini/fomodel .

# 1. generate + histogram — HRG (default model)  -> mass_distribution_hrg.root
./ex_mass_distribution.exe --model hrg

# 2. same for blastwave (needs an ini)           -> mass_distribution_blastwave.root
./ex_mass_distribution.exe --model blastwave --ini ../../ini/blastwave.ini

# 3. plot: overlay the two models -> mass_distribution.pdf
root -b -q 'ex_plot_mass_distribution.C("mass_distribution_hrg.root","mass_distribution_blastwave.root")'
If you'd rather keep the macro's default names, run with --output mass_hrg.root and --output mass_blastwave.root instead, then just root -b -q ex_plot_mass_distribution.C.


Optional — the second function in the same macro, the pT/|p|-sliced plot (-> mass_slices.pdf):
root -b -q -e 'gROOT->ProcessLine(".L ex_plot_mass_distribution.C"); ex_plot_mass_slices("mass_distribution_blastwave.root");'




Purpose

Example program that builds the invariant-mass distribution of correlated final-state particle pairs — pairs of non-decayed particles that descend from the same primordial resonance.
Default pairing: pid1 = 211 (π⁺) with pid2 = -211 (π⁻).
Pipeline (per event, nEvents = 500 default)

GeneratePrimordials() → DecayParticles() on an Event fed by Model + Integrator + ParticleDB (same core pipeline as therm2_events).
Model selectable: Model_HRG (default) or Model_BlastWave (needs --ini).
Correlation logic

Sorts the particle list in place by rooteid (id of the primordial ancestor) so same-ancestor particles are contiguous.
Walks each rooteid group, collects surviving pid1 particles into group1 and pid2 into group2.
Forms every group1 × group2 pair (skipping identical pointers), computes invMass() = √(E² − p²) of the summed 4-momenta.
Outputs (ROOT file, default mass_distribution.root)

hMass — 1D mass spectrum, 100 bins over 0–2 GeV, scaled by 1/nPid1Total → integral = mean number of associated pid2 per pid1.
hMass2D — mass vs. momentum of particle 1 (6 momentum bins, 0–2 GeV).
hPid1Counts — pid1 count vs. momentum (normalization / acceptance reference).
Prints "correlated yield (pairs per pid1 particle)" with error via IntegralAndError.
Momentum variable

BlastWave → transverse pT (radial flow breaks boost invariance).
HRG → full |p| (isotropic source).
