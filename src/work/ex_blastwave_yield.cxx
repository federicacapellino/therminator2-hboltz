// clang-format off
// ex_blastwave_yield.cxx — yield test for the blast-wave Monte Carlo
//
// Checks that Integrator::Integrate agrees with the semi-analytic yield formula
// derived in fix.tex (Section "A Simple Yield Test").
//
// Physics
// -------
// For a cylindrical blast-wave freeze-out at proper time tau, spatial rapidity
// range DeltaEta, and transverse radius rhoMax, with a uniform transverse
// velocity v_T (so u^0 = 1/sqrt(1-v_T^2) is constant over the disk), the
// Cooper-Frye multiplicity reduces to
//
//   N_i = V * u^0 * n_i(T, mu_i)
//
// where the freeze-out proper volume is
//
//   V = tau * DeltaEta * pi * rhoMax^2
//
// and n_i is the isotropic thermal number density in the fluid rest frame:
//
//   n_i = (2s+1)/(2 pi^2) * integral_0^inf p^2 dp / (exp((E-mu)/T) +/- 1)
//
// with E = sqrt(m^2 + p^2), + for fermions (FD), - for bosons (BE).
// The chemical potential is mu_i = B_i*mu_B + I3_i*mu_I + S_i*mu_S + C_i*mu_C.
//
// The u^0 factor arises because the Cooper-Frye surface element contracted with
// a constant 4-velocity u^mu on a constant-tau hypersurface produces a factor
// of u^0 relative to the comoving volume element.  See fix.tex for the
// derivation.
//
// Test strategy
// -------------
// The semi-analytic n_i is evaluated by ROOT::Math::GaussIntegrator::IntegralUp,
// which integrates over [0, inf) directly — no change of variables needed.
// This is deliberately independent of the MC sampler's own p -> zeta/(1-zeta)
// substitution, so agreement between the two validates the full MC pipeline:
// the hypercube volume, the importance-sampling jacobian, and the Cooper-Frye
// integrand implementation in Model_BlastWave.
//
// Particle names follow the SHARE convention XX####yyy (e.g. Ka0492plu for K+).
// Use share/particles.data to look up names; do not guess from PDG masses.
//
// Usage:
//   ex_blastwave_yield.exe --ini <blastwave.ini> [--share <dir>]
//                          [--samples <N>] [--tol <eps>] [--seed <N>]
// clang-format on

#include "CLI11.hpp"
#include "Integrator.h"
#include "Model_BlastWave.h"
#include "ParticleDB.h"
#include "ParticleType.h"
#include "THGlobal.h"
#include <Math/Functor.h>
#include <Math/GaussIntegrator.h>
#include <TRandom.h>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

// thermalDensity — isotropic equilibrium number density in natural units (GeV^3)
//
//   n_i = (2s+1)/(2pi^2) * integral_0^inf p^2 dp / (exp((E-mu)/T) + stat)
//
// stat = +1 for half-integer spin (Fermi-Dirac),
// stat = -1 for integer spin     (Bose-Einstein).
//
// GaussIntegrator::IntegralUp maps [0,inf) internally; the integrand is smooth
// and exponentially suppressed for large p, so convergence is rapid.
double thermalDensity(ParticleType* pt, const Thermodynamics* thermo) {
  const double T    = thermo->GetTemperature();
  const double mu   = thermo->GetChemicalPotential(pt);
  const double m    = pt->GetMass();
  const double spin = pt->GetSpin();
  const double stat = ((spin - static_cast<int>(spin)) < 0.01) ? -1.0 : +1.0;
  const double g    = 2.0 * spin + 1.0;

  auto integrand = [=](double p) -> double {
    const double E = std::sqrt(m * m + p * p);
    return p * p / (std::exp((E - mu) / T) + stat);
  };

  ROOT::Math::Functor1D f(integrand);
  ROOT::Math::GaussIntegrator gi;
  gi.SetFunction(f);
  gi.SetRelTolerance(1.e-8);
  gi.SetAbsTolerance(1.e-14);

  return g / (2.0 * M_PI * M_PI) * gi.IntegralUp(0.0);
}

int main(int argc, char** argv) {
  std::string iniFile  = "";
  std::string shareDir = "";
  int    nSamples      = 2000000;
  double tolerance     = 0.01;   // acceptable fractional difference
  int    seed          = 42;

  CLI::App app{"ex_blastwave_yield: compare MC yield to semi-analytic formula"};
  app.add_option("--ini",     iniFile,   "blast-wave ini file")->required();
  app.add_option("--share",   shareDir,  "path to SHARE particle database");
  app.add_option("--samples", nSamples,  "MC samples per species")->default_val(2000000);
  app.add_option("--tol",     tolerance, "max allowed fractional difference")->default_val(0.01);
  app.add_option("--seed",    seed,      "Seed: set to 0 for a random seed")->default_val(42);
  CLI11_PARSE(app, argc, argv);

  gRandom->SetSeed(seed);

  // Build model and integrator
  Model_BlastWave model(iniFile.c_str());
  auto db = std::make_unique<ParticleDB>(shareDir);
  Integrator integ(nSamples, &model);

  // Extract geometry and flow from model
  const Thermodynamics* thermo = model.GetThermodynamics();
  const double vt     = model.GetVt();
  const double u0     = 1.0 / std::sqrt(1.0 - vt * vt);
  const double tau    = model.GetTau();              // GeV^-1
  const double deta   = model.GetRapSRange();
  const double rhomax = model.GetRhoMax();           // GeV^-1
  const double volume = tau * deta * M_PI * rhomax * rhomax;  // GeV^-3

  printf("\nBlast-wave parameters:\n");
  printf("  T      = %.1f MeV\n", thermo->GetTemperature() * 1000.0);
  printf("  v_T    = %.4f c    u^0 = %.4f\n", vt, u0);
  printf("  tau    = %.2f fm\n", tau    * kHbarC);
  printf("  rhoMax = %.2f fm\n", rhomax * kHbarC);
  printf("  Deta   = %.2f\n",    deta);
  printf("  V*u0   = %.4e GeV^-3\n", volume * u0);

  // Test species
  const std::vector<std::string> species = {"pi0139plu", "Ka0492plu", "pr0938plu"};

  printf("\n%-12s  %10s  %10s  %8s  %s\n",
         "particle", "MC yield", "semi-anal.", "ratio", "status");
  printf("%s\n", std::string(60, '-').c_str());

  int nFail = 0;
  for (const auto& name : species) {
    if (!db->ExistsParticleType(name.c_str())) {
      printf("%-12s  not found in database\n", name.c_str());
      nFail++;
      continue;
    }
    ParticleType* pt = db->GetParticleType(name.c_str());

    // Monte Carlo yield
    double maxInt   = 0.0;
    double mc_yield = integ.Integrate(pt, maxInt);

    // Semi-analytic yield
    double n_semi       = thermalDensity(pt, thermo);
    double semi_yield   = volume * u0 * n_semi;

    double ratio  = mc_yield / semi_yield;
    bool   passed = std::abs(ratio - 1.0) < tolerance;
    if (!passed) nFail++;

    printf("%-12s  %10.5e  %10.5e  %8.6e  %s\n",
           name.c_str(), mc_yield, semi_yield, ratio,
           passed ? "OK" : "FAIL");
  }

  printf("\n%s\n", nFail == 0 ? "All tests passed." : "SOME TESTS FAILED.");
  return nFail;
}
