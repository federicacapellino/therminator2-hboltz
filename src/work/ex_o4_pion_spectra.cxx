// t2.cxx -- pion spectra from HRG with O(4) critical enhancement
//
// Generates pion momentum spectra using the Therminator 2 HRG framework and
// estimates the effect of a modified pion dispersion relation near the O(4)
// chiral critical point (see o4_dispersion below).
//
// Therminator 2 generates primordial pions with the vacuum Bose-Einstein
// distribution fvac(p). To obtain spectra corresponding to a modified
// dispersion relation, primordial pions are reweighted by f(p)/fvac(p), where
// f(p) uses the softened O(4) dispersion curve. Stimulated emission from short-
// lived resonance decays is estimated by an additional factor of (1 + f(p)),
// accounting for the enhanced probability of emission into an already-occupied
// state.
//
// This file also serves as an example of the Therminator 2 event loop:
// generating events, iterating the particle list, and identifying primordial
// vs. decay particles and whether a particle has already decayed.

#include "CLI11.hpp"
#include "Event.h"
#include "Integrator.h"
#include "Model_HRG.h"
#include "ParticleDB.h"
#include "ParticleType.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TH1D.h>
#include <TRandom.h>
#include <cstdio>
#include <iostream>

//! o4_dispersion: modified pion dispersion relation near the O(4) chiral
//! critical point. Parameterizes the softening of the pion mode via a
//! momentum-dependent pole mass and chiral velocity.
class o4_dispersion {
private:
  double fT;       // temperature (GeV)
  double fLambda;  // momentum cutoff (GeV)
  double fMvacuum; // vacuum pion mass (GeV)
  double fRatio;   // chiral condensate ratio sigma/sigma_0

public:
  o4_dispersion(double lambda_by_T = M_PI, double ratio = 0.5);

  double f(double p);      // modified Bose-Einstein distribution
  double fvac(double p);   // vacuum Bose-Einstein distribution
  double eofp(double p);   // modified dispersion E(p) in GeV
  double m2pole(double p); // modified pole mass squared (GeV^2)
  double v2(double p);     // chiral velocity squared

  // writes p, v2, m2pole, E(p), f, fvac vs p to a text file
  void plot(const std::string& filename = "pidispersion.out");
};

o4_dispersion::o4_dispersion(double lambda_by_T, double ratio)
    : fT(0.155), fLambda(lambda_by_T * 0.155), fMvacuum(0.1396), fRatio(ratio) {}

double o4_dispersion::eofp(double p) {
  return sqrt(v2(p) * p * p + m2pole(p));
}

double o4_dispersion::fvac(double p) {
  double x = sqrt(fMvacuum * fMvacuum + p * p) / fT;
  return exp(-x) / (1. - exp(-x));
}

double o4_dispersion::f(double p) {
  double x = eofp(p) / fT;
  return exp(-x) / (1. - exp(-x));
}

double o4_dispersion::m2pole(double p) {
  double mv2 = fMvacuum * fMvacuum;
  double m02 = mv2 * fRatio;
  double x = p * p / (fLambda * fLambda);
  return mv2 - (mv2 - m02) / (1 + x / 2. + x * x);
}

double o4_dispersion::v2(double p) {
  double v02 = fRatio * fRatio;
  double x = p * p / (fLambda * fLambda);
  return 1. - (1. - v02) / (1 + x / 2. + x * x);
}

void o4_dispersion::plot(const std::string& filename) {
  int np = 100;
  double pmin = 0.;
  double pmax = 1.5;
  double dp = (pmax - pmin) / np;
  FILE* fp = fopen(filename.c_str(), "w");
  for (int ip = 0; ip < np; ip++) {
    double p = pmin + ip * dp;
    fprintf(fp, "%15.5e %15.5e %15.5e %15.5e %15.5e %15.5e\n",
            p, v2(p), m2pole(p), eofp(p), f(p), fvac(p));
  }
  fclose(fp);
}

////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
  std::string shareDir   = "";
  std::string outputFile = "o4_pion_spectra.root";
  int nEvents            = 500;
  int nSamples           = 1000000;
  int seed               = 0;

  CLI::App app{"t2: HRG thermal model pion spectra with O(4) critical enhancement"};
  app.add_option("--share",   shareDir,   "path to SHARE particle database (default: $THERMINATOR_SHARE)");
  app.add_option("--output",  outputFile, "output ROOT file")->default_str("o4_pion_spectra.root");
  app.add_option("--events",  nEvents,    "number of events to generate")->default_val(500);
  app.add_option("--samples", nSamples,   "Monte Carlo integration samples")->default_val(1000000);
  app.add_option("--seed",    seed,       "random seed")->default_val(0);
  CLI11_PARSE(app, argc, argv);

  gRandom->SetSeed(seed);

  o4_dispersion o4;
  // o4.plot(); // uncomment to write dispersion curve to pidispersion.out

  auto tPartDB = std::make_unique<ParticleDB>(shareDir);
  Model_HRG hrg;
  auto tInteg = std::make_unique<Integrator>(nSamples, &hrg);
  tInteg->SetMultiplicities(tPartDB.get());
  auto tEvent = std::make_unique<Event>(tPartDB.get(), tInteg.get());

  TFile fout(outputFile.c_str(), "RECREATE", "O(4) pion spectra");
  TH1::SetDefaultSumw2();

  // Direct pions (primordial, fathereid == -1)
  TH1D* hpi           = new TH1D("pi_dNdp",           "dN/dp", 40, 0.0001, 2.);
  TH1D* hpi_stimulated = new TH1D("pi_stimulated_dNdp", "dN/dp", 40, 0.0001, 2.);
  TH1D* hpi_critical   = new TH1D("pi_critical_dNdp",   "dN/dp", 40, 0.0001, 2.);

  // All pions (primordial + decay daughters)
  TH1D* hpi_all           = new TH1D("pi_all_dNdp",           "dN/dp", 40, 0., 2.);
  TH1D* hpi_all_stimulated = new TH1D("pi_all_stimulated_dNdp", "dN/dp", 40, 0., 2.);
  TH1D* hpi_all_critical   = new TH1D("pi_all_critical_dNdp",   "dN/dp", 40, 0., 2.);
  TH1D* hpi_all_weak       = new TH1D("pi_all_weak_dNdp",       "dN/dp", 40, 0., 2.);

  // Pions from decays, split by parent lifetime (tcut = 2.5 fm/c)
  const double tcut = 2.5;
  TH1D* hpi_decayed       = new TH1D("pi_decayed_dNdp",       "dN/dp", 40, 0., 2.);
  TH1D* hpi_decayed_stimulated = new TH1D("pi_decayed_stimulated_dNdp", "dN/dp", 40, 0., 2.);
  TH1D* hpi_decayed_critical   = new TH1D("pi_decayed_critical_dNdp",   "dN/dp", 40, 0., 2.);
  TH1D* hpi_decayed_short = new TH1D("pi_decayed_short_25_dNdp", "dN/dp", 40, 0., 2.);
  TH1D* hpi_decayed_long  = new TH1D("pi_decayed_long_25_dNdp",  "dN/dp", 40, 0., 2.);
  TH1D* hpi_decayed_short_stimulated =
      new TH1D("pi_decayed_short_stimulated_25_dNdp", "dN/dp", 40, 0., 2.);
  TH1D* hpi_decayed_short_critical =
      new TH1D("pi_decayed_short_critical_25_dNdp", "dN/dp", 40, 0., 2.);

  TH1D* htimes = new TH1D("times", "dNdt", 40, 0.0001, 10.);

  for (int ievent = 0; ievent < nEvents; ievent++) {
    tEvent->Reset();
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();

    for (auto& part : tEvent->GetParticleList()) {
      if (ievent == 0) {
        std::cout << part.MakeTEXTEntry() << std::endl;
      }

      const int pid      = part.GetParticleType()->GetPDGCode();
      const bool isPion  = (std::abs(pid) == 211);
      const bool decayed = part.GetDecayed();

      if (isPion && !decayed) {
        hpi_all_weak->Fill(part.GetP());
      }

      if (decayed) {
        continue;
      }

      const double time     = part.t * kHbarC;
      const double p        = part.GetP();
      const bool primordial = (part.fathereid == -1);

      if (time > 10000.) {
        continue;
      }

      if (isPion) {
        hpi_all->Fill(p);
      }

      if (isPion && primordial) {
        hpi->Fill(p);
        hpi_stimulated->Fill(p);
        hpi_critical->Fill(p, o4.f(p) / o4.fvac(p));
        hpi_all_stimulated->Fill(p);
        hpi_all_critical->Fill(p, o4.f(p) / o4.fvac(p));
      }

      if (isPion && !primordial) {
        hpi_decayed->Fill(p);
        htimes->Fill(time);
        if (time < tcut) {
          hpi_decayed_short->Fill(p);
          hpi_decayed_short_stimulated->Fill(p, 1. + o4.fvac(p));
          hpi_decayed_short_critical->Fill(p, 1. + o4.f(p));
        } else {
          hpi_decayed_long->Fill(p);
        }
      }
    }
  }

  // Normalize stimulated/critical short-decay histograms to vacuum integral
  double IVac = hpi_decayed_short->Integral();
  double IStm = hpi_decayed_short_stimulated->Integral();
  double ICrt = hpi_decayed_short_critical->Integral();
  hpi_decayed_short_stimulated->Scale(IVac / IStm);
  hpi_decayed_short_critical->Scale(IVac / ICrt);

  // Build total decay contributions
  hpi_decayed_critical->Add(hpi_decayed_short_critical);
  hpi_decayed_critical->Add(hpi_decayed_long);
  hpi_all_critical->Add(hpi_decayed_critical);

  hpi_decayed_stimulated->Add(hpi_decayed_long);
  hpi_decayed_stimulated->Add(hpi_decayed_short_stimulated);
  hpi_all_stimulated->Add(hpi_decayed_stimulated);

  fout.Write();
  return 0;
}
