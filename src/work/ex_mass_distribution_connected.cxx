// ex_mass_distribution_connected.cxx — correlated pair invariant mass distribution
//
// Generates thermal events and histograms the invariant mass of final-state
// particle pairs that share the same primordial ancestor (rooteid).  Two
// non-decayed particles are "correlated" if they were both produced in the
// decay chain of the same primordial resonance.
//
// After all events the histogram is normalized by the total number of type-1
// particles produced.  The integral of the normalized histogram is therefore
// the mean number of associated type-2 particles per type-1 particle.
//
// Usage:
//   ex_mass_distribution_connected.exe [--model hrg|blastwave] [--ini <file>]
//                            [--share <dir>] [--output <file>]
//                            [--events <N>] [--samples <N>] [--seed <N>]
//                            [--pid1 <pdg>] [--pid2 <pdg>]

#include "CLI11.hpp"
#include "Event.h"
#include "Integrator.h"
#include "Model.h"
#include "Model_BlastWave.h"
#include "Model_HRG.h"
#include "ParticleCoor.h"
#include "ParticleDB.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TRandom.h>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <vector>

// Rapidity-cut second-pass parameters (defaults; override with --ymax / --deltay).
static constexpr double kYmaxDefault   = 2.0;  // symmetric window is [-Ymax/2, +Ymax/2]
static constexpr double kDeltaYDefault = 1.0;  // minimum |y1 - y2|

static double invMass(const ParticleCoor& a, const ParticleCoor& b) {
  double e  = a.e  + b.e;
  double px = a.px + b.px;
  double py = a.py + b.py;
  double pz = a.pz + b.pz;
  double m2 = e*e - px*px - py*py - pz*pz;
  return (m2 > 0.) ? std::sqrt(m2) : 0.;
}

// Momentum-space rapidity, with the same high-rapidity guard used by the other
// src/work examples (ex_balance_function.cxx, ex_connected_correlations_2d.cxx).
static double getRapidity(const Particle& p) {
  double e = p.e, pz = p.pz;
  if (e <= std::abs(pz)) return (pz >= 0) ? 100.0 : -100.0;
  return 0.5 * std::log((e + pz) / (e - pz));
}

// Replace p's rapidity with a uniform draw in [-yHalf, +yHalf], keeping pT, phi
// and mass fixed (mT-based redraw; cf. randomizeRapidity in ex_mixed_events.cxx).
static double randomizeRapidity(ParticleCoor& p, double yHalf) {
  double pT   = std::sqrt(p.px * p.px + p.py * p.py);
  double mT   = std::sqrt(p.mass * p.mass + pT * pT);
  double yNew = gRandom->Uniform(0., yHalf);
  p.e  = mT * std::cosh(yNew);
  p.pz = mT * std::sinh(yNew);
  return yNew;
}

int main(int argc, char** argv) {
  std::string shareDir   = "";
  std::string outputFile = "";
  std::string modelName  = "hrg";
  std::string iniFile    = "";
  int nEvents  = 500;
  int nSamples = 1000000;
  int seed     = 0;
  int pid1     = 211;   // pi+
  int pid2     = -211;  // pi-
  double Ymax   = kYmaxDefault;
  double DeltaY = kDeltaYDefault;

  CLI::App app{"ex_mass_distribution_connected: correlated pair invariant mass distribution"};
  app.add_option("--share",   shareDir,   "path to SHARE particle database");
  auto* outOpt = app.add_option("--output", outputFile,
      "output ROOT file (default: mass_con/mass_con_<model>_<ymax>_<deltay>_<events>.root)");
  app.add_option("--model",   modelName,  "model: hrg (default) or blastwave")->default_str("hrg");
  app.add_option("--ini",     iniFile,    "ini file (required for blastwave model)");
  app.add_option("--events",  nEvents,    "number of events to generate")->default_val(500);
  app.add_option("--samples", nSamples,   "Monte Carlo integration samples")->default_val(1000000);
  app.add_option("--seed",    seed,       "random seed")->default_val(0);
  app.add_option("--pid1",    pid1,       "PDG code of particle type 1")->default_val(211);
  app.add_option("--pid2",    pid2,       "PDG code of particle type 2")->default_val(-211);
  app.add_option("--ymax",    Ymax,       "rapidity-cut pass: window is [-Ymax/2, +Ymax/2]")
      ->default_val(kYmaxDefault);
  app.add_option("--deltay",  DeltaY,     "rapidity-cut pass: minimum |y1 - y2|")
      ->default_val(kDeltaYDefault);
  CLI11_PARSE(app, argc, argv);

  // pT for blastwave (flow breaks boost invariance), |p| for HRG (isotropic).
  bool useTransverse = (modelName == "blastwave");

  if (outOpt->count() == 0) {
    std::string tag = useTransverse ? "blastwave" : "hrg";
    char buf[256];
    std::snprintf(buf, sizeof(buf), "mass_con/mass_con_%s_%g_%g_%d.root", tag.c_str(), Ymax, DeltaY,
                  nEvents);
    outputFile = buf;
  }

  // Make sure the output directory exists before ROOT opens the file.
  {
    std::filesystem::path outPath(outputFile);
    if (outPath.has_parent_path())
      std::filesystem::create_directories(outPath.parent_path());
  }
  std::string pVarLabel = useTransverse ? "p_{T,1} (GeV)" : "|p_{1}| (GeV)";

  gRandom->SetSeed(seed);

  std::unique_ptr<Model> model;
  if (modelName == "blastwave") {
    if (iniFile.empty()) {
      std::cerr << "Error: --ini is required for the blastwave model\n";
      return 1;
    }
    model = std::make_unique<Model_BlastWave>(iniFile.c_str());
  } else {
    model = std::make_unique<Model_HRG>();
  }

  auto tPartDB = std::make_unique<ParticleDB>(shareDir);
  auto tInteg  = std::make_unique<Integrator>(nSamples, model.get());
  tInteg->SetMultiplicities(tPartDB.get());
  auto tEvent  = std::make_unique<Event>(tPartDB.get(), tInteg.get());

  TFile fout(outputFile.c_str(), "RECREATE", "Correlated pair invariant mass");
  TH1::SetDefaultSumw2();
  TH1D* hMass = new TH1D("hMass",
      "Correlated pair invariant mass;M (GeV);(1/N_{1}) dN_{pairs}/dM",
      100, 0., 2.);

  const int    nPbins = 6;
  const double pMax   = 2.0;
  TH2D* hMass2D = new TH2D("hMass2D",
      (std::string(";") + pVarLabel + ";M (GeV)").c_str(),
      nPbins, 0., pMax, 100, 0., 2.);
  TH1D* hPid1Counts = new TH1D("hPid1Counts",
      (std::string(";") + pVarLabel + ";N_{1}").c_str(),
      nPbins, 0., pMax);

  // Rapidity-cut second pass: raw (un-normalized) pair counts.
  TH1D* hMassRapCut = new TH1D("hMassRapCut",
      "Correlated pair invariant mass (rapidity cut);M (GeV);dN_{pairs}/dM",
      100, 0., 2.);
  TH2D* hMass2DRapCut = new TH2D("hMass2DRapCut",
      (std::string(";") + pVarLabel + ";M (GeV)").c_str(),
      nPbins, 0., pMax, 100, 0., 2.);

  long nPid1Total = 0;

  for (int ievent = 0; ievent < nEvents; ievent++) {
    tEvent->Reset();
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();

    // Sort in place by rooteid so same-ancestor particles are contiguous.
    auto& particles = tEvent->GetParticleList();
    particles.sort([](const Particle& a, const Particle& b) {
      return a.rooteid < b.rooteid;
    });

    // Walk through rooteid groups and form all type1 x type2 pairs.
    auto it = particles.begin();
    while (it != particles.end()) {
      Int_t currentRooteid = it->rooteid;
      std::vector<const Particle*> group1, group2;

      while (it != particles.end() && it->rooteid == currentRooteid) {
        if (!it->GetDecayed()) {
          if (it->pid == pid1) {
            group1.push_back(&*it);
            ++nPid1Total;
            double pVar = useTransverse
                ? std::sqrt(it->px * it->px + it->py * it->py)
                : std::sqrt(it->px * it->px + it->py * it->py + it->pz * it->pz);
            hPid1Counts->Fill(pVar);
          }
          if (it->pid == pid2) {
            group2.push_back(&*it);
          }
        }
        ++it;
      }

      for (const auto* p1 : group1) {
        for (const auto* p2 : group2) {
          if (p1 != p2) {
            double m = invMass(*p1, *p2);
            hMass->Fill(m);
            double pVar1 = useTransverse
                ? std::sqrt(p1->px * p1->px + p1->py * p1->py)
                : std::sqrt(p1->px * p1->px + p1->py * p1->py + p1->pz * p1->pz);
            hMass2D->Fill(pVar1, m);
          }
        }
      }

      // Second pass: same correlated pairs, restricted to a symmetric
      // mid-rapidity window with a minimum rapidity separation.
      for (const auto* p1 : group1) {
        for (const auto* p2 : group2) {
          if (p1 == p2) continue;
          double y1 = getRapidity(*p1);
          double y2 = getRapidity(*p2);
          if (y1 < -Ymax / 2. || y1 > Ymax / 2.) continue;
          if (y2 < -Ymax / 2. || y2 > Ymax / 2.) continue;
          if (std::abs(y1 - y2) <= DeltaY) continue;

          // Scramble each partner's rapidity uniformly in the window before
          // forming the pair, so hMassRapCut carries no rapidity correlation.
          ParticleCoor q1 = *p1;
          ParticleCoor q2 = *p2;
          randomizeRapidity(q1, Ymax);
          randomizeRapidity(q2, Ymax);

          double m = invMass(q1, q2);
          hMassRapCut->Fill(m);
          double pVar1 = useTransverse
              ? std::sqrt(q1.px * q1.px + q1.py * q1.py)
              : std::sqrt(q1.px * q1.px + q1.py * q1.py + q1.pz * q1.pz);
          hMass2DRapCut->Fill(pVar1, m);
        }
      }
    }
  }

  if (nPid1Total > 0) {
    hMass->Scale(1.0 / nPid1Total);
  }

   if (nPid1Total > 0) {
    hMassRapCut->Scale(Ymax * Ymax / (Ymax - DeltaY) / (Ymax - DeltaY) / nPid1Total);
  }

  // Connected part: full correlated distribution minus the rapidity-cut estimate.
  TH1D* hMassConn = static_cast<TH1D*>(hMass->Clone("hMassConn"));
  hMassConn->SetTitle("Connected pair invariant mass (hMass - hMassRapCut);M (GeV);dN_{pairs}/dM");
  hMassConn->Add(hMassRapCut, -1.0);

  double yieldErr = 0.;
  double yield = hMass->IntegralAndError(0, hMass->GetNbinsX() + 1, yieldErr);
  std::cout << "Correlated yield (pairs per pid1 particle): "
            << yield << " +/- " << yieldErr << "\n";

  double cutErr = 0.;
  double cutCount = hMassRapCut->IntegralAndError(0, hMassRapCut->GetNbinsX() + 1, cutErr);
  std::cout << "Rapidity-cut raw pair count (Ymax=" << Ymax << ", DeltaY=" << DeltaY
            << "): " << cutCount << " +/- " << cutErr << "\n";

  // Generation parameters, for the plotting macro to annotate the figure.
  fout.cd();
  TNamed runModel("model", modelName.c_str());
  runModel.Write();
  TParameter<int>("nEvents", nEvents).Write();
  TParameter<double>("Ymax", Ymax).Write();
  TParameter<double>("DeltaY", DeltaY).Write();

  fout.Write();
  return 0;
}
