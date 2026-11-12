// ex_phi_distribution_connected.cxx — correlated pair azimuthal-difference distribution
//
// Generates thermal events and histograms Delta phi = phi1 - phi2, wrapped to
// [-pi/2, 3pi/2), of final-state
// particle pairs that share the same primordial ancestor (rooteid).  Two
// non-decayed particles are "correlated" if they were both produced in the
// decay chain of the same primordial resonance.
//
// After all events the histogram is normalized by the total number of type-1
// particles produced.  The integral of the normalized histogram is therefore
// the mean number of associated type-2 particles per type-1 particle.
//
// Usage:
//   ex_phi_distribution_connected.exe [--model hrg|blastwave|blastwaveaniso] [--ini <file>]
//                            [--share <dir>] [--output <file>]
//                            [--events <N>] [--samples <N>] [--seed <N>]
//                            [--pid1 <pdg>] [--pid2 <pdg>]

#include "CLI11.hpp"
#include "Event.h"
#include "Integrator.h"
#include "Model.h"
#include "Model_BlastWave.h"
#include "Model_BlastWaveAniso.h"
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

// Azimuthal difference phi(a) - phi(b), wrapped to [-pi/2, 3pi/2).
static double deltaPhi(const ParticleCoor& a, const ParticleCoor& b) {
  double d = std::atan2(a.py, a.px) - std::atan2(b.py, b.px);
  while (d < -M_PI / 2.) d += 2. * M_PI;
  while (d >= 3. * M_PI / 2.) d -= 2. * M_PI;
  return d;
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

  CLI::App app{"ex_phi_distribution_connected: correlated pair Delta phi distribution"};
  app.add_option("--share",   shareDir,   "path to SHARE particle database");
  auto* outOpt = app.add_option("--output", outputFile,
      "output ROOT file (default: phi_con/phi_con_<model>_<ymax>_<deltay>_<events>.root)");
  app.add_option("--model",   modelName,  "model: hrg (default), blastwave or blastwaveaniso")
      ->default_str("hrg")
      ->check(CLI::IsMember({"hrg", "blastwave", "blastwaveaniso"}));
  app.add_option("--ini",     iniFile,    "ini file (required for the blast-wave models)");
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

  // pT for the blast-wave models (flow breaks boost invariance), |p| for HRG (isotropic).
  bool useTransverse = (modelName != "hrg");

  if (outOpt->count() == 0) {
    std::string tag = modelName;
    char buf[256];
    std::snprintf(buf, sizeof(buf), "phi_con/phi_con_%s_%g_%g_%d.root", tag.c_str(), Ymax, DeltaY,
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
  if (modelName == "blastwave" || modelName == "blastwaveaniso") {
    if (iniFile.empty()) {
      std::cerr << "Error: --ini is required for the " << modelName << " model\n";
      return 1;
    }
    if (modelName == "blastwave")
      model = std::make_unique<Model_BlastWave>(iniFile.c_str());
    else
      model = std::make_unique<Model_BlastWaveAniso>(iniFile.c_str());
  } else {
    model = std::make_unique<Model_HRG>();
  }

  auto tPartDB = std::make_unique<ParticleDB>(shareDir);
  auto tInteg  = std::make_unique<Integrator>(nSamples, model.get());
  tInteg->SetMultiplicities(tPartDB.get());
  auto tEvent  = std::make_unique<Event>(tPartDB.get(), tInteg.get());

  TFile fout(outputFile.c_str(), "RECREATE", "Correlated pair Delta phi");
  TH1::SetDefaultSumw2();
  TH1D* hMass = new TH1D("hMass",
      "Correlated pair #Delta#phi;#Delta#phi (rad);(1/N_{1}) dN_{pairs}/d#Delta#phi",
      100, -M_PI / 2., 3. * M_PI / 2.);

  const int    nPbins = 6;
  const double pMax   = 2.0;
  TH2D* hMass2D = new TH2D("hMass2D",
      (std::string(";") + pVarLabel + ";#Delta#phi (rad)").c_str(),
      nPbins, 0., pMax, 100, -M_PI / 2., 3. * M_PI / 2.);
  TH1D* hPid1Counts = new TH1D("hPid1Counts",
      (std::string(";") + pVarLabel + ";N_{1}").c_str(),
      nPbins, 0., pMax);

  // Rapidity-cut second pass: raw (un-normalized) pair counts.
  TH1D* hMassRapCut = new TH1D("hMassRapCut",
      "Correlated pair #Delta#phi (rapidity cut);#Delta#phi (rad);dN_{pairs}/d#Delta#phi",
      100, -M_PI / 2., 3. * M_PI / 2.);
  TH2D* hMass2DRapCut = new TH2D("hMass2DRapCut",
      (std::string(";") + pVarLabel + ";#Delta#phi (rad)").c_str(),
      nPbins, 0., pMax, 100, -M_PI / 2., 3. * M_PI / 2.);

  // Window-normalized variants: per pid1 particle inside |y| < Ymax/2 instead of per
  // pid1 over all rapidities, so models with different rapidity coverage compare fairly.
  TH1D* hMassWin = new TH1D("hMassWin",
      "Correlated pair #Delta#phi, trigger in window (per N_{1} in |y|<Y_{max}/2);#Delta#phi (rad);(1/N_{1}^{win}) dN_{pairs}/d#Delta#phi",
      100, -M_PI / 2., 3. * M_PI / 2.);

  long nPid1Total = 0;
  long nPid1Win   = 0;

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
            if (std::abs(getRapidity(*it)) <= Ymax / 2.) ++nPid1Win;
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
            double dphi = deltaPhi(*p1, *p2);
            hMass->Fill(dphi);
            if (std::abs(getRapidity(*p1)) <= Ymax / 2.) hMassWin->Fill(dphi);
            double pVar1 = useTransverse
                ? std::sqrt(p1->px * p1->px + p1->py * p1->py)
                : std::sqrt(p1->px * p1->px + p1->py * p1->py + p1->pz * p1->pz);
            hMass2D->Fill(pVar1, dphi);
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

          double dphi = deltaPhi(q1, q2);
          hMassRapCut->Fill(dphi);
          double pVar1 = useTransverse
              ? std::sqrt(q1.px * q1.px + q1.py * q1.py)
              : std::sqrt(q1.px * q1.px + q1.py * q1.py + q1.pz * q1.pz);
          hMass2DRapCut->Fill(pVar1, dphi);
        }
      }
    }
  }

  if (nPid1Total > 0) {
    hMass->Scale(1.0 / nPid1Total);
  }

  // Same raw rapidity-cut pairs, normalized per pid1 inside the window (clone before
  // hMassRapCut is scaled).
  TH1D* hMassRapCutWin = static_cast<TH1D*>(hMassRapCut->Clone("hMassRapCutWin"));
  hMassRapCutWin->SetTitle(
      "Correlated pair #Delta#phi (rapidity cut) (per N_{1} in |y|<Y_{max}/2);#Delta#phi (rad);(1/N_{1}^{win}) dN_{pairs}/d#Delta#phi");

   if (nPid1Total > 0) {
    hMassRapCut->Scale(Ymax * Ymax / (Ymax - DeltaY) / (Ymax - DeltaY) / nPid1Total);
  }

  if (nPid1Win > 0) {
    hMassWin->Scale(1.0 / nPid1Win);
    hMassRapCutWin->Scale(Ymax * Ymax / (Ymax - DeltaY) / (Ymax - DeltaY) / nPid1Win);
  }

  // Connected part: full correlated distribution minus the rapidity-cut estimate.
  TH1D* hMassConn = static_cast<TH1D*>(hMass->Clone("hMassConn"));
  hMassConn->SetTitle("Connected pair #Delta#phi (hMass - hMassRapCut);#Delta#phi (rad);dN_{pairs}/d#Delta#phi");
  hMassConn->Add(hMassRapCut, -1.0);

  double yieldErr = 0.;
  double yield = hMass->IntegralAndError(0, hMass->GetNbinsX() + 1, yieldErr);
  std::cout << "Correlated yield (pairs per pid1 particle): "
            << yield << " +/- " << yieldErr << "\n";

  double cutErr = 0.;
  double cutCount = hMassRapCut->IntegralAndError(0, hMassRapCut->GetNbinsX() + 1, cutErr);
  std::cout << "Rapidity-cut raw pair count (Ymax=" << Ymax << ", DeltaY=" << DeltaY
            << "): " << cutCount << " +/- " << cutErr << "\n";

  // Window-normalized connected part: hMassWin - hMassRapCutWin.
  TH1D* hMassConnWin = static_cast<TH1D*>(hMassWin->Clone("hMassConnWin"));
  hMassConnWin->SetTitle(
      "Connected pair #Delta#phi (hMassWin - hMassRapCutWin) (per N_{1} in |y|<Y_{max}/2);#Delta#phi (rad);(1/N_{1}^{win}) dN_{pairs}/d#Delta#phi");
  hMassConnWin->Add(hMassRapCutWin, -1.0);

  std::cout << "pid1 in window |y| <= " << Ymax / 2. << ": " << nPid1Win << " of " << nPid1Total
            << "\n";
  double winErr = 0., winCutErr = 0.;
  double winYield = hMassWin->IntegralAndError(0, hMassWin->GetNbinsX() + 1, winErr);
  double winCut   = hMassRapCutWin->IntegralAndError(0, hMassRapCutWin->GetNbinsX() + 1, winCutErr);
  std::cout << "Window-normalized correlated yield: " << winYield << " +/- " << winErr << "\n";
  std::cout << "Window-normalized rapidity-cut count: " << winCut << " +/- " << winCutErr << "\n";

  // Generation parameters, for the plotting macro to annotate the figure.
  fout.cd();
  TNamed runModel("model", modelName.c_str());
  runModel.Write();
  TParameter<int>("nEvents", nEvents).Write();
  TParameter<double>("Ymax", Ymax).Write();
  TParameter<double>("DeltaY", DeltaY).Write();
  TParameter<double>("nPid1", nPid1Total).Write();
  TParameter<double>("nPid1Win", nPid1Win).Write();

  fout.Write();
  return 0;
}
