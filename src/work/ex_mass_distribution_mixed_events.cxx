// ex_mass_distribution_mixed_events.cxx — same-event vs mixed-event invariant mass
//
// Generates thermal events and fills two invariant mass histograms:
//
//   hMass   — same-event (pid1, pid2) pairs that share the same primordial
//             ancestor (rooteid).  Two non-decayed particles are "correlated" if
//             they were both produced in the decay chain of the same primordial
//             resonance.  hMass is normalized by the total number of pid1
//             particles produced; its integral is the mean number of associated
//             pid2 particles per pid1 particle.
//
//   hMassME — mixed-event background.  Events are generated two at a time: pid1
//             particles are collected from event A, pid2 particles from event B,
//             and an A-pid1 x B-pid2 pair is histogrammed only when the two
//             partners carry the same rooteid.  This mirrors hMass's
//             same-ancestor combinatorics, but with the partners drawn from
//             independent events so no real correlation survives.  It is
//             normalized by the number of event-A pid1 particles.
//
// Note: events are consumed in pairs, so an odd --events is rounded down to the
// nearest even number.
//
// Usage:
//   ex_mass_distribution_mixed_events.exe [--model hrg|blastwave] [--ini <file>]
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
#include <list>
#include <unordered_map>
#include <vector>

static double invMass(const ParticleCoor& a, const ParticleCoor& b) {
  double e  = a.e  + b.e;
  double px = a.px + b.px;
  double py = a.py + b.py;
  double pz = a.pz + b.pz;
  double m2 = e*e - px*px - py*py - pz*pz;
  return (m2 > 0.) ? std::sqrt(m2) : 0.;
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

  CLI::App app{"ex_mass_distribution_mixed_events: same-event vs mixed-event invariant mass"};
  app.add_option("--share",   shareDir,   "path to SHARE particle database");
  auto* outOpt = app.add_option("--output", outputFile,
      "output ROOT file (default: mass_me/mass_me_<model>_<events>.root)");
  app.add_option("--model",   modelName,  "model: hrg (default) or blastwave")->default_str("hrg");
  app.add_option("--ini",     iniFile,    "ini file (required for blastwave model)");
  app.add_option("--events",  nEvents,    "number of events to generate")->default_val(500);
  app.add_option("--samples", nSamples,   "Monte Carlo integration samples")->default_val(1000000);
  app.add_option("--seed",    seed,       "random seed")->default_val(0);
  app.add_option("--pid1",    pid1,       "PDG code of particle type 1")->default_val(211);
  app.add_option("--pid2",    pid2,       "PDG code of particle type 2")->default_val(-211);
  CLI11_PARSE(app, argc, argv);

  // pT for blastwave (flow breaks boost invariance), |p| for HRG (isotropic).
  bool useTransverse = (modelName == "blastwave");

  if (outOpt->count() == 0) {
    std::string tag = useTransverse ? "blastwave" : "hrg";
    char buf[256];
    std::snprintf(buf, sizeof(buf), "mass_me/mass_me_%s_%d.root", tag.c_str(), nEvents);
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

  TFile fout(outputFile.c_str(), "RECREATE", "Same-event vs mixed-event invariant mass");
  TH1::SetDefaultSumw2();
  TH1D* hMass = new TH1D("hMass",
      "Correlated pair invariant mass;M (GeV);(1/N_{1}) dN_{pairs}/dM",
      100, 0., 2.);

  TH1D* hMassME = new TH1D("hMassME",
      "Mixed-event same-rooteid pair invariant mass;M (GeV);(1/N_{1}^{ME}) dN_{pairs}/dM",
      100, 0., 2.);

  const int    nPbins = 6;
  const double pMax   = 2.0;
  TH2D* hMass2D = new TH2D("hMass2D",
      (std::string(";") + pVarLabel + ";M (GeV)").c_str(),
      nPbins, 0., pMax, 100, 0., 2.);
  TH1D* hPid1Counts = new TH1D("hPid1Counts",
      (std::string(";") + pVarLabel + ";N_{1}").c_str(),
      nPbins, 0., pMax);

  long nPid1Total = 0;
  long nPid1ME    = 0;

  // Same-event correlated-pair fill: sort the particle list by rooteid so
  // same-ancestor particles are contiguous, then form every type1 x type2 pair
  // within each rooteid group.  Applied to every generated event.
  auto fillCorrelated = [&](std::list<Particle>& particles) {
    particles.sort([](const Particle& a, const Particle& b) {
      return a.rooteid < b.rooteid;
    });

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
    }
  };

  // Generate events two at a time: event A supplies pid1, event B supplies pid2
  // for the mixed-event histogram.  Both events also feed the same-event fill,
  // so hMass accumulates over all nEvents generated events.
  for (int ievent = 0; ievent + 1 < nEvents; ievent += 2) {
    // --- event A ---
    tEvent->Reset();
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();
    fillCorrelated(tEvent->GetParticleList());

    // Copy event A's non-decayed pid1 particles before tEvent->Reset() below
    // invalidates the list (Particle -> ParticleCoor slice-copy).
    std::vector<ParticleCoor> group1;
    for (const auto& part : tEvent->GetParticleList()) {
      if (!part.GetDecayed() && part.pid == pid1)
        group1.push_back(part);
    }

    // --- event B ---
    tEvent->Reset();
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();
    fillCorrelated(tEvent->GetParticleList());

    std::vector<ParticleCoor> group2;
    for (const auto& part : tEvent->GetParticleList()) {
      if (!part.GetDecayed() && part.pid == pid2)
        group2.push_back(part);
    }

    // --- mixed-event pairs: event-A pid1 x event-B pid2, same rooteid only ---
    std::unordered_map<Int_t, std::vector<const ParticleCoor*>> group2ByRoot;
    for (const auto& p2 : group2) {
      group2ByRoot[p2.rooteid].push_back(&p2);
    }
    for (const auto& p1 : group1) {
      auto mit = group2ByRoot.find(p1.rooteid);
      if (mit == group2ByRoot.end()) continue;
      for (const auto* p2 : mit->second) {
        hMassME->Fill(invMass(p1, *p2));
      }
    }
    nPid1ME += static_cast<long>(group1.size());
  }

  if (nPid1Total > 0) {
    hMass->Scale(1.0 / nPid1Total);
  }
  if (nPid1ME > 0) {
    hMassME->Scale(1.0 / nPid1ME);
  }

  double yieldErr = 0.;
  double yield = hMass->IntegralAndError(0, hMass->GetNbinsX() + 1, yieldErr);
  std::cout << "Correlated yield (pairs per pid1 particle): "
            << yield << " +/- " << yieldErr << "\n";

  double meErr = 0.;
  double meYield = hMassME->IntegralAndError(0, hMassME->GetNbinsX() + 1, meErr);
  std::cout << "Mixed-event (same rooteid) yield (pairs per pid1 particle): "
            << meYield << " +/- " << meErr << "\n";

  // Generation parameters, for the plotting macro to annotate the figure.
  fout.cd();
  TNamed runModel("model", modelName.c_str());
  runModel.Write();
  TParameter<int>("nEvents", nEvents).Write();

  fout.Write();
  return 0;
}
