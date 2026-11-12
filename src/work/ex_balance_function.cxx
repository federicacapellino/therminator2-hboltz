// ex_balance_function.cxx — Charge Balance Function B(Delta y)
//
// Calculates the charge balance function as defined in STAR paper (arXiv:1005.2307):
//   B(Delta y) = 0.5 * [ (N_{+-} - N_{++})/N_+ + (N_{-+} - N_{--})/N_- ]
//
// Usage:
//   ex_balance_function.exe [--model hrg|blastwave|blastwaveaniso] [--ini <file>]
//                           [--share <dir>] [--output <file>]
//                           [--events <N>] [--samples <N>] [--seed <N>]
//                           [--ymax <val>] [--pid1 <pdg>] [--pid2 <pdg>]

#include "CLI11.hpp"
#include "Event.h"
#include "Integrator.h"
#include "Model.h"
#include "Model_BlastWave.h"
#include "Model_BlastWaveAniso.h"
#include "Model_HRG.h"
#include "ParticleDB.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TH1D.h>
#include <TRandom.h>
#include <cmath>
#include <iostream>
#include <vector>

static double getRapidity(const Particle& p) {
  double e  = p.e;
  double pz = p.pz;
  // Guard against numerical issues at high rapidity
  if (e <= std::abs(pz)) return (pz >= 0) ? 100.0 : -100.0;
  return 0.5 * std::log((e + pz) / (e - pz));
}

int main(int argc, char** argv) {
  std::string shareDir   = "";
  std::string outputFile = "balance_function.root";
  std::string modelName  = "hrg";
  std::string iniFile    = "";
  int nEvents  = 500;
  int nSamples = 1000000;
  int seed     = 0;
  double yMax  = 0.5;   // Mid-rapidity window: |y| < yMax
  int pid1     = 211;   // pi+
  int pid2     = -211;  // pi-

  CLI::App app{"ex_balance_function: charge balance function B(Delta y)"};
  app.add_option("--share",   shareDir,   "path to SHARE particle database");
  app.add_option("--output",  outputFile, "output ROOT file")->default_str("balance_function.root");
  app.add_option("--model",   modelName,  "model: hrg (default), blastwave or blastwaveaniso")
      ->default_str("hrg")
      ->check(CLI::IsMember({"hrg", "blastwave", "blastwaveaniso"}));
  app.add_option("--ini",     iniFile,    "ini file (required for the blast-wave models)");
  app.add_option("--events",  nEvents,    "number of events to generate")->default_val(500);
  app.add_option("--samples", nSamples,   "Monte Carlo integration samples")->default_val(1000000);
  app.add_option("--seed",    seed,       "random seed")->default_val(0);
  app.add_option("--ymax",    yMax,       "rapidity window cut |y| < yMax")->default_val(0.5);
  app.add_option("--pid1",    pid1,       "PDG code of positive species (e.g. 211 for pi+)")->default_val(211);
  app.add_option("--pid2",    pid2,       "PDG code of negative species (e.g. -211 for pi-)")->default_val(-211);
  CLI11_PARSE(app, argc, argv);

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

  TFile fout(outputFile.c_str(), "RECREATE", "Charge Balance Function");
  TH1::SetDefaultSumw2();

  // Histograms for numerator terms
  const int nBins = 50;
  const double dyMax = 2.0 * yMax; // Max possible delta_y within |y| < yMax window
  TH1D* hN_pm = new TH1D("hN_pm", "N_{+-}(#Delta y);#Delta y;Pairs", nBins, 0., dyMax);
  TH1D* hN_mp = new TH1D("hN_mp", "N_{-+}(#Delta y);#Delta y;Pairs", nBins, 0., dyMax);
  TH1D* hN_pp = new TH1D("hN_pp", "N_{++}(#Delta y);#Delta y;Pairs", nBins, 0., dyMax);
  TH1D* hN_mm = new TH1D("hN_mm", "N_{--}(#Delta y);#Delta y;Pairs", nBins, 0., dyMax);

  long nPlusTotal  = 0;
  long nMinusTotal = 0;

  for (int ievent = 0; ievent < nEvents; ievent++) {
    tEvent->Reset();
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();

    std::vector<const Particle*> posList;
    std::vector<const Particle*> negList;

    // 1. Filter final-state particles inside rapidity window |y| < yMax
    for (const auto& p : tEvent->GetParticleList()) {
      if (!p.GetDecayed()) {
        double y = getRapidity(p);
        if (std::abs(y) <= yMax) {
          if (p.pid == pid1) {
            posList.push_back(&p);
            ++nPlusTotal;
          } else if (p.pid == pid2) {
            negList.push_back(&p);
            ++nMinusTotal;
          }
        }
      }
    }

    // 2. Form Unlike-Sign Pairs: N_{+-} and N_{-+}
    for (const auto* p1 : posList) {
      for (const auto* p2 : negList) {
        double dy = std::abs(getRapidity(*p1) - getRapidity(*p2));
        hN_pm->Fill(dy); // Trigger = (+), Assoc = (-)
        hN_mp->Fill(dy); // Trigger = (-), Assoc = (+)
      }
    }

    // 3. Form Like-Sign Pairs: N_{++}
    for (size_t i = 0; i < posList.size(); ++i) {
      for (size_t j = i + 1; j < posList.size(); ++j) {
        double dy = std::abs(getRapidity(*posList[i]) - getRapidity(*posList[j]));
        // Each pair (i,j) contributes twice (i as trigger and j as trigger)
        hN_pp->Fill(dy, 2.0);
      }
    }

    // 4. Form Like-Sign Pairs: N_{--}
    for (size_t i = 0; i < negList.size(); ++i) {
      for (size_t j = i + 1; j < negList.size(); ++j) {
        double dy = std::abs(getRapidity(*negList[i]) - getRapidity(*negList[j]));
        hN_mm->Fill(dy, 2.0);
      }
    }
  }

  // 5. Construct Balance Function B(Delta y)
  TH1D* hBF = new TH1D("hBF", 
      "Balance Function B(#Delta y);#Delta y;B(#Delta y)", 
      nBins, 0., dyMax);

  if (nPlusTotal > 0 && nMinusTotal > 0) {
    // Term 1: (N_{+-} - N_{++}) / N_+
    TH1D* hTerm1 = (TH1D*) hN_pm->Clone("hTerm1");
    hTerm1->Add(hN_pp, -1.0);
    hTerm1->Scale(1.0 / nPlusTotal);

    // Term 2: (N_{-+} - N_{--}) / N_-
    TH1D* hTerm2 = (TH1D*) hN_mp->Clone("hTerm2");
    hTerm2->Add(hN_mm, -1.0);
    hTerm2->Scale(1.0 / nMinusTotal);

    // B(Delta y) = 0.5 * (Term1 + Term2)
    hBF->Add(hTerm1, hTerm2, 0.5, 0.5);

    delete hTerm1;
    delete hTerm2;
  }

  std::cout << "Processed " << nEvents << " events.\n";
  std::cout << "N_pos = " << nPlusTotal << ", N_neg = " << nMinusTotal << "\n";
  std::cout << "Integral B(Delta y) = " << hBF->Integral("width") << "\n";

  fout.Write();
  return 0;
}