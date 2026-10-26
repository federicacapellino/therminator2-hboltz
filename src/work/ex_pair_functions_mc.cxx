// ex_pair_functions_mc.cxx — HRG Monte Carlo for the pair-function comparisons
//
// See README_pair_functions.md, Sec. 4. For every ordered pair of distinct final particles of a
// measured species pair (spec Sec. 2) that share an ancestor, fills (M, kappa, cos theta*) at
// three exclusive ancestry levels:
//   L0  same father, father primordial (what the table contains)
//   L1x same father, father not primordial (L1 = L0 + L1x)
//   L2x same root ancestor only, i.e. cascades (L2 = L1 + L2x)
// Histograms hold raw counts; hEvents holds the number of events, so files from several seeds
// can be merged with hadd. The comparison applies the 1/(N_ev V) normalization.
// Also prints test T0.1 (weak switch). With --spectra (milestone M6, test T6.4) it also fills, for
// all particles in the list (decayed or not), the per-event count and count^2 of every species
// (hSpeciesN, hSpeciesN2, bins = ParticleDB index) and dN/dp of a few species (hSpec_<pdg>).
//
// Usage:
//   ex_pair_functions_mc.exe [--weak on|off] [--share <dir>] [--output <file>]
//                            [--events <N>] [--samples <N>] [--seed <N>] [--spectra]

#include "CLI11.hpp"
#include "Event.h"
#include "Integrator.h"
#include "Model_HRG.h"
#include "PairFunctions.h"
#include "ParticleDB.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TProfile2D.h>
#include <TRandom.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <map>
#include <unordered_map>

static constexpr int kNLevels = 3;
static const char* const kLevelTag[kNLevels] = {"L0", "L1x", "L2x"};
static constexpr int kNCos = 10;

// Measured species, for the single-particle counts.
static constexpr int kMeasuredPDG[] = {211, -211, 321, -321, 2212, -2212};

// Species with a dN/dp histogram (--spectra): pi+, K+, p, rho0, K*0, omega, Lambda, Delta++.
static constexpr int kSpectraPDG[] = {211, 321, 2212, 113, 313, 223, 3122, 2224};
static constexpr int kNSpectra     = sizeof(kSpectraPDG) / sizeof(kSpectraPDG[0]);

int main(int argc, char** argv) {
  std::string shareDir   = "";
  std::string outputFile = "pair_mc.root";
  std::string weak       = "on";
  int nEvents  = 100;
  int nSamples = 1000000;
  int seed     = 0;
  bool spectra = false;

  CLI::App app{"ex_pair_functions_mc: HRG pair densities at three ancestry levels"};
  app.add_option("--share",   shareDir,   "path to SHARE particle database");
  app.add_option("--output",  outputFile, "output ROOT file")->default_val("pair_mc.root");
  app.add_option("--weak",    weak,       "weak decays: on (default) or off")
      ->default_str("on")
      ->check(CLI::IsMember({"on", "off"}));
  app.add_option("--events",  nEvents,    "number of events to generate")->default_val(100);
  app.add_option("--samples", nSamples,   "Monte Carlo integration samples")->default_val(1000000);
  app.add_option("--seed",    seed,       "random seed")->default_val(0);
  app.add_flag("--spectra",   spectra,    "also fill species counts and momentum spectra (T6.4)");
  CLI11_PARSE(app, argc, argv);

  bool weakOn = (weak == "on");

  {
    std::filesystem::path outPath(outputFile);
    if (outPath.has_parent_path())
      std::filesystem::create_directories(outPath.parent_path());
  }

  gRandom->SetSeed(seed);

  auto model   = std::make_unique<Model_HRG>();
  auto tPartDB = std::make_unique<ParticleDB>(shareDir);
  int nStable  = SetWeakDecays(tPartDB.get(), weakOn);
  if (!weakOn && nStable != kNWeakTypes) {
    std::cerr << "Error: --weak off changed " << nStable << " particle types, expected "
              << kNWeakTypes << "\n";
    return 1;
  }
  auto tInteg = std::make_unique<Integrator>(nSamples, model.get());
  tInteg->SetMultiplicities(tPartDB.get());
  auto tEvent = std::make_unique<Event>(tPartDB.get(), tInteg.get());

  const double T = model->GetTemperature();
  const double V = model->GetVolume();
  std::cout << "HRG: T = " << T << " GeV, V = " << V << " GeV^-3 ("
            << V * kHbarC * kHbarC * kHbarC << " fm^3), weak decays " << weak << "\n";

  // ---- Histograms -------------------------------------------------------------------------------
  TFile fout(outputFile.c_str(), "RECREATE", "Pair functions MC");
  TH1D* hEvents = new TH1D("hEvents", "number of events", 1, 0., 1.);
  TH1D* hSingles = new TH1D("hSingles", "final measured particles;species;N", 6, 0., 6.);
  for (int s = 0; s < 6; s++)
    hSingles->GetXaxis()->SetBinLabel(s + 1, std::to_string(kMeasuredPDG[s]).c_str());

  std::vector<double> kEdges = MCKappaEdges();
  std::vector<double> cEdges(kNCos + 1);
  for (int c = 0; c <= kNCos; c++) cEdges[c] = -1. + 2. * c / kNCos;
  TH1D* hM[kNLevels][kNPairs];
  TH3D* hMKC[kNLevels][kNPairs];
  TProfile2D* hP1[kNLevels][kNPairs];
  TProfile2D* hP2[kNLevels][kNPairs];
  for (int p = 0; p < kNPairs; p++) {
    auto masses = PairMasses(tPartDB.get(), p);
    std::vector<double> mEdges = MassEdges(masses.first + masses.second);
    const int nM = mEdges.size() - 1, nK = kEdges.size() - 1;
    for (int l = 0; l < kNLevels; l++) {
      std::string sfx = std::string(kLevelTag[l]) + "_" + kPairTag[p];
      hM[l][p] = new TH1D(("hM_" + sfx).c_str(), (sfx + ";M (GeV);pairs").c_str(), nM,
                          mEdges.data());
      hMKC[l][p] = new TH3D(("hMKC_" + sfx).c_str(),
                            (sfx + ";M (GeV);#kappa (GeV);cos#theta*").c_str(), nM, mEdges.data(),
                            nK, kEdges.data(), kNCos, cEdges.data());
      hP1[l][p] = new TProfile2D(("hP1_" + sfx).c_str(),
                                 (sfx + " <P_{1}(cos#theta*)>;M (GeV);#kappa (GeV)").c_str(), nM,
                                 mEdges.data(), nK, kEdges.data());
      hP2[l][p] = new TProfile2D(("hP2_" + sfx).c_str(),
                                 (sfx + " <P_{2}(cos#theta*)>;M (GeV);#kappa (GeV)").c_str(), nM,
                                 mEdges.data(), nK, kEdges.data());
    }
  }

  // Species counts and spectra (--spectra), all particles in the list, decayed or not.
  const int nTypes = tPartDB->GetParticleTypeCount();
  std::unordered_map<const ParticleType*, int> typeIndex;
  for (int i = 0; i < nTypes; i++) typeIndex[tPartDB->GetParticleType(i)] = i;
  TH1D* hSpeciesN  = nullptr;
  TH1D* hSpeciesN2 = nullptr;
  TH1D* hSpec[kNSpectra] = {nullptr};
  std::vector<long> evCount(nTypes);
  if (spectra) {
    hSpeciesN  = new TH1D("hSpeciesN", "sum over events of the count;ParticleDB index;N", nTypes,
                          0., nTypes);
    hSpeciesN2 = new TH1D("hSpeciesN2", "sum over events of count^{2};ParticleDB index;N^{2}",
                          nTypes, 0., nTypes);
    for (int i = 0; i < nTypes; i++) {
      hSpeciesN->GetXaxis()->SetBinLabel(i + 1, tPartDB->GetParticleType(i)->GetName());
      hSpeciesN2->GetXaxis()->SetBinLabel(i + 1, tPartDB->GetParticleType(i)->GetName());
    }
    for (int s = 0; s < kNSpectra; s++) {
      std::string tag = std::to_string(kSpectraPDG[s]);
      hSpec[s] = new TH1D(("hSpec_" + tag).c_str(), ("PDG " + tag + ";p (GeV);N").c_str(), 300,
                          0., 3.);
    }
  }

  // T0.1 bookkeeping, per weak species (signed PDG code): final and decayed counts.
  std::map<int, long> nFinal, nDecayed;
  for (int code : kWeakPDG) {
    nFinal[code] = nFinal[-code] = 0;
    nDecayed[code] = nDecayed[-code] = 0;
  }
  long nWeakDaughters = 0; // final particles whose father is a weak species
  long nPairs[kNLevels] = {0, 0, 0};

  std::unordered_map<int, const Particle*> byEid;
  std::unordered_map<int, std::vector<const Particle*>> byRoot;

  for (int ievent = 0; ievent < nEvents; ievent++) {
    tEvent->Reset();
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();
    hEvents->Fill(0.5);

    byEid.clear();
    byRoot.clear();
    if (spectra) std::fill(evCount.begin(), evCount.end(), 0);
    for (const Particle& p : tEvent->GetParticleList()) {
      byEid[p.eid] = &p;
      if (spectra) {
        ++evCount[typeIndex.at(p.GetParticleType())];
        for (int s = 0; s < kNSpectra; s++)
          if (p.pid == kSpectraPDG[s])
            hSpec[s]->Fill(std::sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz));
      }
      if (IsWeakSpecies(p.pid)) {
        if (p.GetDecayed())
          ++nDecayed[p.pid];
        else
          ++nFinal[p.pid];
      }
      if (p.GetDecayed()) continue;
      if (p.fathereid != -1 && IsWeakSpecies(p.fatherpid)) ++nWeakDaughters;
      for (int s = 0; s < 6; s++) {
        if (p.pid == kMeasuredPDG[s]) {
          hSingles->Fill(s + 0.5);
          byRoot[p.rooteid].push_back(&p);
        }
      }
    }

    if (spectra) {
      for (int i = 0; i < nTypes; i++) {
        if (evCount[i] == 0) continue;
        hSpeciesN->AddBinContent(i + 1, evCount[i]);
        hSpeciesN2->AddBinContent(i + 1, static_cast<double>(evCount[i]) * evCount[i]);
      }
    }

    for (auto& kv : byRoot) {
      const std::vector<const Particle*>& group = kv.second;
      for (const Particle* p1 : group) {
        for (const Particle* p2 : group) {
          if (p1 == p2) continue;
          int pair = PairIndex(p1->pid, p2->pid);
          if (pair < 0) continue;
          int level = 2;
          if (p1->fathereid != -1 && p1->fathereid == p2->fathereid)
            level = (byEid.at(p1->fathereid)->fathereid == -1) ? 0 : 1;

          TLorentzVector v1(p1->px, p1->py, p1->pz, p1->e);
          TLorentzVector v2(p2->px, p2->py, p2->pz, p2->e);
          TLorentzVector K = v1 + v2;
          double kappa = K.P();
          if (kappa < 1e-6) continue;
          double M = K.M();
          v1.Boost(-K.BoostVector());
          double cth = v1.Vect().Dot(K.Vect()) / (v1.P() * kappa);

          ++nPairs[level];
          hM[level][pair]->Fill(M);
          hMKC[level][pair]->Fill(M, kappa, cth);
          hP1[level][pair]->Fill(M, kappa, cth);
          hP2[level][pair]->Fill(M, kappa, 0.5 * (3. * cth * cth - 1.));
        }
      }
    }
  }
  std::printf("Pairs: L0 %ld, L1x %ld, L2x %ld (all measured ordered pairs, %d events)\n",
              nPairs[0], nPairs[1], nPairs[2], nEvents);

  // ---- T0.1 report ------------------------------------------------------------------------------
  TH1D* hWeakFinal = new TH1D("hWeakFinal", "Final weakly decaying particles;species;N", 2 * kNWeakPDG,
                              0., 2. * kNWeakPDG);
  bool pass = true;
  int bin   = 1;
  std::printf("T0.1  %8s %10s %10s\n", "PDG", "final", "decayed");
  for (int code : kWeakPDG) {
    for (int sgn : {+1, -1}) {
      int pdg = sgn * code;
      std::printf("T0.1  %8d %10ld %10ld\n", pdg, nFinal[pdg], nDecayed[pdg]);
      hWeakFinal->SetBinContent(bin, nFinal[pdg]);
      hWeakFinal->GetXaxis()->SetBinLabel(bin++, std::to_string(pdg).c_str());
      if (nFinal[pdg] == 0 || nDecayed[pdg] != 0) pass = false;
    }
  }
  std::printf("T0.1  final particles with a weakly decaying father: %ld\n", nWeakDaughters);
  if (nWeakDaughters != 0) pass = false;
  if (!weakOn)
    std::cout << "T0.1  " << (pass ? "PASS" : "FAIL") << "\n";
  else
    std::cout << "T0.1  (weak on: reference run, test applies to --weak off)\n";

  fout.cd();
  TNamed("weak", weak.c_str()).Write();
  TNamed("dbhash", ParticleDBHash(tPartDB.get()).c_str()).Write();
  TNamed("modelhash", model->GetHash()).Write();
  TParameter<double>("T", T).Write();
  TParameter<double>("V", V).Write();
  fout.Write();
  return (!weakOn && !pass) ? 2 : 0;
}
