// ex_mass_pt_distribution.cxx — 2D (pT1, M) distribution of pairs from a primordial parent
//
// Generates thermal events and histograms the invariant mass of final-state
// particle pairs versus the transverse momentum of particle 1.  Only pairs whose
// members are both direct daughters of the same primordial parent resonance
// (default: omega(782), PDG 223) are kept; parents produced in the decay of a
// heavier resonance (feed-down) are excluded.
//
// A particle passes the filter if fatherpid == parent && rootpid == parent: its
// direct parent is the resonance and that resonance is primordial.  In that case
// fathereid == rooteid, so grouping by rooteid guarantees both daughters come
// from the same parent.
//
// With --parent 0 the parent filter is disabled: every pair whose members share
// the same primordial ancestor (rooteid) is kept, whatever its species and
// however long the decay chain.
//
// After all events the histograms are normalized by the total number of type-1
// particles produced (all of them, not only those passing the filter).  The
// integral is therefore the mean number of parent-correlated type-2 particles
// per type-1 particle.
//
// hMeanPt1VsM is a TProfile of the mean pT of particle 1 versus the pair
// invariant mass M, averaged over all parent-correlated pairs in all events:
// <pT1>(M) = sum(pT1) / N_pairs in each M bin.  It is not rescaled.
//
// hV2PtM (TProfile2D over pT1, M) and hV2VsM (TProfile over M, integrated over pT1) hold the
// elliptic flow of the pair, v2 = <cos 2(phi_pair - Psi_RP)>, with phi_pair the azimuth of
// p1 + p2.  The reaction plane is the x-z plane in every model, so Psi_RP = 0 exactly.  For
// like-sign pairs (pid1 == pid2) every pair enters twice, once per ordering; both orderings have
// the same phi_pair and M, so hV2VsM is unaffected, while hV2PtM gets one entry at each member's pT.
// Like hMeanPt1VsM, these profiles are not rescaled.
//
// Usage:
//   ex_mass_pt_distribution.exe [--model hrg|blastwave|blastwaveaniso] [--ini <file>]
//                               [--share <dir>] [--output <file>]
//                               [--events <N>] [--samples <N>] [--seed <N>]
//                               [--pid1 <pdg>] [--pid2 <pdg>] [--parent <pdg>]

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
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom.h>
#include <TString.h>
#include <cmath>
#include <iostream>
#include <vector>

static double invMass(const Particle& a, const Particle& b) {
  double e  = a.e  + b.e;
  double px = a.px + b.px;
  double py = a.py + b.py;
  double pz = a.pz + b.pz;
  double m2 = e*e - px*px - py*py - pz*pz;
  return (m2 > 0.) ? std::sqrt(m2) : 0.;
}

static double transverseMomentum(const Particle& p) {
  return std::sqrt(p.px * p.px + p.py * p.py);
}

int main(int argc, char** argv) {
  std::string shareDir   = "";
  std::string outputFile = "";
  std::string modelName  = "hrg";
  std::string iniFile    = "";
  int nEvents   = 500;
  int nSamples  = 1000000;
  int seed      = 0;
  int pid1      = 211;   // pi+
  int pid2      = -211;  // pi-
  int parentPid = 223;   // omega(782)

  CLI::App app{"ex_mass_pt_distribution: 2D (pT1, M) of pairs from a primordial parent"};
  app.add_option("--share",   shareDir,   "path to SHARE particle database");
  auto* outOpt = app.add_option("--output", outputFile,
      "output ROOT file (default: mass_pt_<model>.root)");
  app.add_option("--model",   modelName,  "model: hrg (default), blastwave or blastwaveaniso")
      ->default_str("hrg")
      ->check(CLI::IsMember({"hrg", "blastwave", "blastwaveaniso"}));
  app.add_option("--ini",     iniFile,    "ini file (required for the blast-wave models)");
  app.add_option("--events",  nEvents,    "number of events to generate")->default_val(500);
  app.add_option("--samples", nSamples,   "Monte Carlo integration samples")->default_val(1000000);
  app.add_option("--seed",    seed,       "random seed")->default_val(0);
  app.add_option("--pid1",    pid1,       "PDG code of particle type 1")->default_val(211);
  app.add_option("--pid2",    pid2,       "PDG code of particle type 2")->default_val(-211);
  app.add_option("--parent",  parentPid,  "PDG code of the primordial parent (0: any ancestor)")->default_val(223);
  CLI11_PARSE(app, argc, argv);

  if (outOpt->count() == 0) {
    std::string tag = modelName;
    outputFile = "mass_pt_" + tag + ".root";
  }

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

  TFile fout(outputFile.c_str(), "RECREATE", "Parent-correlated pair invariant mass vs pT1");
  TH1::SetDefaultSumw2();

  const int    nPtBins = 20;
  const double ptMax   = 2.0;
  const int    nMBins  = 100;
  const double mMax    = 2.0;

  const TString sel = (parentPid == 0) ? "Ancestor-correlated" : "Parent-correlated";
  TH1D* hMass = new TH1D("hMass",
      sel + " pair invariant mass;M (GeV);(1/N_{1}) dN_{pairs}",
      nMBins, 0., mMax);
  TH2D* hMassPt = new TH2D("hMassPt",
      sel + " pairs;p_{T,1} (GeV);M (GeV);(1/N_{1}) dN_{pairs}",
      nPtBins, 0., ptMax, nMBins, 0., mMax);
  TProfile* hMeanPt1VsM = new TProfile("hMeanPt1VsM",
      sel + " pairs;M (GeV);#LTp_{T,1}#GT (GeV)",
      nMBins, 0., mMax);
  TProfile2D* hV2PtM = new TProfile2D("hV2PtM",
      sel + " pairs;p_{T,1} (GeV);M (GeV);v_{2}^{pair}",
      nPtBins, 0., ptMax, nMBins, 0., mMax);
  TProfile* hV2VsM = new TProfile("hV2VsM",
      sel + " pairs;M (GeV);v_{2}^{pair}",
      nMBins, 0., mMax);
  TH1D* hPid1Counts = new TH1D("hPid1Counts",
      ";p_{T,1} (GeV);N_{1}",
      nPtBins, 0., ptMax);

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

    // Walk through rooteid groups and form all type1 x type2 pairs among
    // direct daughters of a primordial parent.
    auto it = particles.begin();
    while (it != particles.end()) {
      Int_t currentRooteid = it->rooteid;
      std::vector<const Particle*> group1, group2;

      while (it != particles.end() && it->rooteid == currentRooteid) {
        if (!it->GetDecayed()) {
          bool fromParent = parentPid == 0 ||
                            (it->fatherpid == parentPid && it->rootpid == parentPid);
          if (it->pid == pid1) {
            ++nPid1Total;
            hPid1Counts->Fill(transverseMomentum(*it));
            if (fromParent) group1.push_back(&*it);
          }
          if (it->pid == pid2 && fromParent) {
            group2.push_back(&*it);
          }
        }
        ++it;
      }

      for (const auto* p1 : group1) {
        for (const auto* p2 : group2) {
          if (p1 != p2) {
            double m   = invMass(*p1, *p2);
            double pt1 = transverseMomentum(*p1);
            hMass->Fill(m);
            hMassPt->Fill(pt1, m);
            hMeanPt1VsM->Fill(m, pt1);
            double cos2phi = std::cos(2.0 * std::atan2(p1->py + p2->py, p1->px + p2->px));
            hV2PtM->Fill(pt1, m, cos2phi);
            hV2VsM->Fill(m, cos2phi);
          }
        }
      }
    }
  }

  if (nPid1Total > 0) {
    hMass->Scale(1.0 / nPid1Total);
    hMassPt->Scale(1.0 / nPid1Total);
  }

  double yieldErr = 0.;
  double yield = hMass->IntegralAndError(0, hMass->GetNbinsX() + 1, yieldErr);
  double yield2DErr = 0.;
  double yield2D = hMassPt->IntegralAndError(0, hMassPt->GetNbinsX() + 1,
                                             0, hMassPt->GetNbinsY() + 1, yield2DErr);
  std::cout << sel << " yield (pairs per pid1 particle), 1D: "
            << yield << " +/- " << yieldErr << "\n";
  std::cout << sel << " yield (pairs per pid1 particle), 2D: "
            << yield2D << " +/- " << yield2DErr << "\n";
  std::cout << "Mean pT of pid1 over all pairs: "
            << hMeanPt1VsM->GetMean(2) << " +/- " << hMeanPt1VsM->GetMeanError(2)
            << " GeV\n";
  std::cout << "Pair v2 integrated over pT1 and M: "
            << hV2VsM->GetMean(2) << " +/- " << hV2VsM->GetMeanError(2) << "\n";

  fout.Write();
  return 0;
}
