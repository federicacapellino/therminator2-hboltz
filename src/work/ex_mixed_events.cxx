// ex_mixed_events.cxx — signal and background invariant mass distributions
//
// Generates thermal events and fills two invariant mass histograms:
//
//   hSignalAndBackground — all same-event (pid1, pid2) pairs within |y| < ycut
//
//   hBackground — uncorrelated background, estimated by one of two methods:
//     randomize (default): same-event pairs with |Δy| ≥ DY, both rapidities
//                          randomized uniformly over [-ycut, ycut], reweighted
//                          by Y²/(Y−DY)²  (see documents/extraction.tex)
//     mixed:               pid2 taken from the previous event; different events
//                          are uncorrelated by construction, no reweighting needed
//
// Both histograms are normalized by the total number of pid1 particles.
// The connected correlation is hSignal = hSignalAndBackground − hBackground.
//
// Usage:
//   ex_mixed_events.exe [--model hrg|blastwave] [--ini <file>]
//                       [--share <dir>] [--output <file>]
//                       [--events <N>] [--samples <N>] [--seed <N>]
//                       [--pid1 <pdg>] [--pid2 <pdg>]
//                       [--ycut <Y>] [--dy <DY>]
//                       [--background randomize|mixed]

#include "CLI11.hpp"
#include "Event.h"
#include "Integrator.h"
#include "Model.h"
#include "Model_BlastWave.h"
#include "Model_HRG.h"
#include "ParticleDB.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TRandom.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

static double invMass(const ParticleCoor& a, const ParticleCoor& b) {
  double e  = a.e  + b.e;
  double px = a.px + b.px;
  double py = a.py + b.py;
  double pz = a.pz + b.pz;
  double m2 = e*e - px*px - py*py - pz*pz;
  return (m2 > 0.) ? std::sqrt(m2) : 0.;
}

// Replace rapidity of p with a uniform random value in [-ycut, ycut],
// keeping pT and phi fixed. Returns the new rapidity.
static double randomizeRapidity(ParticleCoor& p, double ycut) {
  double pT   = std::sqrt(p.px*p.px + p.py*p.py);
  double mT   = std::sqrt(p.mass*p.mass + pT*pT);
  double yNew = gRandom->Uniform(-ycut, ycut);
  p.e  = mT * std::cosh(yNew);
  p.pz = mT * std::sinh(yNew);
  return yNew;
}

int main(int argc, char** argv) {
  std::string shareDir   = "";
  std::string outputFile = "mixed_events.root";
  std::string modelName  = "blastwave";
  std::string iniFile    = "";
  std::string bgMethod   = "randomize";
  int    nEvents  = 500;
  int    nSamples = 1000000;
  int    seed     = 0;
  int    pid1     = 211;    // pi+
  int    pid2     = -211;   // pi-
  double ycut     = 2.0;
  double DY       = 2.0;

  CLI::App app{"ex_mixed_events: signal and background invariant mass distributions"};
  app.add_option("--share",      shareDir,   "path to SHARE particle database");
  app.add_option("--output",     outputFile, "output ROOT file")->default_str("mixed_events.root");
  app.add_option("--model",      modelName,  "model: blastwave (default) or hrg")->default_str("blastwave");
  app.add_option("--ini",        iniFile,    "ini file (required for blastwave model)");
  app.add_option("--events",     nEvents,    "number of events to generate")->default_val(500);
  app.add_option("--samples",    nSamples,   "Monte Carlo integration samples")->default_val(1000000);
  app.add_option("--seed",       seed,       "random seed")->default_val(0);
  app.add_option("--pid1",       pid1,       "PDG code of particle type 1")->default_val(211);
  app.add_option("--pid2",       pid2,       "PDG code of particle type 2")->default_val(-211);
  app.add_option("--ycut",       ycut,       "rapidity acceptance |y| < ycut")->default_val(2.0);
  app.add_option("--dy",         DY,         "rapidity window for randomize background")->default_val(2.0);
  app.add_option("--background", bgMethod,   "background method: randomize (default) or mixed")->default_str("randomize");
  CLI11_PARSE(app, argc, argv);

  if (bgMethod != "randomize" && bgMethod != "mixed") {
    std::cerr << "Error: --background must be 'randomize' or 'mixed'\n";
    return 1;
  }

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

  TFile fout(outputFile.c_str(), "RECREATE", "Signal and background invariant mass");
  TH1::SetDefaultSumw2();

  TH1D* hSignalAndBackground = new TH1D("hSignalAndBackground",
      "Signal + background (same event);M (GeV);(1/N_{1}) dN_{pairs}/dM",
      100, 0., 2.);
  TH1D* hBackground = new TH1D("hBackground",
      "Background;M (GeV);(1/N_{1}) dN_{pairs}/dM",
      100, 0., 2.);

  // Diagnostic histograms — filled over full rapidity range (no ycut)
  TH1D* hRapidityPid1 = new TH1D("hRapidityPid1",
      "Rapidity of pid1;y;dN/dy",
      80, -8, 8);
  TH1D* hRapidityPid2 = new TH1D("hRapidityPid2",
      "Rapidity of pid2;y;dN/dy",
      80, -8, 8);
  TH1D* hSameRoot = new TH1D("hSameRoot",
      "Same rooteid, all #Deltay;M (GeV);(1/N_{1}) dN_{pairs}/dM",
      100, 0., 2.);
  TH1D* hSameRootClose = new TH1D("hSameRootClose",
      "Same rooteid, |#Deltay| < DY;M (GeV);(1/N_{1}) dN_{pairs}/dM",
      100, 0., 2.);
  TH1D* hSameRootFar = new TH1D("hSameRootFar",
      "Same rooteid, |#Deltay| #geq DY;M (GeV);(1/N_{1}) dN_{pairs}/dM",
      100, 0., 2.);

  const double Y      = 2. * ycut;
  const double weight = Y * Y / ((Y - DY) * (Y - DY));

  long nPid1Total = 0;
  TH1D* hSabEvt   = new TH1D("hSabEvt",   "", 100, 0., 2.);  hSabEvt->SetDirectory(nullptr);
  TH1D* hBkgEvt   = new TH1D("hBkgEvt",   "", 100, 0., 2.);  hBkgEvt->SetDirectory(nullptr);
  TProfile* hSigProfile = new TProfile("hSigProfile", "", 100, 0., 2.);
  hSigProfile->SetDirectory(nullptr);
  std::vector<ParticleCoor> v2prev;   // pid2 copies from previous event (mixed mode)

  // Seed v2prev with a burn-in event so that both hSignalAndBackground and
  // hBackground accumulate over exactly nEvents events.
  if (bgMethod == "mixed") {
    tEvent->Reset();
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();
    for (const auto& part : tEvent->GetParticleList()) {
      if (part.GetDecayed()) {
        continue;
      }
      if (part.pid == pid2 && std::abs(part.GetRapidityP()) < ycut) {
        v2prev.push_back(part);
      }
    }
  }

  for (int ievent = 0; ievent < nEvents; ievent++) {
    if (ievent % 100 == 0) {
      std::cout << "Event " << ievent << " / " << nEvents << "\n";
    }
    tEvent->Reset();
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();

    // Build pid1 and pid2 vectors. Rapidity histograms fill over full range;
    // the acceptance vectors apply the ycut.
    std::vector<const Particle*> v1, v2;
    for (const auto& part : tEvent->GetParticleList()) {
      if (part.GetDecayed()) {
        continue;
      }
      if (part.pid == pid1) {
        hRapidityPid1->Fill(part.GetRapidityP());
        ++nPid1Total;
        if (std::abs(part.GetRapidityP()) < ycut) {
          v1.push_back(&part);
        }
      }
      if (part.pid == pid2) {
        hRapidityPid2->Fill(part.GetRapidityP());
        if (std::abs(part.GetRapidityP()) < ycut) {
          v2.push_back(&part);
        }
      }
    }

    // Signal: same-event pid1 × pid2 pairs within acceptance.
    // Sameroot diagnostics fill for all |Δy|; the signal window requires |Δy| < DY.
    for (const auto* pa : v1) {
      for (const auto* pb : v2) {
        double dy = std::abs(pa->GetRapidityP() - pb->GetRapidityP());
        if (pa->rooteid == pb->rooteid) {
          double m = invMass(*pa, *pb);
          hSameRoot->Fill(m);
          if (dy < DY) {
            hSameRootClose->Fill(m);
          } else {
            hSameRootFar->Fill(m);
          }
        }
        if (dy >= DY) {
          continue;
        }
        double mSig = invMass(*pa, *pb);
        hSignalAndBackground->Fill(mSig);
        hSabEvt->Fill(mSig);
      }
    }

    // Background.
    if (bgMethod == "randomize") {
      for (const auto* pa : v1) {
        for (const auto* pb : v2) {
          if (std::abs(pa->GetRapidityP() - pb->GetRapidityP()) < DY) {
            continue;
          }
          ParticleCoor acoor = *pa;
          ParticleCoor bcoor = *pb;
          double yaRand = randomizeRapidity(acoor, ycut);
          double ybRand = randomizeRapidity(bcoor, ycut);
          if (std::abs(yaRand - ybRand) < DY) {
            double mBkg = invMass(acoor, bcoor);
            hBackground->Fill(mBkg, weight);
            hBkgEvt->Fill(mBkg, weight);
          }
        }
      }
    } else {
      // mixed: pair current v1 with previous event's v2, same |Δy| < DY window
      for (const auto* pa : v1) {
        for (const auto& pb : v2prev) {
          if (std::abs(pa->GetRapidityP() - pb.GetRapidityP()) >= DY) {
            continue;
          }
          double mBkg = invMass(*pa, pb);
          hBackground->Fill(mBkg);
          hBkgEvt->Fill(mBkg);
        }
      }
      // save current v2 for next event
      v2prev.clear();
      for (const auto* p : v2) {
        v2prev.push_back(*p);
      }
    }

    // Accumulate per-event signal s_i(M) = SAB_i - Bkg_i into TProfile.
    hSabEvt->Add(hBkgEvt, -1.);
    for (int b = 1; b <= hSabEvt->GetNbinsX(); b++) {
      hSigProfile->Fill(hSabEvt->GetBinCenter(b), hSabEvt->GetBinContent(b));
    }
    hSabEvt->Reset();
    hBkgEvt->Reset();
  }

  delete hSabEvt;
  delete hBkgEvt;

  // hSignal: derive from TProfile (mean and correct error on mean per bin).
  // TProfile normalizes by nEvents; scale by nEvents/nPid1Total to get per-pid1 units.
  TH1D* hSignal = (TH1D*) hSigProfile->ProjectionX("hSignal");
  hSignal->SetTitle("Signal (connected correlation);M (GeV);(1/N_{1}) dN_{pairs}/dM");
  delete hSigProfile;

  if (nPid1Total > 0) {
    hSignalAndBackground->Scale(1.0 / nPid1Total);
    hBackground->Scale(1.0 / nPid1Total);
    hSignal->Scale(double(nEvents) / nPid1Total);
    hSameRoot->Scale(1.0 / nPid1Total);
    hSameRootClose->Scale(1.0 / nPid1Total);
    hSameRootFar->Scale(1.0 / nPid1Total);
  }

  double sabErr = 0., bkgErr = 0., sigErr = 0.;
  double srErr = 0., srcErr = 0., srfErr = 0.;
  double sabYield = hSignalAndBackground->IntegralAndError(0, hSignalAndBackground->GetNbinsX() + 1, sabErr);
  double bkgYield = hBackground->IntegralAndError(0, hBackground->GetNbinsX() + 1, bkgErr);
  double sigYield = hSignal->IntegralAndError(0, hSignal->GetNbinsX() + 1, sigErr);
  double srYield  = hSameRoot->IntegralAndError(0, hSameRoot->GetNbinsX() + 1, srErr);
  double srcYield = hSameRootClose->IntegralAndError(0, hSameRootClose->GetNbinsX() + 1, srcErr);
  double srfYield = hSameRootFar->IntegralAndError(0, hSameRootFar->GetNbinsX() + 1, srfErr);
  std::cout << "Signal+background yield (pairs per pid1): " << sabYield << " +/- " << sabErr << "\n";
  std::cout << "Background yield        (pairs per pid1): " << bkgYield << " +/- " << bkgErr << "\n";
  std::cout << "Signal yield            (pairs per pid1): " << sigYield << " +/- " << sigErr << "\n";
  std::cout << "SameRoot yield          (pairs per pid1): " << srYield  << " +/- " << srErr  << "\n";
  std::cout << "SameRootClose yield     (pairs per pid1): " << srcYield << " +/- " << srcErr << "\n";
  std::cout << "SameRootFar yield       (pairs per pid1): " << srfYield << " +/- " << srfErr << "\n";



  fout.Write();
  return 0;
}
