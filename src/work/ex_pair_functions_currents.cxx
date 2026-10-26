// ex_pair_functions_currents.cxx — single-particle currents of a hadron resonance gas at rest
//
// See hrg_pair_density_spec.md, Sec. 9.2 and milestone M6. Builds S_c = A_c p + B_c eps u for every
// species in the ParticleDB (thermal + fed-down, recursion in decreasing mass over the THERMINATOR
// decay channels), writes them to a ROOT file and runs
//   T6.1  one step from thermal parents vs the closed-form two-body F_1, F_1 + F_u
//   T6.2  n_c = int d^3p (A_c + B_c) vs the counting of Sec. 9.5; grid-cut dependence
//   T6.3  positivity, S_c future timelike, finite at small k
//   T6.4  (with --mc) species densities and momentum spectra vs the output of
//         ex_pair_functions_mc.exe --spectra (use --norm mc)
// With --norm mc the thermal densities are the Integrator multiplicities that the MC sampled: run
// in the directory that holds the fmultiplicity cache used by the MC.
//
// Usage:
//   ex_pair_functions_currents.exe [--weak on|off] [--stat qs|boltzmann] [--norm quad|mc]
//                                  [--share <dir>] [--output <file>] [--samples <N>]
//                                  [--kcut <x>] [--cutcheck <x>] [--threads <N>] [--mc <file>]

#include "CLI11.hpp"
#include "Integrator.h"
#include "Model_HRG.h"
#include "PairFunctions.h"
#include "ParticleDB.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <memory>

static const char* PassFail(bool ok) { return ok ? "PASS" : "FAIL"; }

static double Seconds(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// Species with a momentum-spectrum test (same list as ex_pair_functions_mc.cxx --spectra).
static constexpr int kSpectraPDG[] = {211, 321, 2212, 113, 313, 223, 3122, 2224};

int main(int argc, char** argv) {
  std::string shareDir   = "";
  std::string outputFile = "";
  std::string weak       = "on";
  std::string stat       = "qs";
  std::string norm       = "quad";
  std::string mcFile     = "";
  int nSamples    = 1000000;
  double kcut     = 100.;
  double cutCheck = 40.;
  int nThreads    = 8;

  CLI::App app{"ex_pair_functions_currents: single-particle currents A_c, B_c of a resonance gas"};
  app.add_option("--share",    shareDir,   "path to SHARE particle database");
  app.add_option("--output",   outputFile, "output ROOT file (default: currents_<weak>.root)");
  app.add_option("--weak",     weak,       "weak decays: on (default) or off")
      ->default_str("on")
      ->check(CLI::IsMember({"on", "off"}));
  app.add_option("--stat",     stat,       "statistics: qs (default) or boltzmann")
      ->default_str("qs")
      ->check(CLI::IsMember({"qs", "boltzmann"}));
  app.add_option("--norm",     norm,       "thermal densities: quad (default) or mc (Integrator)")
      ->default_str("quad")
      ->check(CLI::IsMember({"quad", "mc"}));
  app.add_option("--samples",  nSamples,   "Integrator samples (--norm mc, must match the MC)")
      ->default_val(1000000);
  app.add_option("--kcut",     kcut,       "momentum grid cut in units of T (<= 100)")
      ->default_val(100.);
  app.add_option("--cutcheck", cutCheck,   "second grid cut for T6.2 (0: skip)")->default_val(40.);
  app.add_option("--threads",  nThreads,   "threads for the node loop")->default_val(8);
  app.add_option("--mc",       mcFile,     "merged ex_pair_functions_mc.exe --spectra output (T6.4)");
  CLI11_PARSE(app, argc, argv);

  const bool weakOn  = (weak == "on");
  const bool quantum = (stat == "qs");
  if (outputFile.empty()) outputFile = "currents_weak" + weak + ".root";
  {
    std::filesystem::path outPath(outputFile);
    if (outPath.has_parent_path()) std::filesystem::create_directories(outPath.parent_path());
  }
  auto tStart = std::chrono::steady_clock::now();

  auto model   = std::make_unique<Model_HRG>();
  auto tPartDB = std::make_unique<ParticleDB>(shareDir);
  int nStable  = SetWeakDecays(tPartDB.get(), weakOn);
  if (!weakOn && nStable != kNWeakTypes) {
    std::cerr << "Error: --weak off changed " << nStable << " particle types, expected "
              << kNWeakTypes << "\n";
    return 1;
  }
  ParticleDB* db   = tPartDB.get();
  const int nTypes = db->GetParticleTypeCount();
  const double T   = model->GetTemperature();
  const double V   = model->GetVolume();
  const std::string dbHash = ParticleDBHash(db);
  std::printf("HRG: T = %g GeV, V = %g GeV^-3, weak %s, stat %s, norm %s, kcut %g T, "
              "ParticleDB hash %s\n",
              T, V, weak.c_str(), stat.c_str(), norm.c_str(), kcut, dbHash.c_str());

  // ---- Thermal normalization ------------------------------------------------------------------
  std::vector<double> thScale(nTypes, 1.);
  if (norm == "mc") {
    Integrator integ(nSamples, model.get());
    integ.SetMultiplicities(db);
    for (int c = 0; c < nTypes; c++) {
      ParticleType* type = db->GetParticleType(c);
      double nq = NumberDensity(type->GetMass(), MakeParentStat(type, T, quantum));
      thScale[c] = type->GetMultiplicity() / V / nq;
    }
  }

  // ---- Currents ----------------------------------------------------------------------------------
  auto t0 = std::chrono::steady_clock::now();
  CurrentSet cs = BuildCurrents(db, T, quantum, thScale, kcut, 16, nThreads);
  const std::vector<SpeciesCurrent>& cur = cs.cur;
  int nTwo = 0, nThree = 0;
  for (const FeedSlot& fs : cs.slots) (fs.mX.size() == 1 ? nTwo : nThree)++;
  std::printf("Currents: %d species, %zu nodes per species, %zu feed slots (%d two-body, %d "
              "three-body), built in %.1f s\n",
              nTypes, cur[0].k.size(), cs.slots.size(), nTwo, nThree, Seconds(t0));

  // ---- T6.1: one step from thermal parents vs closed-form two-body F_1, F_1 + F_u ----------------
  {
    std::vector<SpeciesCurrent> thermal(nTypes);
    for (int c = 0; c < nTypes; c++) {
      SpeciesCurrent& s = thermal[c];
      s.m = cur[c].m; s.st = cur[c].st; s.thScale = cur[c].thScale;
      s.edges = cur[c].edges; s.nPerPanel = cur[c].nPerPanel; s.hasFeed = false;
    }
    double worstA = 0., worstAB = 0.;
    std::string nameA, nameAB;
    long nCmp = 0, nSkip = 0;
    for (const FeedSlot& fs : cs.slots) {
      if (fs.mX.size() != 1) continue;
      const SpeciesCurrent& par = thermal[fs.parent];
      const SpeciesCurrent& dau = cur[fs.daughter];
      if (par.thScale == 0.) continue;
      double ma = par.m, mb = dau.m, mX = fs.mX[0];
      double ps = PStar(ma, mb, mX), Es = (ma * ma + mb * mb - mX * mX) / (2. * ma);
      for (double k : dau.k) {
        double eps = std::hypot(mb, k);
        double Em  = ma * (Es * Es + k * k) / (eps * Es + k * ps);
        // The recursion cuts the parent at E_max; compare where the cut is below e^-50.
        if (Em > par.Emax() - 50. * T) {
          ++nSkip;
          continue;
        }
        double a, ab, ca, cab;
        TwoBodyFeed(k, mb, par, mX, fs.B, T, a, ab);
        TwoBodyClosedForm(k, mb, ma, mX, fs.B, par.st, par.thScale, ca, cab);
        if (cab <= 0.) continue;
        ++nCmp;
        // A changes sign: relative to max(|A|, A + B).
        double rA  = std::abs(a - ca) / std::max(std::abs(ca), cab);
        double rAB = std::abs(ab / cab - 1.);
        char buf[200];
        std::snprintf(buf, sizeof(buf), "%s -> %s, k = %.3e", cur[fs.parent].name.c_str(),
                      dau.name.c_str(), k);
        if (rA > worstA) { worstA = rA; nameA = buf; }
        if (rAB > worstAB) { worstAB = rAB; nameAB = buf; }
      }
    }
    std::printf("T6.1  %d two-body slots x nodes: %ld comparisons (%ld nodes beyond E_max - 50 T "
                "skipped)\n", nTwo, nCmp, nSkip);
    std::printf("T6.1  A_c   vs F_1      : worst rel. diff %.2e (%s)  %s\n", worstA, nameA.c_str(),
                PassFail(worstA < 1e-10));
    std::printf("T6.1  A_c+B_c vs F_1+F_u: worst rel. diff %.2e (%s)  %s\n", worstAB,
                nameAB.c_str(), PassFail(worstAB < 1e-10));
  }

  // ---- T6.2: grid density vs counting; grid cut --------------------------------------------------
  {
    double worst = 0.;
    std::string worstName;
    for (const SpeciesCurrent& s : cur) {
      double r = (s.nCount > 0.) ? std::abs(s.nGrid / s.nCount - 1.) : std::abs(s.nGrid);
      if (r > worst) { worst = r; worstName = s.name; }
    }
    std::printf("T6.2  n_grid vs counting, %d species: worst rel. diff %.2e (%s)  %s\n", nTypes,
                worst, worstName.c_str(), PassFail(worst < 1e-6));
    if (cutCheck > 0.) {
      auto t1 = std::chrono::steady_clock::now();
      CurrentSet cs2 = BuildCurrents(db, T, quantum, thScale, cutCheck, 16, nThreads);
      double worstCut = 0., worstCount = 0.;
      std::string nameCut, nameCount;
      for (int c = 0; c < nTypes; c++) {
        const SpeciesCurrent& s = cs2.cur[c];
        if (cur[c].nGrid <= 0.) continue;
        double r  = std::abs(s.nGrid / cur[c].nGrid - 1.);
        double r2 = std::abs(s.nGrid / s.nCount - 1.);
        if (r > worstCut) { worstCut = r; nameCut = s.name; }
        if (r2 > worstCount) { worstCount = r2; nameCount = s.name; }
      }
      std::printf("T6.2  grid cut %g T vs %g T: max |dn/n| %.2e (%s); at %g T n_grid vs counting "
                  "%.2e (%s)  %s  [%.1f s]\n",
                  cutCheck, kcut, worstCut, nameCut.c_str(), cutCheck, worstCount,
                  nameCount.c_str(), PassFail(worstCut < 1e-6), Seconds(t1));
    }
    std::printf("T6.2  %-12s %8s %14s %14s %14s %8s\n", "species", "pdg", "n_th", "n_count",
                "n_grid", "n/n_th");
    for (const char* nm : {"pi0139plu", "pi0139min", "pi0135zer", "Ka0492plu", "Ka0492zer",
                           "pr0938plu", "rho770zer", "Ka0892zer", "om0782zer", "eta0547zer",
                           "Lm1115zer", "Sg1192zer", "Dl1232plp", "gam000zer"}) {
      if (!db->ExistsParticleType(nm)) continue;
      const SpeciesCurrent& s = cur[db->GetParticleTypeIndex(nm)];
      std::printf("T6.2  %-12s %8d %14.7e %14.7e %14.7e %8.4f\n", nm, s.pdg, s.nth, s.nCount,
                  s.nGrid, s.nth > 0. ? s.nCount / s.nth : 0.);
    }
  }

  // ---- T6.3: positivity, timelike, finite at small k --------------------------------------------
  {
    long nNodes = 0, nNeg = 0, nSpace = 0, nNonFinite = 0, nNegA = 0;
    long gNodes = 0, gNeg = 0, gSpace = 0;
    double worstMargin = 1e300;
    std::string marginName;
    for (const SpeciesCurrent& s : cur) {
      bool photon = (s.thScale == 0. && s.name.find("gam000zer") != std::string::npos);
      for (size_t i = 0; i < s.k.size(); i++) {
        double eps = std::hypot(s.m, s.k[i]);
        double a = s.A[i] + s.ThermalA(eps), ab = s.AB[i] + s.ThermalA(eps);
        if (ab == 0. && a == 0.) continue; // beyond the support (e.g. underflow)
        bool bad = !std::isfinite(a) || !std::isfinite(ab);
        bool neg = !(ab > 0.), space = !(eps * ab > s.k[i] * std::abs(a));
        if (photon) {
          ++gNodes; gNeg += neg; gSpace += space;
          continue;
        }
        ++nNodes; nNonFinite += bad; nNeg += neg; nSpace += space; nNegA += (a < 0.);
        double margin = (eps * ab - s.k[i] * std::abs(a)) / (eps * ab);
        if (margin < worstMargin) { worstMargin = margin; marginName = s.name; }
      }
    }
    std::printf("T6.3  hadrons: %ld nodes: %ld non-finite, %ld with A+B<=0, %ld spacelike, %ld "
                "with A<0 (allowed); min (eps(A+B)-k|A|)/(eps(A+B)) = %.3e (%s)  %s\n",
                nNodes, nNonFinite, nNeg, nSpace, nNegA, worstMargin, marginName.c_str(),
                PassFail(nNonFinite == 0 && nNeg == 0 && nSpace == 0));
    std::printf("T6.3  photons: %ld nodes: %ld with A+B<=0, %ld spacelike (m_gamma = %.3g GeV in "
                "particles.data)\n",
                gNodes, gNeg, gSpace, cur[db->GetParticleTypeIndex("gam000zer")].m);
    int nSmall = 0, nSmallOK = 0;
    double worstJump = 0.;
    std::string jumpName;
    for (int b = 0; b < nTypes; b++) {
      if (cs.feeds[b].empty()) continue;
      const SpeciesCurrent& s = cur[b];
      double a6, ab6, a8, ab8;
      FeedAt(1e-6, s, cs.feeds[b], cur, T, a6, ab6);
      FeedAt(1e-8, s, cs.feeds[b], cur, T, a8, ab8);
      double th = s.ThermalA(s.m);
      ++nSmall;
      bool ok = std::isfinite(a6) && std::isfinite(ab6) && std::isfinite(a8) && std::isfinite(ab8) &&
                ab8 + th > 0. && ab6 + th > 0.;
      nSmallOK += ok;
      // Smoothness: k = 1e-8 vs 1e-6 (the currents are even in k, difference O(k^2)).
      double scale = std::max(std::abs(ab8), 1e-300);
      double jump  = std::max(std::abs(a6 - a8), std::abs(ab6 - ab8)) / scale;
      if (s.name.find("gam000zer") == std::string::npos && jump > worstJump) {
        worstJump = jump;
        jumpName  = s.name;
      }
    }
    std::printf("T6.3  k = 1e-6, 1e-8: %d/%d fed species finite with A+B > 0; max rel. change of "
                "the feed part between the two (hadrons) %.2e (%s)  %s\n",
                nSmallOK, nSmall, worstJump, jumpName.c_str(),
                PassFail(nSmallOK == nSmall && worstJump < 1e-4));
  }

  // ---- Output -------------------------------------------------------------------------------------
  TFile fout(outputFile.c_str(), "RECREATE", "Single-particle currents of a resonance gas at rest");
  {
    TTree* tree = new TTree("currents", "A_c, A_c+B_c (thermal + feed-down) on the momentum grid");
    char name[64];
    int pdg;
    double m, nth, nCount, nGrid;
    std::vector<double> k, w, A, AB;
    tree->Branch("name", name, "name/C");
    tree->Branch("pdg", &pdg, "pdg/I");
    tree->Branch("mass", &m, "mass/D");
    tree->Branch("n_th", &nth, "n_th/D");
    tree->Branch("n_count", &nCount, "n_count/D");
    tree->Branch("n_grid", &nGrid, "n_grid/D");
    tree->Branch("k", &k);
    tree->Branch("w", &w);
    tree->Branch("A", &A);
    tree->Branch("AB", &AB);
    for (const SpeciesCurrent& s : cur) {
      std::snprintf(name, sizeof(name), "%s", s.name.c_str());
      pdg = s.pdg; m = s.m; nth = s.nth; nCount = s.nCount; nGrid = s.nGrid;
      k = s.k; w = s.w; A = s.A; AB = s.AB;
      for (size_t i = 0; i < k.size(); i++) {
        double th = s.ThermalA(std::hypot(s.m, k[i]));
        A[i] += th;
        AB[i] += th;
      }
      tree->Fill();
    }
    tree->Write();
    delete tree;
  }

  // ---- T6.4: MC species densities and momentum spectra ------------------------------------------
  if (!mcFile.empty()) {
    if (norm != "mc")
      std::printf("T6.4  warning: --norm %s; the MC samples the Integrator multiplicities (use "
                  "--norm mc)\n", norm.c_str());
    TFile fmc(mcFile.c_str());
    TH1D* hEv = dynamic_cast<TH1D*>(fmc.Get("hEvents"));
    TH1D* hN  = dynamic_cast<TH1D*>(fmc.Get("hSpeciesN"));
    TH1D* hN2 = dynamic_cast<TH1D*>(fmc.Get("hSpeciesN2"));
    TNamed* mcWeak = dynamic_cast<TNamed*>(fmc.Get("weak"));
    if (!hEv || !hN || !hN2) {
      std::cerr << "Error: " << mcFile << " has no --spectra histograms\n";
      return 1;
    }
    if (mcWeak && weak != mcWeak->GetTitle()) {
      std::cerr << "Error: MC file has weak " << mcWeak->GetTitle() << ", this run " << weak << "\n";
      return 1;
    }
    const double nEv = hEv->GetBinContent(1);
    double sumP = 0., sumP2 = 0., worst = 0.;
    int nP = 0, nOver4 = 0;
    std::string worstName;
    fout.cd();
    TH1D* hPull = new TH1D("hPullSpecies", "species densities;(n_{MC}-n_{c})/#sigma;species", 50,
                           -5., 5.);
    for (int c = 0; c < nTypes; c++) {
      const SpeciesCurrent& s = cur[c];
      double S = hN->GetBinContent(c + 1), S2 = hN2->GetBinContent(c + 1);
      if (s.nCount <= 0. && S == 0.) continue;
      double mean = S / nEv;
      double var  = std::max(S2 / nEv - mean * mean, 0.);
      double sig  = (S > 0.) ? std::sqrt(var / nEv) : std::sqrt(s.nCount * V / nEv);
      if (sig <= 0.) continue;
      double pull = (mean - s.nCount * V) / sig;
      hPull->Fill(pull);
      ++nP;
      sumP += pull;
      sumP2 += pull * pull;
      nOver4 += (std::abs(pull) > 4.);
      if (std::abs(pull) > std::abs(worst)) {
        worst = pull;
        char buf[200];
        std::snprintf(buf, sizeof(buf), "%s: MC %.6e, n_c V %.6e per event", s.name.c_str(), mean,
                      s.nCount * V);
        worstName = buf;
      }
    }
    double pMean = sumP / nP, pRms = std::sqrt(sumP2 / nP);
    std::printf("T6.4  species densities, %.0f events, %d species: pull mean %.3f, rms %.3f, "
                "%d with |pull|>4, max |pull| %.2f (%s)  %s\n",
                nEv, nP, pMean, pRms, nOver4, std::abs(worst), worstName.c_str(),
                PassFail(std::abs(worst) < 5. && pRms < 1.3));
    for (const char* nm : {"pi0139plu", "Ka0492plu", "pr0938plu", "rho770zer", "Ka0892zer",
                           "om0782zer", "Lm1115zer", "Dl1232plp", "Ka0492zer", "gam000zer"}) {
      int c = db->GetParticleTypeIndex(nm);
      double S = hN->GetBinContent(c + 1), S2 = hN2->GetBinContent(c + 1);
      double mean = S / nEv, sig = std::sqrt(std::max(S2 / nEv - mean * mean, 0.) / nEv);
      std::printf("T6.4  %-10s MC %12.5f +- %8.5f   n_c V %12.5f  (n_th V %12.5f)  pull %6.2f\n", nm,
                  mean, sig, cur[c].nCount * V, cur[c].nth * V,
                  sig > 0. ? (mean - cur[c].nCount * V) / sig : 0.);
    }

    bool allOK = true;
    const GaussLegendre& gl = GLRule(8);
    for (int pdgS : kSpectraPDG) {
      std::string tag = std::to_string(pdgS);
      TH1D* hmc = dynamic_cast<TH1D*>(fmc.Get(("hSpec_" + tag).c_str()));
      int c = -1;
      for (int i = 0; i < nTypes; i++)
        if (cur[i].pdg == pdgS) c = i;
      if (!hmc || c < 0) {
        std::printf("T6.4  dN/dp %s: missing\n", tag.c_str());
        allOK = false;
        continue;
      }
      const SpeciesCurrent& s = cur[c];
      fout.cd();
      TH1D* hmcOut = static_cast<TH1D*>(hmc->Clone(("hSpecMC_" + tag).c_str()));
      TH1D* hpred  = static_cast<TH1D*>(hmc->Clone(("hSpecPred_" + tag).c_str()));
      hpred->Reset();
      double chi2 = 0.;
      int ndf = 0;
      for (int ib = 1; ib <= hmc->GetNbinsX(); ib++) {
        double p1 = hmc->GetXaxis()->GetBinLowEdge(ib), p2 = hmc->GetXaxis()->GetBinUpEdge(ib);
        double hp = 0.5 * (p2 - p1), cp = 0.5 * (p2 + p1), sum = 0.;
        for (int i = 0; i < 8; i++) {
          double p = cp + hp * gl.x[i];
          double a, ab;
          s.AtEnergy(std::hypot(s.m, p), a, ab);
          sum += gl.w[i] * hp * 4. * M_PI * p * p * ab;
        }
        double expct = nEv * V * sum;
        hpred->SetBinContent(ib, expct);
        hpred->SetBinError(ib, 0.);
        hmcOut->SetBinError(ib, std::sqrt(hmc->GetBinContent(ib)));
        if (expct < 5.) continue;
        double d = hmc->GetBinContent(ib) - expct;
        chi2 += d * d / expct;
        ++ndf;
      }
      double prob = TMath::Prob(chi2, ndf);
      allOK = allOK && prob > 1e-3;
      std::printf("T6.4  dN/dp %-10s (pdg %5d): MC %.0f, predicted %.1f, chi2/ndf = %.1f/%d, p = "
                  "%.3g  %s\n",
                  s.name.c_str(), pdgS, hmc->Integral(), hpred->Integral(), chi2, ndf, prob,
                  PassFail(prob > 1e-3));
      hmcOut->Write();
      hpred->Write();
    }
    std::printf("T6.4  momentum spectra, all species p > 1e-3  %s\n", PassFail(allOK));
    fout.cd();
    hPull->Write();
    TParameter<double>("nEvents", nEv).Write();
  }

  fout.cd();
  TParameter<double>("T", T).Write();
  TParameter<double>("V", V).Write();
  TParameter<double>("kcut", kcut).Write();
  TNamed("weak", weak.c_str()).Write();
  TNamed("stat", stat.c_str()).Write();
  TNamed("norm", norm.c_str()).Write();
  TNamed("dbhash", dbHash.c_str()).Write();
  fout.Close();
  std::printf("Done in %.1f s, output %s\n", Seconds(tStart), outputFile.c_str());
  return 0;
}
