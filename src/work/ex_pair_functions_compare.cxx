// ex_pair_functions_compare.cxx — MC at levels L0, L1, L2 versus the pair-function table
//
// See README_pair_functions.md, milestones M2 and M3. Reads the (hadd-merged) output of
// ex_pair_functions_mc.exe, rebuilds the table contributions from the same ParticleDB and
// weak setting, integrates the predictions over the MC bins and runs
//   T2.1  rho(M), L0 vs table, chi2/ndf per ordered pair (pole bins listed)
//   T2.2  dP/dkappa at fixed M (pole bins and continuum windows), shape chi2
//   T2.3  Legendre moments a1, a2 at L0 consistent with 0 (pulls)
//   M3    rho(M) at the cumulative levels L1 = L0 + L1x and L2 = L1 + L2x, with the table split by
//         parent; dP/dkappa at fixed M for L0, L1, L2 and the table, divided by the Boltzmann
//         reference; a1(kappa), a2(kappa) and raw cos theta* at fixed M; the L2 table of a1, a2
//         and C_Q/(C_K + C_u) in (M, eps_K) cells
//   T3.1  a1, a2 at L1 consistent with 0 (pulls)
// With --level L1 (milestone M7) the table is the L1 table (parent currents, spec Sec. 9.3) and the
// MC is the cumulative L1 level; only the rho(M) and dP/dkappa tests run, labelled T7.4rho and
// T7.4kap (same criteria as T2.1, T2.2). With --level L2 --table <file> (milestone M8) the chain
// part of the L2 table (ex_pair_functions_table.exe --level L2) is added to the L1 prediction, its
// kernel variance is added to the chi2, the MC is the cumulative L2 level, and the tests are labelled
// T8.3rho and T8.3kap.
// Predictions are normalized to the Integrator multiplicities that the MC sampled (--norm mc,
// default) or to the quadrature densities of the table (--norm quad). Run in the directory that
// holds the fmultiplicity cache used by the MC. The plots are made by ex_plot_pair_functions.C.
//
// Usage:
//   ex_pair_functions_compare.exe --mc <file> [--share <dir>] [--output <file>]
//                                 [--norm mc|quad] [--samples <N>] [--minexp <x>]
//                                 [--level L0|L1|L2] [--table <L2 table file>]

#include "CLI11.hpp"
#include "Integrator.h"
#include "Model_HRG.h"
#include "PairFunctions.h"
#include "ParticleDB.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <memory>

static const char* PassFail(bool ok) { return ok ? "PASS" : "FAIL"; }

// Cumulative ancestry levels L0, L1 = L0 + L1x, L2 = L1 + L2x, built from the exclusive MC levels.
static constexpr int kNLevels = 3;
static const char* const kLevel[kNLevels]   = {"L0", "L1", "L2"};
static const char* const kMCLevel[kNLevels] = {"L0", "L1x", "L2x"};

// chi2 of counts n against expectations mu (variance varFactor * mu). Bins with mu < minExp are
// pooled into one extra bin, used if its total expectation reaches minExp.
// For like-sign pairs every pair of particles is filled twice, (i,j) and (j,i), with the same M
// and kappa (and the same P2, opposite P1): entries come in identical doublets, so the variance
// of a bin count is 2 mu and a profile's error is too small by sqrt(2).
struct Chi2 {
  double chi2 = 0.;
  int ndf     = 0;
};

// extraVar (optional): variance of the expectation itself (L2 kernel statistics), added in
// quadrature.
static Chi2 PoissonChi2(const std::vector<double>& n, const std::vector<double>& mu, double minExp,
                        int constraints = 0, double varFactor = 1.,
                        const std::vector<double>* extraVar = nullptr) {
  Chi2 r;
  double nPool = 0., muPool = 0., varPool = 0.;
  for (size_t i = 0; i < n.size(); i++) {
    double ev = extraVar ? (*extraVar)[i] : 0.;
    if (mu[i] >= minExp) {
      r.chi2 += (n[i] - mu[i]) * (n[i] - mu[i]) / (varFactor * mu[i] + ev);
      ++r.ndf;
    } else {
      nPool += n[i];
      muPool += mu[i];
      varPool += ev;
    }
  }
  if (muPool >= minExp) {
    r.chi2 += (nPool - muPool) * (nPool - muPool) / (varFactor * muPool + varPool);
    ++r.ndf;
  }
  r.ndf -= constraints;
  return r;
}

// Density (per unit volume) of one contribution in the bins 1..n of the M axis; index n holds
// M above the last edge.
static std::vector<double> PredictM(const Contribution& co, const TAxis* ax) {
  const int nM = ax->GetNbins();
  std::vector<double> mu(nM + 1, 0.);
  if (co.third < 0) {
    int b = ax->FindFixBin(co.ma);
    mu[std::min(b, nM + 1) - 1] = co.count * co.B * co.na;
  } else {
    for (int b = 1; b <= nM; b++)
      mu[b - 1] = ContDensityM(co, ax->GetBinLowEdge(b), ax->GetBinUpEdge(b));
    mu[nM] = ContDensityM(co, ax->GetXmax(), 1e9);
  }
  return mu;
}

// Density (per unit volume) of one contribution with M in [M1, M2), in the bins of the kappa axis.
static std::vector<double> PredictKappa(const Contribution& co, double M1, double M2,
                                        const TAxis* ay) {
  std::vector<double> mu(ay->GetNbins(), 0.);
  bool isPole = (co.third < 0);
  if (isPole && (co.ma < M1 || co.ma >= M2)) return mu;
  for (int k = 1; k <= ay->GetNbins(); k++) {
    double k1 = ay->GetBinLowEdge(k), k2 = ay->GetBinUpEdge(k);
    mu[k - 1] = isPole ? PoleDensityKappa(co, k1, k2) : ContDensityMK(co, M1, M2, k1, k2);
  }
  return mu;
}

// Sum of the exclusive MC levels up to `level` (cumulative L0, L1, L2), detached from any file.
template <class H>
static std::unique_ptr<H> Cumulative(TFile& f, const char* prefix, const std::string& tag,
                                     int level) {
  std::unique_ptr<H> h(static_cast<H*>(f.Get<H>(Form("%s_L0_%s", prefix, tag.c_str()))
                                           ->Clone(Form("%s_%s_%s_cum", prefix, kLevel[level],
                                                        tag.c_str()))));
  h->SetDirectory(nullptr);
  for (int l = 1; l <= level; l++)
    h->Add(f.Get<H>(Form("%s_%s_%s", prefix, kMCLevel[l], tag.c_str())));
  return h;
}

// Pulls a_l / sigma(a_l) over the (M, kappa) bins of a profile with at least minEntries entries,
// split into all bins, bins containing a pole and continuum-only bins. Prints one line per
// non-empty split, headed by `head`; returns the verdict for "all" (true if no bins). With
// info == true the lines carry no verdict.
static bool PullTest(const std::string& head, TProfile2D* prof, int l, bool likeSign,
                     const std::vector<double>& poleMasses, double minEntries, bool info = false) {
  std::vector<double> pulls[3]; // all, pole bins, continuum-only bins
  for (int ix = 1; ix <= prof->GetNbinsX(); ix++) {
    double lo = prof->GetXaxis()->GetBinLowEdge(ix), hi = prof->GetXaxis()->GetBinUpEdge(ix);
    bool hasPole = false;
    for (double m : poleMasses) hasPole = hasPole || (m >= lo && m < hi);
    for (int iy = 1; iy <= prof->GetNbinsY(); iy++) {
      int b = prof->GetBin(ix, iy);
      if (prof->GetBinEntries(b) < minEntries || prof->GetBinError(b) <= 0.) continue;
      double pull = prof->GetBinContent(b) / prof->GetBinError(b) / (likeSign ? std::sqrt(2.) : 1.);
      pulls[0].push_back(pull);
      pulls[hasPole ? 1 : 2].push_back(pull);
    }
  }
  const char* names[3] = {"all", "pole", "continuum"};
  bool okAll = true;
  for (int s = 0; s < 3; s++) {
    const std::vector<double>& v = pulls[s];
    if (v.empty()) {
      if (s == 0)
        std::printf("%s a%d  %-10s no bins with >= %.0f entries\n", head.c_str(), l, names[s],
                    minEntries);
      continue;
    }
    double sum = 0., sum2 = 0.;
    for (double x : v) {
      sum += x;
      sum2 += x * x;
    }
    double nb = v.size(), mean = sum / nb, width = std::sqrt(sum2 / nb - mean * mean);
    double tolM = 3. / std::sqrt(nb), tolW = 3. / std::sqrt(2. * nb);
    bool ok = std::abs(mean) < tolM && std::abs(width - 1.) < tolW;
    if (s == 0) okAll = ok;
    std::printf("%s a%d  %-10s %6.0f %+8.3f %8.3f  (%.3f, %.3f)  %s\n", head.c_str(), l, names[s],
                nb, mean, width, tolM, tolW, info ? "(info)" : PassFail(ok));
  }
  return okAll;
}

// Profile-to-histogram: mean and error per kappa bin of a TProfile; errors scaled by errScale.
static TH1D* ProfileToHist(TProfile* prof, const std::string& name, double errScale) {
  TH1D* h = prof->ProjectionX(name.c_str());
  for (int k = 1; k <= h->GetNbinsX(); k++) h->SetBinError(k, h->GetBinError(k) * errScale);
  return h;
}

int main(int argc, char** argv) {
  std::string mcFile     = "";
  std::string shareDir   = "";
  std::string outputFile = "";
  std::string norm       = "mc";
  std::string level      = "L0";
  std::string tableFile  = "";
  int nSamples  = 1000000;
  double minExp = 5.;

  CLI::App app{"ex_pair_functions_compare: MC (L0, L1, L2) versus the pair-function table"};
  app.add_option("--mc",      mcFile,     "merged output of ex_pair_functions_mc.exe")->required();
  app.add_option("--share",   shareDir,   "path to SHARE particle database");
  app.add_option("--output",  outputFile, "output ROOT file (default: compare_<weak>.root)");
  app.add_option("--norm",    norm,       "parent densities: mc (Integrator, default) or quad")
      ->default_str("mc")
      ->check(CLI::IsMember({"mc", "quad"}));
  app.add_option("--samples", nSamples,   "Integrator samples (must match the MC)")
      ->default_val(1000000);
  app.add_option("--minexp",  minExp,     "minimum expected count per chi2 bin")->default_val(5.);
  app.add_option("--level",   level,      "table and MC level: L0 (default), L1 or L2")
      ->default_str("L0")
      ->check(CLI::IsMember({"L0", "L1", "L2"}));
  app.add_option("--table",   tableFile,  "L2: output of ex_pair_functions_table.exe --level L2");
  CLI11_PARSE(app, argc, argv);

  TFile fin(mcFile.c_str());
  if (fin.IsZombie()) return 1;
  std::string weak   = fin.Get<TNamed>("weak")->GetTitle();
  std::string mcHash = fin.Get<TNamed>("dbhash")->GetTitle();
  const double nEv   = fin.Get<TH1D>("hEvents")->GetBinContent(1);
  const bool l2 = (level == "L2");
  const bool l1 = (level == "L1") || l2; // L2 prediction = L1 table + chain part
  const int mcLevel = l2 ? 2 : (l1 ? 1 : 0);
  if (outputFile.empty())
    outputFile = (l1 ? "compare_" + level + "_weak" : std::string("compare_weak")) + weak + ".root";
  const char* tRho = l2 ? "T8.3rho" : (l1 ? "T7.4rho" : "T2.1");
  const char* tKap = l2 ? "T8.3kap" : (l1 ? "T7.4kap" : "T2.2");
  const char* mcLev = kLevel[mcLevel];

  std::unique_ptr<TFile> ftab;
  if (l2) {
    ftab.reset(TFile::Open(tableFile.c_str()));
    if (!ftab || ftab->IsZombie() || !ftab->Get<TNamed>("level") ||
        std::string(ftab->Get<TNamed>("level")->GetTitle()) != "L2") {
      std::cerr << "Error: --level L2 needs --table <output of ex_pair_functions_table.exe "
                   "--level L2>\n";
      return 1;
    }
    for (const char* key : {"weak", "norm"}) {
      std::string want = std::string(key) == "weak" ? weak : norm;
      if (want != ftab->Get<TNamed>(key)->GetTitle()) {
        std::cerr << "Error: table file has " << key << " " << ftab->Get<TNamed>(key)->GetTitle()
                  << ", this run " << want << "\n";
        return 1;
      }
    }
  }

  auto model   = std::make_unique<Model_HRG>();
  auto tPartDB = std::make_unique<ParticleDB>(shareDir);
  ParticleDB* db = tPartDB.get();
  SetWeakDecays(db, weak == "on");
  if (ParticleDBHash(db) != mcHash) {
    std::cerr << "Error: ParticleDB hash " << ParticleDBHash(db) << " differs from the MC file ("
              << mcHash << ")\n";
    return 1;
  }
  if (l2 && mcHash != ftab->Get<TNamed>("dbhash")->GetTitle()) {
    std::cerr << "Error: ParticleDB hash of the table file differs from the MC file\n";
    return 1;
  }
  Integrator integ(nSamples, model.get());
  if (!integ.SetMultiplicities(db))
    std::cerr << "WARNING: multiplicity cache recomputed, not the one used by the MC\n";

  const double T = model->GetTemperature();
  const double V = model->GetVolume();
  ContributionSet cs = BuildContributions(db, T, true);
  // Per-contribution normalization: Integrator density / quadrature density.
  std::vector<double> scale(cs.contribs.size(), 1.);
  // L1: parent currents at the same normalization; the scale is then inside n_c and the currents.
  CurrentSet curs;
  if (l1) {
    std::vector<double> thScale(db->GetParticleTypeCount(), 1.);
    if (norm == "mc")
      for (int c = 0; c < db->GetParticleTypeCount(); c++)
        thScale[c] = db->GetParticleType(c)->GetMultiplicity() / V / cs.nquad[c];
    curs = BuildCurrents(db, T, true, thScale);
    ToL1(cs, db, curs);
  }
  if (norm == "mc" && !l1)
    for (size_t i = 0; i < cs.contribs.size(); i++) {
      const Contribution& co = cs.contribs[i];
      scale[i] = db->GetParticleType(co.parent)->GetMultiplicity() / V / co.na;
    }
  const double NV = nEv * V;
  std::printf("MC file %s: %.0f events, weak %s, norm %s, level %s, T = %g GeV, V = %g GeV^-3\n",
              mcFile.c_str(), nEv, weak.c_str(), norm.c_str(), level.c_str(), T, V);

  TFile fout(outputFile.c_str(), "RECREATE", "MC L0, L1, L2 versus table");
  fout.cd();
  TNamed("weak", weak.c_str()).Write();
  bool allT21 = true, allT22 = true, allT23 = true, allT31 = true;

  // ---- T2.1: rho(M) ---------------------------------------------------------------------------
  std::printf("%s  %-9s %12s %12s %8s %5s %8s %9s\n", tRho, "pair",
              (std::string("N_") + mcLev).c_str(), "N_pred", "chi2", "ndf", "chi2/ndf", "p");
  std::vector<std::string> poleLines;
  std::vector<std::vector<double>> muMPair(kNPairs); // table counts per M bin, for M3
  for (int p = 0; p < kNPairs; p++) {
    std::string tag = kPairTag[p];
    auto hCum = Cumulative<TH1D>(fin, "hM", tag, mcLevel);
    TH1D* h   = hCum.get();
    const int nM = h->GetNbinsX();
    // bins 1..nM and the overflow (M > 2 GeV)
    std::vector<double> n(nM + 1), mu(nM + 1, 0.);
    for (int b = 1; b <= nM + 1; b++) n[b - 1] = h->GetBinContent(b);
    // Per-parent table decomposition (M3), indexed by ParticleDB index.
    std::map<int, std::vector<double>> muParent;
    for (size_t i = 0; i < cs.contribs.size(); i++) {
      const Contribution& co = cs.contribs[i];
      if (co.pair != p) continue;
      std::vector<double> m = PredictM(co, h->GetXaxis());
      std::vector<double>& mp = muParent[co.parent];
      mp.resize(nM + 1, 0.);
      for (int b = 0; b <= nM; b++) {
        mu[b] += NV * scale[i] * m[b];
        mp[b] += scale[i] * m[b];
      }
    }
    // L2: chain part (M bins and the overflow) and its kernel variance
    std::vector<double> extraVar(nM + 1, 0.);
    if (l2) {
      TH1D* hc  = ftab->Get<TH1D>(("rho_chain_" + tag).c_str());
      TH1D* hcv = ftab->Get<TH1D>(("rho_chain_var_" + tag).c_str());
      TH1D* hChain = static_cast<TH1D*>(h->Clone(("rho_pred_chain_" + tag).c_str()));
      hChain->Reset();
      for (int b = 0; b <= nM; b++) {
        mu[b] += NV * hc->GetBinContent(b + 1);
        extraVar[b] = NV * NV * hcv->GetBinContent(b + 1);
        hChain->SetBinContent(b + 1, hc->GetBinContent(b + 1));
      }
      fout.cd();
      hChain->Write();
    }
    muMPair[p] = mu;
    double nTot = 0., muTot = 0.;
    for (int b = 0; b <= nM; b++) {
      nTot += n[b];
      muTot += mu[b];
    }
    TH1D* hPred = static_cast<TH1D*>(h->Clone(("rho_pred_" + tag).c_str()));
    TH1D* hMC   = static_cast<TH1D*>(h->Clone(("rho_" + std::string(mcLev) + "_" + tag).c_str()));
    hPred->Reset();
    for (int b = 1; b <= nM + 1; b++) hPred->SetBinContent(b, mu[b - 1] / NV);
    hMC->Sumw2();
    hMC->Scale(1. / NV);
    fout.cd();
    hPred->Write();
    hMC->Write();

    // M3: per-parent decomposition, the 7 largest parents and the rest as "other".
    std::vector<std::pair<double, int>> order;
    for (auto& kv : muParent) {
      double s = 0.;
      for (double x : kv.second) s += x;
      order.push_back({s, kv.first});
    }
    std::sort(order.rbegin(), order.rend());
    std::string list;
    TH1D* hOther = nullptr;
    for (size_t k = 0; k < order.size(); k++) {
      const std::vector<double>& mp = muParent[order[k].second];
      TH1D* hp;
      if (k < 7) {
        std::string pname = db->GetParticleType(order[k].second)->GetName();
        hp = static_cast<TH1D*>(h->Clone(("rho_pred_" + tag + "_" + pname).c_str()));
        hp->Reset();
        list += (list.empty() ? "" : ",") + pname;
      } else {
        if (!hOther) {
          hOther = static_cast<TH1D*>(h->Clone(("rho_pred_" + tag + "_other").c_str()));
          hOther->Reset();
          list += (list.empty() ? "" : ",") + std::string("other");
        }
        hp = hOther;
      }
      for (int b = 1; b <= nM + 1; b++) hp->AddBinContent(b, mp[b - 1]);
      if (k < 7) hp->Write();
    }
    if (hOther) hOther->Write();
    TNamed(("decomp_" + tag).c_str(), list.c_str()).Write();

    if (muTot == 0.) {
      bool ok = (nTot == 0.);
      allT21 = allT21 && ok;
      std::printf("%s  %-9s %12.0f %12.1f   (table empty: %s must be empty)  %s\n", tRho,
                  tag.c_str(), nTot, muTot, mcLev, PassFail(ok));
      continue;
    }
    const double varFactor = (kPairPDG[p][0] == kPairPDG[p][1]) ? 2. : 1.;
    Chi2 c = PoissonChi2(n, mu, minExp, 0, varFactor, &extraVar);
    double prob = TMath::Prob(c.chi2, c.ndf);
    bool ok = prob > 1e-3;
    allT21 = allT21 && ok;
    std::printf("%s  %-9s %12.0f %12.1f %8.1f %5d %8.3f %9.3g  %s\n", tRho, tag.c_str(), nTot, muTot,
                c.chi2, c.ndf, c.chi2 / c.ndf, prob, PassFail(ok));
    for (const Contribution& co : cs.contribs) {
      if (co.pair != p || co.third >= 0 || co.ma >= 2.0) continue;
      int b = h->FindBin(co.ma);
      double dEdge = std::min(co.ma - h->GetXaxis()->GetBinLowEdge(b),
                              h->GetXaxis()->GetBinUpEdge(b) - co.ma);
      char buf[200];
      std::snprintf(buf, sizeof(buf),
                    "%s  pole %-9s %-10s M %.5f bin [%.4f,%.4f] (edge %.1e): N %9.0f  pred %11.1f"
                    "  pull %+6.2f",
                    tRho, tag.c_str(), db->GetParticleType(co.parent)->GetName(), co.ma,
                    h->GetXaxis()->GetBinLowEdge(b), h->GetXaxis()->GetBinUpEdge(b), dEdge,
                    n[b - 1], mu[b - 1], (n[b - 1] - mu[b - 1]) / std::sqrt(varFactor * mu[b - 1]));
      poleLines.push_back(buf);
    }
  }
  for (const std::string& s : poleLines) std::printf("%s\n", s.c_str());
  std::printf("%s  all ordered pairs  %s\n", tRho, PassFail(allT21));

  // ---- M3: rho(M) at the cumulative levels ----------------------------------------------------
  // Like-sign errors include the doublet factor; unlike-sign errors are Poisson (pairs sharing a
  // particle at L1 and L2 make them a lower bound).
  if (!l1) {
    std::printf("M3rho %-9s %-3s %12s %12s %10s %10s\n", "pair", "lev", "N", "N_table", "N/N_L0",
                "N/N_table");
    for (int p = 0; p < kNPairs; p++) {
      std::string tag = kPairTag[p];
      const double errScale = (kPairPDG[p][0] == kPairPDG[p][1]) ? std::sqrt(2.) : 1.;
      double muTot = 0.;
      for (double x : muMPair[p]) muTot += x;
      double nL0 = 0.;
      for (int l = 0; l < kNLevels; l++) {
        auto hc = Cumulative<TH1D>(fin, "hM", tag, l);
        double nTot = hc->Integral(0, hc->GetNbinsX() + 1);
        if (l == 0) nL0 = nTot;
        std::printf("M3rho %-9s %-3s %12.0f %12.1f %10.4f %10.4f\n", tag.c_str(), kLevel[l], nTot,
                    muTot, nL0 > 0. ? nTot / nL0 : 0., muTot > 0. ? nTot / muTot : 0.);
        if (l == 0) continue; // rho_L0_ written by T2.1
        TH1D* hMC = static_cast<TH1D*>(hc->Clone(Form("rho_%s_%s", kLevel[l], tag.c_str())));
        hMC->Sumw2();
        for (int b = 0; b <= hMC->GetNbinsX() + 1; b++)
          hMC->SetBinError(b, hMC->GetBinError(b) * errScale);
        hMC->Scale(1. / NV);
        fout.cd();
        hMC->Write();
      }
    }
  }

  // ---- T2.2 and M3: dP/dkappa at fixed M --------------------------------------------------------
  struct TestMass {
    int pdg1;                // first species (|pdg|) of the pairs it applies to
    const char* poleParent;  // pole parent name, or nullptr for a continuum window
    double M;                // continuum window center
  };
  const TestMass tests[] = {{211, "Ka0492zer", 0.},  {211, "rho770zer", 0.}, {211, "om0782zer", 0.},
                            {211, "f00980zer", 0.},  {211, nullptr, 0.30},     {211, nullptr, 0.45},
                            {211, nullptr, 0.60},    {321, "ph1020zer", 0.},   {321, nullptr, 1.10},
                            {321, nullptr, 1.30}};
  // kappa bands of the raw cos theta* histograms
  const double cthBand[2][2] = {{0.10, 0.30}, {0.50, 0.80}};
  const double minM3 = 500.; // minimum number of pairs (or predicted pairs) for an M3 curve
  std::printf("%s  %-9s %-14s %-17s %10s %10s %7s %7s %4s %8s %9s\n", tKap, "pair", "mass",
              "M window", (std::string("N_") + mcLev).c_str(), "N_pred", "N/pred", "chi2", "ndf", "chi2/ndf", "p");
  std::vector<std::string> m3Lines;
  for (int p = 0; p < kNPairs; p++) {
    std::string tag = kPairTag[p];
    const bool likeSign    = (kPairPDG[p][0] == kPairPDG[p][1]);
    const bool unlikeSign  = (kPairPDG[p][0] == -kPairPDG[p][1]);
    const double errScale  = likeSign ? std::sqrt(2.) : 1.;
    auto h3Cum = Cumulative<TH3D>(fin, "hMKC", tag, mcLevel);
    TH3D* h3   = h3Cum.get();
    const TAxis* ax = h3->GetXaxis();
    const TAxis* ay = h3->GetYaxis();
    const int nK = ay->GetNbins();
    std::unique_ptr<TH3D> h3c[kNLevels];
    std::unique_ptr<TProfile2D> hPc[kNLevels][2];
    for (int l = 0; l < kNLevels; l++) {
      h3c[l]    = Cumulative<TH3D>(fin, "hMKC", tag, l);
      hPc[l][0] = Cumulative<TProfile2D>(fin, "hP1", tag, l);
      hPc[l][1] = Cumulative<TProfile2D>(fin, "hP2", tag, l);
    }
    for (const TestMass& tm : tests) {
      if (std::abs(kPairPDG[p][0]) != tm.pdg1 || std::abs(kPairPDG[p][1]) != tm.pdg1) continue;
      int b1, b2;
      double poleMass = 0.;
      std::string label;
      if (tm.poleParent) {
        int a = db->GetParticleTypeIndex(tm.poleParent);
        poleMass = db->GetParticleType(a)->GetMass();
        b1 = b2 = ax->FindFixBin(poleMass);
        label = std::string("pole ") + tm.poleParent;
      } else {
        int b = ax->FindFixBin(tm.M);
        b1 = b - 2;
        b2 = b + 2;
        char buf[32];
        std::snprintf(buf, sizeof(buf), "cont %.2f", tm.M);
        label = buf;
      }
      double M1 = ax->GetBinLowEdge(b1), M2 = ax->GetBinUpEdge(b2);
      std::vector<double> n(nK), mu(nK, 0.);
      std::string hkName = "kappa_" + std::string(mcLev) + "_" + tag + "_" + label;
      TH1D* hk = h3->ProjectionY(hkName.c_str(), b1, b2, 1, h3->GetZaxis()->GetNbins());
      for (int k = 1; k <= nK; k++) n[k - 1] = hk->GetBinContent(k);
      for (size_t i = 0; i < cs.contribs.size(); i++) {
        const Contribution& co = cs.contribs[i];
        if (co.pair != p) continue;
        std::vector<double> m = PredictKappa(co, M1, M2, ay);
        for (int k = 0; k < nK; k++) mu[k] += NV * scale[i] * m[k];
      }
      // L2: chain part in the window's M bins (same M and kappa binning as the MC)
      std::vector<double> extraK(nK, 0.);
      if (l2) {
        TH2D* hc  = ftab->Get<TH2D>(("chainMK_" + tag).c_str());
        TH2D* hcv = ftab->Get<TH2D>(("chainMK_var_" + tag).c_str());
        for (int b = b1; b <= b2; b++)
          for (int k = 1; k <= nK; k++) {
            mu[k - 1] += NV * hc->GetBinContent(b, k);
            extraK[k - 1] += NV * NV * hcv->GetBinContent(b, k);
          }
      }
      double nTot = 0., muTot = 0.;
      for (int k = 0; k < nK; k++) {
        nTot += n[k];
        muTot += mu[k];
      }

      // Boltzmann reference in the kappa bins, averaged over the M bins of the window with the
      // weights w (pole window: at the pole mass).
      auto reference = [&](const std::vector<double>& w) {
        std::vector<double> ref(nK, 0.);
        double wTot = 0.;
        for (int b = b1; b <= b2; b++) {
          double Mb = tm.poleParent ? poleMass : ax->GetBinCenter(b);
          for (int k = 1; k <= nK; k++)
            ref[k - 1] += w[b - b1] *
                          BoltzmannRefKappa(Mb, ay->GetBinLowEdge(k), ay->GetBinUpEdge(k), T);
          wTot += w[b - b1];
        }
        for (double& r : ref) r /= wTot;
        return ref;
      };
      // Writes dP/dkappa (per GeV) and its ratio to the reference for counts c with errors e.
      auto writeCurve = [&](const std::string& lev, const std::vector<double>& c,
                            const std::vector<double>& e, const std::vector<double>& ref) {
        double tot = 0.;
        for (double x : c) tot += x;
        TH1D* hd = static_cast<TH1D*>(hk->Clone(("dPdk_" + lev + "_" + tag + "_" + label).c_str()));
        TH1D* hr =
            static_cast<TH1D*>(hk->Clone(("dPdkRatio_" + lev + "_" + tag + "_" + label).c_str()));
        hd->Reset();
        hr->Reset();
        for (int k = 1; k <= nK; k++) {
          double w = ay->GetBinWidth(k);
          hd->SetBinContent(k, c[k - 1] / tot / w);
          hd->SetBinError(k, e[k - 1] / tot / w);
          if (ref[k - 1] > 0.) {
            hr->SetBinContent(k, c[k - 1] / tot / ref[k - 1]);
            hr->SetBinError(k, e[k - 1] / tot / ref[k - 1]);
          }
        }
        fout.cd();
        hd->Write();
        hr->Write();
        return hr;
      };

      // T2.2 (unlike-sign pairs only); at L2 also like-sign pairs in the continuum windows, with the
      // doublet variance 2 mu
      const bool tested = unlikeSign || (l2 && likeSign && !tm.poleParent);
      if (tested && muTot < 500.) {
        std::printf("%s  %-9s %-14s [%.4f,%.4f]  skipped (%.0f predicted pairs)\n", tKap,
                    tag.c_str(), label.c_str(), M1, M2, muTot);
      } else if (tested) {
        // Shape: expectations rescaled to the MC total (one constraint).
        std::vector<double> muShape(nK);
        std::vector<double> extraShape(nK);
        for (int k = 0; k < nK; k++) {
          muShape[k]    = mu[k] * nTot / muTot;
          extraShape[k] = extraK[k] * (nTot / muTot) * (nTot / muTot);
        }
        Chi2 c = PoissonChi2(n, muShape, minExp, 1, likeSign ? 2. : 1., &extraShape);
        double prob = TMath::Prob(c.chi2, c.ndf);
        bool ok = prob > 1e-3;
        allT22 = allT22 && ok;
        std::printf("%s  %-9s %-14s [%.4f,%.4f] %10.0f %10.1f %7.4f %7.1f %4d %8.3f %9.3g  %s\n",
                    tKap, tag.c_str(), label.c_str(), M1, M2, nTot, muTot, nTot / muTot, c.chi2, c.ndf,
                    c.chi2 / c.ndf, prob, PassFail(ok));
      }

      // M3: table curve and its ratio to the reference (weights: table counts per M bin).
      if (muTot >= minM3) {
        std::vector<double> w(b2 - b1 + 1);
        for (int b = b1; b <= b2; b++) w[b - b1] = muMPair[p][b - 1];
        std::vector<double> ref = reference(w);
        TH1D* hr = writeCurve("pred", mu, std::vector<double>(nK, 0.), ref);
        char buf[200];
        std::snprintf(buf, sizeof(buf),
                      "M3ref %-9s %-14s table/Boltzmann: kappa 0-50 MeV %.4f, 200-250 MeV %.4f, "
                      "500-550 MeV %.4f, 1.0-1.1 GeV %.4f",
                      tag.c_str(), label.c_str(), hr->GetBinContent(1), hr->GetBinContent(5),
                      hr->GetBinContent(11), hr->GetBinContent(21));
        m3Lines.push_back(buf);
      }
      // M3: MC levels
      for (int l = 0; l < kNLevels; l++) {
        TH1D* hkl = h3c[l]->ProjectionY(Form("kappaLev_%s_%s_%s", kLevel[l], tag.c_str(),
                                             label.c_str()),
                                        b1, b2, 1, h3c[l]->GetZaxis()->GetNbins());
        std::vector<double> c(nK), e(nK);
        double tot = 0.;
        for (int k = 1; k <= nK; k++) {
          c[k - 1] = hkl->GetBinContent(k);
          e[k - 1] = std::sqrt(c[k - 1]) * errScale;
          tot += c[k - 1];
        }
        delete hkl;
        if (tot < minM3) continue;
        std::vector<double> w(b2 - b1 + 1);
        TH1D* hml = h3c[l]->ProjectionX("_px_tmp", 1, nK, 1, h3c[l]->GetZaxis()->GetNbins());
        for (int b = b1; b <= b2; b++) w[b - b1] = hml->GetBinContent(b);
        delete hml;
        writeCurve(kLevel[l], c, e, reference(w));

        // a1(kappa), a2(kappa) at this M
        for (int il = 0; il < 2; il++) {
          if (il == 0 && likeSign) continue; // antisymmetric by construction
          std::unique_ptr<TProfile> py(hPc[l][il]->ProfileY("_pfy_tmp", b1, b2));
          TH1D* ha = ProfileToHist(py.get(), Form("a%d_%s_%s_%s", il + 1, kLevel[l], tag.c_str(),
                                                  label.c_str()),
                                   errScale);
          fout.cd();
          ha->Write();
        }
        // raw cos theta*
        for (int ib = 0; ib < 2; ib++) {
          int k1 = ay->FindFixBin(cthBand[ib][0] + 1e-9), k2 = ay->FindFixBin(cthBand[ib][1] - 1e-9);
          TH1D* hc = h3c[l]->ProjectionZ(Form("cth_%s_%s_%s_k%d", kLevel[l], tag.c_str(),
                                              label.c_str(), ib),
                                         b1, b2, k1, k2);
          hc->SetTitle(Form("%s %s %s, #kappa in [%.2f,%.2f] GeV", kLevel[l], tag.c_str(),
                            label.c_str(), cthBand[ib][0], cthBand[ib][1]));
          fout.cd();
          hc->Write();
        }
      }
    }
  }
  std::printf("%s  all tested masses  %s\n", tKap, PassFail(allT22));
  for (const std::string& s : m3Lines) std::printf("%s\n", s.c_str());

  if (!l1) {
    // ---- T2.3: a1, a2 at L0; T3.1: a1, a2 at L1 (and L1x alone); L2 for information ----------
    const double minEntries = 50.;
    std::printf("T2.3  %-9s %3s %-10s %6s %8s %8s  %s\n", "pair", "a_l", "bins", "n", "mean",
                "width", "(tolerance mean, width)");
    std::vector<std::string> t31Head;
    for (int p = 0; p < kNPairs; p++) {
      std::string tag = kPairTag[p];
      bool likeSign = (kPairPDG[p][0] == kPairPDG[p][1]);
      // M bins containing a pole of this pair
      std::vector<double> poleMasses;
      for (const Contribution& co : cs.contribs)
        if (co.pair == p && co.third < 0) poleMasses.push_back(co.ma);
      for (int l = 1; l <= 2; l++) {
        if (l == 1 && likeSign) continue; // antisymmetric by construction
        TProfile2D* prof = fin.Get<TProfile2D>(Form("hP%d_L0_%s", l, tag.c_str()));
        char head[64];
        std::snprintf(head, sizeof(head), "T2.3  %-9s", tag.c_str());
        allT23 = PullTest(head, prof, l, likeSign, poleMasses, minEntries) && allT23;
      }
    }
    std::printf("T2.3  all ordered pairs  %s\n", PassFail(allT23));

    std::printf("T3.1  %-9s %-3s %3s %-10s %6s %8s %8s  %s\n", "pair", "lev", "a_l", "bins", "n",
                "mean", "width", "(tolerance mean, width)");
    for (int p = 0; p < kNPairs; p++) {
      std::string tag = kPairTag[p];
      bool likeSign = (kPairPDG[p][0] == kPairPDG[p][1]);
      std::vector<double> poleMasses;
      for (const Contribution& co : cs.contribs)
        if (co.pair == p && co.third < 0) poleMasses.push_back(co.ma);
      for (int l = 1; l <= 2; l++) {
        if (l == 1 && likeSign) continue;
        const std::string pre = "hP" + std::to_string(l);
        auto profL1 = Cumulative<TProfile2D>(fin, pre.c_str(), tag, 1);
        TProfile2D* profL1x = fin.Get<TProfile2D>((pre + "_L1x_" + tag).c_str());
        auto profL2 = Cumulative<TProfile2D>(fin, pre.c_str(), tag, 2);
        char head[64];
        std::snprintf(head, sizeof(head), "T3.1  %-9s %-3s", tag.c_str(), "L1");
        allT31 = PullTest(head, profL1.get(), l, likeSign, poleMasses, minEntries) && allT31;
        std::snprintf(head, sizeof(head), "T3.1  %-9s %-3s", tag.c_str(), "L1x");
        allT31 = PullTest(head, profL1x, l, likeSign, poleMasses, minEntries) && allT31;
        std::snprintf(head, sizeof(head), "L2pul %-9s %-3s", tag.c_str(), "L2");
        PullTest(head, profL2.get(), l, likeSign, poleMasses, minEntries, true);
      }
    }
    std::printf("T3.1  all ordered pairs  %s\n", PassFail(allT31));

    // ---- M3: a1, a2 and C_Q/(C_K + C_u) in (M, eps_K) cells, levels L1 and L2 --------------------
    // Cell averages from the profile sums; C_Q/(C_K + C_u) = a1 / <2 kappa p*/(3 M eps_K)>, the
    // factor averaged over the profile bins of the cell with their entries as weights (equal masses).
    const std::vector<double> kBands = {0.0, 0.2, 0.4, 0.7, 1.0, 2.0};
    std::printf("M3tab %-9s %-3s %-13s %-11s %6s %6s %10s %10s %10s %10s %10s %10s %11s %11s\n",
                "pair", "lev", "M band", "kappa band", "<M>", "<epsK>", "N", "a1", "err", "a2",
                "err", "factor", "CQ/(CK+Cu)", "err");
    for (int p = 0; p < kNPairs; p++) {
      std::string tag = kPairTag[p];
      bool likeSign = (kPairPDG[p][0] == kPairPDG[p][1]);
      const double errScale = likeSign ? std::sqrt(2.) : 1.;
      if (p == PairIndex(-211, 211) || p == PairIndex(-321, 321)) continue; // mirror pairs
      std::vector<double> mBands;
      if (std::abs(kPairPDG[p][0]) == 211)
        mBands = {0.279, 0.40, 0.50, 0.70, 0.85, 1.20, 2.0};
      else if (std::abs(kPairPDG[p][0]) == 321)
        mBands = {0.987, 1.03, 1.20, 1.50, 2.0};
      else
        continue;
      const double m1 = PairMasses(db, p).first;
      for (int l = 1; l <= 2; l++) {
        auto pr1 = Cumulative<TProfile2D>(fin, "hP1", tag, l);
        auto pr2 = Cumulative<TProfile2D>(fin, "hP2", tag, l);
        const int nMB = mBands.size() - 1, nKB = kBands.size() - 1;
        TH2D* hA[2];
        for (int il = 0; il < 2; il++)
          hA[il] = new TH2D(Form("a%d_cells_%s_%s", il + 1, kLevel[l], tag.c_str()),
                            Form("a_{%d} %s %s;M (GeV);#kappa (GeV)", il + 1, kLevel[l], tag.c_str()),
                            nMB, mBands.data(), nKB, kBands.data());
        TH2D* hR = new TH2D(Form("CQratio_cells_%s_%s", kLevel[l], tag.c_str()),
                            Form("C_{Q}/(C_{K}+C_{u}) %s %s;M (GeV);#kappa (GeV)", kLevel[l],
                                 tag.c_str()),
                            nMB, mBands.data(), nKB, kBands.data());
        TH2D* hE = new TH2D(Form("epsK_cells_%s_%s", kLevel[l], tag.c_str()),
                            Form("<#epsilon_{K}> %s %s;M (GeV);#kappa (GeV)", kLevel[l], tag.c_str()),
                            nMB, mBands.data(), nKB, kBands.data());
        for (int im = 0; im < nMB; im++) {
          for (int ik = 0; ik < nKB; ik++) {
            double N = 0., sy[2] = {0., 0.}, syy[2] = {0., 0.}, sM = 0., sE = 0., sF = 0.;
            for (int ix = 1; ix <= pr1->GetNbinsX(); ix++) {
              double Mc = pr1->GetXaxis()->GetBinCenter(ix);
              if (Mc < mBands[im] || Mc >= mBands[im + 1]) continue;
              for (int iy = 1; iy <= pr1->GetNbinsY(); iy++) {
                double kc = pr1->GetYaxis()->GetBinCenter(iy);
                if (kc < kBands[ik] || kc >= kBands[ik + 1]) continue;
                int b = pr1->GetBin(ix, iy);
                double nb = pr1->GetBinEntries(b);
                if (nb <= 0.) continue;
                N += nb;
                // profile sums: sum y = mean * entries; sum y^2 is kept in fSumw2
                sy[0] += pr1->GetBinContent(b) * nb;
                syy[0] += pr1->GetSumw2()->At(b);
                sy[1] += pr2->GetBinContent(b) * nb;
                syy[1] += pr2->GetSumw2()->At(b);
                double eK = std::hypot(Mc, kc);
                sM += nb * Mc;
                sE += nb * eK;
                sF += nb * 2. * kc * PStar(Mc, m1, m1) / (3. * Mc * eK);
              }
            }
            if (N < minEntries) continue;
            double a[2], ea[2];
            for (int il = 0; il < 2; il++) {
              a[il]  = sy[il] / N;
              ea[il] = std::sqrt(std::max(syy[il] / N - a[il] * a[il], 0.) / N) * errScale;
              hA[il]->SetBinContent(im + 1, ik + 1, a[il]);
              hA[il]->SetBinError(im + 1, ik + 1, ea[il]);
            }
            double fac = sF / N;
            hE->SetBinContent(im + 1, ik + 1, sE / N);
            if (!likeSign) {
              hR->SetBinContent(im + 1, ik + 1, a[0] / fac);
              hR->SetBinError(im + 1, ik + 1, ea[0] / fac);
            }
            char mb[32], kb[32];
            std::snprintf(mb, sizeof(mb), "[%.3f,%.2f)", mBands[im], mBands[im + 1]);
            std::snprintf(kb, sizeof(kb), "[%.1f,%.1f)", kBands[ik], kBands[ik + 1]);
            if (likeSign)
              std::printf("M3tab %-9s %-3s %-13s %-11s %6.3f %6.3f %10.0f %10s %10s %+10.5f %10.5f "
                          "%10.5f %11s %11s\n",
                          tag.c_str(), kLevel[l], mb, kb, sM / N, sE / N, N, "-", "-", a[1], ea[1],
                          fac, "-", "-");
            else
              std::printf("M3tab %-9s %-3s %-13s %-11s %6.3f %6.3f %10.0f %+10.5f %10.5f %+10.5f "
                          "%10.5f %10.5f %+11.4f %11.4f\n",
                          tag.c_str(), kLevel[l], mb, kb, sM / N, sE / N, N, a[0], ea[0], a[1],
                          ea[1], fac, a[0] / fac, ea[0] / fac);
          }
        }
        fout.cd();
        hA[0]->Write();
        hA[1]->Write();
        hR->Write();
        hE->Write();
      }
    }
  }

  fout.Close();
  return 0;
}
