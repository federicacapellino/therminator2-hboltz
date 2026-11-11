// ex_pair_functions_harmonics.cxx — spherical-harmonic anisotropy test of the pair MC (M9)
//
// See hrg_pair_density_spec.md, Sec. 10. Without --mc only T9.0 runs. With --mc (the hadd-merged
// output of ex_pair_functions_mc.exe --harmonics) it runs
//   T9.0  the 24 real Y_lm (l <= 4) are orthonormal and have parity (-1)^l
//   T9.1  a10, a20 equal sqrt(3/4pi) <P1>, sqrt(5/4pi) <P2> of the reference M2 file (--ref) bin by
//         bin, exclusive levels L0, L1x, L2x (same seeds = same events); statistical if the
//         entries differ
//   T9.2  m != 0 (c and s), l <= 4 (even l only for like-sign), levels L0, L1, L2: pulls
//   T9.3  odd l, m = 0, unlike-sign pairs, levels L0, L1, L2: pulls
//   T9.4  even l, m = 0, levels L0, L1: pulls
// and reports a20, a40 at the cumulative levels versus kappa at the M3 test masses (M9kap) and in
// the M3 (M, eps_K) cells (M9tab). Pull criterion as T2.3: bins with >= 50 entries,
// |mean| < 3/sqrt(n), |width - 1| < 3/sqrt(2n); like-sign errors of even-l moments are scaled by
// sqrt(2) (doublets). Each test is run per (pair, level, function) and pooled per (level, function)
// and per level over the pairs; pooling leaves out the mirror pairs (pi-pi+, K-K+, pbar p), which
// are the (+-) pairs with n -> -n. A test passes if every pooled test passes and the number of
// failed individual tests is consistent with the false-alarm rate of the criterion.
//
// Usage:
//   ex_pair_functions_harmonics.exe [--mc <file>] [--ref <file>] [--share <dir>] [--output <file>]
//                                   [--minentries <N>]
// The profile errors are estimated from the spread in the bin, so pulls of bins with few entries
// are Student-t like, width about sqrt(N/(N-3)) (+2.5% at N = 50). Pooled over 1e5 bins this
// exceeds the width tolerance; --minentries 500 keeps the bias below it (M9 report).

#include "CLI11.hpp"
#include "PairFunctions.h"
#include "ParticleDB.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TNamed.h>
#include <TProfile2D.h>
#include <TRandom3.h>
#include <cstdio>
#include <iostream>
#include <memory>

static const char* PassFail(bool ok) { return ok ? "PASS" : "FAIL"; }

static constexpr int kNLevels = 3;
static const char* const kLevel[kNLevels]   = {"L0", "L1", "L2"};
static const char* const kMCLevel[kNLevels] = {"L0", "L1x", "L2x"};

// Sum of the exclusive MC levels up to `level` (cumulative L0, L1, L2), detached from any file.
static std::unique_ptr<TProfile2D> Cumulative(TFile& f, const std::string& prefix,
                                              const std::string& tag, int level) {
  std::unique_ptr<TProfile2D> h(static_cast<TProfile2D*>(
      f.Get<TProfile2D>(Form("%s_L0_%s", prefix.c_str(), tag.c_str()))
          ->Clone(Form("%s_%s_%s_cum", prefix.c_str(), kLevel[level], tag.c_str()))));
  h->SetDirectory(nullptr);
  for (int l = 1; l <= level; l++)
    h->Add(f.Get<TProfile2D>(Form("%s_%s_%s", prefix.c_str(), kMCLevel[l], tag.c_str())));
  return h;
}

static bool IsLikeSign(int p) { return kPairPDG[p][0] == kPairPDG[p][1]; }
static bool IsMirror(int p) { return kPairPDG[p][0] < 0 && kPairPDG[p][1] > 0; }

// Pulls a / sigma(a) over the bins with at least minEntries entries; errors scaled by errScale.
static std::vector<double> Pulls(const TProfile2D* prof, double errScale, double minEntries) {
  std::vector<double> v;
  for (int ix = 1; ix <= prof->GetNbinsX(); ix++)
    for (int iy = 1; iy <= prof->GetNbinsY(); iy++) {
      int b = prof->GetBin(ix, iy);
      if (prof->GetBinEntries(b) < minEntries || prof->GetBinError(b) <= 0.) continue;
      v.push_back(prof->GetBinContent(b) / (prof->GetBinError(b) * errScale));
    }
  return v;
}

struct PullStat {
  double n = 0., mean = 0., width = 0., tolM = 0., tolW = 0.;
  bool ok = true;
};

static PullStat Stat(const std::vector<double>& v) {
  PullStat s;
  s.n = v.size();
  if (v.empty()) return s;
  double sum = 0., sum2 = 0.;
  for (double x : v) {
    sum += x;
    sum2 += x * x;
  }
  s.mean  = sum / s.n;
  s.width = std::sqrt(std::max(sum2 / s.n - s.mean * s.mean, 0.));
  s.tolM  = 3. / std::sqrt(s.n);
  s.tolW  = 3. / std::sqrt(2. * s.n);
  s.ok    = std::abs(s.mean) < s.tolM && std::abs(s.width - 1.) < s.tolW;
  return s;
}

// False-alarm probability of the pull criterion for n standard normal pulls. The sample mean and
// width are independent, mean ~ N(0, 1/n) and n width^2 ~ chi2(n - 1).
static double FalseAlarm(int n) {
  double pMean = 1. - TMath::Erfc(3. / std::sqrt(2.));
  if (n < 2) return 1. - pMean;
  double tolW = 3. / std::sqrt(2. * n);
  double lo = std::max(1. - tolW, 0.), hi = 1. + tolW;
  double pWidth = TMath::Prob(n * lo * lo, n - 1) - TMath::Prob(n * hi * hi, n - 1);
  return 1. - pMean * pWidth;
}

// Poisson probability of at least k events for mean mu.
static double PoissonTail(int k, double mu) {
  double term = std::exp(-mu), cdf = 0.;
  for (int i = 0; i < k; i++) {
    cdf += term;
    term *= mu / (i + 1);
  }
  return std::max(0., 1. - cdf);
}

// One pull test of the M9 family: accumulates individual and pooled results.
struct TestFamily {
  std::string name;
  int nTests = 0, nFail = 0;
  double expFail = 0.;
  bool pooledOk = true;
  std::map<std::string, std::vector<double>> pooled; // key: level and function, or level only
};

static void Individual(TestFamily& fam, const std::string& head, const std::vector<double>& v) {
  if (v.empty()) return;
  PullStat s = Stat(v);
  ++fam.nTests;
  fam.nFail += !s.ok;
  fam.expFail += FalseAlarm(v.size());
  std::printf("%si %s %6.0f %+8.3f %8.3f  (%.3f, %.3f)  %s\n", fam.name.c_str(), head.c_str(), s.n,
              s.mean, s.width, s.tolM, s.tolW, PassFail(s.ok));
}

static bool Finish(TestFamily& fam) {
  for (auto& kv : fam.pooled) {
    PullStat s = Stat(kv.second);
    if (kv.second.empty()) continue;
    fam.pooledOk = fam.pooledOk && s.ok;
    std::printf("%sp %-24s %8.0f %+8.3f %8.3f  (%.3f, %.3f)  %s\n", fam.name.c_str(),
                kv.first.c_str(), s.n, s.mean, s.width, s.tolM, s.tolW, PassFail(s.ok));
  }
  double p  = PoissonTail(fam.nFail, fam.expFail);
  bool ok   = fam.pooledOk && p > 1e-3;
  std::printf("%s  individual tests %d, failed %d, expected %.2f (Poisson p = %.3g); pooled %s"
              "  ->  %s\n",
              fam.name.c_str(), fam.nTests, fam.nFail, fam.expFail, p, PassFail(fam.pooledOk),
              PassFail(ok));
  return ok;
}

// ---- T9.0 -----------------------------------------------------------------------------------------
static bool TestT90() {
  const GaussLegendre& gl = GLRule(32);
  const int nPhi = 64;
  std::vector<double> G(kNHarm * kNHarm, 0.);
  double y[kNHarm];
  for (size_t i = 0; i < gl.x.size(); i++)
    for (int j = 0; j < nPhi; j++) {
      RealYlm(gl.x[i], 2. * M_PI * j / nPhi, y);
      double w = gl.w[i] * 2. * M_PI / nPhi;
      for (int a = 0; a < kNHarm; a++)
        for (int b = 0; b < kNHarm; b++) G[a * kNHarm + b] += w * y[a] * y[b];
    }
  double worstOrtho = 0.;
  for (int a = 0; a < kNHarm; a++)
    for (int b = 0; b < kNHarm; b++)
      worstOrtho = std::max(worstOrtho, std::abs(G[a * kNHarm + b] - (a == b ? 1. : 0.)));

  TRandom3 rng(9);
  double worstParity = 0., worstPl = 0., ym[kNHarm];
  for (int t = 0; t < 10000; t++) {
    double c = rng.Uniform(-1., 1.), phi = rng.Uniform(0., 2. * M_PI);
    RealYlm(c, phi, y);
    RealYlm(-c, phi + M_PI, ym);
    for (int i = 0; i < kNHarm; i++)
      worstParity = std::max(worstParity, std::abs(ym[i] - (HarmL(i) % 2 ? -1. : 1.) * y[i]));
    const double P[5] = {1., c, 0.5 * (3 * c * c - 1), 0.5 * (5 * c * c * c - 3 * c),
                         (35 * c * c * c * c - 30 * c * c + 3) / 8.};
    for (int l = 1; l <= kHarmLMax; l++)
      worstPl = std::max(worstPl,
                         std::abs(y[l * l - 1] - std::sqrt((2. * l + 1.) / (4. * M_PI)) * P[l]));
  }
  bool ok = worstOrtho < 1e-12 && worstParity < 1e-12 && worstPl < 1e-12;
  std::printf("T9.0  functions:");
  for (int i = 0; i < kNHarm; i++) std::printf(" %s", HarmTag(i).c_str());
  std::printf("\nT9.0  orthonormality max |<Y_a,Y_b> - delta_ab| = %.2e (32 x %d quadrature)\n",
              worstOrtho, nPhi);
  std::printf("T9.0  parity max |Y(-n) - (-1)^l Y(n)| = %.2e (1e4 random directions)\n", worstParity);
  std::printf("T9.0  Y_l0 vs sqrt((2l+1)/4pi) P_l: max diff %.2e\n", worstPl);
  std::printf("T9.0  %s\n", PassFail(ok));
  return ok;
}

// ---- T9.1 -----------------------------------------------------------------------------------------
static bool TestT91(TFile& fin, TFile& fref, double minEntries) {
  bool exact = true; // same entries everywhere
  double worstRel = 0.;
  std::string worstWhere;
  std::vector<double> statPulls;
  long nBins = 0;
  for (int l = 0; l < kNLevels; l++)
    for (int p = 0; p < kNPairs; p++)
      for (int il = 1; il <= 2; il++) {
        std::string sfx = std::string(kMCLevel[l]) + "_" + kPairTag[p];
        auto* hy = fin.Get<TProfile2D>(Form("hY%d0_%s", il, sfx.c_str()));
        auto* hp = fref.Get<TProfile2D>(Form("hP%d_%s", il, sfx.c_str()));
        const double c = std::sqrt((2. * il + 1.) / (4. * M_PI));
        for (int b = 0; b < hy->GetNcells(); b++) {
          double ny = hy->GetBinEntries(b), np = hp->GetBinEntries(b);
          if (ny != np) exact = false;
          if (ny <= 0. || np <= 0.) continue;
          ++nBins;
          double a = hy->GetBinContent(b), r = c * hp->GetBinContent(b);
          double ea = hy->GetBinError(b), er = c * hp->GetBinError(b);
          double denom = std::max({std::abs(a), std::abs(r), ea});
          if (denom > 0. && std::abs(a - r) / denom > worstRel) {
            worstRel   = std::abs(a - r) / denom;
            worstWhere = Form("a%d0 %s bin %d (%.0f entries)", il, sfx.c_str(), b, ny);
          }
          if (ny >= minEntries && np >= minEntries && ea > 0. && er > 0.)
            statPulls.push_back((a - r) / std::hypot(ea, er));
        }
      }
  std::printf("T9.1  %ld non-empty bins (a10, a20; L0, L1x, L2x; 12 pairs); entries %s\n", nBins,
              exact ? "identical in every bin" : "DIFFER");
  if (exact) {
    bool ok = worstRel < 1e-10;
    std::printf("T9.1  max |a - sqrt((2l+1)/4pi)<P_l>| / max(|a|, sigma) = %.2e at %s  %s\n",
                worstRel, worstWhere.c_str(), PassFail(ok));
    return ok;
  }
  PullStat s = Stat(statPulls);
  std::printf("T9.1  events not reproduced: statistical comparison (independent errors) over %.0f "
              "bins, pull mean %+.3f width %.3f (%.3f, %.3f)  %s\n",
              s.n, s.mean, s.width, s.tolM, s.tolW, PassFail(s.ok));
  return s.ok;
}

// ---- Report: a20, a40 ----------------------------------------------------------------------------
// Mean and error of a profile over the M bins [b1, b2] and kappa bin iy (sums of the profile).
struct Moment {
  double N = 0., a = 0., e = 0.;
};
static Moment Combine(const TProfile2D* pr, int b1, int b2, int iy1, int iy2, double errScale) {
  Moment m;
  double sy = 0., syy = 0.;
  for (int ix = b1; ix <= b2; ix++)
    for (int iy = iy1; iy <= iy2; iy++) {
      int b     = pr->GetBin(ix, iy);
      double nb = pr->GetBinEntries(b);
      if (nb <= 0.) continue;
      m.N += nb;
      sy += pr->GetBinContent(b) * nb;
      syy += pr->GetSumw2()->At(b);
    }
  if (m.N <= 0.) return m;
  m.a = sy / m.N;
  m.e = std::sqrt(std::max(syy / m.N - m.a * m.a, 0.) / m.N) * errScale;
  return m;
}

int main(int argc, char** argv) {
  std::string mcFile     = "";
  std::string refFile    = "";
  std::string shareDir   = "";
  std::string outputFile = "";
  double minEntries      = 50.;

  CLI::App app{"ex_pair_functions_harmonics: spherical-harmonic anisotropy tests (M9)"};
  app.add_option("--mc",     mcFile,     "merged output of ex_pair_functions_mc.exe --harmonics");
  app.add_option("--ref",    refFile,    "M2 output (hP1, hP2) of the same seeds, for T9.1");
  app.add_option("--share",  shareDir,   "path to SHARE particle database");
  app.add_option("--output", outputFile, "output ROOT file (default: harmonics_<weak>.root)");
  app.add_option("--minentries", minEntries, "minimum entries of a tested bin")->default_val(50.);
  CLI11_PARSE(app, argc, argv);

  bool okT90 = TestT90();
  if (mcFile.empty()) return okT90 ? 0 : 2;

  TFile fin(mcFile.c_str(), "READ");
  if (fin.IsZombie() || !fin.Get<TProfile2D>("hY10_L0_pippim")) {
    std::cerr << "Error: " << mcFile << " has no harmonic profiles\n";
    return 1;
  }
  std::string weak = fin.Get<TNamed>("weak")->GetTitle();
  double nEv       = fin.Get<TH1D>("hEvents")->GetBinContent(1);
  std::printf("MC file %s: weak %s, %.0f events\n", mcFile.c_str(), weak.c_str(), nEv);
  if (outputFile.empty()) outputFile = "harmonics_weak" + weak + ".root";

  auto db = std::make_unique<ParticleDB>(shareDir);

  bool okT91 = true;
  if (!refFile.empty()) {
    TFile fref(refFile.c_str(), "READ");
    std::string refWeak = fref.Get<TNamed>("weak")->GetTitle();
    std::printf("T9.1  reference %s: weak %s, %.0f events\n", refFile.c_str(), refWeak.c_str(),
                fref.Get<TH1D>("hEvents")->GetBinContent(1));
    okT91 = (refWeak == weak) && TestT91(fin, fref, minEntries);
  } else {
    std::printf("T9.1  skipped (no --ref)\n");
  }

  // ---- T9.2 - T9.4 ------------------------------------------------------------------------------
  TestFamily t92{"T9.2"}, t93{"T9.3"}, t94{"T9.4"};
  std::printf("T9.xi %-9s %-3s %-4s %6s %8s %8s  %s\n", "pair", "lev", "Y", "bins", "mean", "width",
              "(tolerance mean, width)");
  TFile fout(outputFile.c_str(), "RECREATE");
  TH1D* hPull[3] = {new TH1D("pulls_T92", "T9.2 pulls;a/#sigma;bins", 100, -5., 5.),
                    new TH1D("pulls_T93", "T9.3 pulls;a/#sigma;bins", 100, -5., 5.),
                    new TH1D("pulls_T94", "T9.4 pulls;a/#sigma;bins", 100, -5., 5.)};
  for (int l = 0; l < kNLevels; l++) {
    for (int p = 0; p < kNPairs; p++) {
      const bool likeSign = IsLikeSign(p);
      for (int i = 0; i < kNHarm; i++) {
        const int L = HarmL(i), m = HarmM(i);
        TestFamily* fam = nullptr;
        if (m != 0 && (!likeSign || L % 2 == 0))
          fam = &t92;
        else if (m == 0 && L % 2 == 1 && !likeSign)
          fam = &t93;
        else if (m == 0 && L % 2 == 0 && l <= 1)
          fam = &t94;
        if (!fam) continue;
        auto prof = Cumulative(fin, "hY" + HarmTag(i), kPairTag[p], l);
        double errScale       = (likeSign && L % 2 == 0) ? std::sqrt(2.) : 1.;
        std::vector<double> v = Pulls(prof.get(), errScale, minEntries);
        char head[64];
        std::snprintf(head, sizeof(head), "%-9s %-3s %-4s", kPairTag[p], kLevel[l],
                      HarmTag(i).c_str());
        Individual(*fam, head, v);
        if (IsMirror(p)) continue;
        auto& byFun = fam->pooled[std::string(kLevel[l]) + " Y" + HarmTag(i)];
        auto& byLev = fam->pooled[std::string(kLevel[l]) + " all"];
        byFun.insert(byFun.end(), v.begin(), v.end());
        byLev.insert(byLev.end(), v.begin(), v.end());
        for (double x : v) hPull[fam == &t92 ? 0 : (fam == &t93 ? 1 : 2)]->Fill(x);
      }
    }
  }
  bool okT92 = Finish(t92);
  bool okT93 = Finish(t93);
  bool okT94 = Finish(t94);

  // ---- Report: a20, a40 versus kappa at the M3 test masses ----------------------------------------
  struct TestMass {
    int pdg1;
    const char* poleParent;
    double M;
  };
  const TestMass tests[] = {{211, "Ka0492zer", 0.}, {211, "rho770zer", 0.}, {211, "om0782zer", 0.},
                            {211, "f00980zer", 0.}, {211, nullptr, 0.30},   {211, nullptr, 0.45},
                            {211, nullptr, 0.60},   {321, "ph1020zer", 0.}, {321, nullptr, 1.10},
                            {321, nullptr, 1.30}};
  const int reportPairs[] = {PairIndex(211, -211), PairIndex(211, 211), PairIndex(321, -321)};
  const int harmIdx[2]    = {3, 15}; // Y20, Y40
  std::printf("M9kap %-9s %-3s %-16s %-12s %10s %10s %9s %10s %9s %9s %8s\n", "pair", "lev", "mass",
              "kappa", "N", "a20", "err", "a40", "err", "a40/a20", "err");
  for (int p : reportPairs) {
    const std::string tag = kPairTag[p];
    const double errScale = IsLikeSign(p) ? std::sqrt(2.) : 1.;
    for (int l = 0; l < kNLevels; l++) {
      std::unique_ptr<TProfile2D> pr[2] = {Cumulative(fin, "hY20", tag, l),
                                           Cumulative(fin, "hY40", tag, l)};
      const TAxis* ax = pr[0]->GetXaxis();
      const TAxis* ay = pr[0]->GetYaxis();
      for (const TestMass& tm : tests) {
        if (std::abs(kPairPDG[p][0]) != tm.pdg1) continue;
        int b1, b2;
        std::string label;
        if (tm.poleParent) {
          double mass = db->GetParticleType(db->GetParticleTypeIndex(tm.poleParent))->GetMass();
          b1 = b2 = ax->FindFixBin(mass);
          label   = std::string("pole ") + tm.poleParent;
        } else {
          int b = ax->FindFixBin(tm.M);
          b1    = b - 2;
          b2    = b + 2;
          label = Form("cont %.2f", tm.M);
        }
        TH1D* h[2];
        for (int j = 0; j < 2; j++) {
          h[j] = new TH1D(Form("a%s_%s_%s_%s", HarmTag(harmIdx[j]).c_str(), kLevel[l], tag.c_str(),
                               label.c_str()),
                          Form("a_{%s} %s %s %s;#kappa (GeV);a", HarmTag(harmIdx[j]).c_str(),
                               kLevel[l], tag.c_str(), label.c_str()),
                          ay->GetNbins(), ay->GetXbins()->GetArray());
        }
        for (int iy = 1; iy <= ay->GetNbins(); iy++) {
          Moment m2 = Combine(pr[0].get(), b1, b2, iy, iy, errScale);
          Moment m4 = Combine(pr[1].get(), b1, b2, iy, iy, errScale);
          if (m2.N < minEntries) continue;
          h[0]->SetBinContent(iy, m2.a);
          h[0]->SetBinError(iy, m2.e);
          h[1]->SetBinContent(iy, m4.a);
          h[1]->SetBinError(iy, m4.e);
          if (l != 2) continue;
          double r = m4.a / m2.a;
          double er = std::abs(r) * std::hypot(m4.e / m4.a, m2.e / m2.a);
          std::printf("M9kap %-9s %-3s %-16s [%.2f,%.2f) %10.0f %+10.5f %9.5f %+10.5f %9.5f %+9.3f "
                      "%8.3f\n",
                      tag.c_str(), kLevel[l], label.c_str(), ay->GetBinLowEdge(iy),
                      ay->GetBinUpEdge(iy), m2.N, m2.a, m2.e, m4.a, m4.e, r, er);
        }
        fout.cd();
        h[0]->Write();
        h[1]->Write();
      }
    }
  }

  // ---- Report: (M, eps_K) cells as in M3 -----------------------------------------------------------
  const std::vector<double> kBands = {0.0, 0.2, 0.4, 0.7, 1.0, 2.0};
  std::printf("M9tab %-9s %-3s %-13s %-11s %6s %6s %10s %10s %9s %10s %9s %9s %8s\n", "pair", "lev",
              "M band", "kappa band", "<M>", "<epsK>", "N", "a20", "err", "a40", "err", "a40/a20",
              "err");
  for (int p : reportPairs) {
    const std::string tag = kPairTag[p];
    const double errScale = IsLikeSign(p) ? std::sqrt(2.) : 1.;
    std::vector<double> mBands;
    if (std::abs(kPairPDG[p][0]) == 211)
      mBands = {0.279, 0.40, 0.50, 0.70, 0.85, 1.20, 2.0};
    else
      mBands = {0.987, 1.03, 1.20, 1.50, 2.0};
    const int nMB = mBands.size() - 1, nKB = kBands.size() - 1;
    for (int l = 1; l <= 2; l++) {
      std::unique_ptr<TProfile2D> pr[2] = {Cumulative(fin, "hY20", tag, l),
                                           Cumulative(fin, "hY40", tag, l)};
      TH2D* hA[2];
      for (int j = 0; j < 2; j++)
        hA[j] = new TH2D(Form("a%s_cells_%s_%s", HarmTag(harmIdx[j]).c_str(), kLevel[l], tag.c_str()),
                         Form("a_{%s} %s %s;M (GeV);#kappa (GeV)", HarmTag(harmIdx[j]).c_str(),
                              kLevel[l], tag.c_str()),
                         nMB, mBands.data(), nKB, kBands.data());
      const TAxis* ax = pr[0]->GetXaxis();
      const TAxis* ay = pr[0]->GetYaxis();
      for (int im = 0; im < nMB; im++)
        for (int ik = 0; ik < nKB; ik++) {
          // profile bins whose centres lie in the cell (as M3)
          int bx1 = ax->GetNbins() + 1, bx2 = 0, by1 = ay->GetNbins() + 1, by2 = 0;
          for (int ix = 1; ix <= ax->GetNbins(); ix++) {
            double c = ax->GetBinCenter(ix);
            if (c >= mBands[im] && c < mBands[im + 1]) {
              bx1 = std::min(bx1, ix);
              bx2 = std::max(bx2, ix);
            }
          }
          for (int iy = 1; iy <= ay->GetNbins(); iy++) {
            double c = ay->GetBinCenter(iy);
            if (c >= kBands[ik] && c < kBands[ik + 1]) {
              by1 = std::min(by1, iy);
              by2 = std::max(by2, iy);
            }
          }
          Moment m2 = Combine(pr[0].get(), bx1, bx2, by1, by2, errScale);
          Moment m4 = Combine(pr[1].get(), bx1, bx2, by1, by2, errScale);
          if (m2.N < minEntries) continue;
          double sM = 0., sE = 0., sN = 0.;
          for (int ix = bx1; ix <= bx2; ix++)
            for (int iy = by1; iy <= by2; iy++) {
              double nb = pr[0]->GetBinEntries(pr[0]->GetBin(ix, iy));
              double Mc = ax->GetBinCenter(ix), kc = ay->GetBinCenter(iy);
              sN += nb;
              sM += nb * Mc;
              sE += nb * std::hypot(Mc, kc);
            }
          hA[0]->SetBinContent(im + 1, ik + 1, m2.a);
          hA[0]->SetBinError(im + 1, ik + 1, m2.e);
          hA[1]->SetBinContent(im + 1, ik + 1, m4.a);
          hA[1]->SetBinError(im + 1, ik + 1, m4.e);
          double r  = m4.a / m2.a;
          double er = std::abs(r) * std::hypot(m4.e / m4.a, m2.e / m2.a);
          std::printf("M9tab %-9s %-3s [%.3f,%.2f) [%.1f,%.1f) %6.3f %6.3f %10.0f %+10.5f %9.5f "
                      "%+10.5f %9.5f %+9.3f %8.3f\n",
                      tag.c_str(), kLevel[l], mBands[im], mBands[im + 1], kBands[ik],
                      kBands[ik + 1], sM / sN, sE / sN, m2.N, m2.a, m2.e, m4.a, m4.e, r, er);
        }
      fout.cd();
      hA[0]->Write();
      hA[1]->Write();
    }
  }

  bool all = okT90 && okT91 && okT92 && okT93 && okT94;
  std::printf("M9  weak %s: T9.0 %s, T9.1 %s, T9.2 %s, T9.3 %s, T9.4 %s  ->  %s\n", weak.c_str(),
              PassFail(okT90), PassFail(okT91), PassFail(okT92), PassFail(okT93), PassFail(okT94),
              PassFail(all));
  fout.cd();
  TNamed("weak", weak.c_str()).Write();
  fout.Write();
  fout.Close();
  return all ? 0 : 2;
}
