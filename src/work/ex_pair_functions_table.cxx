// ex_pair_functions_table.cxx — pair-function table for a hadron resonance gas at rest
//
// See README_pair_functions.md, Sec. 3 and milestone M1. Tabulates C_K and C_u for the 12
// ordered measured pairs from the THERMINATOR decay channels (one-step decays of thermal
// parents), writes them to a ROOT file and runs tests T1.1-T1.6.
// With --level L1 (spec Sec. 9.3, milestone M7) the thermal parent distribution is replaced by the
// full single-particle current of the parent (BuildCurrents, M6) and the tests are T7.1-T7.3;
// --norm mc normalizes the thermal densities to the Integrator multiplicities (as the MC).
// With --level L2 (spec Sec. 9.4, milestone M8) the L1 table is complemented by the chain part:
// per parent, a kernel of cascade pairs from Monte Carlo decays at rest, convolved with the parent
// current (PFC master formula); tests T8.0-T8.2 and T8.4.
//
// Usage:
//   ex_pair_functions_table.exe [--level L0|L1|L2] [--weak on|off] [--stat qs|boltzmann]
//                               [--norm quad|mc] [--share <dir>] [--output <file>]
//                               [--samples <N>] [--varsamples <N>] [--seed <N>]
//                               [--kratio <x>] [--kmin <N>] [--kmax <N>] [--mcevents <N>]
//                               [--threads <N>]

#include "CLI11.hpp"
#include "Integrator.h"
#include "Model_HRG.h"
#include "PairFunctions.h"
#include "ParticleCoor.h"
#include "ParticleDB.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TRandom.h>
#include <TTree.h>
#include <TVectorD.h>
#include <array>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <map>
#include <thread>
#include <tuple>

static double RelDiff(double a, double b) { return (b != 0.) ? std::abs(a / b - 1.) : std::abs(a); }

static const char* PassFail(bool ok) { return ok ? "PASS" : "FAIL"; }

int main(int argc, char** argv) {
  std::string shareDir   = "";
  std::string outputFile = "pair_table.root";
  std::string weak       = "on";
  std::string stat       = "qs";
  std::string level      = "L0";
  std::string norm       = "quad";
  int nSamples    = 1000000;
  int nVarSamples = 200000;
  int seed        = 0;
  double kRatio   = 4.;
  double kMin     = 2e4;
  double kMax     = 2e7;
  double mcEvents = 2e5;
  int nThreads    = 8;

  CLI::App app{"ex_pair_functions_table: pair functions C_K, C_u of a resonance gas at rest"};
  app.add_option("--share",      shareDir,    "path to SHARE particle database");
  app.add_option("--output",     outputFile,  "output ROOT file")->default_val("pair_table.root");
  app.add_option("--weak",       weak,        "weak decays: on (default) or off")
      ->default_str("on")
      ->check(CLI::IsMember({"on", "off"}));
  app.add_option("--stat",       stat,        "parent statistics: qs (default) or boltzmann")
      ->default_str("qs")
      ->check(CLI::IsMember({"qs", "boltzmann"}));
  app.add_option("--level",      level,       "L0 (thermal parents), L1 (currents), L2 (+ chains)")
      ->default_str("L0")
      ->check(CLI::IsMember({"L0", "L1", "L2"}));
  app.add_option("--kratio",     kRatio,      "L2: kernel decays per parent / parents in the MC")
      ->default_val(4.);
  app.add_option("--kmin",       kMin,        "L2: minimum kernel decays per parent")->default_val(2e4);
  app.add_option("--kmax",       kMax,        "L2: maximum kernel decays per parent")->default_val(2e7);
  app.add_option("--mcevents",   mcEvents,    "L2: events of the MC the table is compared with")
      ->default_val(2e5);
  app.add_option("--threads",    nThreads,    "L2: threads for the convolution")->default_val(8);
  app.add_option("--norm",       norm,        "L1 thermal densities: quad (default) or mc")
      ->default_str("quad")
      ->check(CLI::IsMember({"quad", "mc"}));
  app.add_option("--samples",    nSamples,    "Integrator samples (T1.5, --norm mc)")
      ->default_val(1000000);
  app.add_option("--varsamples", nVarSamples, "samples for the T1.5 error estimate")
      ->default_val(200000);
  app.add_option("--seed",       seed,        "random seed")->default_val(0);
  CLI11_PARSE(app, argc, argv);

  const bool weakOn  = (weak == "on");
  const bool quantum = (stat == "qs");
  const bool l2      = (level == "L2");
  const bool l1      = (level == "L1") || l2; // L2 = L1 table + chain part
  auto tStart = std::chrono::steady_clock::now();

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
  ParticleDB* db = tPartDB.get();
  const double T = model->GetTemperature();
  const double V = model->GetVolume();
  const std::string dbHash = ParticleDBHash(db);
  std::cout << "HRG: T = " << T << " GeV, V = " << V << " GeV^-3, weak " << weak << ", stat "
            << stat << ", level " << level << (l1 ? ", norm " + norm : "") << ", ParticleDB hash "
            << dbHash << "\n";

  auto mass = [&](int idx) { return static_cast<double>(db->GetParticleType(idx)->GetMass()); };
  auto pdg  = [&](int idx) { return db->GetParticleType(idx)->GetPDGCode(); };
  auto name = [&](int idx) { return std::string(db->GetParticleType(idx)->GetName()); };

  // ---- Channel loop (spec Sec. 3.1) and one-step counting (Sec. 6.1) --------------------------
  ContributionSet cs = BuildContributions(db, T, quantum);
  const ContributionSet cs0 = cs; // L0 contributions (T7.1)

  // ---- L1: parent currents (spec Sec. 9.2-9.3) ---------------------------------------------------
  CurrentSet curs, cursTh;
  if (l1) {
    std::vector<double> thScale(db->GetParticleTypeCount(), 1.);
    if (norm == "mc") {
      Integrator integ(nSamples, model.get());
      integ.SetMultiplicities(db);
      for (int c = 0; c < db->GetParticleTypeCount(); c++)
        thScale[c] = db->GetParticleType(c)->GetMultiplicity() / V / cs.nquad[c];
    }
    auto t0 = std::chrono::steady_clock::now();
    curs = BuildCurrents(db, T, quantum, thScale);
    std::printf("Currents built in %.1f s\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    ToL1(cs, db, curs);
    // Thermal-only currents, quadrature normalization (T7.1).
    cursTh.cur.resize(curs.cur.size());
    for (size_t c = 0; c < curs.cur.size(); c++) {
      SpeciesCurrent& s = cursTh.cur[c];
      s.m = curs.cur[c].m; s.st = curs.cur[c].st; s.thScale = 1.;
      s.edges = curs.cur[c].edges; s.nPerPanel = curs.cur[c].nPerPanel; s.hasFeed = false;
    }
  }
  const std::vector<Contribution>& contribs = cs.contribs;
  const std::vector<ParentStat>& pstat     = cs.pstat;
  const std::vector<double>& nquad         = cs.nquad;
  const double* c1step                     = cs.c1step;
  const int nSkippedGamma                  = cs.nSkippedGamma;
  int nPoles = 0, nCont = 0;
  for (auto& co : contribs) (co.third < 0 ? nPoles : nCont)++;
  std::cout << "Channel loop: " << contribs.size() << " (pair, channel) terms: " << nPoles
            << " poles, " << nCont << " continua; " << nSkippedGamma
            << " parents with channels skipped (Gamma < 0)\n";

  // ---- Tables and output ----------------------------------------------------------------------
  TFile fout(outputFile.c_str(), "RECREATE", "Pair functions of a resonance gas at rest");
  std::vector<double> kEdges = KappaEdges();
  const int nK = kEdges.size() - 1;
  TVectorD kappaCenters(nK);
  for (int k = 0; k < nK; k++) kappaCenters[k] = 0.5 * (kEdges[k] + kEdges[k + 1]);

  long nCells = 0, nNegCK = 0, nNegCu = 0, nNegSum = 0, nNonFinite = 0, nSpacelike = 0;
  double gridIntegral[kNPairs] = {0.};
  std::array<std::vector<double>, kNPairs> contCK, contCu; // L1 continuum on the grid (T8.4)
  for (int p = 0; p < kNPairs; p++) {
    std::string tag = kPairTag[p];
    double m1 = 0., m2 = 0.;
    for (int idx = 0; idx < db->GetParticleTypeCount(); idx++) {
      if (pdg(idx) == kPairPDG[p][0]) m1 = mass(idx);
      if (pdg(idx) == kPairPDG[p][1]) m2 = mass(idx);
    }
    std::vector<double> mEdges = MassEdges(m1 + m2);
    const int nM = mEdges.size() - 1;
    contCK[p].assign(nM * nK, 0.);
    contCu[p].assign(nM * nK, 0.);

    fout.cd();
    TH2D* hCK = new TH2D(("CK_cont_" + tag).c_str(),
                         ("C_{K} continuum " + tag + ";M (GeV);#kappa (GeV)").c_str(), nM,
                         mEdges.data(), nK, kEdges.data());
    TH2D* hCu = new TH2D(("Cu_cont_" + tag).c_str(),
                         ("C_{u} continuum " + tag + ";M (GeV);#kappa (GeV)").c_str(), nM,
                         mEdges.data(), nK, kEdges.data());
    TH1D* hRho = new TH1D(("rho_cont_" + tag).c_str(),
                          ("#rho continuum " + tag + ";M (GeV);#rho (GeV^{2})").c_str(), nM,
                          mEdges.data());

    TTree* poles = new TTree(("poles_" + tag).c_str(), ("two-body poles " + tag).c_str());
    char parentName[64];
    int parentPDG, mult;
    double tMa, tAoverMa, tBR, tNa;
    std::vector<double> fa(nK), abc(nK); // L1: fa holds A_c, abc A_c + B_c
    poles->Branch("parent", parentName, "parent/C");
    poles->Branch("pdg", &parentPDG, "pdg/I");
    poles->Branch("ma", &tMa, "ma/D");
    poles->Branch("AoverMa", &tAoverMa, "AoverMa/D");
    poles->Branch("BR", &tBR, "BR/D");
    poles->Branch("mult", &mult, "mult/I");
    poles->Branch("na", &tNa, "na/D");
    poles->Branch("fa", &fa);
    if (l1) poles->Branch("abc", &abc);

    TDirectory* chDir = fout.mkdir(("channels/" + tag).c_str(), "", true);
    for (const Contribution& co : contribs) {
      if (co.pair != p) continue;
      if (co.third < 0) {
        std::snprintf(parentName, sizeof(parentName), "%s", name(co.parent).c_str());
        parentPDG = pdg(co.parent);
        tMa = co.ma; tAoverMa = co.A / co.ma; tBR = co.B; mult = co.count; tNa = co.na;
        for (int k = 0; k < nK; k++) {
          if (l1)
            co.cur->AtEnergy(std::hypot(co.ma, kappaCenters[k]), fa[k], abc[k]);
          else
            fa[k] = co.st.f(std::hypot(co.ma, kappaCenters[k]));
        }
        poles->Fill();
        if (co.ma < 2.0) gridIntegral[p] += co.count * co.B * co.na;
        continue;
      }
      chDir->cd();
      TH2F* hCKch = new TH2F(("CK_" + co.name).c_str(), ("C_{K} " + co.name).c_str(), nM,
                             mEdges.data(), nK, kEdges.data());
      TH2F* hCuch = new TH2F(("Cu_" + co.name).c_str(), ("C_{u} " + co.name).c_str(), nM,
                             mEdges.data(), nK, kEdges.data());
      TH1D* hRhoch = new TH1D(("rho_" + co.name).c_str(), ("#rho " + co.name).c_str(), nM,
                              mEdges.data());
      double lo = co.m1 + co.m2, hi = co.ma - co.m3;
      for (int im = 1; im <= nM; im++) {
        double M = hCK->GetXaxis()->GetBinCenter(im);
        if (M <= lo || M >= hi) continue;
        // rho_ch(M) from the Dalitz projection (identical to the kappa integral by T1.2).
        double rho = Rho3Dalitz(M, co.m1, co.m2, co.ma, co.m3, co.count * co.B, co.na, co.phi3);
        hRhoch->SetBinContent(im, rho);
        hRho->AddBinContent(im, rho);
        double sumK = 0.;
        for (int ik = 1; ik <= nK; ik++) {
          double kap = kappaCenters[ik - 1];
          double eK  = std::hypot(M, kap);
          CPair cc = l1 ? CKCu3Current(M, kap, co.ma, co.m3, co.D, *co.cur)
                        : CKCu3(M, kap, co.ma, co.m3, co.D, co.st);
          ++nCells;
          if (l1) nSpacelike += !(eK * (cc.CK + cc.Cu) > kap * std::abs(cc.CK));
          if (!std::isfinite(cc.CK) || !std::isfinite(cc.Cu)) ++nNonFinite;
          nNegCK += (cc.CK < 0.);
          nNegCu += (cc.Cu < 0.);
          nNegSum += (cc.CK + cc.Cu < 0.);
          hCKch->SetBinContent(im, ik, cc.CK);
          hCuch->SetBinContent(im, ik, cc.Cu);
          hCK->AddBinContent(hCK->GetBin(im, ik), cc.CK);
          hCu->AddBinContent(hCu->GetBin(im, ik), cc.Cu);
          contCK[p][(im - 1) * nK + ik - 1] += cc.CK;
          contCu[p][(im - 1) * nK + ik - 1] += cc.Cu;
          sumK += kap * kap * (cc.CK + cc.Cu) * (kEdges[ik] - kEdges[ik - 1]);
        }
        gridIntegral[p] += 16. * M_PI * M_PI * PStar(M, co.m1, co.m2) * sumK *
                           hCK->GetXaxis()->GetBinWidth(im);
      }
      hCKch->Write();
      hCuch->Write();
      hRhoch->Write();
      delete hCKch;
      delete hCuch;
      delete hRhoch;
    }
    fout.cd();
    hCK->Write();
    hCu->Write();
    hRho->Write();
    poles->Write();
    delete hCK;
    delete hCu;
    delete hRho;
    delete poles;
  }

  if (!l1) {
    // ---- T1.1: Phi3 from the mass integral and from the Dalitz area ----------------------------
    {
      double worst = 0.;
      std::string worstName;
      int nChan = 0;
      for (int a = 0; a < db->GetParticleTypeCount(); a++) {
        ParticleType* type = db->GetParticleType(a);
        if (!DecaysInEvent(type)) continue;
        for (int c = 0; c <= type->GetTable()->GetChannelCount(); c++) {
          const DecayChannel* ch = type->GetTable()->GetDecayChannel(c);
          if (!ch->Is3Particle()) continue;
          ++nChan;
          std::vector<int> d = ChannelDaughters(ch);
          int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
          for (auto& pm : perm) {
            double m1 = mass(d[pm[0]]), m2 = mass(d[pm[1]]), m3 = mass(d[pm[2]]);
            double r = RelDiff(Phi3Mass(mass(a), m1, m2, m3), Phi3Dalitz(mass(a), m1, m2, m3));
            if (r > worst) {
              worst = r;
              worstName = name(a) + " -> " + name(d[pm[0]]) + " " + name(d[pm[1]]) + " " +
                          name(d[pm[2]]);
            }
          }
        }
      }
      std::printf("T1.1  Phi3 mass integral vs Dalitz area, %d three-body channels x 6 orderings: "
                  "worst rel. diff %.2e (%s)  %s\n",
                  nChan, worst, worstName.c_str(), PassFail(worst < 1e-8));
    }

    // ---- T1.2: rho(M) from C_K + C_u vs Dalitz projection --------------------------------------
    {
      double worst = 0.;
      std::string worstName;
      int nChan = 0, nPoints = 0;
      for (int a = 0; a < db->GetParticleTypeCount(); a++) {
        ParticleType* type = db->GetParticleType(a);
        if (!DecaysInEvent(type)) continue;
        for (int c = 0; c <= type->GetTable()->GetChannelCount(); c++) {
          const DecayChannel* ch = type->GetTable()->GetDecayChannel(c);
          if (!ch->Is3Particle()) continue;
          ++nChan;
          std::vector<int> d = ChannelDaughters(ch);
          double B = EffectiveBR(type->GetTable(), c), ma = mass(a);
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              if (i == j) continue;
              double m1 = mass(d[i]), m2 = mass(d[j]), m3 = mass(d[3 - i - j]);
              double phi3 = Phi3Mass(ma, m1, m2, m3);
              double D    = B * pstat[a].g / (std::pow(2. * M_PI, 3) * phi3);
              double lo = m1 + m2, hi = ma - m3;
              for (int q = 0; q < 12; q++) {
                double M = lo + (hi - lo) * (q + 0.5) / 12.;
                double r = RelDiff(Rho3(M, m1, m2, ma, m3, D, pstat[a]),
                                   Rho3Dalitz(M, m1, m2, ma, m3, B, nquad[a], phi3));
                ++nPoints;
                if (r > worst) {
                  worst = r;
                  char buf[256];
                  std::snprintf(buf, sizeof(buf), "%s -> %s %s + %s, M = %.4f", name(a).c_str(),
                                name(d[i]).c_str(), name(d[j]).c_str(), name(d[3 - i - j]).c_str(),
                                M);
                  worstName = buf;
                }
              }
            }
          }
        }
      }
      std::printf("T1.2  rho(M) from C_K+C_u vs Dalitz projection, %d channels, %d points: "
                  "worst rel. diff %.2e (%s)  %s\n",
                  nChan, nPoints, worst, worstName.c_str(), PassFail(worst < 1e-6));
    }

    // ---- T1.3: every pole carries B n_a ----------------------------------------------------------
    {
      double worst = 0.;
      std::string worstName;
      int nChan = 0;
      for (int a = 0; a < db->GetParticleTypeCount(); a++) {
        ParticleType* type = db->GetParticleType(a);
        if (!DecaysInEvent(type)) continue;
        for (int c = 0; c <= type->GetTable()->GetChannelCount(); c++) {
          const DecayChannel* ch = type->GetTable()->GetDecayChannel(c);
          if (ch->Is3Particle()) continue;
          ++nChan;
          double B = EffectiveBR(type->GetTable(), c), ma = mass(a);
          double m1 = mass(ch->GetParticle1()), m2 = mass(ch->GetParticle2());
          double A  = PoleA(ma, m1, m2, B, pstat[a].g);
          double r  = RelDiff(PoleRhoIntegral(ma, m1, m2, A, pstat[a]), B * nquad[a]);
          if (r > worst) {
            worst = r;
            worstName = name(a) + " -> " + name(ch->GetParticle1()) + " " + name(ch->GetParticle2());
          }
        }
      }
      std::printf("T1.3  pole integral vs B n_a, %d two-body channels: worst rel. diff %.2e (%s)  %s\n",
                  nChan, worst, worstName.c_str(), PassFail(worst < 1e-6));
    }

    // ---- T1.4: int dM rho_ss'(M) = c^1step_ss' ----------------------------------------------------
    {
      bool ok = true;
      std::printf("T1.4  %-9s %16s %16s %10s   %16s %8s\n", "pair", "int dM rho", "c^1step",
                  "rel.diff", "grid (M<2 GeV)", "grid/c");
      for (int p = 0; p < kNPairs; p++) {
        double lhs = 0.;
        for (const Contribution& co : contribs) {
          if (co.pair != p) continue;
          if (co.third < 0)
            lhs += PoleRhoIntegral(co.ma, co.m1, co.m2, co.A, co.st);
          else
            lhs += IntegrateCos(
                [&](double M) { return Rho3(M, co.m1, co.m2, co.ma, co.m3, co.D, co.st); },
                co.m1 + co.m2, co.ma - co.m3, 128);
        }
        double r = RelDiff(lhs, c1step[p]);
        if (c1step[p] == 0. && lhs == 0.) r = 0.;
        ok = ok && (r < 1e-6);
        std::printf("T1.4  %-9s %16.9e %16.9e %10.2e   %16.9e %8.5f\n", kPairTag[p], lhs, c1step[p],
                    r, gridIntegral[p], c1step[p] > 0. ? gridIntegral[p] / c1step[p] : 0.);
      }
      std::printf("T1.4  int dM rho = c^1step for all 12 ordered pairs  %s\n", PassFail(ok));
    }

    // ---- T1.5: n_a from quadrature vs Integrator ---------------------------------------------------
    if (quantum) {
      Integrator integ(nSamples, model.get());
      integ.SetMultiplicities(db);
      double sumP = 0., sumP2 = 0., worst = 0.;
      std::string worstName;
      int n = 0;
      ParticleCoor coor;
      for (int a = 0; a < db->GetParticleTypeCount(); a++) {
        ParticleType* type = db->GetParticleType(a);
        // Per-sample spread of the Integrator estimate HyperCube * mean(GetIntegrand).
        double s = 0., s2 = 0.;
        for (int i = 0; i < nVarSamples; i++) {
          double v = model->GetIntegrand(type, coor);
          s += v;
          s2 += v * v;
        }
        double mean = s / nVarSamples, var = s2 / nVarSamples - mean * mean;
        double sigma = model->GetHyperCubeVolume() * std::sqrt(var / nSamples) / V;
        double pull  = (type->GetMultiplicity() / V - nquad[a]) / sigma;
        sumP += pull;
        sumP2 += pull * pull;
        ++n;
        if (std::abs(pull) > std::abs(worst)) {
          worst = pull;
          char buf[160];
          std::snprintf(buf, sizeof(buf), "%s: n_quad %.6e, n_MC %.6e, sigma %.2e", name(a).c_str(),
                        nquad[a], type->GetMultiplicity() / V, sigma);
          worstName = buf;
        }
      }
      double pMean = sumP / n, pRms = std::sqrt(sumP2 / n);
      std::printf("T1.5  (n_MC - n_quad)/sigma over %d types (Integrator %d samples, error from %d "
                  "samples): mean %.3f, rms %.3f, max |pull| %.2f (%s)  %s\n",
                  n, nSamples, nVarSamples, pMean, pRms, std::abs(worst), worstName.c_str(),
                  PassFail(std::abs(worst) < 5. && pRms > 0.7 && pRms < 1.3));
      std::printf("T1.5  note: the fmultiplicity cache in the run directory must have been made with "
                  "--samples %d\n", nSamples);
    } else {
      std::printf("T1.5  skipped (--stat boltzmann; THERMINATOR uses quantum statistics)\n");
    }

    // ---- T1.6: positivity, finiteness, kappa -> 0 limit --------------------------------------------
    {
      double worst = 0.;
      std::string worstName;
      int nSmallFinite = 0, nSmall = 0;
      for (const Contribution& co : contribs) {
        if (co.third < 0) continue;
        double lo = co.m1 + co.m2, hi = co.ma - co.m3;
        for (int q = 0; q < 5; q++) {
          double M  = lo + (hi - lo) * (q + 0.5) / 5.;
          CPair c   = CKCu3(M, 1e-6, co.ma, co.m3, co.D, co.st);
          CPair c0  = CKCu3Limit0(M, co.ma, co.m3, co.D, co.st);
          CPair c8  = CKCu3(M, 1e-8, co.ma, co.m3, co.D, co.st);
          ++nSmall;
          nSmallFinite += (std::isfinite(c8.CK) && std::isfinite(c8.Cu));
          double r  = std::max(RelDiff(c.CK, c0.CK), RelDiff(c.Cu, c0.Cu));
          if (r > worst) {
            worst = r;
            worstName = co.name + " (" + kPairTag[co.pair] + ")";
          }
        }
      }
      // C_u >= 0 follows from f decreasing; C_K + C_u is the rest-frame density (> 0). C_K alone
      // can be negative near threshold at small kappa (also in the Appendix A prototype).
      std::printf("T1.6  grid: %ld cells: %ld non-finite, %ld with C_u<0, %ld with C_K+C_u<0  %s\n",
                  nCells, nNonFinite, nNegCu, nNegSum,
                  PassFail(nNonFinite == 0 && nNegCu == 0 && nNegSum == 0));
      std::printf("T1.6  grid: %ld cells with C_K<0 (spec requires C_K>=0)  %s\n", nNegCK,
                  PassFail(nNegCK == 0));
      std::printf("T1.6  kappa = 1e-8: %d/%d evaluations finite  %s\n", nSmallFinite, nSmall,
                  PassFail(nSmallFinite == nSmall));
      std::printf("T1.6  kappa = 1e-6 vs analytic kappa -> 0 limit, %d continua x 5 masses: worst rel. "
                  "diff %.2e (%s)  %s\n",
                  nCont, worst, worstName.c_str(), PassFail(worst < 1e-6));
    }

    // ---- Appendix A reference values (omega -> pi+ pi- pi0, Boltzmann, B = g = 1) -----------------
    {
      double ma = 0.78266, m1 = 0.13957, m2 = 0.13957, m3 = 0.13498;
      ParentStat boltz{1., 0., 0.155};
      double P3a = Phi3Mass(ma, m1, m2, m3), P3b = Phi3Dalitz(ma, m1, m2, m3);
      double D = 1. / (std::pow(2. * M_PI, 3) * P3a);
      double na = NumberDensity(ma, boltz);
      std::printf("REF   Phi3_mass %.17g Phi3_dalitz %.17g\n", P3a, P3b);
      std::printf("REF   n_a %.17g\n", na);
      for (double M : {0.30, 0.40, 0.50, 0.60})
        std::printf("REF   rho M %.2f %.17g dalitz %.17g\n", M, Rho3(M, m1, m2, ma, m3, D, boltz),
                    Rho3Dalitz(M, m1, m2, ma, m3, 1., na, P3a));
      for (double k : {1e-4, 1e-2, 0.5}) {
        CPair c = CKCu3(0.4, k, ma, m3, D, boltz);
        std::printf("REF   CK_Cu M 0.4 kappa %g %.17g %.17g\n", k, c.CK, c.Cu);
      }
    }

  } else if (!l2) {
    // ---- T7.1: thermal-only currents reproduce the M1 table ----------------------------------
    {
      double worstC = 0., worstP = 0.;
      std::string nameC, nameP;
      long nCmp = 0;
      for (size_t i = 0; i < cs0.contribs.size(); i++) {
        const Contribution& c0 = cs0.contribs[i];
        const SpeciesCurrent& th = cursTh.cur[c0.parent];
        auto mm = PairMasses(db, c0.pair);
        std::vector<double> mEdges = MassEdges(mm.first + mm.second);
        if (c0.third < 0) {
          // pole: (A/m_c) A_c vs (A/m_a) f, and (A+B)_c = A_c
          double AL1 = c0.count * c0.B / Phi2(c0.ma, c0.m1, c0.m2);
          for (int k = 0; k < nK; k++) {
            double E = std::hypot(c0.ma, kappaCenters[k]), a, ab;
            th.AtEnergy(E, a, ab);
            double ref = c0.A / c0.ma * c0.st.f(E);
            double r   = std::max(RelDiff(AL1 / c0.ma * a, ref), RelDiff(AL1 / c0.ma * ab, ref));
            ++nCmp;
            if (r > worstP) { worstP = r; nameP = c0.name; }
          }
          continue;
        }
        double DL1 = c0.count * c0.B / c0.phi3;
        double lo = c0.m1 + c0.m2, hi = c0.ma - c0.m3;
        for (size_t im = 0; im + 1 < mEdges.size(); im++) {
          double M = 0.5 * (mEdges[im] + mEdges[im + 1]);
          if (M <= lo || M >= hi) continue;
          for (int k = 0; k < nK; k++) {
            CPair r0 = CKCu3(M, kappaCenters[k], c0.ma, c0.m3, c0.D, c0.st);
            CPair r1 = CKCu3Current(M, kappaCenters[k], c0.ma, c0.m3, DL1, th);
            double scale = std::max(std::abs(r0.CK), r0.CK + r0.Cu); // C_K changes sign
            double r = std::max(std::abs(r1.CK - r0.CK), std::abs(r1.Cu - r0.Cu)) / scale;
            ++nCmp;
            if (r > worstC) {
              worstC = r;
              char buf[200];
              std::snprintf(buf, sizeof(buf), "%s (%s), M = %.4f, kappa = %.4f", c0.name.c_str(),
                            kPairTag[c0.pair], M, kappaCenters[k]);
              nameC = buf;
            }
          }
        }
      }
      std::printf("T7.1  thermal-only currents vs M1 table, %ld grid values: continua worst rel. "
                  "diff %.2e (%s), poles %.2e (%s)  %s\n",
                  nCmp, worstC, nameC.c_str(), worstP, nameP.c_str(),
                  PassFail(worstC < 1e-12 && worstP < 1e-12));
    }

    // ---- T7.2: rho^L1 vs Dalitz projection with n_c; int dM rho = c^L1; poles carry B n_c ------
    {
      double worst = 0.;
      std::string worstName;
      int nPoints = 0;
      for (const Contribution& co : contribs) {
        if (co.third < 0) continue;
        double lo = co.m1 + co.m2, hi = co.ma - co.m3;
        for (int q = 0; q < 12; q++) {
          double M = lo + (hi - lo) * (q + 0.5) / 12.;
          double r = RelDiff(Rho3Current(M, co),
                             Rho3Dalitz(M, co.m1, co.m2, co.ma, co.m3, co.count * co.B, co.na,
                                        co.phi3));
          ++nPoints;
          if (r > worst) {
            worst = r;
            char buf[200];
            std::snprintf(buf, sizeof(buf), "%s (%s), M = %.4f", co.name.c_str(),
                          kPairTag[co.pair], M);
            worstName = buf;
          }
        }
      }
      std::printf("T7.2  rho(M) from C_K+C_u vs Dalitz projection with n_c, %d points: worst rel. "
                  "diff %.2e (%s)  %s\n",
                  nPoints, worst, worstName.c_str(), PassFail(worst < 1e-6));

      double worstPole = 0.;
      std::string namePole;
      for (const Contribution& co : contribs) {
        if (co.third >= 0) continue;
        double I = 16. * M_PI * M_PI * PStar(co.ma, co.m1, co.m2) * co.A / co.ma *
                   IntegratePanels([&](double k) { return PoleKappaIntegrandL1(co, k); },
                                   ThermalPanels(T));
        double r = RelDiff(I, co.count * co.B * co.na);
        if (r > worstPole) { worstPole = r; namePole = co.name; }
      }
      std::printf("T7.2  every pole carries count B n_c: worst rel. diff %.2e (%s)  %s\n",
                  worstPole, namePole.c_str(), PassFail(worstPole < 1e-6));

      bool ok = true;
      std::printf("T7.2  %-9s %16s %16s %10s   %16s %8s\n", "pair", "int dM rho", "c^L1",
                  "rel.diff", "grid (M<2 GeV)", "grid/c");
      for (int p = 0; p < kNPairs; p++) {
        // One thread per contribution of this pair.
        std::vector<const Contribution*> mine;
        for (const Contribution& co : contribs)
          if (co.pair == p) mine.push_back(&co);
        std::vector<double> parts(mine.size(), 0.);
        std::vector<std::thread> pool;
        for (size_t i = 0; i < mine.size(); i++) {
          pool.emplace_back([&, i]() {
            const Contribution& co = *mine[i];
            if (co.third < 0)
              parts[i] = 16. * M_PI * M_PI * PStar(co.ma, co.m1, co.m2) * co.A / co.ma *
                         IntegratePanels([&](double k) { return PoleKappaIntegrandL1(co, k); },
                                         ThermalPanels(T));
            else
              parts[i] = IntegrateCos([&](double M) { return Rho3Current(M, co); },
                                      co.m1 + co.m2, co.ma - co.m3, 64);
          });
        }
        for (auto& th : pool) th.join();
        double lhs = 0.;
        for (double x : parts) lhs += x;
        double r = RelDiff(lhs, c1step[p]);
        if (c1step[p] == 0. && lhs == 0.) r = 0.;
        ok = ok && (r < 1e-6);
        std::printf("T7.2  %-9s %16.9e %16.9e %10.2e   %16.9e %8.5f\n", kPairTag[p], lhs,
                    c1step[p], r, gridIntegral[p], c1step[p] > 0. ? gridIntegral[p] / c1step[p] : 0.);
      }
      std::printf("T7.2  int dM rho^L1 = c^L1 for all 12 ordered pairs  %s\n", PassFail(ok));
      std::printf("T7.2  %-9s %16s %16s %8s\n", "pair", "c^L0 (quad)", "c^L1", "L1/L0");
      for (int p = 0; p < kNPairs; p++)
        std::printf("T7.2  %-9s %16.9e %16.9e %8.4f\n", kPairTag[p], cs0.c1step[p], c1step[p],
                    cs0.c1step[p] > 0. ? c1step[p] / cs0.c1step[p] : 0.);
    }

    // ---- T7.3: positivity, timelike, finite as kappa -> 0 --------------------------------------
    {
      std::printf("T7.3  grid: %ld cells: %ld non-finite, %ld with C_K+C_u<0, %ld spacelike "
                  "(eps_K(C_K+C_u) <= kappa|C_K|)  %s\n",
                  nCells, nNonFinite, nNegSum, nSpacelike,
                  PassFail(nNonFinite == 0 && nNegSum == 0 && nSpacelike == 0));
      std::printf("T7.3  grid: %ld cells with C_u<0, %ld with C_K<0 (allowed, as at L0)  %s\n",
                  nNegCu, nNegCK, PassFail(nNegCu == 0));
      int nSmall = 0, nFinite = 0;
      double worst = 0.;
      std::string worstName;
      for (const Contribution& co : contribs) {
        if (co.third < 0) continue;
        double lo = co.m1 + co.m2, hi = co.ma - co.m3;
        for (int q = 0; q < 5; q++) {
          double M = lo + (hi - lo) * (q + 0.5) / 5.;
          CPair c8 = CKCu3Current(M, 1e-8, co.ma, co.m3, co.D, *co.cur);
          CPair c6 = CKCu3Current(M, 1e-6, co.ma, co.m3, co.D, *co.cur);
          CPair c4 = CKCu3Current(M, 1e-4, co.ma, co.m3, co.D, *co.cur);
          ++nSmall;
          nFinite += std::isfinite(c8.CK) && std::isfinite(c8.Cu) && c8.CK + c8.Cu > 0.;
          // the functions are even in kappa: change O(kappa^2) between 1e-6 and 1e-4
          double scale = std::max(std::abs(c4.CK), c4.CK + c4.Cu);
          double r = std::max(std::abs(c6.CK - c4.CK), std::abs(c6.Cu - c4.Cu)) / scale;
          if (r > worst) { worst = r; worstName = co.name + " (" + kPairTag[co.pair] + ")"; }
        }
      }
      std::printf("T7.3  kappa = 1e-8: %d/%d evaluations finite with C_K+C_u>0; kappa = 1e-6 vs "
                  "1e-4: worst rel. change %.2e (%s)  %s\n",
                  nFinite, nSmall, worst, worstName.c_str(),
                  PassFail(nFinite == nSmall && worst < 1e-4));
    }
  }

  if (l2) {
    // ---- L2: chain kernels and convolution (spec Sec. 9.4) -------------------------------------
    const ChainCounting cnt = CountChains(db);
    std::vector<std::vector<double>> mEdgesP(kNPairs);
    std::vector<std::pair<double, double>> pairMasses(kNPairs);
    for (int p = 0; p < kNPairs; p++) {
      pairMasses[p] = PairMasses(db, p);
      mEdgesP[p] = MassEdges(pairMasses[p].first + pairMasses[p].second);
    }
    std::vector<double> kcEdges;
    for (int i = 0; i <= 150; i++) kcEdges.push_back(0.01 * i); // k_c: 0 - 1.5 GeV
    std::vector<double> kTable(kappaCenters.GetMatrixArray(), kappaCenters.GetMatrixArray() + nK);
    ChainTable chain;
    chain.Init(mEdgesP, kTable, MCKappaEdges());

    fout.cd();
    TTree* tk = new TTree("kernels", "chain kernels per parent");
    char kName[64];
    long kN;
    double kNmc, kSec, kChain[kNPairs], kAll[kNPairs], kAllExp[kNPairs], kOver[kNPairs];
    tk->Branch("parent", kName, "parent/C");
    tk->Branch("N", &kN, "N/L");
    tk->Branch("Nmc", &kNmc, "Nmc/D");
    tk->Branch("seconds", &kSec, "seconds/D");
    tk->Branch("chain", kChain, Form("chain[%d]/D", kNPairs));
    tk->Branch("all", kAll, Form("all[%d]/D", kNPairs));
    tk->Branch("allExpected", kAllExp, Form("allExpected[%d]/D", kNPairs));
    tk->Branch("overflow", kOver, Form("overflow[%d]/D", kNPairs));

    std::vector<double> pulls81;
    double worst81 = 0.;
    std::string worst81Name;
    double chainSum[kNPairs] = {0.}, chainVar[kNPairs] = {0.};
    long totalDecays = 0;
    int nParents = 0;
    auto tChain = std::chrono::steady_clock::now();
    for (int c : curs.order) {
      ParticleType* type = db->GetParticleType(c);
      if (!DecaysInEvent(type)) continue;
      double exp = 0.;
      for (int p = 0; p < kNPairs; p++) exp += cnt.allLCA[c][p] - cnt.oneStep[c][p];
      if (exp < 1e-12) continue;
      const SpeciesCurrent& sc = curs.cur[c];
      double nmc = mcEvents * V * sc.nCount;
      long N = static_cast<long>(std::min(std::max(kRatio * nmc, kMin), kMax));
      auto t0 = std::chrono::steady_clock::now();
      ChainKernel kern = GenerateKernel(db, c, N, mEdgesP, kcEdges);
      ConvolveChain(kern, sc, T, pairMasses, chain, nThreads);
      kSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      ++nParents;
      totalDecays += N;
      for (int p = 0; p < kNPairs; p++) {
        double mAll = kern.sAll[p] / N, mCh = kern.sChain[p] / N;
        double vAll = std::max(kern.s2All[p] / N - mAll * mAll, 0.) / N;
        double vCh  = std::max(kern.s2Chain[p] / N - mCh * mCh, 0.) / N;
        chainSum[p] += sc.nCount * mCh;
        chainVar[p] += sc.nCount * sc.nCount * vCh;
        kChain[p] = mCh; kAll[p] = mAll; kAllExp[p] = cnt.allLCA[c][p];
        kOver[p] = kern.mOverflow[p] / N;
        if (cnt.allLCA[c][p] <= 0.) continue;
        double sig  = std::sqrt(vAll);
        double pull = (sig > 0.) ? (mAll - cnt.allLCA[c][p]) / sig : 0.;
        pulls81.push_back(pull);
        if (std::abs(pull) > std::abs(worst81)) {
          worst81 = pull;
          char buf[200];
          std::snprintf(buf, sizeof(buf), "%s %s: %.6e vs %.6e", sc.name.c_str(), kPairTag[p],
                        mAll, cnt.allLCA[c][p]);
          worst81Name = buf;
        }
      }
      std::snprintf(kName, sizeof(kName), "%s", sc.name.c_str());
      kN = N;
      kNmc = nmc;
      tk->Fill();
      std::printf("KERNEL %-10s N %9ld (%.2f x MC) chain pairs/decay: pi+pi- %.4e pi+pi+ %.4e "
                  "K+K- %.4e  [%.1f s]%s\n",
                  sc.name.c_str(), N, N / nmc, kChain[1], kChain[0], kChain[5], kSec,
                  kern.nNonFinite ? Form("  (%ld non-finite final particles skipped)",
                                         kern.nNonFinite)
                                  : "");
    }
    std::printf("Chain kernels: %d parents, %ld decays in total, %.1f s\n", nParents, totalDecays,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - tChain).count());

    // ---- T8.0: master formula with delta kernels reproduces the L1 functions -------------------
    {
      double worstC = 0., worstP = 0.;
      std::string nameC, nameP;
      std::map<int, CurrentPrimitives> prims;
      for (const Contribution& co : contribs) {
        if (!prims.count(co.parent)) prims[co.parent].Build(*co.cur, T);
        const CurrentPrimitives& prim = prims[co.parent];
        const double mc = co.ma;
        if (co.third < 0) {
          // pole: h = count B delta(M - m_c) delta(k_c)
          double ps = PStar(mc, co.m1, co.m2);
          for (double kap : {0.01, 0.3, 1.5}) {
            double IA, IAB;
            ChainT(mc, kap, mc, 0., *co.cur, prim, T, IA, IAB);
            double ck = 1. / (4. * M_PI * ps) * 0.5 * co.count * co.B * IA;
            double cku = 1. / (4. * M_PI * ps) * 0.5 * co.count * co.B * IAB;
            double a, ab;
            co.cur->AtEnergy(std::hypot(mc, kap), a, ab);
            double r = std::max(RelDiff(ck, co.A / mc * a), RelDiff(cku, co.A / mc * ab));
            if (r > worstP) { worstP = r; nameP = co.name; }
          }
          continue;
        }
        // three-body: h = count B M Phi2(M; m1, m2) Phi2(m_c; M, m3)/Phi3 delta(k_c - p3)
        double lo = co.m1 + co.m2, hi = co.ma - co.m3;
        for (int q = 0; q < 3; q++) {
          double M  = lo + (hi - lo) * (q + 0.5) / 3.;
          double ps = PStar(M, co.m1, co.m2), kc = PStar(mc, M, co.m3);
          double h  = co.count * co.B * M * Phi2(M, co.m1, co.m2) * Phi2(mc, M, co.m3) / co.phi3;
          double EP = mc / M * std::hypot(M, kc), pP = mc / M * kc;
          for (double kap : {0.01, 0.3, 1.5}) {
            double IA, IAB;
            ChainT(M, kap, EP, pP, *co.cur, prim, T, IA, IAB);
            double pref = 1. / (4. * M_PI * ps) * mc * mc / (2. * M * M) * h;
            CPair ref = CKCu3Current(M, kap, co.ma, co.m3, co.D, *co.cur, 256);
            double scale = std::max(std::abs(ref.CK), ref.CK + ref.Cu);
            double r = std::max(std::abs(pref * IA - ref.CK), std::abs(pref * IAB - ref.CK - ref.Cu)) /
                       scale;
            if (r > worstC) {
              worstC = r;
              char buf[200];
              std::snprintf(buf, sizeof(buf), "%s, M = %.4f, kappa = %.2f", co.name.c_str(), M, kap);
              nameC = buf;
            }
          }
        }
      }
      std::printf("T8.0  master formula with delta kernels vs L1 functions: continua worst rel. diff "
                  "%.2e (%s), poles %.2e (%s)  %s\n",
                  worstC, nameC.c_str(), worstP, nameP.c_str(),
                  PassFail(worstC < 1e-6 && worstP < 1e-6));
      std::printf("T8.0  convolution normalization, sum_kappa rho cells = n_c pairs/N per M row: "
                  "worst rel. diff %.2e (%s)  %s\n",
                  chain.worstNorm, chain.worstNormName.c_str(), PassFail(chain.worstNorm < 1e-6));
    }

    // ---- T8.1: all-LCA pairs per decay vs counting -----------------------------------------------
    {
      double sum = 0., sum2 = 0.;
      for (double x : pulls81) { sum += x; sum2 += x * x; }
      double n = pulls81.size(), mean = sum / n, rms = std::sqrt(sum2 / n);
      std::printf("T8.1  all-LCA pairs per decay vs sum_ch B sum_{i!=j} nbar nbar, %.0f (parent, "
                  "pair) terms: pull mean %.3f, rms %.3f, max |pull| %.2f (%s)  %s\n",
                  n, mean, rms, std::abs(worst81), worst81Name.c_str(),
                  PassFail(std::abs(worst81) < 5. && rms > 0.7 && rms < 1.3));
    }

    // ---- T8.2: int dM rho^L2 = c^L2 --------------------------------------------------------------
    {
      bool ok = true;
      std::printf("T8.2  %-9s %16s %16s %16s %16s %8s %8s\n", "pair", "c^L1", "chain", "c^L2",
                  "L1+chain", "pull", "L2/L1");
      for (int p = 0; p < kNPairs; p++) {
        double cL2 = 0.;
        for (int c = 0; c < db->GetParticleTypeCount(); c++)
          cL2 += curs.cur[c].nCount * cnt.allLCA[c][p];
        double lhs = c1step[p] + chainSum[p], sig = std::sqrt(chainVar[p]);
        double pull = sig > 0. ? (lhs - cL2) / sig : (std::abs(lhs - cL2) <= 1e-12 * cL2 ? 0. : 99.);
        ok = ok && std::abs(pull) < 4.;
        std::printf("T8.2  %-9s %16.9e %16.9e %16.9e %16.9e %+8.2f %8.3f\n", kPairTag[p], c1step[p],
                    chainSum[p], cL2, lhs, pull, c1step[p] > 0. ? cL2 / c1step[p] : 0.);
      }
      std::printf("T8.2  int dM rho^L2 = c^L2 within the kernel statistics  %s\n", PassFail(ok));
    }

    // ---- T8.4: L2 table positive and timelike on the grid ----------------------------------------
    {
      long n = 0, nNeg = 0, nSpace = 0, nBad = 0;
      for (int p = 0; p < kNPairs; p++) {
        const int nM = mEdgesP[p].size() - 1;
        for (int im = 0; im < nM; im++) {
          double M = 0.5 * (mEdgesP[p][im] + mEdgesP[p][im + 1]);
          for (int ik = 0; ik < nK; ik++) {
            double ck  = contCK[p][im * nK + ik] + chain.CK[p][im * nK + ik];
            double sum = contCK[p][im * nK + ik] + contCu[p][im * nK + ik] +
                         chain.CKCu[p][im * nK + ik];
            if (ck == 0. && sum == 0.) continue;
            ++n;
            double kap = kTable[ik], eK = std::hypot(M, kap);
            nBad += !std::isfinite(ck) || !std::isfinite(sum);
            nNeg += !(sum > 0.);
            nSpace += !(eK * sum > kap * std::abs(ck));
          }
        }
      }
      std::printf("T8.4  L2 continuum + chain grid: %ld cells: %ld non-finite, %ld with C_K+C_u<=0, "
                  "%ld spacelike  %s\n",
                  n, nBad, nNeg, nSpace, PassFail(nBad == 0 && nNeg == 0 && nSpace == 0));
    }

    // ---- Output -----------------------------------------------------------------------------------
    fout.cd();
    tk->Write();
    std::vector<double> kmcEdges = MCKappaEdges();
    const int nKm = kmcEdges.size() - 1;
    for (int p = 0; p < kNPairs; p++) {
      std::string tag = kPairTag[p];
      const std::vector<double>& me = mEdgesP[p];
      const int nM = me.size() - 1;
      TH2D hCK(("CKbar_chain_" + tag).c_str(), ("#bar{C}_{K} chain " + tag + ";M;#kappa").c_str(),
               nM, me.data(), nK, kEdges.data());
      TH2D hCu(("Cubar_chain_" + tag).c_str(), ("#bar{C}_{u} chain " + tag + ";M;#kappa").c_str(),
               nM, me.data(), nK, kEdges.data());
      TH1D hRho(("rho_chain_" + tag).c_str(), ("#rho chain " + tag + ";M (GeV)").c_str(), nM,
                me.data());
      TH1D hRhoV(("rho_chain_var_" + tag).c_str(), ("variance of rho_chain " + tag).c_str(), nM,
                 me.data());
      TH2D hMK(("chainMK_" + tag).c_str(), ("chain density " + tag + ";M;#kappa").c_str(), nM,
               me.data(), nKm, kmcEdges.data());
      TH2D hMKV(("chainMK_var_" + tag).c_str(), ("variance of chainMK " + tag).c_str(), nM,
                me.data(), nKm, kmcEdges.data());
      for (int im = 0; im < nM; im++) {
        for (int ik = 0; ik < nK; ik++) {
          hCK.SetBinContent(im + 1, ik + 1, chain.CK[p][im * nK + ik]);
          hCu.SetBinContent(im + 1, ik + 1, chain.CKCu[p][im * nK + ik] - chain.CK[p][im * nK + ik]);
        }
        for (int b = 0; b < nKm; b++) {
          hMK.SetBinContent(im + 1, b + 1, chain.mk[p][im * nKm + b]);
          hMKV.SetBinContent(im + 1, b + 1, chain.mkVar[p][im * nKm + b]);
        }
      }
      for (int im = 0; im <= nM; im++) {
        hRho.SetBinContent(im + 1, chain.rho[p][im]);
        hRhoV.SetBinContent(im + 1, chain.rhoVar[p][im]);
      }
      hCK.Write(); hCu.Write(); hRho.Write(); hRhoV.Write(); hMK.Write(); hMKV.Write();
    }
    TParameter<double>("kratio", kRatio).Write();
    TParameter<double>("mcevents", mcEvents).Write();
  }

  // ---- Report: dropped channels and BR sums -------------------------------------------------------
  std::string dropped, brBelow;
  for (const TString& s : db->GetDroppedChannels()) dropped += std::string(s.Data()) + "\n";
  int nBelow = 0, nAbove = 0;
  std::printf("REPORT dropped channels (%zu):\n%s", db->GetDroppedChannels().size(),
              dropped.c_str());
  for (int a = 0; a < db->GetParticleTypeCount(); a++) {
    ParticleType* type = db->GetParticleType(a);
    if (!DecaysInEvent(type)) continue;
    double sum = 0.;
    for (int c = 0; c <= type->GetTable()->GetChannelCount(); c++)
      sum += type->GetTable()->GetDecayChannel(c)->GetBranchingRatio();
    if (std::abs(sum - 1.) > 1e-4) {
      char buf[128];
      std::snprintf(buf, sizeof(buf), "%s %.6f\n", name(a).c_str(), sum);
      brBelow += buf;
      (sum < 1. ? nBelow : nAbove)++;
    }
  }
  std::printf("REPORT parents with BR sum != 1 (|sum-1| > 1e-4): %d below, %d above:\n%s", nBelow,
              nAbove, brBelow.c_str());

  fout.cd();
  TParameter<double>("T", T).Write();
  TParameter<double>("V", V).Write();
  TNamed("weak", weak.c_str()).Write();
  TNamed("stat", stat.c_str()).Write();
  TNamed("level", level.c_str()).Write();
  if (l1) TNamed("norm", norm.c_str()).Write();
  TNamed("dbhash", dbHash.c_str()).Write();
  TNamed("dropped_channels", dropped.c_str()).Write();
  TNamed("br_sum_not_1", brBelow.c_str()).Write();
  kappaCenters.Write("kappa_centers");
  fout.Close();

  double secs =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();
  std::printf("Done in %.1f s, output %s\n", secs, outputFile.c_str());
  return 0;
}
