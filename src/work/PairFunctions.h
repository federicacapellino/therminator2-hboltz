// PairFunctions.h — shared helpers for the pair-function table and MC programs
//
// See README_pair_functions.md. Gas at rest, constant matrix elements, zero widths, GeV units.
//
// The --weak switch is implemented by setting Gamma < 0 on the
// weakly decaying species in the ParticleDB: Event::DecayParticles() never decays a particle
// with Gamma < 0, and the table code skips such parents. The decay tables are left untouched.

#ifndef _TH2_PAIR_FUNCTIONS_H_
#define _TH2_PAIR_FUNCTIONS_H_

#include "Crc32.h"
#include "Particle.h"
#include "ParticleDB.h"
#include "ParticleDecayer.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

// K0, Lambda, Sigma+, Sigma-, Xi0, Xi-, Omega- (antiparticles: negative codes).
static constexpr int kWeakPDG[] = {311, 3122, 3222, 3112, 3322, 3312, 3334};
static constexpr int kNWeakPDG  = sizeof(kWeakPDG) / sizeof(kWeakPDG[0]);
// Number of ParticleDB entries the switch must touch (species and antiparticles).
static constexpr int kNWeakTypes = 2 * kNWeakPDG;

inline bool IsWeakSpecies(int pdg) {
  for (int code : kWeakPDG)
    if (std::abs(pdg) == code) return true;
  return false;
}

// With weakOn == false, makes every weakly decaying species stable. Returns the number of
// ParticleDB entries changed (0 for weakOn == true, kNWeakTypes otherwise).
inline int SetWeakDecays(ParticleDB* aDB, bool weakOn) {
  if (weakOn) return 0;
  int nChanged = 0;
  for (int i = 0; i < aDB->GetParticleTypeCount(); i++) {
    ParticleType* type = aDB->GetParticleType(i);
    if (IsWeakSpecies(type->GetPDGCode())) {
      type->SetGamma(-1.0);
      ++nChanged;
    }
  }
  return nChanged;
}

// ---------------------------------------------------------------------------------------------
// Measured species and ordered pairs (spec Sec. 2)

static constexpr int kNPairs = 12;
static constexpr int kPairPDG[kNPairs][2] = {
    {211, 211},   {211, -211},   {-211, 211},   {-211, -211},  {321, 321},     {321, -321},
    {-321, 321},  {-321, -321},  {2212, 2212},  {2212, -2212}, {-2212, 2212},  {-2212, -2212}};
static const char* const kPairTag[kNPairs] = {"pippip", "pippim", "pimpip", "pimpim",
                                              "KpKp",   "KpKm",   "KmKp",   "KmKm",
                                              "pp",     "ppbar",  "pbarp",  "pbarpbar"};

// Index of the ordered pair (pdg1, pdg2) in kPairPDG, or -1.
inline int PairIndex(int pdg1, int pdg2) {
  for (int i = 0; i < kNPairs; i++)
    if (kPairPDG[i][0] == pdg1 && kPairPDG[i][1] == pdg2) return i;
  return -1;
}

// CRC32 of the parsed database: types, masses, widths, spins, PDG codes, channels and BRs.
inline std::string ParticleDBHash(ParticleDB* aDB) {
  std::ostringstream oss;
  oss.precision(9);
  for (int i = 0; i < aDB->GetParticleTypeCount(); i++) {
    ParticleType* type = aDB->GetParticleType(i);
    oss << type->GetName() << ' ' << type->GetMass() << ' ' << type->GetGamma() << ' '
        << type->GetSpin() << ' ' << type->GetPDGCode() << ';';
    DecayTable* table = type->GetTable();
    for (int c = 0; table && c <= table->GetChannelCount(); c++) {
      const DecayChannel* ch = table->GetDecayChannel(c);
      oss << ch->GetParticle1() << ',' << ch->GetParticle2() << ',' << ch->GetParticle3() << ','
          << ch->GetBranchingRatio() << ';';
    }
  }
  std::string str = oss.str();
  Crc32 hash;
  hash.Update(str.data(), str.size());
  hash.Finish();
  return hash.GetValueHex();
}

// ---------------------------------------------------------------------------------------------
// Quadrature

struct GaussLegendre {
  std::vector<double> x, w; // nodes and weights on [-1, 1]
};

// n-point Gauss-Legendre rule (Newton iteration on P_n), cached per n.
inline const GaussLegendre& GLRule(int n) {
  static std::map<int, GaussLegendre> cache;
  auto it = cache.find(n);
  if (it != cache.end()) return it->second;
  GaussLegendre& rule = cache[n];
  rule.x.resize(n);
  rule.w.resize(n);
  for (int i = 0; i < (n + 1) / 2; i++) {
    double z = std::cos(M_PI * (i + 0.75) / (n + 0.5));
    double dp = 0.;
    for (int iter = 0; iter < 100; iter++) {
      double p0 = 1., p1 = z;
      for (int k = 2; k <= n; k++) {
        double p2 = ((2. * k - 1.) * z * p1 - (k - 1.) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (z * p1 - p0) / (z * z - 1.);
      double dz = p1 / dp;
      z -= dz;
      if (std::abs(dz) < 1e-16) break;
    }
    rule.x[i] = -z;
    rule.x[n - 1 - i] = z;
    rule.w[i] = rule.w[n - 1 - i] = 2. / ((1. - z * z) * dp * dp);
  }
  return rule;
}

// Composite Gauss-Legendre: integral of f over consecutive panels [edges[k], edges[k+1]].
template <class F>
double IntegratePanels(F&& f, const std::vector<double>& edges, int nPerPanel = 16) {
  const GaussLegendre& gl = GLRule(nPerPanel);
  double sum = 0.;
  for (size_t k = 0; k + 1 < edges.size(); k++) {
    double h = 0.5 * (edges[k + 1] - edges[k]), c = 0.5 * (edges[k + 1] + edges[k]);
    for (int i = 0; i < nPerPanel; i++) sum += gl.w[i] * h * f(c + h * gl.x[i]);
  }
  return sum;
}

// Integral of f over [lo, hi] with x = lo + (hi - lo)(1 - cos th)/2: removes square-root
// behaviour at both endpoints (two-body thresholds), giving exponential convergence.
template <class F> double IntegrateCos(F&& f, double lo, double hi, int n = 256) {
  const GaussLegendre& gl = GLRule(n);
  double sum = 0.;
  for (int i = 0; i < n; i++) {
    double th = 0.5 * M_PI * (gl.x[i] + 1.);
    double x  = lo + 0.5 * (hi - lo) * (1. - std::cos(th));
    sum += gl.w[i] * f(x) * 0.5 * (hi - lo) * std::sin(th);
  }
  return 0.5 * M_PI * sum;
}

// Momentum-like integration panels on [0, ~100 T]: fine near 0, wider in the tail.
inline std::vector<double> ThermalPanels(double T) {
  std::vector<double> edges;
  for (double x = 0.; x < 10.; x += 0.5) edges.push_back(x * T);
  for (double x = 10.; x < 40.; x += 1.) edges.push_back(x * T);
  for (double x = 40.; x <= 100. + 1e-9; x += 5.) edges.push_back(x * T);
  return edges;
}

// ---------------------------------------------------------------------------------------------
// Kinematics and phase space (LCS Sec. 2)

inline double Lambda(double x, double y, double z) {
  return x * x + y * y + z * z - 2. * x * y - 2. * y * z - 2. * z * x;
}

// Two-body momentum in the rest frame of M.
inline double PStar(double M, double m1, double m2) {
  double l = Lambda(M * M, m1 * m1, m2 * m2);
  return (l > 0.) ? std::sqrt(l) / (2. * M) : 0.;
}

inline double Phi2(double M, double m1, double m2) { return 4. * M_PI * PStar(M, m1, m2) / M; }

// Phi3(ma) = int M dM Phi2(ma; M, m3) Phi2(M; m1, m2).
inline double Phi3Mass(double ma, double m1, double m2, double m3, int n = 256) {
  return IntegrateCos([&](double M) { return M * Phi2(ma, M, m3) * Phi2(M, m1, m2); }, m1 + m2,
                      ma - m3, n);
}

// Phi3(ma) = 2 pi^2 A_D / ma^2, A_D = Dalitz area in (s12, s23).
inline double Phi3Dalitz(double ma, double m1, double m2, double m3, int n = 256) {
  auto width = [&](double s12) { // (s23)_max - (s23)_min = 4 p2 p3 in the 12 rest frame
    double r  = std::sqrt(s12);
    double E2 = (s12 - m1 * m1 + m2 * m2) / (2. * r);
    double E3 = (ma * ma - s12 - m3 * m3) / (2. * r);
    return 4. * std::sqrt(std::max(E2 * E2 - m2 * m2, 0.)) *
           std::sqrt(std::max(E3 * E3 - m3 * m3, 0.));
  };
  double AD = IntegrateCos(width, (m1 + m2) * (m1 + m2), (ma - m3) * (ma - m3), n);
  return 2. * M_PI * M_PI * AD / (ma * ma);
}

// ---------------------------------------------------------------------------------------------
// Parent statistics: f(eps) = 1/(exp(eps/T) + eta), eta = -1 (BE), +1 (FD), 0 (Boltzmann).

struct ParentStat {
  double g;   // degeneracy 2J+1
  double eta;
  double T;
  double f(double eps) const {
    return (eta == 0.) ? std::exp(-eps / T) : 1. / (std::exp(eps / T) + eta);
  }
  double df(double eps) const { // d f / d eps
    if (eta == 0.) return -std::exp(-eps / T) / T;
    double e = std::exp(eps / T), d = e + eta;
    return -e / (T * d * d);
  }
};

// Same statistics as Model_HRG::GetIntegrand: BE for integer spin, FD otherwise.
inline ParentStat MakeParentStat(const ParticleType* type, double T, bool quantum) {
  double spin = type->GetSpin();
  double eta  = ((spin - static_cast<int>(spin)) < 0.01) ? -1. : +1.;
  return {2. * spin + 1., quantum ? eta : 0., T};
}

// n_a = g/(2 pi^2) int k^2 f(sqrt(m^2 + k^2)) dk.
inline double NumberDensity(double m, const ParentStat& st) {
  return st.g / (2. * M_PI * M_PI) *
         IntegratePanels([&](double k) { return k * k * st.f(std::hypot(m, k)); },
                         ThermalPanels(st.T));
}

// ---------------------------------------------------------------------------------------------
// Pair functions (LCS Eqs. 13-18, spec Sec. 3.1). Gas at rest: eps_K = E_K, kappa = |K|.

struct CPair {
  double CK, Cu;
};

// Three-body continuum, a -> 1 2 3, at pair mass M and kappa; D = B g/((2 pi)^3 Phi3).
// The overall g of f is in D, so st.g is not used here.
inline CPair CKCu3(double M, double kappa, double ma, double m3, double D, const ParentStat& st,
                   int nt = 64) {
  const GaussLegendre& gl = GLRule(nt);
  double eK = std::hypot(M, kappa);
  double EP = (ma * ma + M * M - m3 * m3) / (2. * M);
  double p3 = PStar(M, ma, m3);
  double I0 = 0., I1 = 0.;
  for (int i = 0; i < nt; i++) {
    double fP = st.f((eK * EP + kappa * p3 * gl.x[i]) / M);
    I0 += gl.w[i] * fP;
    I1 += gl.w[i] * gl.x[i] * fP;
  }
  double CK = 2. * M_PI * D * p3 / (kappa * M * M) * (kappa * EP * I0 + eK * p3 * I1);
  double Cu = -2. * M_PI * D * p3 * p3 / (kappa * eK) * I1;
  return {CK, Cu};
}

// Analytic kappa -> 0 limit of CKCu3.
inline CPair CKCu3Limit0(double M, double ma, double m3, double D, const ParentStat& st) {
  double EP = (ma * ma + M * M - m3 * m3) / (2. * M);
  double p3 = PStar(M, ma, m3);
  double f0 = st.f(EP), df0 = st.df(EP);
  return {2. * M_PI * D * p3 / (M * M) * (2. * EP * f0 + 2. / 3. * p3 * p3 * df0),
          -4. * M_PI / 3. * D * p3 * p3 * p3 * df0 / (M * M)};
}

// rho(M) = 16 pi^2 p*(M; m1, m2) int kappa^2 (CK + Cu) dkappa, three-body channel.
inline double Rho3(double M, double m1, double m2, double ma, double m3, double D,
                   const ParentStat& st) {
  return 16. * M_PI * M_PI * PStar(M, m1, m2) *
         IntegratePanels(
             [&](double k) {
               CPair c = CKCu3(M, k, ma, m3, D, st);
               return k * k * (c.CK + c.Cu);
             },
             ThermalPanels(st.T));
}

// Dalitz projection (IM Eq. 9): B n_a M Phi2(M; m1, m2) Phi2(ma; M, m3) / Phi3(ma).
inline double Rho3Dalitz(double M, double m1, double m2, double ma, double m3, double B,
                         double na, double phi3) {
  return B * na * M * Phi2(M, m1, m2) * Phi2(ma, M, m3) / phi3;
}

// Two-body pole: C_K = (A/ma) f(eps_K) delta(M - ma), A = B/Phi2(ma; m1, m2) g/(2 pi)^3.
inline double PoleA(double ma, double m1, double m2, double B, double g) {
  return B / Phi2(ma, m1, m2) * g / std::pow(2. * M_PI, 3);
}

// int dM rho of a pole: 16 pi^2 p*(ma) int kappa^2 (A/ma) f(sqrt(ma^2 + kappa^2)) dkappa.
inline double PoleRhoIntegral(double ma, double m1, double m2, double A, const ParentStat& st) {
  return 16. * M_PI * M_PI * PStar(ma, m1, m2) * A / ma *
         IntegratePanels([&](double k) { return k * k * st.f(std::hypot(ma, k)); },
                         ThermalPanels(st.T));
}

// ---------------------------------------------------------------------------------------------
// Channel loop (spec Sec. 3.1) over the parsed ParticleDB

// Parents that decay in Event::DecayParticles(): Gamma >= 0 and at least one channel.
inline bool DecaysInEvent(const ParticleType* type) {
  return type->GetGamma() >= 0. && type->GetTable() && type->GetTable()->GetChannelCount() >= 0;
}

inline std::vector<int> ChannelDaughters(const DecayChannel* ch) {
  std::vector<int> d = {ch->GetParticle1(), ch->GetParticle2()};
  if (ch->Is3Particle()) d.push_back(ch->GetParticle3());
  return d;
}

// Branching ratio of channel c as ParticleDecayer realizes it: DecayTable::ChooseDecayChannelOrNot
// draws u in [0,1] against the cumulative BRs, so where they sum to more than 1 (the K*(1430)s,
// 1.004) the channels beyond 1 are cut. Equal to GetBranchingRatio() whenever the sum is <= 1.
inline double EffectiveBR(const DecayTable* table, int c) {
  double before = 0.;
  for (int i = 0; i < c; i++) before += table->GetDecayChannel(i)->GetBranchingRatio();
  double B = table->GetDecayChannel(c)->GetBranchingRatio();
  if (before + B <= 1.) return B;
  return std::max(0., 1. - std::min(before, 1.));
}

struct SpeciesCurrent;

// One (ordered pair, parent, channel, third daughter) term of the channel loop.
// `count` is the number of ordered daughter slots (i, j) giving the same term.
struct Contribution {
  int pair, parent, channel, third; // third = -1 for two-body channels
  int count;
  double B;                         // branching ratio as stored after parsing
  double ma, m1, m2, m3;
  ParentStat st;
  double na;                        // parent density from quadrature
  double phi3;                      // three-body only
  double D;                         // three-body: count B g / ((2 pi)^3 Phi3)
  double A;                         // two-body: count B g / ((2 pi)^3 Phi2)
  std::string name;                 // <parent>_<d1>_<d2>[_<d3>]
  // L1 (spec Sec. 9.3, set by ToL1): parent current instead of the thermal f; then na = n_c,
  // D = count B / Phi3 and A = count B / Phi2 (the g/(2 pi)^3 is inside the current).
  const SpeciesCurrent* cur = nullptr;
};

struct ContributionSet {
  std::vector<Contribution> contribs;
  std::vector<ParentStat> pstat; // per ParticleDB index
  std::vector<double> nquad;     // n_a from quadrature, per ParticleDB index
  double c1step[kNPairs];        // one-step counting, sum_a n_a sum_ch B n_s (n_s' - delta)
                                 // (after ToL1: c^L1 with n_a -> n_c)
  int nSkippedGamma;             // parents with channels but Gamma < 0 (--weak off)
};

inline ContributionSet BuildContributions(ParticleDB* db, double T, bool quantum) {
  ContributionSet cs;
  const int nTypes = db->GetParticleTypeCount();
  auto mass = [&](int idx) { return static_cast<double>(db->GetParticleType(idx)->GetMass()); };
  auto pdg  = [&](int idx) { return db->GetParticleType(idx)->GetPDGCode(); };
  auto name = [&](int idx) { return std::string(db->GetParticleType(idx)->GetName()); };

  cs.pstat.resize(nTypes);
  cs.nquad.resize(nTypes);
  for (int a = 0; a < nTypes; a++) {
    cs.pstat[a] = MakeParentStat(db->GetParticleType(a), T, quantum);
    cs.nquad[a] = NumberDensity(mass(a), cs.pstat[a]);
  }
  for (double& c : cs.c1step) c = 0.;
  cs.nSkippedGamma = 0;

  std::map<std::tuple<int, int, int, int>, Contribution> contribMap;
  for (int a = 0; a < nTypes; a++) {
    ParticleType* type = db->GetParticleType(a);
    if (type->GetGamma() < 0. && type->GetTable() && type->GetTable()->GetChannelCount() >= 0)
      ++cs.nSkippedGamma;
    if (!DecaysInEvent(type)) continue;
    DecayTable* table = type->GetTable();
    for (int c = 0; c <= table->GetChannelCount(); c++) {
      const DecayChannel* ch = table->GetDecayChannel(c);
      std::vector<int> d = ChannelDaughters(ch);
      double B = EffectiveBR(table, c);

      // Counting: n_s (n_s' - delta_ss'), no slot loop.
      for (int p = 0; p < kNPairs; p++) {
        int ns = 0, nsp = 0;
        for (int k : d) {
          ns += (pdg(k) == kPairPDG[p][0]);
          nsp += (pdg(k) == kPairPDG[p][1]);
        }
        cs.c1step[p] += cs.nquad[a] * B * ns * (nsp - (kPairPDG[p][0] == kPairPDG[p][1] ? 1 : 0));
      }

      for (size_t i = 0; i < d.size(); i++) {
        for (size_t j = 0; j < d.size(); j++) {
          if (i == j) continue;
          int p = PairIndex(pdg(d[i]), pdg(d[j]));
          if (p < 0) continue;
          int third = (d.size() == 3) ? d[3 - i - j] : -1;
          auto key  = std::make_tuple(p, a, c, third);
          auto it   = contribMap.find(key);
          if (it != contribMap.end()) {
            ++it->second.count;
            continue;
          }
          Contribution co{};
          co.pair = p; co.parent = a; co.channel = c; co.third = third;
          co.count = 1; co.B = B;
          co.ma = mass(a); co.m1 = mass(d[i]); co.m2 = mass(d[j]);
          co.m3 = (third >= 0) ? mass(third) : 0.;
          co.st = cs.pstat[a];
          co.na = cs.nquad[a];
          co.name = name(a);
          for (int k : d) co.name += "_" + name(k);
          contribMap[key] = co;
        }
      }
    }
  }
  for (auto& kv : contribMap) {
    Contribution co = kv.second;
    if (co.third >= 0) {
      co.phi3 = Phi3Mass(co.ma, co.m1, co.m2, co.m3);
      co.D    = co.count * co.B * co.st.g / (std::pow(2. * M_PI, 3) * co.phi3);
    } else {
      co.A = co.count * PoleA(co.ma, co.m1, co.m2, co.B, co.st.g);
    }
    cs.contribs.push_back(co);
  }
  return cs;
}

// ---------------------------------------------------------------------------------------------
// Binning (spec Sec. 3.2) and bin-integrated predictions (per unit volume)

// M: from threshold to 2 GeV, 1 MeV for the first 30 MeV, then 4 MeV.
inline std::vector<double> MassEdges(double threshold) {
  std::vector<double> edges;
  for (int i = 0; i < 30; i++) edges.push_back(threshold + 0.001 * i);
  for (double e = threshold + 0.030; e < 2.0 - 1e-9; e += 0.004) edges.push_back(e);
  edges.push_back(2.0);
  return edges;
}

// Table kappa grid: 0 - 0.1 GeV by 5 MeV, 0.1 - 1 by 10 MeV, 1 - 3 by 25 MeV.
inline std::vector<double> KappaEdges() {
  std::vector<double> edges;
  for (int i = 0; i < 20; i++) edges.push_back(0.005 * i);
  for (int i = 0; i < 90; i++) edges.push_back(0.1 + 0.010 * i);
  for (int i = 0; i <= 80; i++) edges.push_back(1.0 + 0.025 * i);
  return edges;
}

// MC kappa grid (coarser, a subset of the table edges): 0 - 1 GeV by 50 MeV, 1 - 2 by 100 MeV,
// then one bin 2 - 3 GeV.
inline std::vector<double> MCKappaEdges() {
  std::vector<double> edges;
  for (int i = 0; i < 20; i++) edges.push_back(0.05 * i);
  for (int i = 0; i <= 10; i++) edges.push_back(1.0 + 0.1 * i);
  edges.push_back(3.0);
  return edges;
}

// Pair masses of the ordered pair p (first and second particle).
inline std::pair<double, double> PairMasses(ParticleDB* db, int p) {
  double m1 = 0., m2 = 0.;
  for (int idx = 0; idx < db->GetParticleTypeCount(); idx++) {
    if (db->GetParticleType(idx)->GetPDGCode() == kPairPDG[p][0])
      m1 = db->GetParticleType(idx)->GetMass();
    if (db->GetParticleType(idx)->GetPDGCode() == kPairPDG[p][1])
      m2 = db->GetParticleType(idx)->GetMass();
  }
  return {m1, m2};
}

// L1 integrands (defined after the currents section): kappa^2 (A+B)_c(eps_K) of a pole and
// kappa^2 (C_K + C_u) of a continuum.
inline double PoleKappaIntegrandL1(const Contribution& co, double k);
inline double ContKappaIntegrandL1(const Contribution& co, double M, double k);

// Pole: pair density with kappa in [k1, k2] (all of it at M = ma).
inline double PoleDensityKappa(const Contribution& co, double k1, double k2) {
  if (co.cur)
    return 16. * M_PI * M_PI * PStar(co.ma, co.m1, co.m2) * co.A / co.ma *
           IntegratePanels([&](double k) { return PoleKappaIntegrandL1(co, k); }, {k1, k2});
  return 16. * M_PI * M_PI * PStar(co.ma, co.m1, co.m2) * co.A / co.ma *
         IntegratePanels([&](double k) { return k * k * co.st.f(std::hypot(co.ma, k)); },
                         {k1, k2});
}

// Continuum: pair density with M in [M1, M2] (all kappa), from the Dalitz projection (T1.2).
inline double ContDensityM(const Contribution& co, double M1, double M2) {
  double lo = std::max(M1, co.m1 + co.m2), hi = std::min(M2, co.ma - co.m3);
  if (hi <= lo) return 0.;
  return IntegrateCos(
      [&](double M) {
        return Rho3Dalitz(M, co.m1, co.m2, co.ma, co.m3, co.count * co.B, co.na, co.phi3);
      },
      lo, hi, 32);
}

// Continuum: pair density with M in [M1, M2] and kappa in [k1, k2].
inline double ContDensityMK(const Contribution& co, double M1, double M2, double k1, double k2) {
  double lo = std::max(M1, co.m1 + co.m2), hi = std::min(M2, co.ma - co.m3);
  if (hi <= lo) return 0.;
  return IntegrateCos(
      [&](double M) {
        return 16. * M_PI * M_PI * PStar(M, co.m1, co.m2) *
               IntegratePanels(
                   [&](double k) {
                     if (co.cur) return ContKappaIntegrandL1(co, M, k);
                     CPair c = CKCu3(M, k, co.ma, co.m3, co.D, co.st);
                     return k * k * (c.CK + c.Cu);
                   },
                   {k1, k2});
      },
      lo, hi, 16);
}

// Boltzmann reference of a particle of mass M at rest in the fluid (spec Sec. 3.3): fraction of
// dP_th/dkappa = kappa^2 exp(-eps_K/T) / (M^2 T K_2(M/T)) with kappa in [k1, k2]. The normalization
// int_0^inf kappa^2 exp(-eps_K/T) dkappa = M^2 T K_2(M/T) is done by quadrature; the common factor
// exp(-M/T) is dropped from both integrals.
inline double BoltzmannRefKappa(double M, double k1, double k2, double T) {
  auto f = [&](double k) { return k * k * std::exp(-(std::hypot(M, k) - M) / T); };
  return IntegratePanels(f, {k1, k2}) / IntegratePanels(f, ThermalPanels(T));
}

// ---------------------------------------------------------------------------------------------
// Single-particle currents (spec Sec. 9.2): S_c = A_c p + B_c eps u, gas at rest, per unit volume
// and d^3p (A_c + B_c is the rest-frame density). Thermal part analytic, feed-down part on a grid.

// Barycentric weights of the n-point Gauss-Legendre nodes (for Lagrange interpolation).
inline const std::vector<double>& GLBaryWeights(int n) {
  static std::map<int, std::vector<double>> cache;
  auto it = cache.find(n);
  if (it != cache.end()) return it->second;
  const GaussLegendre& gl = GLRule(n);
  std::vector<double>& bw = cache[n];
  bw.resize(n);
  for (int j = 0; j < n; j++)
    bw[j] = ((j % 2) ? -1. : 1.) * std::sqrt((1. - gl.x[j] * gl.x[j]) * gl.w[j]);
  return bw;
}

// GLRule and GLBaryWeights fill static caches on first use, which is not thread-safe: call this
// before starting threads, so that the threads only read the caches.
inline void WarmQuadratureCaches(int nPerPanel = 16) {
  for (int n : {2, 6, 8, 16, 32, 64, 128, 256, nPerPanel}) GLRule(n);
  GLBaryWeights(nPerPanel);
}

struct SpeciesCurrent {
  std::string name;
  int pdg;
  double m;
  ParentStat st;
  double thScale;            // thermal normalization (n_th / n_quad); 0 if not primordial
  int nPerPanel;
  std::vector<double> edges; // k panel edges, [0, kmax]
  std::vector<double> k, w;  // nodes and weights (dk)
  std::vector<double> A, AB; // feed-down part of A_c and A_c + B_c on the nodes
  bool hasFeed = false;
  double nth = 0., nCount = 0., nGrid = 0.;

  double kmax() const { return edges.back(); }
  double Emax() const { return std::hypot(m, kmax()); }
  double ThermalA(double eps) const {
    return (thScale == 0.) ? 0. : thScale * st.g / std::pow(2. * M_PI, 3) * st.f(eps);
  }
  // Feed-down part at momentum q (Lagrange interpolation on the panel's GL nodes; 0 beyond kmax).
  void Feed(double q, double& a, double& ab) const {
    a = ab = 0.;
    if (!hasFeed || q < 0. || q > kmax()) return;
    int ip = std::upper_bound(edges.begin(), edges.end(), q) - edges.begin() - 1;
    ip = std::min(std::max(ip, 0), static_cast<int>(edges.size()) - 2);
    const GaussLegendre& gl = GLRule(nPerPanel);
    const std::vector<double>& bw = GLBaryWeights(nPerPanel);
    double h = 0.5 * (edges[ip + 1] - edges[ip]), c = 0.5 * (edges[ip + 1] + edges[ip]);
    double t = (q - c) / h;
    double num = 0., numAB = 0., den = 0.;
    const int off = ip * nPerPanel;
    for (int j = 0; j < nPerPanel; j++) {
      double d = t - gl.x[j];
      if (d == 0.) {
        a  = A[off + j];
        ab = AB[off + j];
        return;
      }
      double r = bw[j] / d;
      num += r * A[off + j];
      numAB += r * AB[off + j];
      den += r;
    }
    a  = num / den;
    ab = numAB / den;
  }
  // Full current (A_c, A_c + B_c) at energy E >= m.
  void AtEnergy(double E, double& a, double& ab) const {
    double q = std::sqrt(std::max((E - m) * (E + m), 0.));
    Feed(q, a, ab);
    double th = ThermalA(E);
    a += th;
    ab += th;
  }
};

// One daughter slot b of a decay of parent a: two-body (partner mass mX[0], weight 1) or
// three-body (effective partner masses mX with normalized phase-space weights wX).
struct FeedSlot {
  int parent, daughter;
  double B;
  std::vector<double> mX, wX;
};

// Contribution of a -> b X (fixed mX, branching B) to (A_b, A_b + B_b) at daughter momentum k
// (LCS Eq. twobodyF with f_a replaced by the parent current), written in the parent energy E:
//   A_b     += B ma/(2 p* k^3) int dE (eps E - ma E*) A_a(E),
//   A_b+B_b += B ma/(2 p* k eps) int dE E (A_a + B_a)(E),
// over [E-, min(E+, E_max of the parent grid)]. E is the linear map of the spec's w; integrating
// in E lets the panels follow the thermal scale T at any k.
inline void TwoBodyFeed(double k, double mb, const SpeciesCurrent& par, double mX, double B,
                        double T, double& dA, double& dAB) {
  dA = dAB = 0.;
  const double ma = par.m;
  const double ps = PStar(ma, mb, mX);
  if (ps <= 0.) return;
  const double Es  = (ma * ma + mb * mb - mX * mX) / (2. * ma);
  const double eps = std::hypot(mb, k);
  const double Em  = ma * (Es * Es + k * k) / (eps * Es + k * ps);
  const double Ep  = ma * (eps * Es + k * ps) / (mb * mb);
  // eps E - ma E* = eps y + ma E* k^2/mb^2, y = E - Ec, Ec = (E- + E+)/2, half-width h:
  // the y form avoids the cancellation at small k, the direct form the one for wide windows.
  // In the y form the window is [E-, E- + 2h] exactly (E+ - E- loses the symmetry in y at small k).
  const double h     = ma * k * ps / (mb * mb);
  const bool  yForm  = eps * h < ma * Es;
  const double ck2   = ma * Es * k * k / (mb * mb);
  const double oMax  = std::min(yForm ? 2. * h : Ep - Em, par.Emax() - Em);
  if (oMax <= 0.) return;
  static const double kOff[] = {0., 0.5, 1., 2., 4., 8., 16., 32., 64., 100.};
  const GaussLegendre& gl = GLRule(16);
  double IA = 0., IAB = 0.;
  for (int ip = 0; ip + 1 < 10; ip++) {
    double o1 = kOff[ip] * T, o2 = kOff[ip + 1] * T;
    if (o1 >= oMax) break;
    if (o2 > oMax) o2 = oMax;
    double hw = 0.5 * (o2 - o1), cw = 0.5 * (o2 + o1);
    for (int i = 0; i < 16; i++) {
      double o = cw + hw * gl.x[i];
      double E = Em + o;
      double a, ab;
      par.AtEnergy(E, a, ab);
      double fac = yForm ? eps * (o - h) + ck2 : eps * E - ma * Es;
      IA += gl.w[i] * hw * fac * a;
      IAB += gl.w[i] * hw * E * ab;
    }
  }
  dA  = B * ma / (2. * ps * k * k * k) * IA;
  dAB = B * ma / (2. * ps * k * eps) * IAB;
}

// Sum of all feed slots of daughter b at momentum k.
inline void FeedAt(double k, const SpeciesCurrent& dau, const std::vector<const FeedSlot*>& slots,
                   const std::vector<SpeciesCurrent>& cur, double T, double& a, double& ab) {
  a = ab = 0.;
  for (const FeedSlot* s : slots) {
    for (size_t q = 0; q < s->mX.size(); q++) {
      double dA, dAB;
      TwoBodyFeed(k, dau.m, cur[s->parent], s->mX[q], s->B * s->wX[q], T, dA, dAB);
      a += dA;
      ab += dAB;
    }
  }
}

inline bool IsPhoton(const ParticleType* type) {
  return std::string(type->GetName()).find("gam000zer") != std::string::npos;
}

struct CurrentSet {
  std::vector<SpeciesCurrent> cur;               // per ParticleDB index
  std::vector<FeedSlot> slots;
  std::vector<std::vector<const FeedSlot*>> feeds; // per daughter index
  std::vector<int> order;                          // decreasing mass
  double T;
};

// thScale per species: 1 (quadrature densities) or n_th/n_quad. Photons are never primordial.
// kcut: grid cut in units of T; nPerPanel: GL nodes per k panel; nThreads for the node loop.
inline CurrentSet BuildCurrents(ParticleDB* db, double T, bool quantum,
                                const std::vector<double>& thScale, double kcut = 100.,
                                int nPerPanel = 16, int nThreads = 8) {
  CurrentSet cs;
  cs.T = T;
  WarmQuadratureCaches(nPerPanel);
  const int nTypes = db->GetParticleTypeCount();
  auto mass = [&](int idx) { return static_cast<double>(db->GetParticleType(idx)->GetMass()); };

  std::vector<double> edges;
  for (double e : ThermalPanels(T))
    if (e < kcut * T - 1e-12) edges.push_back(e);
  edges.push_back(kcut * T);
  const GaussLegendre& gl = GLRule(nPerPanel);

  cs.cur.resize(nTypes);
  for (int c = 0; c < nTypes; c++) {
    ParticleType* type = db->GetParticleType(c);
    SpeciesCurrent& s = cs.cur[c];
    s.name = type->GetName();
    s.pdg = type->GetPDGCode();
    s.m = mass(c);
    s.st = MakeParentStat(type, T, quantum);
    s.thScale = IsPhoton(type) ? 0. : thScale[c];
    s.nPerPanel = nPerPanel;
    s.edges = edges;
    for (size_t ip = 0; ip + 1 < edges.size(); ip++) {
      double h = 0.5 * (edges[ip + 1] - edges[ip]), cc = 0.5 * (edges[ip + 1] + edges[ip]);
      for (int i = 0; i < nPerPanel; i++) {
        s.k.push_back(cc + h * gl.x[i]);
        s.w.push_back(h * gl.w[i]);
      }
    }
    s.A.assign(s.k.size(), 0.);
    s.AB.assign(s.k.size(), 0.);
    s.nth = (s.thScale == 0.) ? 0. : s.thScale * NumberDensity(s.m, s.st);
  }

  // Feed slots: every daughter slot of every channel of every decaying parent.
  for (int a = 0; a < nTypes; a++) {
    ParticleType* type = db->GetParticleType(a);
    if (!DecaysInEvent(type)) continue;
    for (int c = 0; c <= type->GetTable()->GetChannelCount(); c++) {
      const DecayChannel* ch = type->GetTable()->GetDecayChannel(c);
      std::vector<int> d = ChannelDaughters(ch);
      for (size_t i = 0; i < d.size(); i++) {
        FeedSlot fs;
        fs.parent = a;
        fs.daughter = d[i];
        fs.B = EffectiveBR(type->GetTable(), c);
        if (d.size() == 2) {
          fs.mX = {mass(d[1 - i])};
          fs.wX = {1.};
        } else {
          // Effective partner mass of the other two daughters, weight p*(ma; mb, mt) p*(mt; mc, md)
          // (three-body phase space with constant matrix element), cos-mapped GL nodes.
          double ma = mass(a), mb = mass(d[i]), mc = mass(d[(i + 1) % 3]), md = mass(d[(i + 2) % 3]);
          double lo = mc + md, hi = ma - mb;
          const int nm = 32;
          const GaussLegendre& gm = GLRule(nm);
          double sum = 0.;
          for (int q = 0; q < nm; q++) {
            double th = 0.5 * M_PI * (gm.x[q] + 1.);
            double mt = lo + 0.5 * (hi - lo) * (1. - std::cos(th));
            double wt = gm.w[q] * std::sin(th) * PStar(ma, mb, mt) * PStar(mt, mc, md);
            fs.mX.push_back(mt);
            fs.wX.push_back(wt);
            sum += wt;
          }
          for (double& wt : fs.wX) wt /= sum;
        }
        cs.slots.push_back(fs);
      }
    }
  }
  cs.feeds.resize(nTypes);
  for (const FeedSlot& fs : cs.slots) cs.feeds[fs.daughter].push_back(&fs);

  cs.order.resize(nTypes);
  for (int c = 0; c < nTypes; c++) cs.order[c] = c;
  std::stable_sort(cs.order.begin(), cs.order.end(),
                   [&](int x, int y) { return cs.cur[x].m > cs.cur[y].m; });

  // Recursion in decreasing mass: when b is reached all its (heavier) parents are complete.
  // Counting (Sec. 9.5) in the same order: n_b = n_b^th + sum_a n_a sum_ch B n_b^ch.
  for (int b : cs.order) {
    SpeciesCurrent& s = cs.cur[b];
    s.nCount = s.nth;
    for (const FeedSlot* fs : cs.feeds[b]) s.nCount += cs.cur[fs->parent].nCount * fs->B;
    if (!cs.feeds[b].empty()) {
      const int nNodes = s.k.size();
      auto work = [&](int t0) {
        for (int i = t0; i < nNodes; i += nThreads)
          FeedAt(s.k[i], s, cs.feeds[b], cs.cur, T, s.A[i], s.AB[i]);
      };
      std::vector<std::thread> pool;
      for (int t = 0; t < nThreads; t++) pool.emplace_back(work, t);
      for (auto& th : pool) th.join();
      s.hasFeed = true;
    }
    s.nGrid = 0.;
    for (size_t i = 0; i < s.k.size(); i++)
      s.nGrid += 4. * M_PI * s.k[i] * s.k[i] * s.w[i] * (s.AB[i] + s.ThermalA(std::hypot(s.m, s.k[i])));
  }
  return cs;
}

// Closed form of LCS Eq. (twobodyF) for a thermal parent (f expanded as sum_n c_n exp(-n E/T)):
// (F_1, F_1 + F_u) of a -> b X at daughter momentum k, parent normalization thScale g/(2 pi)^3.
inline void TwoBodyClosedForm(double k, double mb, double ma, double mX, double B,
                              const ParentStat& st, double thScale, double& F1, double& F1u) {
  const double ps  = PStar(ma, mb, mX);
  const double Es  = (ma * ma + mb * mb - mX * mX) / (2. * ma);
  const double eps = std::hypot(mb, k);
  const double Em  = ma * (Es * Es + k * k) / (eps * Es + k * ps);
  const double Ep  = ma * (eps * Es + k * ps) / (mb * mb);
  const double h   = ma * k * ps / (mb * mb);
  const double ck2 = ma * Es * k * k / (mb * mb);
  double IA = 0., IAB = 0.;
  for (int n = 1; n < 400; n++) {
    double cn   = (st.eta == 0.) ? (n == 1 ? 1. : 0.) : ((st.eta < 0. || n % 2) ? 1. : -1.);
    double beta = n / st.T;
    double em = std::exp(-beta * Em), ep = std::exp(-beta * Ep);
    double J0 = (em - ep) / beta;                                       // int e^{-bE}
    double JE = (Em / beta + 1. / (beta * beta)) * em - (Ep / beta + 1. / (beta * beta)) * ep;
    double x  = beta * h;
    double termA;
    if (x < 0.5) {
      // int y e^{-b(Ec+y)} dy over [-h, h] = e^{-b Ec} 2 (sinh x - x cosh x)/b^2,
      // sinh x - x cosh x = -sum_j 2j x^{2j+1}/(2j+1)!
      double s = 0., term = x; // x^{2j+1}/(2j+1)!
      for (int j = 1; j < 30; j++) {
        term *= x * x / ((2. * j) * (2. * j + 1.));
        s -= 2. * j * term;
      }
      double eEc = std::exp(-beta * (Em + h));
      double Jy  = eEc * 2. * s / (beta * beta);
      termA = eps * Jy + ck2 * eEc * 2. * std::sinh(x) / beta;
    } else {
      termA = eps * JE - ma * Es * J0;
    }
    IA += cn * termA;
    IAB += cn * JE;
    if (cn != 0. && std::exp(-(n - 1.) * Em / st.T) < 1e-18) break;
    if (st.eta == 0.) break;
  }
  const double norm = thScale * st.g / std::pow(2. * M_PI, 3);
  F1  = B * ma / (2. * ps * k * k * k) * norm * IA;
  F1u = B * ma / (2. * ps * k * eps) * norm * IAB;
}

// ---------------------------------------------------------------------------------------------
// L1 pair functions (spec Sec. 9.3): the thermal f_a replaced by the parent current S_c.

// Three-body continuum with parent current c, D = count B / Phi3. Written as the L0 arithmetic of
// CKCu3 for the A_c part plus the B_c term of C_u, so a thermal-only current reproduces CKCu3:
//   C_K = 2 pi D p3/(kappa M^2) (kappa E_P I0[A] + eps_K p3 I1[A]),
//   C_u = -2 pi D p3^2/(kappa eps_K) I1[A] + 2 pi D p3/(M eps_K) int dt eps_P B_c(eps_P).
inline CPair CKCu3Current(double M, double kappa, double ma, double m3, double D,
                          const SpeciesCurrent& c, int nt = 64) {
  const GaussLegendre& gl = GLRule(nt);
  double eK = std::hypot(M, kappa);
  double EP = (ma * ma + M * M - m3 * m3) / (2. * M);
  double p3 = PStar(M, ma, m3);
  double I0 = 0., I1 = 0., IB = 0.;
  for (int i = 0; i < nt; i++) {
    double eP = (eK * EP + kappa * p3 * gl.x[i]) / M;
    double a, ab;
    c.AtEnergy(eP, a, ab);
    I0 += gl.w[i] * a;
    I1 += gl.w[i] * gl.x[i] * a;
    IB += gl.w[i] * eP * (ab - a);
  }
  double CK = 2. * M_PI * D * p3 / (kappa * M * M) * (kappa * EP * I0 + eK * p3 * I1);
  double Cu = -2. * M_PI * D * p3 * p3 / (kappa * eK) * I1 + 2. * M_PI * D * p3 / (M * eK) * IB;
  return {CK, Cu};
}

inline double PoleKappaIntegrandL1(const Contribution& co, double k) {
  double a, ab;
  co.cur->AtEnergy(std::hypot(co.ma, k), a, ab);
  return k * k * ab;
}

inline double ContKappaIntegrandL1(const Contribution& co, double M, double k) {
  CPair c = CKCu3Current(M, k, co.ma, co.m3, co.D, *co.cur);
  return k * k * (c.CK + c.Cu);
}

// rho(M) of an L1 continuum from the kappa integral (T7.2).
inline double Rho3Current(double M, const Contribution& co) {
  return 16. * M_PI * M_PI * PStar(M, co.m1, co.m2) *
         IntegratePanels([&](double k) { return ContKappaIntegrandL1(co, M, k); },
                         ThermalPanels(co.st.T));
}

// Turns the L0 contributions into L1 ones (parent current, n_c) and c1step into c^L1.
inline void ToL1(ContributionSet& cs, ParticleDB* db, const CurrentSet& cur) {
  for (Contribution& co : cs.contribs) {
    co.cur = &cur.cur[co.parent];
    co.na  = cur.cur[co.parent].nCount;
    if (co.third >= 0)
      co.D = co.count * co.B / co.phi3;
    else
      co.A = co.count * co.B / Phi2(co.ma, co.m1, co.m2);
  }
  for (double& c : cs.c1step) c = 0.;
  for (int a = 0; a < db->GetParticleTypeCount(); a++) {
    ParticleType* type = db->GetParticleType(a);
    if (!DecaysInEvent(type)) continue;
    for (int c = 0; c <= type->GetTable()->GetChannelCount(); c++) {
      const DecayChannel* ch = type->GetTable()->GetDecayChannel(c);
      std::vector<int> d = ChannelDaughters(ch);
      for (int p = 0; p < kNPairs; p++) {
        int ns = 0, nsp = 0;
        for (int k : d) {
          ns += (db->GetParticleType(k)->GetPDGCode() == kPairPDG[p][0]);
          nsp += (db->GetParticleType(k)->GetPDGCode() == kPairPDG[p][1]);
        }
        cs.c1step[p] += cur.cur[a].nCount * EffectiveBR(type->GetTable(), c) * ns *
                        (nsp - (kPairPDG[p][0] == kPairPDG[p][1] ? 1 : 0));
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------
// Chain kernels and the L2 table (spec Sec. 9.4-9.5, PFC Sec. 4-5)

// Measured species (index 0..5) of the ordered pairs.
static constexpr int kNMeasured = 6;
static constexpr int kMeasPDG[kNMeasured] = {211, -211, 321, -321, 2212, -2212};
inline int MeasuredIndex(int pdg) {
  for (int i = 0; i < kNMeasured; i++)
    if (kMeasPDG[i] == pdg) return i;
  return -1;
}

// Counting of Sec. 9.5. nbar[d][s]: mean number of final measured s per particle d (decayed with
// the parsed BRs, surviving with 1 - sum B). Per parent c and ordered pair p, per decay:
// allLCA = sum_ch B sum_{i!=j} nbar_{d_i -> s} nbar_{d_j -> s'} (pairs with LCA at the vertex),
// oneStep = sum_ch B n_s (n_s' - delta) (both members direct daughters); chain = allLCA - oneStep.
struct ChainCounting {
  std::vector<std::array<double, kNMeasured>> nbar;
  std::vector<std::array<double, kNPairs>> allLCA, oneStep;
};

inline ChainCounting CountChains(ParticleDB* db) {
  const int nTypes = db->GetParticleTypeCount();
  ChainCounting cc;
  cc.nbar.assign(nTypes, {});
  cc.allLCA.assign(nTypes, {});
  cc.oneStep.assign(nTypes, {});
  std::vector<int> order(nTypes);
  for (int c = 0; c < nTypes; c++) order[c] = c;
  std::stable_sort(order.begin(), order.end(), [&](int x, int y) {
    return db->GetParticleType(x)->GetMass() < db->GetParticleType(y)->GetMass();
  });
  for (int d : order) {
    ParticleType* type = db->GetParticleType(d);
    std::array<double, kNMeasured>& nb = cc.nbar[d];
    nb.fill(0.);
    double surv = 1.;
    if (DecaysInEvent(type)) {
      for (int c = 0; c <= type->GetTable()->GetChannelCount(); c++) {
        const DecayChannel* ch = type->GetTable()->GetDecayChannel(c);
        double B = EffectiveBR(type->GetTable(), c);
        surv -= B;
        std::vector<int> dd = ChannelDaughters(ch);
        for (int k : dd)
          for (int s = 0; s < kNMeasured; s++) nb[s] += B * cc.nbar[k][s];
        for (int p = 0; p < kNPairs; p++) {
          int s1 = MeasuredIndex(kPairPDG[p][0]), s2 = MeasuredIndex(kPairPDG[p][1]);
          int n1 = 0, n2 = 0;
          for (size_t i = 0; i < dd.size(); i++) {
            n1 += (db->GetParticleType(dd[i])->GetPDGCode() == kPairPDG[p][0]);
            n2 += (db->GetParticleType(dd[i])->GetPDGCode() == kPairPDG[p][1]);
            for (size_t j = 0; j < dd.size(); j++)
              if (i != j) cc.allLCA[d][p] += B * cc.nbar[dd[i]][s1] * cc.nbar[dd[j]][s2];
          }
          cc.oneStep[d][p] += B * n1 * (n2 - (s1 == s2 ? 1 : 0));
        }
      }
    }
    int sm = MeasuredIndex(type->GetPDGCode());
    if (sm >= 0) nb[sm] += surv;
  }
  return cc;
}

// Kernel of one parent c: ordered measured pairs per particle c at rest, from N decays.
// counts[p][im * (nKc + 1) + ik]: chain pairs (different first-generation branches, not both direct
// daughters) in M bin im of the pair's MassEdges and k_c bin ik (ik = nKc: k_c overflow); counts with
// M above the last edge go to mOverflow[p]. sAll/s2All, sChain/s2Chain: sums over trees of the
// per-tree number of all-LCA and chain pairs and of its square.
struct ChainKernel {
  int parent = -1;
  long N = 0;
  std::vector<double> kcEdges;
  std::vector<std::vector<double>> mEdges; // per pair
  std::vector<std::vector<double>> counts; // per pair
  std::array<double, kNPairs> mOverflow{}, sAll{}, s2All{}, sChain{}, s2Chain{};
  long nNonFinite = 0; // final measured particles with non-finite momentum (skipped)
  int nKc() const { return static_cast<int>(kcEdges.size()) - 1; }
};

inline ChainKernel GenerateKernel(ParticleDB* db, int c, long N,
                                  const std::vector<std::vector<double>>& mEdges,
                                  const std::vector<double>& kcEdges, long batch = 10000) {
  ChainKernel kern;
  kern.parent = c;
  kern.N = N;
  kern.kcEdges = kcEdges;
  kern.mEdges = mEdges;
  const int nKc = kern.nKc();
  for (int p = 0; p < kNPairs; p++)
    kern.counts.emplace_back((mEdges[p].size() - 1) * (nKc + 1), 0.);
  ParticleType* type = db->GetParticleType(c);
  const double m = type->GetMass();

  struct Final {
    int s;          // measured index
    int branch;     // eid of the first-generation daughter it descends from
    bool direct;    // direct daughter of c
    double e, px, py, pz;
  };
  std::vector<const Particle*> byEid;
  std::vector<std::vector<Final>> trees;
  for (long done = 0; done < N; done += batch) {
    const long nb = std::min(batch, N - done);
    Particle::ZeroEID();
    std::list<Particle> particles;
    // "At rest" with a tiny momentum: ParticleDecayer's boost has (gamma - 1)/v^2 = 0/0 at v = 0.
    // The shift of k_c is below 1e-7 GeV.
    const double pz = 1e-7 * m;
    for (long i = 0; i < nb; i++) {
      Particle part(type);
      part.SetParticlePX(std::hypot(m, pz), 0., 0., pz, 0., 0., 0., 0.);
      particles.push_back(part);
    }
    // As Event::DecayParticles: the end of the list moves while daughters are appended.
    ParticleDecayer decayer(db, &particles);
    for (auto it = particles.begin(); it != particles.end(); ++it)
      if (DecaysInEvent(it->GetParticleType())) decayer.DecayParticle(&(*it));

    int maxEid = 0;
    for (const Particle& q : particles) maxEid = std::max(maxEid, q.eid);
    byEid.assign(maxEid + 1, nullptr);
    for (const Particle& q : particles) byEid[q.eid] = &q;
    trees.assign(nb, {});
    for (const Particle& q : particles) {
      if (q.GetDecayed() || q.fathereid == -1) continue;
      int s = MeasuredIndex(q.pid);
      if (s < 0) continue;
      if (!std::isfinite(q.e) || !std::isfinite(q.px) || !std::isfinite(q.py) ||
          !std::isfinite(q.pz)) {
        ++kern.nNonFinite;
        continue;
      }
      const Particle* a = &q;
      while (a->fathereid != a->rooteid) a = byEid[a->fathereid];
      // roots are the first nb particles created after ZeroEID: rooteid in 0..nb-1
      trees[q.rooteid].push_back({s, a->eid, q.fathereid == q.rooteid, q.e, q.px, q.py, q.pz});
    }
    for (const std::vector<Final>& t : trees) {
      std::array<double, kNPairs> nAll{}, nChain{};
      for (const Final& f1 : t) {
        for (const Final& f2 : t) {
          if (&f1 == &f2 || f1.branch == f2.branch) continue;
          int p = PairIndex(kMeasPDG[f1.s], kMeasPDG[f2.s]);
          if (p < 0) continue;
          nAll[p] += 1.;
          if (f1.direct && f2.direct) continue;
          nChain[p] += 1.;
          double E = f1.e + f2.e, px = f1.px + f2.px, py = f1.py + f2.py, pz = f1.pz + f2.pz;
          double kc = std::sqrt(px * px + py * py + pz * pz);
          double M  = std::sqrt(std::max(E * E - kc * kc, 0.));
          const std::vector<double>& me = mEdges[p];
          if (M >= me.back()) {
            kern.mOverflow[p] += 1.;
            continue;
          }
          int im = std::upper_bound(me.begin(), me.end(), M) - me.begin() - 1;
          if (im < 0) im = 0; // rounding at threshold
          int ik = std::upper_bound(kcEdges.begin(), kcEdges.end(), kc) - kcEdges.begin() - 1;
          if (ik >= nKc) ik = nKc;
          kern.counts[p][im * (nKc + 1) + ik] += 1.;
        }
      }
      for (int p = 0; p < kNPairs; p++) {
        kern.sAll[p] += nAll[p];
        kern.s2All[p] += nAll[p] * nAll[p];
        kern.sChain[p] += nChain[p];
        kern.s2Chain[p] += nChain[p] * nChain[p];
      }
    }
  }
  return kern;
}

// Complementary primitives of the current of one parent in the fluid-frame energy, integrated from
// e to infinity: G0 = int A, G1 = int e A, F1 = int e (A+B), on a uniform grid from m to E_max (the
// current is cut there), evaluated by cubic Hermite interpolation with the exact derivatives
// -A, -e A, -e (A+B). Integrating from above keeps differences deep in the thermal tail accurate.
struct CurrentPrimitives {
  double e0 = 0., h = 0.;
  std::vector<double> A, AB, G0, G1, F1;

  void Build(const SpeciesCurrent& c, double T, double step = 0.025) {
    e0 = c.m;
    h  = step * T;
    const int n = static_cast<int>((c.Emax() - e0) / h) + 2;
    A.resize(n); AB.resize(n); G0.assign(n, 0.); G1.assign(n, 0.); F1.assign(n, 0.);
    for (int i = 0; i < n; i++) c.AtEnergy(e0 + i * h, A[i], AB[i]);
    const GaussLegendre& gl = GLRule(6);
    for (int i = 0; i + 1 < n; i++) {
      double g0 = 0., g1 = 0., f1 = 0.;
      for (int q = 0; q < 6; q++) {
        double e = e0 + (i + 0.5 * (1. + gl.x[q])) * h, a, ab;
        c.AtEnergy(e, a, ab);
        g0 += 0.5 * h * gl.w[q] * a;
        g1 += 0.5 * h * gl.w[q] * e * a;
        f1 += 0.5 * h * gl.w[q] * e * ab;
      }
      G0[i] = g0;
      G1[i] = g1;
      F1[i] = f1;
    }
    for (int i = n - 2; i >= 0; i--) { // accumulate from above
      G0[i] += G0[i + 1];
      G1[i] += G1[i + 1];
      F1[i] += F1[i + 1];
    }
  }
  void Eval(double e, double& g0, double& g1, double& f1) const {
    const int n = A.size();
    double x = (e - e0) / h;
    if (x <= 0.) { g0 = G0[0]; g1 = G1[0]; f1 = F1[0]; return; }
    if (x >= n - 1) { g0 = g1 = f1 = 0.; return; }
    int i = static_cast<int>(x);
    double t = x - i, t2 = t * t, t3 = t2 * t;
    double h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + t, h01 = -2 * t3 + 3 * t2,
           h11 = t3 - t2;
    double ea = e0 + i * h, eb = ea + h;
    g0 = h00 * G0[i] + h01 * G0[i + 1] - h * (h10 * A[i] + h11 * A[i + 1]);
    g1 = h00 * G1[i] + h01 * G1[i + 1] - h * (h10 * ea * A[i] + h11 * eb * A[i + 1]);
    f1 = h00 * F1[i] + h01 * F1[i + 1] - h * (h10 * ea * AB[i] + h11 * eb * AB[i + 1]);
  }
};

// t integrals of the master formula at pair mass M, kappa, for parent energy/momentum EP, pP in the
// pair rest frame: IA = int dt A(eps_P)(kappa EP + eps_K pP t)/(M kappa),
// IAB = int dt (A+B)(eps_P) eps_P/eps_K, eps_P = (eps_K EP + kappa pP t)/M. Wide windows use the
// primitives (eps = linear in t), narrow ones 16-pt Gauss-Legendre in t.
inline void ChainT(double M, double kappa, double EP, double pP, const SpeciesCurrent& c,
                   const CurrentPrimitives& prim, double T, double& IA, double& IAB) {
  const double eK = std::hypot(M, kappa);
  const double w  = 2. * kappa * pP / M;
  if (w > T) {
    double em = (eK * EP - kappa * pP) / M, ep = (eK * EP + kappa * pP) / M;
    double a0, a1, af, b0, b1, bf;
    prim.Eval(em, a0, a1, af);
    prim.Eval(ep, b0, b1, bf);
    // complementary primitives: int_{em}^{ep} = G(em) - G(ep)
    IA  = M / (kappa * kappa * kappa * pP) * (eK * (a1 - b1) - M * EP * (a0 - b0));
    IAB = M / (kappa * pP * eK) * (af - bf);
    return;
  }
  const GaussLegendre& gl = GLRule(16);
  IA = IAB = 0.;
  for (int i = 0; i < 16; i++) {
    double t = gl.x[i], eP = (eK * EP + kappa * pP * t) / M, a, ab;
    c.AtEnergy(eP, a, ab);
    IA += gl.w[i] * a * (kappa * EP + eK * pP * t) / (M * kappa);
    IAB += gl.w[i] * ab * eP / eK;
  }
}

// Chain contribution of all parents, per ordered pair, accumulated by ConvolveChain:
// CK, CKCu: C_K and C_K + C_u at the table (M bin centre, kappa centre) points, row-major in M;
// rho, rhoVar: density per M bin (index nM = above the last edge) from rho = n_c int dk_c h_c;
// mk, mkVar: density in (M bin, MC kappa bin) cells.
struct ChainTable {
  std::vector<double> kTable, kMCEdges;
  std::vector<std::vector<double>> mEdges;
  std::array<std::vector<double>, kNPairs> CK, CKCu, rho, rhoVar, mk, mkVar;
  double worstNorm = 0.; // max |sum_kappa mk / (n_c counts/N) - 1| over parents and M rows
                         // (kappa integrated to 15 GeV for this check)
  std::string worstNormName;

  void Init(const std::vector<std::vector<double>>& me, const std::vector<double>& kt,
            const std::vector<double>& kmc) {
    mEdges = me; kTable = kt; kMCEdges = kmc;
    for (int p = 0; p < kNPairs; p++) {
      const size_t nM = me[p].size() - 1;
      CK[p].assign(nM * kt.size(), 0.);
      CKCu[p].assign(nM * kt.size(), 0.);
      rho[p].assign(nM + 1, 0.);
      rhoVar[p].assign(nM + 1, 0.);
      mk[p].assign(nM * (kmc.size() - 1), 0.);
      mkVar[p].assign(nM * (kmc.size() - 1), 0.);
    }
  }
};

// Master formula (PFC Sec. 4) for one parent: the current carries n_c, the kernel the pairs per c.
// The t integrals are shared by the pairs of one species (same M bins); rows are split over threads.
inline void ConvolveChain(const ChainKernel& kern, const SpeciesCurrent& cur, double T,
                          const std::vector<std::pair<double, double>>& pairMasses,
                          ChainTable& out, int nThreads = 8) {
  WarmQuadratureCaches(cur.nPerPanel);
  CurrentPrimitives prim;
  prim.Build(cur, T);
  const double mc = cur.m, nc = cur.nCount;
  const int nKc = kern.nKc();
  const int nKt = out.kTable.size(), nKm = out.kMCEdges.size() - 1;
  const GaussLegendre& gk = GLRule(2);
  const GaussLegendre& gm = GLRule(6);
  // MC kappa quadrature nodes and weights (dkappa), per MC bin
  std::vector<double> mcK, mcW;
  for (int b = 0; b < nKm; b++) {
    double h = 0.5 * (out.kMCEdges[b + 1] - out.kMCEdges[b]),
           c = 0.5 * (out.kMCEdges[b + 1] + out.kMCEdges[b]);
    for (int q = 0; q < 6; q++) {
      mcK.push_back(c + h * gm.x[q]);
      mcW.push_back(h * gm.w[q]);
    }
  }
  // kappa beyond the MC range, for the normalization check only
  std::vector<double> tailK, tailW;
  {
    const GaussLegendre& gt = GLRule(16);
    double lo = out.kMCEdges.back(), hi = 15.;
    for (int q = 0; q < 16; q++) {
      tailK.push_back(0.5 * (lo + hi) + 0.5 * (hi - lo) * gt.x[q]);
      tailW.push_back(0.5 * (hi - lo) * gt.w[q]);
    }
  }
  std::string worstName;
  double worst = 0.;
  for (int fam = 0; fam < 3; fam++) {
    std::vector<int> pairs;
    for (int p = 4 * fam; p < 4 * fam + 4; p++)
      if (kern.sChain[p] > 0.) pairs.push_back(p);
    if (pairs.empty()) continue;
    const std::vector<double>& me = kern.mEdges[pairs[0]];
    const int nM = me.size() - 1;
    const double m1 = pairMasses[pairs[0]].first, m2 = pairMasses[pairs[0]].second;
    std::vector<double> rowWorst(nM, 0.);
    auto work = [&](int t0) {
      std::vector<double> IAt(nKt), IABt(nKt), IABm(mcK.size());
      for (int im = t0; im < nM; im += nThreads) {
        const double M = 0.5 * (me[im] + me[im + 1]), dM = me[im + 1] - me[im];
        const double ps = PStar(M, m1, m2);
        if (ps <= 0.) continue;
        double rowCounts = 0., rowMk = 0., rowTail = 0.;
        for (int ik = 0; ik < nKc; ik++) {
          bool any = false;
          for (int p : pairs) any = any || kern.counts[p][im * (nKc + 1) + ik] > 0.;
          if (!any) continue;
          // per-count weights of this k_c bin in the table points and MC cells
          std::vector<double> wCK(nKt, 0.), wCKCu(nKt, 0.), wMk(nKm, 0.);
          double wTail = 0.;
          for (int g = 0; g < 2; g++) {
            double k1 = kern.kcEdges[ik], k2 = kern.kcEdges[ik + 1];
            double kc = 0.5 * (k1 + k2) + 0.5 * (k2 - k1) * gk.x[g];
            double EP = mc / M * std::hypot(M, kc), pP = mc / M * kc;
            double pref = 1. / (4. * M_PI * ps) * mc * mc / (2. * M * M) * 0.5 * gk.w[g] /
                          (kern.N * dM);
            for (int ik2 = 0; ik2 < nKt; ik2++) {
              double IA, IAB;
              ChainT(M, out.kTable[ik2], EP, pP, cur, prim, T, IA, IAB);
              wCK[ik2] += pref * IA;
              wCKCu[ik2] += pref * IAB;
            }
            for (size_t q = 0; q < mcK.size(); q++) {
              double IA, IAB;
              ChainT(M, mcK[q], EP, pP, cur, prim, T, IA, IAB);
              wMk[q / 6] += dM * 16. * M_PI * M_PI * ps * mcW[q] * mcK[q] * mcK[q] * pref * IAB;
            }
            for (size_t q = 0; q < tailK.size(); q++) {
              double IA, IAB;
              ChainT(M, tailK[q], EP, pP, cur, prim, T, IA, IAB);
              wTail += dM * 16. * M_PI * M_PI * ps * tailW[q] * tailK[q] * tailK[q] * pref * IAB;
            }
          }
          for (int p : pairs) {
            double cnt = kern.counts[p][im * (nKc + 1) + ik];
            if (cnt <= 0.) continue;
            double infl = kern.s2Chain[p] / kern.sChain[p];
            for (int ik2 = 0; ik2 < nKt; ik2++) {
              out.CK[p][im * nKt + ik2] += wCK[ik2] * cnt;
              out.CKCu[p][im * nKt + ik2] += wCKCu[ik2] * cnt;
            }
            for (int b = 0; b < nKm; b++) {
              out.mk[p][im * nKm + b] += wMk[b] * cnt;
              out.mkVar[p][im * nKm + b] += wMk[b] * wMk[b] * cnt * infl;
            }
            if (p == pairs[0]) {
              rowCounts += cnt;
              rowTail += wTail * cnt;
              for (int b = 0; b < nKm; b++) rowMk += wMk[b] * cnt;
            }
          }
        }
        // normalization: sum over kappa of the cells (plus the kappa > 3 GeV tail, which the cells
        // do not hold) = n_c (pairs in the row) / N
        if (rowCounts > 0.)
          rowWorst[im] = std::abs((rowMk + rowTail) / (nc * rowCounts / kern.N) - 1.);
      }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < nThreads; t++) pool.emplace_back(work, t);
    for (auto& th : pool) th.join();
    for (int im = 0; im < nM; im++)
      if (rowWorst[im] > worst) {
        worst = rowWorst[im];
        char buf[64];
        std::snprintf(buf, sizeof(buf), "M = %.4f", 0.5 * (me[im] + me[im + 1]));
        worstName = buf;
      }
    // rho per M bin from the normalization identity (k_c overflow included)
    for (int p : pairs) {
      double infl = kern.s2Chain[p] / kern.sChain[p];
      for (int im = 0; im < nM; im++) {
        double cnt = 0.;
        for (int ik = 0; ik <= nKc; ik++) cnt += kern.counts[p][im * (nKc + 1) + ik];
        out.rho[p][im] += nc * cnt / kern.N;
        out.rhoVar[p][im] += (nc / kern.N) * (nc / kern.N) * cnt * infl;
      }
      out.rho[p][nM] += nc * kern.mOverflow[p] / kern.N;
      out.rhoVar[p][nM] += (nc / kern.N) * (nc / kern.N) * kern.mOverflow[p] * infl;
    }
  }
  if (worst > out.worstNorm) {
    out.worstNorm = worst;
    out.worstNormName = cur.name + ", " + worstName;
  }
}

#endif
