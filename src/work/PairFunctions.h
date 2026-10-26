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
#include "ParticleDB.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
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
      double B = ch->GetBranchingRatio();

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
        fs.B = ch->GetBranchingRatio();
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
        cs.c1step[p] += cur.cur[a].nCount * ch->GetBranchingRatio() * ns *
                        (nsp - (kPairPDG[p][0] == kPairPDG[p][1] ? 1 : 0));
      }
    }
  }
}

#endif
