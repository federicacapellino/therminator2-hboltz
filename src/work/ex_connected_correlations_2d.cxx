// ex_connected_correlations_2d.cxx — Connected two-particle correlations, following Eq. (2) of
// Short_range_hadron_correlations.pdf: C_ss',conn = <2-body> - <uncorrelated product>.
//
// The full differential object is a 3D histogram per term (Delta y, Delta pT, Delta phi); the
// "uncorrelated" term is estimated per Section II (Eqs. 5-8) by randomizing only the rapidity of
// each particle in a pair over the full accepted window W = 2*Y (particles are accepted with
// |y| <= Y), reweighting by W^2/(W-delta)^2, and keeping the pair's real Delta pT/Delta phi. Both
// 2D projections (Delta y vs Delta pT, Delta y vs Delta phi) are then derived from the same 3D
// cube via TH3::Project3D ("yx"/"zx" — Project3D puts its *second* letter on the new X axis), so
// they are consistent slices of one object rather than independently filled histograms. Pure 1D
// marginals (Delta y, Delta pT, Delta phi alone) are likewise derived via the single-letter form
// ("x"/"y"/"z" — unambiguous, unlike the 2-letter form).
//
// The pair invariant mass M^2 = (E1+E2)^2 - (p1+p2)^2 is NOT a function of (Delta y, Delta pT,
// Delta phi) alone, so it is filled independently from each pair's real 4-momenta rather than
// projected from the 3D cube. Its "uncorrelated" term reconstructs each particle's E/pz from its
// (rapidity-invariant) transverse mass and the same randomized rapidity used for the other terms,
// keeping pT/phi/mass at their real values, so the M^2 background is consistent with the
// Delta y/Delta pT/Delta phi background (Eqs. 5-8) rather than an independent recipe.
//
// eps1 = p1.u, eps2 = p2.u is the comoving energy of each paired particle in the local fluid rest
// frame -- the same quantity that enters the model's thermal Boltzmann/Bose/Fermi exponent
// (Model_HRG.cxx's and Model_BlastWave.cxx's "PdotU"). For HRG, u=(1,0,0,0) everywhere (static
// box), so eps is just the lab energy. For BlastWave, u^mu depends on the emission/decay point, so
// eps is reconstructed from the particle's stored 4-position (spatial rapidity/azimuthal angle) and
// the model's transverse flow velocity (Model_BlastWave::GetVt()), mirroring the model's own p.u
// formula. This is computed for ALL final-state particles (primordial and resonance-decay
// daughters alike, same population as everything else here) -- for BlastWave decay daughters this
// evaluates the flow field at the decay vertex rather than the strict freeze-out point, since decay
// products are displaced from freeze-out by the resonance's lifetime (see ParticleDecayer.cxx).
// eps is a single-particle quantity, not a pair difference, so unlike the observables above it has
// no "uncorrelated"/"connected" counterpart -- just one histogram per role (eps1/eps2) per
// category, filled once per real pair occurrence exactly like the other observables here.
//
// Direct PDG ID filtering (no unsafe ParticleDB charge lookups during event iteration).
// Usage (generate in-process):
//   ex_connected_correlations_2d.exe [--model hrg|blastwave] [--ini <file>]
//                                    [--events <N>] [--ymax <Y>] [--delta <val>]
// Usage (read pre-generated events from therm2_events.exe instead of regenerating them --
// avoids paying the Monte-Carlo generation cost again for a re-analysis):
//   ex_connected_correlations_2d.exe --input-dir events/ [--events <N>] [--ymax <Y>] [--delta <val>]

#include "CLI11.hpp"
#include "Event.h"
#include "Integrator.h"
#include "Model.h"
#include "Model_BlastWave.h"
#include "Model_HRG.h"
#include "ParticleCoor.h"
#include "ParticleDB.h"
#include "StructModel.h"
#include "THGlobal.h"
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TRandom.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>
#include <memory>
#include <string>

struct CachedParticle {
  double y;
  double pt;
  double phi;
  double px, py, pz, e; // raw 4-momentum, for invariant mass
  double mt;            // sqrt(e*e - pz*pz) == sqrt(mass*mass + pt*pt); rapidity-invariant, so it
                         // lets the uncorrelated (randomized-y) term reconstruct e/pz consistently
  double eps;           // comoving energy p.u in the local fluid rest frame (== e for HRG, since
                         // u=(1,0,0,0) everywhere; flow-boosted for BlastWave)
  int pid;
};

struct SpeciesPairConfig {
  std::string name;
  int pid1; // e.g.  211 for pi+; 0 for all charged positive
  int pid2; // e.g. -211 for pi-; 0 for all charged negative
};

// Takes ParticleCoor (rather than Particle) so it works equally on a live Particle from in-process
// generation and on a raw ParticleCoor read back from a saved event*.root file (Particle is-a
// ParticleCoor, so live Particle objects still bind here without any change at call sites).
static double getRapidity(const ParticleCoor& p) {
  double e  = p.e;
  double pz = p.pz;
  if (e <= std::abs(pz)) return (pz >= 0) ? 100.0 : -100.0;
  return 0.5 * std::log((e + pz) / (e - pz));
}

// Spatial rapidity of the particle's emission/decay point (same form as getRapidity(), but on
// position rather than momentum) — used to reconstruct the BlastWave fluid velocity u^mu at that
// point.
static double getSpatialRapidity(const ParticleCoor& p) {
  double t = p.t;
  double z = p.z;
  if (t <= std::abs(z)) return (z >= 0) ? 100.0 : -100.0;
  return 0.5 * std::log((t + z) / (t - z));
}

// Only used by the 3D (Delta y, Delta pT, Delta phi) histogram fills, currently commented out in
// processSpeciesPair2D() — re-enable together with those.
// static double getDeltaPhi(double phi1, double phi2) {
//   double dphi = std::abs(phi1 - phi2);
//   while (dphi > M_PI) dphi -= 2.0 * M_PI;
//   return std::abs(dphi);
// }

void processSpeciesPair2D(const SpeciesPairConfig& cfg,
                          const std::vector<std::vector<CachedParticle>>& eventStore,
                          double Y, double delta,
                          int nBinsY, int nBinsPt, int nBinsPhi, int nBinsM2, int nBinsEps,
                          double ptMax, double m2Max, double epsMax) {

  std::string filename = "corr2d_" + cfg.name + ".root";
  TFile fout(filename.c_str(), "RECREATE");
  TH1::SetDefaultSumw2();

  double dyMax = 2.0 * Y;
  // Full width of the accepted |y| <= Y rapidity window — this is the "Y" of Eqs. 5-8 in the
  // paper (there, particles are produced over a window of width Y; here that width is 2*Y since
  // particles are accepted with |y| <= Y). Must be used consistently for the reweight factor and
  // the randomization support below, or the phase-space correction desyncs from the real acceptance.
  const double W = dyMax;
  const double pairWeight = (W * W) / ((W - delta) * (W - delta));

  // --- Full 3D (Delta y, Delta pT, Delta phi) histogram machinery commented out for now: only the
  // invariant-mass and eps1/eps2 projections are needed at the moment, and skipping the TH3D fills
  // (pricier per-call than TH1D) + Project3D walks cuts a meaningful chunk of the per-pair cost.
  // Re-enable by uncommenting this block, the h3_*->Fill(...) calls below, subtract3D, and
  // project2D/project1D.
  //
  // auto make3DSet = [&](const std::string& tag, TH3D*& h2body, TH3D*& huncorr, TH3D*& hconn) {
  //   std::string name2 = "h3_2body_" + tag;
  //   std::string nameU = "h3_uncorr_" + tag;
  //   std::string nameC = "h3_conn_" + tag;
  //
  //   std::string title = ";#Delta y;#Delta p_{T} [GeV/c];#Delta#phi [rad]";
  //
  //   h2body  = new TH3D(name2.c_str(), ("Physical 2-Body " + title).c_str(),
  //                       nBinsY, 0., dyMax, nBinsPt, 0., ptMax, nBinsPhi, 0., M_PI);
  //   huncorr = new TH3D(nameU.c_str(), ("Uncorrelated " + title).c_str(),
  //                       nBinsY, 0., dyMax, nBinsPt, 0., ptMax, nBinsPhi, 0., M_PI);
  //   hconn   = new TH3D(nameC.c_str(), ("Connected C_{conn} " + title).c_str(),
  //                       nBinsY, 0., dyMax, nBinsPt, 0., ptMax, nBinsPhi, 0., M_PI);
  // };
  //
  // TH3D *h3_2body_ss1, *h3_uncorr_ss1, *h3_conn_ss1; // opposite charge: species s, s'
  // TH3D *h3_2body_s,   *h3_uncorr_s,   *h3_conn_s;   // same charge: species s, s
  // TH3D *h3_2body_s1,  *h3_uncorr_s1,  *h3_conn_s1;  // same charge: species s', s'
  //
  // make3DSet("ss1", h3_2body_ss1, h3_uncorr_ss1, h3_conn_ss1);
  // make3DSet("s",   h3_2body_s,   h3_uncorr_s,   h3_conn_s);
  // make3DSet("s1",  h3_2body_s1,  h3_uncorr_s1,  h3_conn_s1);

  // Invariant mass M^2 = (E1+E2)^2 - (p1+p2)^2 is not a function of (Delta y, Delta pT, Delta phi)
  // alone, so it is not a projection of the 3D cube above; it needs its own fill from each pair's
  // real 4-momenta, in parallel with the 3D fills below.
  auto make1DM2Set = [&](const std::string& tag, TH1D*& h2body, TH1D*& huncorr, TH1D*& hconn) {
    std::string name2 = "h1_2body_" + tag + "_M2";
    std::string nameU = "h1_uncorr_" + tag + "_M2";
    std::string nameC = "h1_conn_" + tag + "_M2";
    std::string title = ";M^{2} [(GeV/c^{2})^{2}]";

    h2body  = new TH1D(name2.c_str(), ("Physical 2-Body " + title).c_str(),  nBinsM2, 0., m2Max);
    huncorr = new TH1D(nameU.c_str(), ("Uncorrelated " + title).c_str(),    nBinsM2, 0., m2Max);
    hconn   = new TH1D(nameC.c_str(), ("Connected C_{conn} " + title).c_str(), nBinsM2, 0., m2Max);
  };

  TH1D *h1_2body_M2_ss1, *h1_uncorr_M2_ss1, *h1_conn_M2_ss1;
  TH1D *h1_2body_M2_s,   *h1_uncorr_M2_s,   *h1_conn_M2_s;
  TH1D *h1_2body_M2_s1,  *h1_uncorr_M2_s1,  *h1_conn_M2_s1;

  make1DM2Set("ss1", h1_2body_M2_ss1, h1_uncorr_M2_ss1, h1_conn_M2_ss1);
  make1DM2Set("s",   h1_2body_M2_s,   h1_uncorr_M2_s,   h1_conn_M2_s);
  make1DM2Set("s1",  h1_2body_M2_s1,  h1_uncorr_M2_s1,  h1_conn_M2_s1);

  // M^2 of a pair from raw 4-momenta: (E1+E2)^2 - (px1+px2)^2 - (py1+py2)^2 - (pz1+pz2)^2
  auto invMass2 = [](double e1, double px1, double py1, double pz1,
                     double e2, double px2, double py2, double pz2) {
    double E  = e1 + e2, Px = px1 + px2, Py = py1 + py2, Pz = pz1 + pz2;
    return E * E - Px * Px - Py * Py - Pz * Pz;
  };

  // eps1 = p1.u, eps2 = p2.u: the comoving energy of each paired particle in the local fluid rest
  // frame (already precomputed per particle at cache time — see CachedParticle::eps). Unlike the
  // other observables this is a single-particle quantity, not a pair difference, so there is no
  // "uncorrelated"/"connected" counterpart — just one plain histogram per role per category.
  auto make1DEpsSet = [&](const std::string& tag, TH1D*& hEps1, TH1D*& hEps2) {
    std::string title = ";#epsilon [GeV]";
    hEps1 = new TH1D(("h1_eps1_" + tag).c_str(), ("Comoving energy #epsilon_{1} " + title).c_str(), nBinsEps, 0., epsMax);
    hEps2 = new TH1D(("h1_eps2_" + tag).c_str(), ("Comoving energy #epsilon_{2} " + title).c_str(), nBinsEps, 0., epsMax);
  };

  TH1D *h1_eps1_ss1, *h1_eps2_ss1;
  TH1D *h1_eps1_s,   *h1_eps2_s;
  TH1D *h1_eps1_s1,  *h1_eps2_s1;

  make1DEpsSet("ss1", h1_eps1_ss1, h1_eps2_ss1);
  make1DEpsSet("s",   h1_eps1_s,   h1_eps2_s);
  make1DEpsSet("s1",  h1_eps1_s1,  h1_eps2_s1);

  // --- Pair Filling Loop ---
  for (const auto& evtParticles : eventStore) {
    std::vector<CachedParticle> posList;
    std::vector<CachedParticle> negList;

    // Filter using direct PDG comparisons (No ParticleDB lookups!)
    for (const auto& p : evtParticles) {
      if (cfg.pid1 != 0) {
        if (p.pid == cfg.pid1)      posList.push_back(p);
        else if (p.pid == cfg.pid2) negList.push_back(p);
      } else {
        // Inclusive charged hadrons
        if (p.pid == 211 || p.pid == 321 || p.pid == 2212) {
          posList.push_back(p);
        } else if (p.pid == -211 || p.pid == -321 || p.pid == -2212) {
          negList.push_back(p);
        }
      }
    }

    // --- Opposite charge pairs (s, s') ---
    for (const auto& p1 : posList) {
      for (const auto& p2 : negList) {
        double dy   = std::abs(p1.y - p2.y);
        // double dpt  = std::abs(p1.pt - p2.pt);       // only needed by the 3D histogram (below)
        // double dphi = getDeltaPhi(p1.phi, p2.phi);   // only needed by the 3D histogram (below)

        if (dy > delta) {
          // h3_2body_ss1->Fill(dy, dpt, dphi);
          double m2 = invMass2(p1.e, p1.px, p1.py, p1.pz, p2.e, p2.px, p2.py, p2.pz);
          h1_2body_M2_ss1->Fill(m2);
          h1_eps1_ss1->Fill(p1.eps);
          h1_eps2_ss1->Fill(p2.eps);
        }

        // Randomize only the rapidity (paper's Step 2); keep the real Delta pT/Delta phi.
        double y1_rand = gRandom->Uniform(0.0, W);
        double y2_rand = gRandom->Uniform(0.0, W);
        double dy_rand = std::abs(y1_rand - y2_rand);
        if (dy_rand > delta) {
          // h3_uncorr_ss1->Fill(dy_rand, dpt, dphi, pairWeight);

          // Reconstruct e/pz consistent with the randomized rapidity, keeping px/py (hence pT,
          // phi, and mass) fixed at their real values, so M^2 reflects the same randomization.
          double e1_rand  = p1.mt * std::cosh(y1_rand);
          double pz1_rand = p1.mt * std::sinh(y1_rand);
          double e2_rand  = p2.mt * std::cosh(y2_rand);
          double pz2_rand = p2.mt * std::sinh(y2_rand);
          double m2_rand  = invMass2(e1_rand, p1.px, p1.py, pz1_rand, e2_rand, p2.px, p2.py, pz2_rand);
          h1_uncorr_M2_ss1->Fill(m2_rand, pairWeight);
        }
      }
    }

    // --- Same charge pairs (s, s) ---
    for (size_t i = 0; i < posList.size(); ++i) {
      for (size_t j = i + 1; j < posList.size(); ++j) {
        double dy   = std::abs(posList[i].y - posList[j].y);
        // double dpt  = std::abs(posList[i].pt - posList[j].pt);          // 3D-histogram only
        // double dphi = getDeltaPhi(posList[i].phi, posList[j].phi);      // 3D-histogram only

        if (dy > delta) {
          // h3_2body_s->Fill(dy, dpt, dphi, 2.0);
          double m2 = invMass2(posList[i].e, posList[i].px, posList[i].py, posList[i].pz,
                                posList[j].e, posList[j].px, posList[j].py, posList[j].pz);
          h1_2body_M2_s->Fill(m2, 2.0);
          h1_eps1_s->Fill(posList[i].eps, 2.0);
          h1_eps2_s->Fill(posList[j].eps, 2.0);
        }

        double y1_rand = gRandom->Uniform(0.0, W);
        double y2_rand = gRandom->Uniform(0.0, W);
        double dy_rand = std::abs(y1_rand - y2_rand);
        if (dy_rand > delta) {
          // h3_uncorr_s->Fill(dy_rand, dpt, dphi, 2.0 * pairWeight);

          double e1_rand  = posList[i].mt * std::cosh(y1_rand);
          double pz1_rand = posList[i].mt * std::sinh(y1_rand);
          double e2_rand  = posList[j].mt * std::cosh(y2_rand);
          double pz2_rand = posList[j].mt * std::sinh(y2_rand);
          double m2_rand  = invMass2(e1_rand, posList[i].px, posList[i].py, pz1_rand,
                                     e2_rand, posList[j].px, posList[j].py, pz2_rand);
          h1_uncorr_M2_s->Fill(m2_rand, 2.0 * pairWeight);
        }
      }
    }

    // --- Same charge pairs (s', s') ---
    for (size_t i = 0; i < negList.size(); ++i) {
      for (size_t j = i + 1; j < negList.size(); ++j) {
        double dy   = std::abs(negList[i].y - negList[j].y);
        // double dpt  = std::abs(negList[i].pt - negList[j].pt);          // 3D-histogram only
        // double dphi = getDeltaPhi(negList[i].phi, negList[j].phi);      // 3D-histogram only

        if (dy > delta) {
          // h3_2body_s1->Fill(dy, dpt, dphi, 2.0);
          double m2 = invMass2(negList[i].e, negList[i].px, negList[i].py, negList[i].pz,
                                negList[j].e, negList[j].px, negList[j].py, negList[j].pz);
          h1_2body_M2_s1->Fill(m2, 2.0);
          h1_eps1_s1->Fill(negList[i].eps, 2.0);
          h1_eps2_s1->Fill(negList[j].eps, 2.0);
        }

        double y1_rand = gRandom->Uniform(0.0, W);
        double y2_rand = gRandom->Uniform(0.0, W);
        double dy_rand = std::abs(y1_rand - y2_rand);
        if (dy_rand > delta) {
          // h3_uncorr_s1->Fill(dy_rand, dpt, dphi, 2.0 * pairWeight);

          double e1_rand  = negList[i].mt * std::cosh(y1_rand);
          double pz1_rand = negList[i].mt * std::sinh(y1_rand);
          double e2_rand  = negList[j].mt * std::cosh(y2_rand);
          double pz2_rand = negList[j].mt * std::sinh(y2_rand);
          double m2_rand  = invMass2(e1_rand, negList[i].px, negList[i].py, pz1_rand,
                                     e2_rand, negList[j].px, negList[j].py, pz2_rand);
          h1_uncorr_M2_s1->Fill(m2_rand, 2.0 * pairWeight);
        }
      }
    }
  }

  // --- Compute C_conn = 2body - uncorr (Eq. 2), bin-wise on the full 3D cube --- (commented out
  // along with the 3D histograms themselves, see note above)
  // auto subtract3D = [](TH3D* h2, TH3D* hu, TH3D* hc) {
  //   hc->Add(h2, hu, 1.0, -1.0);
  // };
  //
  // subtract3D(h3_2body_ss1, h3_uncorr_ss1, h3_conn_ss1);
  // subtract3D(h3_2body_s,   h3_uncorr_s,   h3_conn_s);
  // subtract3D(h3_2body_s1,  h3_uncorr_s1,  h3_conn_s1);

  auto subtract1D = [](TH1D* h2, TH1D* hu, TH1D* hc) {
    hc->Add(h2, hu, 1.0, -1.0);
  };

  subtract1D(h1_2body_M2_ss1, h1_uncorr_M2_ss1, h1_conn_M2_ss1);
  subtract1D(h1_2body_M2_s,   h1_uncorr_M2_s,   h1_conn_M2_s);
  subtract1D(h1_2body_M2_s1,  h1_uncorr_M2_s1,  h1_conn_M2_s1);

  // --- Derive the 2D/1D projections of the 3D cube --- (commented out along with the 3D
  // histograms themselves; see note above make3DSet for how to re-enable)
  // auto project2D = [](TH3D* h3, const std::string& tag) {
  //   TH2D* hPt = (TH2D*) h3->Project3D("yx");
  //   hPt->SetNameTitle(("h2_" + tag + "_pt").c_str(), h3->GetTitle());
  //
  //   TH2D* hPhi = (TH2D*) h3->Project3D("zx");
  //   hPhi->SetNameTitle(("h2_" + tag + "_phi").c_str(), h3->GetTitle());
  // };
  //
  // auto project1D = [](TH3D* h3, const std::string& tag) {
  //   TH1D* hDy = (TH1D*) h3->Project3D("x");
  //   hDy->SetNameTitle(("h1_" + tag + "_dy").c_str(), h3->GetTitle());
  //
  //   TH1D* hPt = (TH1D*) h3->Project3D("y");
  //   hPt->SetNameTitle(("h1_" + tag + "_pt").c_str(), h3->GetTitle());
  //
  //   TH1D* hPhi = (TH1D*) h3->Project3D("z");
  //   hPhi->SetNameTitle(("h1_" + tag + "_phi").c_str(), h3->GetTitle());
  // };
  //
  // project2D(h3_2body_ss1,  "2body_ss1");
  // project2D(h3_uncorr_ss1, "uncorr_ss1");
  // project2D(h3_conn_ss1,   "conn_ss1");
  //
  // project2D(h3_2body_s,  "2body_s");
  // project2D(h3_uncorr_s, "uncorr_s");
  // project2D(h3_conn_s,   "conn_s");
  //
  // project2D(h3_2body_s1,  "2body_s1");
  // project2D(h3_uncorr_s1, "uncorr_s1");
  // project2D(h3_conn_s1,   "conn_s1");
  //
  // project1D(h3_2body_ss1,  "2body_ss1");
  // project1D(h3_uncorr_ss1, "uncorr_ss1");
  // project1D(h3_conn_ss1,   "conn_ss1");
  //
  // project1D(h3_2body_s,  "2body_s");
  // project1D(h3_uncorr_s, "uncorr_s");
  // project1D(h3_conn_s,   "conn_s");
  //
  // project1D(h3_2body_s1,  "2body_s1");
  // project1D(h3_uncorr_s1, "uncorr_s1");
  // project1D(h3_conn_s1,   "conn_s1");

  fout.Write();
  fout.Close();
  std::cout << "  [+] Output written to " << filename << "\n";
}

// Fills currentEvent with the CachedParticle entries for one final-state particle, applying the
// same |y|<=Y acceptance and eps reconstruction used by both the generate and read-from-disk
// paths below. `p` may be a live Particle (in-process) or a raw ParticleCoor read from a tree.
static void cacheParticleIfAccepted(const ParticleCoor& p, bool decayed, double Y, bool isBlastWave,
                                    double vT, std::vector<CachedParticle>& currentEvent) {
  if (decayed) return;
  double y = getRapidity(p);
  if (std::abs(y) > Y) return;

  double pt  = std::sqrt(p.px * p.px + p.py * p.py);
  double phi = std::atan2(p.py, p.px);
  if (phi < 0.0) phi += 2.0 * M_PI;
  double mt  = std::sqrt(p.e * p.e - p.pz * p.pz);

  // eps = p.u, the comoving energy in the local fluid rest frame. HRG: u=(1,0,0,0) everywhere, so
  // eps is just the lab energy. BlastWave: reconstruct the emission/decay point's spatial
  // rapidity/azimuthal angle from the stored position and mirror the p.u formula in
  // Model_BlastWave::GetIntegrand (Model_BlastWave.cxx).
  double eps;
  if (isBlastWave) {
    double rapS   = getSpatialRapidity(p);
    double phiS   = std::atan2(p.y, p.x);
    double gammaT = 1.0 / std::sqrt(1.0 - vT * vT);
    eps = gammaT * (mt * std::cosh(rapS - y) - vT * pt * std::cos(phiS - phi));
  } else {
    eps = p.e;
  }

  // Pure cached kinematic properties — zero ParticleDB lookup overhead
  currentEvent.push_back({y, pt, phi, p.px, p.py, p.pz, p.e, mt, eps, p.pid});
}

// Reads pre-generated events back from event_*.root files written by therm2_events.exe, instead
// of regenerating them in-process -- avoids paying the ~seconds/event Monte-Carlo sampling cost
// again for a re-analysis. Returns false (with a message already printed) on any fatal error.
//
// Position (t,x,y,z) in these files is stored in fm (EventGenerator::SaveAsRoot scales it by
// kHbarC at save time) while momentum (e,px,py,pz) is stored in GeV as generated. This does NOT
// need correcting here: getSpatialRapidity() only ever uses the ratio (t+z)/(t-z), which is
// invariant under an overall rescaling of t and z by the same factor.
static bool readEventsFromDisk(const std::string& inputDir, double Y, int nEvents,
                               std::vector<std::vector<CachedParticle>>& eventStore) {
  namespace fs = std::filesystem;

  std::error_code tEc;
  if (!fs::is_directory(inputDir, tEc)) {
    std::cerr << "Error: --input-dir \"" << inputDir << "\" is not a directory\n";
    return false;
  }

  std::vector<std::string> files;
  for (const auto& entry : fs::directory_iterator(inputDir, tEc)) {
    std::string name = entry.path().filename().string();
    if (name.rfind("event_", 0) == 0 && entry.path().extension() == ".root") {
      files.push_back(entry.path().string());
    }
  }
  std::sort(files.begin(), files.end());
  if (files.empty()) {
    std::cerr << "Error: no event_*.root files found in " << inputDir << "\n";
    return false;
  }

  // Recover the model name and (for BlastWave) the transverse flow velocity from the first file's
  // "parameters" tree (written by EventGenerator::SaveAsRoot / Model::AddParameterBranch) -- there
  // is no live Model object to introspect in this mode, unlike the in-process path.
  bool isBlastWave = false;
  double vT = 0.0;
  {
    std::unique_ptr<TFile> f(TFile::Open(files.front().c_str()));
    if (!f || f->IsZombie()) {
      std::cerr << "Error: cannot open " << files.front() << "\n";
      return false;
    }
    TTree* ptree = (TTree*) f->Get(_PARAMETERS_TREE_);
    if (!ptree) {
      std::cerr << "Error: no \"" << _PARAMETERS_TREE_ << "\" tree in " << files.front() << "\n";
      return false;
    }

    char modelNameBuf[30] = {0};
    ptree->SetBranchAddress(_MODELNAME_BRANCH_, modelNameBuf);
    ptree->GetEntry(0);
    isBlastWave = (std::string(modelNameBuf) == "Blast-Wave"); // literal Model_BlastWave::mName

    if (isBlastWave) {
      Model_t_BlastWave bwPar{};
      ptree->SetBranchAddress(_MODEL_T_BRANCH_, &bwPar);
      ptree->GetEntry(0);
      vT = bwPar.Vt;
    }
  }

  TChain chain(_PARTICLES_TREE_);
  for (const auto& fn : files) chain.Add(fn.c_str());

  ParticleCoor partCoor;
  chain.SetBranchAddress(_PARTICLE_BRANCH_, &partCoor);

  Long64_t nAvailable = chain.GetEntries();
  std::cout << "Reading up to " << nEvents << " event(s) from " << nAvailable
           << " particle entries across " << files.size() << " file(s) in " << inputDir
           << " (model: " << (isBlastWave ? "Blast-Wave" : "other") << ")...\n";

  std::vector<CachedParticle> currentEvent;
  bool inEvent = false;
  UInt_t currentEventId = 0;

  for (Long64_t i = 0; i < nAvailable; i++) {
    chain.GetEntry(i);

    if (inEvent && partCoor.eventid != currentEventId) {
      eventStore.push_back(std::move(currentEvent));
      currentEvent.clear();
      inEvent = false;
      if (static_cast<int>(eventStore.size()) >= nEvents) break;
    }
    if (!inEvent) {
      currentEventId = partCoor.eventid;
      inEvent = true;
    }

    cacheParticleIfAccepted(partCoor, partCoor.decayed, Y, isBlastWave, vT, currentEvent);
  }
  if (inEvent && static_cast<int>(eventStore.size()) < nEvents) {
    eventStore.push_back(std::move(currentEvent));
  }

  return true;
}

int main(int argc, char** argv) {
  std::string shareDir  = "";
  std::string modelName = "hrg";
  std::string iniFile   = "";
  std::string inputDir  = "";
  int nEvents           = 500;
  int nSamples          = 1000000;
  int seed              = 0;
  double Y              = 1.5;
  double delta          = 1.0;
  double ptMax          = 3.0; // GeV/c
  double m2Max          = 6.0; // (GeV/c^2)^2
  double epsMax         = 3.0; // GeV

  CLI::App app{"ex_connected_correlations_2d: 2D Differential Connected Correlations"};
  app.add_option("--share",     shareDir,  "path to SHARE particle database");
  app.add_option("--model",     modelName, "model: hrg (default) or blastwave")->default_str("hrg");
  app.add_option("--ini",       iniFile,   "ini file (required for blastwave model)");
  app.add_option("--input-dir", inputDir,  "read pre-generated event_*.root files from this directory "
                                          "(from therm2_events.exe) instead of generating in-process; "
                                          "mutually exclusive with --model/--ini/--samples");
  app.add_option("--events",  nEvents,   "number of events to generate, or to read (as a cap) when --input-dir is given")->default_val(500);
  app.add_option("--samples", nSamples,  "Monte Carlo integration samples")->default_val(1000000);
  app.add_option("--seed",    seed,      "random seed")->default_val(0);
  app.add_option("--ymax",    Y,         "rapidity window parameter Y (y in [-Y, Y])")->default_val(2.0);
  app.add_option("--delta",   delta,     "minimum rapidity separation cut Delta y > delta")->default_val(0.0);
  app.add_option("--ptmax",   ptMax,     "maximum Delta pT for histogram range")->default_val(3.0);
  app.add_option("--m2max",   m2Max,     "maximum pair invariant mass^2 for histogram range (GeV^2)")->default_val(6.0);
  app.add_option("--epsmax",  epsMax,    "maximum comoving energy eps=p.u for histogram range (GeV)")->default_val(3.0);
  CLI11_PARSE(app, argc, argv);

  if (delta >= Y) {
    std::cerr << "Error: delta (" << delta << ") must be strictly less than Y (" << Y << ")\n";
    return 1;
  }

  bool readMode = !inputDir.empty();
  if (readMode && (app.count("--model") || app.count("--ini") || app.count("--samples"))) {
    std::cerr << "Error: --input-dir cannot be combined with --model/--ini/--samples "
                 "(the model is whatever generated the input files)\n";
    return 1;
  }

  gRandom->SetSeed(seed);

  std::vector<std::vector<CachedParticle>> eventStore;

  if (readMode) {
    if (!readEventsFromDisk(inputDir, Y, nEvents, eventStore)) {
      return 1;
    }
  } else {
    std::unique_ptr<Model> model;
    bool isBlastWave = (modelName == "blastwave");
    if (isBlastWave) {
      if (iniFile.empty()) {
        std::cerr << "Error: --ini is required for the blastwave model\n";
        return 1;
      }
      model = std::make_unique<Model_BlastWave>(iniFile.c_str());
    } else {
      model = std::make_unique<Model_HRG>();
    }

    // Transverse flow velocity used to reconstruct the fluid four-velocity u^mu at each particle's
    // emission/decay point (see getSpatialRapidity() and cacheParticleIfAccepted()). Zero/unused
    // for HRG, whose u^mu = (1,0,0,0) everywhere regardless of position.
    double vT = 0.0;
    if (isBlastWave) {
      vT = dynamic_cast<Model_BlastWave*>(model.get())->GetVt();
    }

    auto tPartDB = std::make_unique<ParticleDB>(shareDir);
    auto tInteg  = std::make_unique<Integrator>(nSamples, model.get());
    tInteg->SetMultiplicities(tPartDB.get());
    auto tEvent  = std::make_unique<Event>(tPartDB.get(), tInteg.get());

    eventStore.reserve(nEvents);

    std::cout << "Generating " << nEvents << " events...\n";

    for (int ievent = 0; ievent < nEvents; ievent++) {
      tEvent->Reset();
      tEvent->GeneratePrimordials();
      tEvent->DecayParticles();

      std::vector<CachedParticle> currentEvent;
      for (const auto& p : tEvent->GetParticleList()) {
        cacheParticleIfAccepted(p, p.GetDecayed(), Y, isBlastWave, vT, currentEvent);
      }
      eventStore.push_back(currentEvent);
    }
  }

  std::vector<SpeciesPairConfig> speciesList = {
    {"pi_pi",   211, -211},
    {"K_K",     321, -321},
    {"p_p",    2212,-2212},
    // {"all_all",   0,    0} // commented out: O(n^2) over all charged hadrons is by far the most
                              // expensive config — re-enable only when the inclusive correlation
                              // is actually needed.
  };

  std::cout << "Building 2D connected correlation matrices...\n";
  for (const auto& cfg : speciesList) {
    processSpeciesPair2D(cfg, eventStore, Y, delta, 40, 30, 32, 50, 50, ptMax, m2Max, epsMax);
  }

  std::cout << "Finished successfully.\n";
  return 0;
}