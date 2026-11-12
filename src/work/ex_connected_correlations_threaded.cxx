// ex_connected_correlations_threaded.cxx — thread-parallel twin of ex_connected_correlations_2d.cxx.
// Same physics/algorithm (connected two-particle correlations, following Eq. (2) of
// Short_range_hadron_correlations.pdf: C_ss',conn = <2-body> - <uncorrelated product>, background
// estimated per Section II Eqs. 5-8 by randomizing pair rapidity over the accepted window). See
// ex_connected_correlations_2d.cxx for the full derivation/physics comments; they are not repeated
// here except where the threading design changes something physics-relevant (the RNG used for the
// "uncorrelated" background).
//
// THREADING DESIGN: event generation and disk-reading (building `eventStore`) stay exactly
// sequential/single-threaded, unchanged from ex_connected_correlations_2d.cxx -- parallelizing
// those would require touching shared library code (ROOT's global gRandom is reached by ~25 call
// sites across Integrator/Model_HRG/Model_BlastWave/ParticleDecayer, and Event::Reset() zeroes a
// process-wide Particle::EIDi counter every event), which is out of scope here.
//
// What IS parallelized is the per-species-pair analysis phase inside processSpeciesPair2D: for
// each species pair, `eventStore` is split into up to --threads contiguous, roughly-equal chunks
// (event-based chunking). Each chunk is filled by its own std::thread into its own private
// histogram set (PerThreadHistos2D) using its own private TRandom3 (never the shared gRandom) for
// the rapidity-randomized background term -- no shared mutable state is touched during the
// parallel section, so no locks are needed in the hot O(events x particles^2) pairing loop. After
// all threads join, the per-thread histogram sets are merged into thread 0's via TH1::Add (safe:
// join() establishes happens-before), the connected (2body - uncorr) subtraction runs once on the
// merged totals, and the result is written out exactly as in the sequential version.
//
// RNG reproducibility: with a fixed --threads N, every run is fully deterministic (chunk
// boundaries are a pure function of eventStore.size() and N; each chunk's TRandom3 is seeded
// deterministically from --seed + species index + chunk index; fill order within a chunk is
// single-threaded; merge order is sequential). Runs at *different* --threads values are NOT
// bit-for-bit identical to each other or to ex_connected_correlations_2d.cxx's gRandom-based
// stream (chunk boundaries shift with N, so identical per-chunk seeds land on different event
// subsets) -- they agree only statistically. The real (non-randomized) "2body" and eps histograms
// have no RNG dependence at all and match ex_connected_correlations_2d.cxx bit-for-bit.
//
// Usage: identical CLI to ex_connected_correlations_2d.cxx, plus:
//   --threads <N>   worker threads per species pair for histogram filling (default: detected core
//                    count)

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
#include <TRandom3.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

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
  int idx;  // stable numeric identity for this species pair, independent of speciesList's vector
            // order -- used to decorrelate per-thread RNG seeds across species (see chunkSeed()).
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

// --- Per-thread histogram set + private RNG for one chunk of eventStore -------------------------
//
// Deliberately named members (mirroring make1DM2Set/make1DEpsSet in ex_connected_correlations_2d.cxx)
// rather than a vector -- keeps the diff against the sequential version easy to follow. Re-enabling
// the 3D (Delta y, Delta pT, Delta phi) histograms (see the commented-out block below, carried over
// unchanged from ex_connected_correlations_2d.cxx) would mean adding the corresponding TH3D members
// here, filling them in fillChunk2D exactly like the 1D members, and adding their Add() calls to
// AddInto() -- TH3D::Add has identical semantics to TH1D::Add (both derive from TH1).
//
// unique_ptr members give RAII cleanup for the N-1 non-master chunks after they're merged and
// discarded, instead of manual delete.
struct PerThreadHistos2D {
  std::unique_ptr<TH1D> h1_2body_M2_ss1, h1_uncorr_M2_ss1;
  std::unique_ptr<TH1D> h1_2body_M2_s,   h1_uncorr_M2_s;
  std::unique_ptr<TH1D> h1_2body_M2_s1,  h1_uncorr_M2_s1;
  std::unique_ptr<TH1D> h1_eps1_ss1, h1_eps2_ss1;
  std::unique_ptr<TH1D> h1_eps1_s,   h1_eps2_s;
  std::unique_ptr<TH1D> h1_eps1_s1,  h1_eps2_s1;
  std::unique_ptr<TRandom3> rng;

  // Adds this chunk's histograms into dst's (bin contents + Sumw2 errors). Called only after every
  // worker thread has join()'d, so no locking is needed -- join() establishes happens-before.
  void AddInto(PerThreadHistos2D& dst) const {
    dst.h1_2body_M2_ss1->Add(h1_2body_M2_ss1.get());
    dst.h1_uncorr_M2_ss1->Add(h1_uncorr_M2_ss1.get());
    dst.h1_2body_M2_s->Add(h1_2body_M2_s.get());
    dst.h1_uncorr_M2_s->Add(h1_uncorr_M2_s.get());
    dst.h1_2body_M2_s1->Add(h1_2body_M2_s1.get());
    dst.h1_uncorr_M2_s1->Add(h1_uncorr_M2_s1.get());
    dst.h1_eps1_ss1->Add(h1_eps1_ss1.get());
    dst.h1_eps2_ss1->Add(h1_eps2_ss1.get());
    dst.h1_eps1_s->Add(h1_eps1_s.get());
    dst.h1_eps2_s->Add(h1_eps2_s.get());
    dst.h1_eps1_s1->Add(h1_eps1_s1.get());
    dst.h1_eps2_s1->Add(h1_eps2_s1.get());
  }
};

// Builds one chunk's private histogram set + private RNG. Relies on TH1::AddDirectory(kFALSE)
// (set once, globally, in main() before any thread runs this) so these `new TH1D(...)` calls never
// race on ROOT's shared, non-thread-safe gDirectory. Histogram names follow the same convention as
// ex_connected_correlations_2d.cxx (species-scoped, not chunk-scoped) -- safe because chunk sets
// are never simultaneously registered in any TDirectory, and only the merged thread-0 set is ever
// written to file.
static PerThreadHistos2D makePerThreadHistos2D(int nBinsM2, double m2Max, int nBinsEps,
                                               double epsMax, unsigned int rngSeed) {
  PerThreadHistos2D h;

  std::string titleM2  = ";M^{2} [(GeV/c^{2})^{2}]";
  h.h1_2body_M2_ss1  = std::make_unique<TH1D>("h1_2body_ss1_M2",  ("Physical 2-Body " + titleM2).c_str(),  nBinsM2, 0., m2Max);
  h.h1_uncorr_M2_ss1 = std::make_unique<TH1D>("h1_uncorr_ss1_M2", ("Uncorrelated " + titleM2).c_str(),     nBinsM2, 0., m2Max);
  h.h1_2body_M2_s    = std::make_unique<TH1D>("h1_2body_s_M2",    ("Physical 2-Body " + titleM2).c_str(),  nBinsM2, 0., m2Max);
  h.h1_uncorr_M2_s   = std::make_unique<TH1D>("h1_uncorr_s_M2",   ("Uncorrelated " + titleM2).c_str(),     nBinsM2, 0., m2Max);
  h.h1_2body_M2_s1   = std::make_unique<TH1D>("h1_2body_s1_M2",   ("Physical 2-Body " + titleM2).c_str(),  nBinsM2, 0., m2Max);
  h.h1_uncorr_M2_s1  = std::make_unique<TH1D>("h1_uncorr_s1_M2",  ("Uncorrelated " + titleM2).c_str(),     nBinsM2, 0., m2Max);

  std::string titleEps = ";#epsilon [GeV]";
  h.h1_eps1_ss1 = std::make_unique<TH1D>("h1_eps1_ss1", ("Comoving energy #epsilon_{1} " + titleEps).c_str(), nBinsEps, 0., epsMax);
  h.h1_eps2_ss1 = std::make_unique<TH1D>("h1_eps2_ss1", ("Comoving energy #epsilon_{2} " + titleEps).c_str(), nBinsEps, 0., epsMax);
  h.h1_eps1_s   = std::make_unique<TH1D>("h1_eps1_s",   ("Comoving energy #epsilon_{1} " + titleEps).c_str(), nBinsEps, 0., epsMax);
  h.h1_eps2_s   = std::make_unique<TH1D>("h1_eps2_s",   ("Comoving energy #epsilon_{2} " + titleEps).c_str(), nBinsEps, 0., epsMax);
  h.h1_eps1_s1  = std::make_unique<TH1D>("h1_eps1_s1",  ("Comoving energy #epsilon_{1} " + titleEps).c_str(), nBinsEps, 0., epsMax);
  h.h1_eps2_s1  = std::make_unique<TH1D>("h1_eps2_s1",  ("Comoving energy #epsilon_{2} " + titleEps).c_str(), nBinsEps, 0., epsMax);

  h.rng = std::make_unique<TRandom3>(rngSeed);
  return h;
}

// M^2 of a pair from raw 4-momenta: (E1+E2)^2 - (px1+px2)^2 - (py1+py2)^2 - (pz1+pz2)^2
static double invMass2(double e1, double px1, double py1, double pz1,
                       double e2, double px2, double py2, double pz2) {
  double E  = e1 + e2, Px = px1 + px2, Py = py1 + py2, Pz = pz1 + pz2;
  return E * E - Px * Px - Py * Py - Pz * Pz;
}

// Fills one contiguous chunk [begin, end) of eventStore into `out`. Touches only `out` (private to
// the calling thread) and read-only shared data (cfg, eventStore) -- safe to run concurrently
// across threads with zero synchronization. Body is identical to the "Pair Filling Loop" in
// ex_connected_correlations_2d.cxx, except it uses out.rng instead of the shared gRandom for the
// rapidity-randomized background term.
static void fillChunk2D(const SpeciesPairConfig& cfg,
                        const std::vector<std::vector<CachedParticle>>& eventStore,
                        size_t begin, size_t end,
                        double delta, double W, double pairWeight,
                        PerThreadHistos2D& out) {
  for (size_t ie = begin; ie < end; ++ie) {
    const auto& evtParticles = eventStore[ie];
    std::vector<CachedParticle> posList;
    std::vector<CachedParticle> negList;

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
        double dy = std::abs(p1.y - p2.y);

        if (dy > delta) {
          double m2 = invMass2(p1.e, p1.px, p1.py, p1.pz, p2.e, p2.px, p2.py, p2.pz);
          out.h1_2body_M2_ss1->Fill(m2);
          out.h1_eps1_ss1->Fill(p1.eps);
          out.h1_eps2_ss1->Fill(p2.eps);
        }

        // Randomize only the rapidity (paper's Step 2); keep the real Delta pT/Delta phi.
        double y1_rand = out.rng->Uniform(0.0, W);
        double y2_rand = out.rng->Uniform(0.0, W);
        double dy_rand = std::abs(y1_rand - y2_rand);
        if (dy_rand > delta) {
          double e1_rand  = p1.mt * std::cosh(y1_rand);
          double pz1_rand = p1.mt * std::sinh(y1_rand);
          double e2_rand  = p2.mt * std::cosh(y2_rand);
          double pz2_rand = p2.mt * std::sinh(y2_rand);
          double m2_rand  = invMass2(e1_rand, p1.px, p1.py, pz1_rand, e2_rand, p2.px, p2.py, pz2_rand);
          out.h1_uncorr_M2_ss1->Fill(m2_rand, pairWeight);
        }
      }
    }

    // --- Same charge pairs (s, s) ---
    for (size_t i = 0; i < posList.size(); ++i) {
      for (size_t j = i + 1; j < posList.size(); ++j) {
        double dy = std::abs(posList[i].y - posList[j].y);

        if (dy > delta) {
          double m2 = invMass2(posList[i].e, posList[i].px, posList[i].py, posList[i].pz,
                                posList[j].e, posList[j].px, posList[j].py, posList[j].pz);
          out.h1_2body_M2_s->Fill(m2, 2.0);
          out.h1_eps1_s->Fill(posList[i].eps, 2.0);
          out.h1_eps2_s->Fill(posList[j].eps, 2.0);
        }

        double y1_rand = out.rng->Uniform(0.0, W);
        double y2_rand = out.rng->Uniform(0.0, W);
        double dy_rand = std::abs(y1_rand - y2_rand);
        if (dy_rand > delta) {
          double e1_rand  = posList[i].mt * std::cosh(y1_rand);
          double pz1_rand = posList[i].mt * std::sinh(y1_rand);
          double e2_rand  = posList[j].mt * std::cosh(y2_rand);
          double pz2_rand = posList[j].mt * std::sinh(y2_rand);
          double m2_rand  = invMass2(e1_rand, posList[i].px, posList[i].py, pz1_rand,
                                     e2_rand, posList[j].px, posList[j].py, pz2_rand);
          out.h1_uncorr_M2_s->Fill(m2_rand, 2.0 * pairWeight);
        }
      }
    }

    // --- Same charge pairs (s', s') ---
    for (size_t i = 0; i < negList.size(); ++i) {
      for (size_t j = i + 1; j < negList.size(); ++j) {
        double dy = std::abs(negList[i].y - negList[j].y);

        if (dy > delta) {
          double m2 = invMass2(negList[i].e, negList[i].px, negList[i].py, negList[i].pz,
                                negList[j].e, negList[j].px, negList[j].py, negList[j].pz);
          out.h1_2body_M2_s1->Fill(m2, 2.0);
          out.h1_eps1_s1->Fill(negList[i].eps, 2.0);
          out.h1_eps2_s1->Fill(negList[j].eps, 2.0);
        }

        double y1_rand = out.rng->Uniform(0.0, W);
        double y2_rand = out.rng->Uniform(0.0, W);
        double dy_rand = std::abs(y1_rand - y2_rand);
        if (dy_rand > delta) {
          double e1_rand  = negList[i].mt * std::cosh(y1_rand);
          double pz1_rand = negList[i].mt * std::sinh(y1_rand);
          double e2_rand  = negList[j].mt * std::cosh(y2_rand);
          double pz2_rand = negList[j].mt * std::sinh(y2_rand);
          double m2_rand  = invMass2(e1_rand, negList[i].px, negList[i].py, pz1_rand,
                                     e2_rand, negList[j].px, negList[j].py, pz2_rand);
          out.h1_uncorr_M2_s1->Fill(m2_rand, 2.0 * pairWeight);
        }
      }
    }
  }
}

void processSpeciesPair2D(const SpeciesPairConfig& cfg,
                          const std::vector<std::vector<CachedParticle>>& eventStore,
                          double Y, double delta,
                          int nBinsY, int nBinsPt, int nBinsPhi, int nBinsM2, int nBinsEps,
                          double ptMax, double m2Max, double epsMax,
                          unsigned int nThreads, unsigned int baseSeed) {

  std::string filename = "corr2d_" + cfg.name + ".root";

  double dyMax = 2.0 * Y;
  // Full width of the accepted |y| <= Y rapidity window — this is the "Y" of Eqs. 5-8 in the
  // paper (there, particles are produced over a window of width Y; here that width is 2*Y since
  // particles are accepted with |y| <= Y). Must be used consistently for the reweight factor and
  // the randomization support below, or the phase-space correction desyncs from the real acceptance.
  const double W = dyMax;
  const double pairWeight = (W * W) / ((W - delta) * (W - delta));

  // --- Chunk eventStore across up to nThreads worker threads (event-based parallelism). Never
  // spawn more threads than events, and never spawn zero threads (an empty eventStore still gets
  // exactly one, trivially-empty chunk, so the rest of the function needs no special-casing).
  unsigned int N = std::max(1u, std::min(nThreads, (unsigned int) eventStore.size()));

  std::vector<size_t> chunkBegin(N), chunkEnd(N);
  {
    size_t total  = eventStore.size();
    size_t base   = total / N;
    size_t rem    = total % N;
    size_t cursor = 0;
    for (unsigned int k = 0; k < N; ++k) {
      size_t sz     = base + (k < rem ? 1 : 0);
      chunkBegin[k] = cursor;
      chunkEnd[k]   = cursor + sz;
      cursor += sz;
    }
  }

  // Per-chunk RNG seed: deterministic for a fixed nThreads (same baseSeed/idx/k always produces
  // the same TRandom3 stream), but NOT reproducible across different --threads values, since chunk
  // boundaries shift with N and so the same stream lands on a different event subset. That's an
  // accepted tradeoff -- solving it (e.g. per-event RNG re-seeding) would be far more expensive
  // than it's worth. The prime multipliers just decorrelate (baseSeed, species, chunk) across the
  // 32-bit seed space; +1u avoids TRandom3's seed==0 "seed from clock" special case.
  auto chunkSeed = [&](unsigned int k) {
    return baseSeed * 1000003u + (unsigned int) cfg.idx * 7919u + k * 104729u + 1u;
  };

  std::vector<PerThreadHistos2D> perThread(N);
  for (unsigned int k = 0; k < N; ++k) {
    perThread[k] = makePerThreadHistos2D(nBinsM2, m2Max, nBinsEps, epsMax, chunkSeed(k));
  }

  // Worker threads handle chunks 1..N-1; the calling thread handles chunk 0 directly (no wasted
  // spawn-then-immediately-join, and it does useful work instead of idling while workers run).
  std::vector<std::thread> workers;
  workers.reserve(N - 1);
  for (unsigned int k = 1; k < N; ++k) {
    workers.emplace_back(fillChunk2D, std::cref(cfg), std::cref(eventStore),
                          chunkBegin[k], chunkEnd[k], delta, W, pairWeight,
                          std::ref(perThread[k]));
  }
  fillChunk2D(cfg, eventStore, chunkBegin[0], chunkEnd[0], delta, W, pairWeight, perThread[0]);
  for (auto& t : workers) t.join();

  // Merge every other chunk's histograms into chunk 0's. Safe with no locking: this runs strictly
  // after every join() above, which establishes happens-before against each worker's writes.
  for (unsigned int k = 1; k < N; ++k) {
    perThread[k].AddInto(perThread[0]);
  }

  // --- Compute C_conn = 2body - uncorr (Eq. 2) on the merged totals -- once, not per-chunk (doing
  // it per-chunk would be wasted work, thrown away at merge time anyway).
  auto subtract1D = [](TH1D* h2, TH1D* hu, TH1D* hc) {
    hc->Add(h2, hu, 1.0, -1.0);
  };

  std::string titleM2 = ";M^{2} [(GeV/c^{2})^{2}]";
  auto h1_conn_M2_ss1 = std::make_unique<TH1D>("h1_conn_ss1_M2", ("Connected C_{conn} " + titleM2).c_str(), nBinsM2, 0., m2Max);
  auto h1_conn_M2_s   = std::make_unique<TH1D>("h1_conn_s_M2",   ("Connected C_{conn} " + titleM2).c_str(), nBinsM2, 0., m2Max);
  auto h1_conn_M2_s1  = std::make_unique<TH1D>("h1_conn_s1_M2",  ("Connected C_{conn} " + titleM2).c_str(), nBinsM2, 0., m2Max);

  subtract1D(perThread[0].h1_2body_M2_ss1.get(), perThread[0].h1_uncorr_M2_ss1.get(), h1_conn_M2_ss1.get());
  subtract1D(perThread[0].h1_2body_M2_s.get(),   perThread[0].h1_uncorr_M2_s.get(),   h1_conn_M2_s.get());
  subtract1D(perThread[0].h1_2body_M2_s1.get(),  perThread[0].h1_uncorr_M2_s1.get(),  h1_conn_M2_s1.get());

  // --- Write output. TH1::AddDirectory(kFALSE) (set once in main(), before any thread constructs
  // a histogram) means none of the histograms above auto-attached to any TDirectory -- write each
  // one explicitly rather than relying on TFile::Write()'s "write everything owned by me" default.
  TFile fout(filename.c_str(), "RECREATE");
  fout.cd();
  perThread[0].h1_2body_M2_ss1->Write();
  perThread[0].h1_uncorr_M2_ss1->Write();
  h1_conn_M2_ss1->Write();
  perThread[0].h1_2body_M2_s->Write();
  perThread[0].h1_uncorr_M2_s->Write();
  h1_conn_M2_s->Write();
  perThread[0].h1_2body_M2_s1->Write();
  perThread[0].h1_uncorr_M2_s1->Write();
  h1_conn_M2_s1->Write();
  perThread[0].h1_eps1_ss1->Write();
  perThread[0].h1_eps2_ss1->Write();
  perThread[0].h1_eps1_s->Write();
  perThread[0].h1_eps2_s->Write();
  perThread[0].h1_eps1_s1->Write();
  perThread[0].h1_eps2_s1->Write();
  fout.Close();
  std::cout << "  [+] Output written to " << filename << " (" << N << " thread"
            << (N == 1 ? "" : "s") << ")\n";
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
// Stays sequential/untouched: its event-boundary detection assumes strictly sequential TChain
// iteration (partCoor.eventid change-detection), and parallelizing it safely would require first
// verifying no single event's particles ever span two event_*.root files -- unverified here, and
// out of scope for this threaded twin, which only parallelizes the analysis phase downstream of
// eventStore. Both ingestion paths (this one and the in-process generation loop in main()) feed
// the same eventStore that processSpeciesPair2D consumes, so both benefit identically either way.
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
  unsigned int hwThreads = std::thread::hardware_concurrency();
  int nThreads           = (hwThreads > 0) ? static_cast<int>(hwThreads) : 1;

  CLI::App app{"ex_connected_correlations_threaded: 2D Differential Connected Correlations (thread-parallel analysis phase)"};
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
  app.add_option("--threads", nThreads,  "worker threads per species pair for histogram filling")->default_val(nThreads);
  CLI11_PARSE(app, argc, argv);

  if (delta >= Y) {
    std::cerr << "Error: delta (" << delta << ") must be strictly less than Y (" << Y << ")\n";
    return 1;
  }
  if (nThreads < 1) {
    std::cerr << "Error: --threads must be >= 1\n";
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
    {"pi_pi",   211, -211, 0},
    {"K_K",     321, -321, 1},
    {"p_p",    2212,-2212, 2},
    // {"all_all",   0,    0, 3} // commented out: O(n^2) over all charged hadrons is by far the
                                  // most expensive config — re-enable only when the inclusive
                                  // correlation is actually needed.
  };

  // Both are process-wide static flags on the TH1 class -- must be set exactly once, single-
  // threaded, before any thread (main or worker) constructs a TH1/TH2/TH3. AddDirectory(kFALSE)
  // is new relative to ex_connected_correlations_2d.cxx: without it, every `new TH1D(...)` inside
  // worker threads would try to register into ROOT's shared, non-thread-safe gDirectory.
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2();

  std::cout << "Building 2D connected correlation matrices (" << nThreads << " thread(s) per species pair)...\n";
  for (const auto& cfg : speciesList) {
    processSpeciesPair2D(cfg, eventStore, Y, delta, 40, 30, 32, 50, 50, ptMax, m2Max, epsMax,
                         static_cast<unsigned int>(nThreads), static_cast<unsigned int>(seed));
  }

  std::cout << "Finished successfully.\n";
  return 0;
}
