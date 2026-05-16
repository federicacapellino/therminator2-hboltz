// clang-format off
/********************************************************************************
*                                                                              *
*             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
*                                                                              *
* Version:                                                                     *
*      Release, 2.0.3, 1 February 2011                                         *
*                                                                              *
* Authors:                                                                     *
*      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
*      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
*      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
*      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
*                                                                              *
* Project homepage:                                                            *
*      http://therminator2.ifj.edu.pl/                                         *
*                                                                              *
* For the detailed description of the program and further references           *
* to the description of the model please refer to                              *
* http://arxiv.org/abs/1102.0273                                               *
*                                                                              *
* This code can be freely used and redistributed. However if you decide to     *
* make modifications to the code, please, inform the authors.                  *
* Any publication of results obtained using this code must include the         *
* reference to arXiv:1102.0273 and the published version of it, when           *
* available.                                                                   *
*                                                                              *
********************************************************************************/
// clang-format on

#include "Configurator.h"
#include "EventGenerator.h"
#include "Model_HRG.h"
#include "Parser.h"
#include "ParticleDB.h"
#include "ParticleDecayer.h"
#include "ParticleType.h"
#include "THGlobal.h"
#include <TH1D.h>
#include <TString.h>
#include <cstdio>
#include <fstream>

#ifndef M_HBARC
#define M_HBARC 0.19732697
#endif
using namespace std;

//! o4_dispersion provides services to compute the modified dispersion curve
//! and modified distribution function.
class o4_dispersion {
public:            // private
  double fT;       //! Temperature
  double fLambda;  //! Cutoff in GeV
  double fMvacuum; //! Vacuum pion mass in GeV
  double fRatio;   //! Ratio of chiral condensate to vacuum one
public:
  o4_dispersion(const double& lambda_by_T = M_PI, const double& ratio = 0.5);
  //! returns the modified distribution function
  double f(const double& p);
  //! returns the vacuum distribution function
  double fvac(const double& p);
  //! returns the modified dispersion curve E(p) in units of GeV
  double eofp(const double& p);
  //! returns the modified pole mass squared in units of GeV**2j
  double m2pole(const double& p);
  //! returns the chiral velocity
  double v2(const double& p);

  //! Makes a plot of the dispersion curve
  void plot(const std::string& filename = "pidispersion.out");
};

o4_dispersion::o4_dispersion(const double& lambda_by_T, const double& ratio)
    : fT(0.155), fLambda(lambda_by_T * 0.155), fMvacuum(0.1396), fRatio(ratio) {}

double o4_dispersion::eofp(const double& p) {
  return sqrt(v2(p) * p * p + m2pole(p));
}

double o4_dispersion::fvac(const double& p) {
  double x = sqrt(fMvacuum * fMvacuum + p * p) / fT;
  return exp(-x) / (1. - exp(-x));
}

double o4_dispersion::f(const double& p) {
  double x = eofp(p) / fT;
  return exp(-x) / (1. - exp(-x));
}

double o4_dispersion::m2pole(const double& p) {
  double mv2 = fMvacuum * fMvacuum;
  double m02 = mv2 * fRatio;
  double x = p * p / (fLambda * fLambda);
  return mv2 - (mv2 - m02) / (1 + x / 2. + x * x);
}

double o4_dispersion::v2(const double& p) {
  double v02 = fRatio * fRatio;
  double x = p * p / (fLambda * fLambda);
  return 1. - (1. - v02) / (1 + x / 2. + x * x);
}

void o4_dispersion::plot(const std::string& filename) {
  int np = 100;
  double pmin = 0.;
  double pmax = 1.5;
  double dp = (pmax - pmin) / (double)np;
  int ip;
  FILE* fp = fopen(filename.c_str(), "w");
  for (ip = 0; ip < np; ip++) {
    double p = pmin + ip * dp;
    fprintf(fp, "%15.5e ", p);
    fprintf(fp, "%15.5e ", v2(p));
    fprintf(fp, "%15.5e ", m2pole(p));
    fprintf(fp, "%15.5e ", eofp(p));
    fprintf(fp, "%15.5e ", f(p));
    fprintf(fp, "%15.5e ", fvac(p));
    fprintf(fp, "\n");
  }
  fclose(fp);
}

////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
  o4_dispersion o4;
  o4.plot();

  auto tMainConfig = std::make_unique<Configurator>("./events.ini");
  tMainConfig->PrintParameters();

  auto tPartDB = std::make_unique<ParticleDB>(tMainConfig->GetParameter("ShareDir", ""));

  // Generate the event multiplicities
  tMainConfig->PrintParameters();
  int tIntegrateSample = tMainConfig->GetParameter("IntegrateSamples").Atoi();

  // Set up the model
  Model_HRG hrg;

  // Initialize the integrator
  auto tInteg = make_unique<Integrator>(tIntegrateSample, &hrg);
  tInteg->SetMultiplicities(tPartDB.get());

  // Create the event structure
  auto tEvent = make_unique<Event>(tPartDB.get(), tInteg.get());

  // Set up root file to store results
  TFile fout("funcs.root", "RECREATE", "a test file");
  TH1::SetDefaultSumw2();

  // Direct pions not modified
  TH1D* hpi = new TH1D("pi_dNdp", "dN/dp", 40, 0.0001, 2.);
  // Direct pions not modified for the simulated case
  TH1D* hpi_stimulated = new TH1D("pi_stimulated_dNdp", "dN/dp", 40, 0.0001, 2.);
  // Direct pions enhanced.
  TH1D* hpi_critical = new TH1D("pi_critical_dNdp", "dN/dp", 40, 0.0001, 2.);

  // All pions not modified
  TH1D* hpi_all = new TH1D("pi_all_dNdp", "dN/dp", 40, 0.0000, 2.);
  // All pions including the stimulated emission
  TH1D* hpi_all_stimulated = new TH1D("pi_all_stimulated_dNdp", "dN/dp", 40, 0.0000, 2.);
  // All pions including the stimulated critical emission
  TH1D* hpi_all_critical = new TH1D("pi_all_critical_dNdp", "dN/dp", 40, 0.0000, 2.);

  // All pions including  weak decays
  TH1D* hpi_all_weak = new TH1D("pi_all_weak_dNdp", "dN/dp", 40, 0.0000, 2.);

  double tcut = 2.5;

  // Vacuum pions from decays
  TH1D* hpi_decayed = new TH1D("pi_decayed_dNdp", "dN/dp", 40, 0.0000, 2.);
  // Pions from stimulated decays
  TH1D* hpi_decayed_stimulated = new TH1D("pi_decayed_stimulated_dNdp", "dN/dp", 40, 0.0000, 2.);
  // Pions from critical decays
  TH1D* hpi_decayed_critical = new TH1D("pi_decayed_critical_dNdp", "dN/dp", 40, 0.0000, 2.);

  // Vacuum pions from short decays
  TH1D* hpi_decayed_short = new TH1D("pi_decayed_short_25_dNdp", "dN/dp", 40, 0.0000, 2.);
  // Vacuum pions from long decays
  TH1D* hpi_decayed_long = new TH1D("pi_decayed_long_25_dNdp", "dN/dp", 40, 0.0000, 2.);

  // Stimulated pions from short lived decays
  TH1D* hpi_decayed_short_stimulated =
      new TH1D("pi_decayed_short_stimulated_25_dNdp", "dN/dp", 40, 0.0000, 2.);
  // Critical pions from short lived decays
  TH1D* hpi_decayed_short_critical =
      new TH1D("pi_decayed_short_critical_25_dNdp", "dN/dp", 40, 0.0000, 2.);

  // Distribution of decay times
  TH1D* htimes = new TH1D("times", "dNdt", 40, 0.0001, 10.);

  int nevents = 500;
  for (int ievent = 0; ievent < nevents; ievent++) {
    tEvent->Reset();
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();

    for (auto ptr = tEvent->GetParticleList().begin(); ptr != tEvent->GetParticleList().end();
         ++ptr) {

      // Write out the first event
      if (ievent == 0) {
        cout << ptr->MakeTEXTEntry() << endl;
      }
      int pid = ptr->GetParticleType()->GetPDGCode();

      // We only keep study the particles that remain
      if (ptr->GetDecayed())
        continue;

      double time = ptr->t * kHbarC;
      double p = ptr->GetP();

      if (abs(pid) == 211) {
        hpi_all_weak->Fill(p);
      }
      if (time > 10000.)
        continue;

      bool primordial = (ptr->fathereid == -1);

      if (abs(pid) == 211) {
        hpi_all->Fill(p);
      }

      // This is a primordial since eid == -1
      if (abs(pid) == 211 && primordial) {
        // Direct pions for the vacuum and stimulated cases
        hpi->Fill(p);
        // Direct pions in the stimulated case
        hpi_stimulated->Fill(p);
        // Direct pions for the critical case. Include the enhancement
        hpi_critical->Fill(p, o4.f(p) / o4.fvac(p));

        // All pions for the simulated case
        hpi_all_stimulated->Fill(p);
        // All pions for the critical case
        hpi_all_critical->Fill(p, o4.f(p) / o4.fvac(p));
      }

      // This is not a primordial, it it is from a decay
      if (abs(pid) == 211 && !primordial) {

        // All decays with vacuum rates
        hpi_decayed->Fill(p);

        // Momentum distribution with short decay times
        if (time < tcut) {
          hpi_decayed_short->Fill(p);
          // Estimate the effect do to stimulated emission with vacuum
          // dispersion curve
          hpi_decayed_short_stimulated->Fill(p, 1. + o4.fvac(p));
          // Estimate the effect do to stimulated emission with vacuum
          // modified curve
          hpi_decayed_short_critical->Fill(p, 1. + o4.f(p));
        } else {
          hpi_decayed_long->Fill(p);
        }
      }
      if (abs(pid) == 211 && ptr->fathereid != -1) {
        htimes->Fill(time);
      }
    }
  }
  // Normalize the histograms so that we don't get more decay pions
  double IVac = hpi_decayed_short->Integral();
  double IStm = hpi_decayed_short_stimulated->Integral();
  double ICrt = hpi_decayed_short_critical->Integral();

  hpi_decayed_short_stimulated->Scale(IVac / IStm);
  hpi_decayed_short_critical->Scale(IVac / ICrt);

  // Compute critical decay component
  hpi_decayed_critical->Add(hpi_decayed_short_critical);
  hpi_decayed_critical->Add(hpi_decayed_long);
  // Add the critical decay component to find the total critical pions
  hpi_all_critical->Add(hpi_decayed_critical);

  // Stimulated decay contributions
  hpi_decayed_stimulated->Add(hpi_decayed_long);
  hpi_decayed_stimulated->Add(hpi_decayed_short_stimulated);
  // Add the decay decay component to find the total decay pions
  hpi_all_stimulated->Add(hpi_decayed_stimulated);

  fout.Write();

  // EventGenerator  *tGenerator = new EventGenerator(tPartDB) ;
  // tGenerator->GenerateEvents();

  return 0;
}
