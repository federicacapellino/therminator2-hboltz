/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2 *
 *                                                                              *
 * Version: * Release, 2.0.3, 1 February 2011 *
 *                                                                              *
 * Authors: * Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl) * Adam Kisiel
 *(kisiel@if.pw.edu.pl)                               * Wojciech Broniowski
 *(Wojciech.Broniowski@ifj.edu.pl)                    * Wojciech Florkowski
 *(Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage: * http://therminator2.ifj.edu.pl/ *
 *                                                                              *
 * For the detailed description of the program and further references * to the
 *description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273 *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to *
 * make modifications to the code, please, inform the authors. * Any publication
 *of results obtained using this code must include the         * reference to
 *arXiv:1102.0273 and the published version of it, when           * available. *
 *                                                                              *
 ********************************************************************************/

#include "Configurator.h"
#include "EventGenerator.h"
#include "Parser.h"
#include "ParticleDB.h"
#include "ParticleDecayer.h"
#include "ParticleType.h"
#include "THGlobal.h"
#include <TH1D.h>
#include <TString.h>
#include <fstream>

extern Configurator *sMainConfig;
extern TString sMainINI;
extern TString sModelINI;
extern TString sHyperXML;
extern TString sEventDIR;
extern TString sTimeStamp;
extern int sModel;
extern int sRandomize;
extern int sIntegrateSample;
extern int sParentPID;

extern void ReadParameters();
extern void ReadSHARE(ParticleDB *aPartDB);
extern void CheckSHARE(ParticleDB *aPartDB);
extern void MessageIntro();
extern void MessageHelp();
extern void MessageVersion();
extern void CopyINIFile();
extern void AddLogEntry(const char *aEntry);

#ifndef M_HBARC
#define M_HBARC 0.19732697
#endif
using namespace std;

//! Returns the phase space density
double phasespace_density_weight(double p, bool critical) {
  double mvac = 0.139;
  double Tpc = 0.155;
  double Lambda = M_PI*Tpc;
  double v0 ;

  if (critical) {

  } else {
    double e = sqrt(p*p + mvac*mvac) ;
  }
  f = 1./(exp(e/Tpc) - 1.) ;
}


//!
int main(int argc, char **argv) {
  sMainINI = "./events.ini";
  if (argc > 1) {
    TString tDummy;
    for (int i = 1; i < argc; i++) {
      tDummy = argv[i];
      if (tDummy.EndsWith(".ini")) {
        sMainINI = tDummy;
      }
    }
  }

  sMainConfig = new Configurator;
  ReadParameters();

  cout << "sEventDIR = " << sEventDIR << endl;
  cout << "sTimeStamp = " << sTimeStamp << endl;
  cout << "sModelINI = " << sModelINI << endl;
  cout << "sModel = " << sModel << endl;
  cout << "sRandomize = " << sRandomize << endl;
  cout << "sIntegrateSample  =" << sIntegrateSample << endl;
  cout << "sParentPID  =" << sParentPID << endl;

  ParticleDB *tPartDB;
  Parser *tParser = new Parser(sMainINI);
  tParser->ReadINI(sMainConfig);

  tPartDB = new ParticleDB();
  ReadSHARE(tPartDB);
  CheckSHARE(tPartDB);

  auto tRandom2 = new TRandom2;

  sMainConfig->PrintParameters();
  sIntegrateSample = sMainConfig->GetParameter("IntegrateSamples").Atoi();
  auto tInteg = new Integrator(sIntegrateSample);
  tInteg->SetMultiplicities(tPartDB);
  auto tEvent = new Event(tPartDB, tInteg);

  TFile fout("funcs.root", "RECREATE", "a test file");
  TH1::SetDefaultSumw2();
  TH1D *hpi = new TH1D("pi_dNdp", "dN/dp", 40, 0.0001, 2.);
  TH1D *hpi_decayed = new TH1D("pi_decayed_dNdp", "dN/dp", 40, 0.0000, 2.);
  TH1D *hpi_all = new TH1D("pi_all_dNdp", "dN/dp", 40, 0.0000, 2.);
  TH1D *hpi_all_weak = new TH1D("pi_all_weak_dNdp", "dN/dp", 40, 0.0000, 2.);

  double tcut = 2.5;
  TH1D *hpi_decayed_short =
      new TH1D("pi_decayed_short_25_dNdp", "dN/dp", 40, 0.0000, 2.);
  TH1D *htimes = new TH1D("times", "dNdt", 40, 0.0001, 10.);

  int nevents = 500;
  for (int ievent = 0; ievent < nevents; ievent++) {
    tEvent->Reset(ievent);
    tEvent->GeneratePrimordials();
    tEvent->DecayParticles();

    for (auto ptr = tEvent->GetParticleList()->begin();
         ptr != tEvent->GetParticleList()->end(); ++ptr) {

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

      double w_primordial = ;
      double w_stimulated = ;
      double w_total = ;

      // This is a primordial since eid == -1
      if (abs(pid) == 211 && primordial) {
        hpi->Fill(p);
        hpi_modified->Fill(p, w_primordial);
        w_total = w_primordial
      }
      // This is not a primordial
      if (abs(pid) == 211 && !primordial) {
        hpi_decayed->Fill(p);
        hpi_decayed->Fill(p, w_stimulated);

        // Momentum distribution with short decay times
        if (time < tcut) {
          hpi_decayed_short->Fill(p);
          hpi_decayed_short->Fill(p, w_stimulated);
          w_total = w_stimulated ;
        }
      }
      if (abs(pid) == 211) {
        hpi_all_weighted->Fill(p, w_total);
      }

      if (abs(pid) == 211 && ptr->fathereid != -1) {
        htimes->Fill(time);
      }
    }
  }
  fout.Write();


  delete tEvent;
  delete tInteg;
  delete tRandom2;
  delete tPartDB;
  delete tParser;
  delete sMainConfig;

  return 0;
}
