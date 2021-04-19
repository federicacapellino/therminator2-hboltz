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

#include <cstdio>
#include <fstream>
#include <TString.h>
#include <TH1D.h>
#include "THGlobal.h"
#include "Configurator.h"
#include "Parser.h"
#include "ParticleDB.h"
#include "ParticleDecayer.h"
#include "ParticleType.h"
#include "EventGenerator.h"
#include "Model_HRG.h"
#include "Model_BlastWave.h"


//! These are some global functions that we will want to get rid of
extern void ReadParameters();
extern void ReadSHARE(ParticleDB* aPartDB);
extern void CheckSHARE(ParticleDB* aPartDB);
extern void AddLogEntry(const char* aEntry);

#ifndef M_HBARC
#define M_HBARC  0.19732697
#endif
using namespace std ;


//! o4_dispersion provides services to compute the modified dispersion curve
//! and modified distribution function.
class o4_dispersion
{
  public: // private
    double fT ;  //! Temperature
    double fLambda ; //!Cutoff in GeV
    double fMvacuum; //! Vacuum pion mass in GeV
    double fRatio ; //! Ratio of chiral condensate to vacuum one
  public:
    o4_dispersion(const double &lambda_by_T=M_PI, const double &ratio=0.5) ;
    //! returns the modified distribution function
    double f(const double &p) ;
    //! returns the vacuum distribution function
    double fvac(const double &p) ;
    //! returns the modified dispersion curve E(p) in units of GeV
    double eofp(const double &p) ;
    //! returns the modified pole mass squared in units of GeV**2j
    double m2pole(const double &p) ;
    //! returns the chiral velocity  
    double v2(const double &p) ;

    //! Makes a plot of the dispersion curve
    void plot(const std::string &filename="pidispersion.out")  ;
} ;

o4_dispersion::o4_dispersion(const double &lambda_by_T, const double &ratio) :
  fT(0.155), fLambda(lambda_by_T*0.155), fMvacuum(0.1396), fRatio(ratio) 
{}

double o4_dispersion::eofp(const double &p)  { 
  return sqrt( v2(p) * p * p + m2pole(p) ) ;
}

double o4_dispersion::fvac(const double &p) {
  double x = sqrt(fMvacuum*fMvacuum + p*p)/fT ;
  return exp(-x)/(1. - exp(-x)) ;
}

double o4_dispersion::f(const double &p) {
  double x = eofp(p)/fT ;
  return exp(-x)/(1. - exp(-x)) ;
}

double o4_dispersion::m2pole(const double &p)  { 
  double mv2 = fMvacuum * fMvacuum ;
  double m02 = mv2 * fRatio ;
  double x = p*p/(fLambda*fLambda) ;
  return mv2 - (mv2 - m02)/(1 + x/2. + x *x) ;
}

double o4_dispersion::v2(const double &p)  { 
  double v02 = fRatio*fRatio ;
  double x = p*p/(fLambda*fLambda) ;
  return 1. - (1.-v02)/(1 + x/2. + x*x) ;
}

void o4_dispersion::plot(const std::string &filename) 
{
   int    np   = 100;
   double pmin = 0.;
   double pmax = 1.5;
   double dp   = (pmax - pmin)/ (double) np ;
   int  ip ;
   FILE *fp = fopen(filename.c_str(), "w") ;
   for (ip = 0 ; ip < np ; ip++) {
     double p = pmin + ip*dp;
     fprintf(fp,"%15.5e ", p) ;
     fprintf(fp,"%15.5e ", v2(p)) ;
     fprintf(fp,"%15.5e ", m2pole(p)) ;
     fprintf(fp,"%15.5e ", eofp(p)) ;
     fprintf(fp,"%15.5e ", f(p)) ;
     fprintf(fp,"%15.5e ", fvac(p)) ;
     fprintf(fp,"\n") ;
   }
   fclose(fp) ;
}

////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  TString tMainINI = "./events.ini" ; 

  o4_dispersion o4 ;
  o4.plot() ;

  // Read in the configuration file
  auto tMainConfig = std::make_unique<Configurator>("./events.ini");
  tMainConfig->PrintParameters() ;

  // Find the particles
  TString tShareDir ;
  try {
    tShareDir = tMainConfig->GetParameter("ShareDir"); 
  } catch (TString tError) {
    PRINT_MESSAGE("\tDid not find SHARE input file location.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  auto tPartDB     = std::make_unique<ParticleDB>(tShareDir);

  // Build the model
  Model_BlastWave bw("./blastwave.ini") ;
  // Initialize the integrator
  int tIntegrateSample = tMainConfig->GetParameter("IntegrateSamples").Atoi() ;
  auto tInteg = make_unique<Integrator>(tIntegrateSample, &bw) ;
  tInteg->SetMultiplicities(tPartDB.get()) ;

  EventGenerator tEvents(tMainConfig.get(), tPartDB.get(), &bw) ; 
  tEvents.GenerateEvents() ;
  
  return 0;
}

