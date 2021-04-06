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

extern Configurator *sMainConfig;
extern TString	sMainINI;
extern TString	sModelINI;
extern TString	sEventDIR;
extern TString	sTimeStamp;
extern int	sModel;
extern int	sRandomize;
extern int	sIntegrateSample;
extern int	sParentPID;


extern void ReadParameters();
extern void ReadSHARE(ParticleDB* aPartDB);
extern void CheckSHARE(ParticleDB* aPartDB);
extern void CopyINIFile();
extern void AddLogEntry(const char* aEntry);

#ifndef M_HBARC
#define M_HBARC  0.19732697
#endif
using namespace std ;


double legendre(const unsigned &l, double x) 
{
  double p = 0.;
  switch(l) {
  case 0:
    p=1.;
    break ;
  case 1:
    p=x ;
    break ;
  case 2:
    p=(3.*x*x - 1.)/2. ;
    break ;
  case 3:
    p=(5.*x*x*x - 3.*x)/2. ;
    break ;
  default:
    std::cout << " Bad selection in switch! " << std::endl;
    p=0. ;
    break ;
  }
  return p * sqrt((2.*l + 1.)) ;
}

double weight(Particle *prt,const double &m, const double &temper, const int &icase) 
{
  // double V = 5. /pow(M_HBARC, 3) ;

  double Spin	= prt->GetParticleType()->GetSpin();
  double Statistics = ( (Spin - static_cast<int>(Spin)) < 0.01 ? -1.0 : +1.0 );
  double e = prt->e ;
  double u = exp(-1.0 * e /temper) ; 
  double p = prt->GetP() ;
  double f = u/(1- Statistics*u) ;
  double phis = (1. + Statistics*f)  ;
  double val = (2. * Spin + 1.) * f * p * p / (kTwoPi3)    ;
  double x = prt->pz/p ;
  switch(icase) { 
    case 0:
      return val ;
      break ;
    case 1:
      return val * legendre(1,x) * p/e ;
      break ;
    case 2:
      return val * legendre(2,x) * p * p * phis;
      break ;
    case 3:
      return val * legendre(3,x) * p * p * p/e * phis ;
      break ;
    case 4:
      return val * legendre(1,x) * 2./5. * p * p * p/e * phis ;
      break;
  }
  return 0 ;

}

int main(int argc, char **argv)
{
  sMainINI = "./events.ini" ; 
  if (argc > 1) {
    TString tDummy;
    for(int i=1; i<argc;i++) {
      tDummy = argv[i];
      if (tDummy.EndsWith(".ini")) {
        sMainINI  = tDummy;
      }
    }
  }

  sMainConfig = new Configurator;
  ReadParameters();

  ParticleDB*	  tPartDB; 
  Parser *tParser = new Parser(sMainINI);
  tParser->ReadINI(sMainConfig);

  tPartDB     = new ParticleDB();
  ReadSHARE(tPartDB);

  std::list<Particle> tPList ;
  ParticleDecayer tPDecay(tPartDB, &tPList) ;

  auto tPType = tPartDB->GetParticleType("om0782zer") ;
  std::string s(tPType->GetName()) ;
  std::cout <<  "Processing decays of " << s << std::endl; 

  TFile fout("funcs.root", "RECREATE",  "a test file") ;
  double pmax = 4.0 ;
  double pmin = 0;
  const int nsamples = 40000000 ;

  TH1::SetDefaultSumw2() ;
  TH1D h0input("om0782zer_dNdp0", "dN/dp l=0", 100, 0., 3.) ;
  TH1D h0output("pi0139plu_from_om0782zer_dNdp0", "dN/dp l=0", 100, 0., 3.) ;

  TH1D h1input("om0782zer_dNdp1", "dN/dp l=1", 100, 0., 3.) ;
  TH1D h1output("pi0139plu_from_om0782zer_dNdp1", "dN/dp l=1", 100, 0., 3.) ;

  TH1D h2input("om0782zer_dNdp2", "dN/dp l=2", 100, 0., 3.) ;
  TH1D h2output("pi0139plu_from_om0782zer_dNdp2", "dN/dp l=2", 100, 0., 3.) ;

  TH1D h3input("om0782zer_dNdp3", "dN/dp l=3", 100, 0., 3.) ;
  TH1D h3output("pi0139plu_from_om0782zer_dNdp3", "dN/dp l=3", 100, 0., 3.) ;

  TH1D h1input_phis("om0782zer_dNdp1_phis", "dN/dp l=1 phi3", 100, 0., 3.) ;
  TH1D h1output_phis("pi0139plu_from_om0782zer_dNdp1_phis", "dN/dp l=1 phis", 100, 0., 3.) ;


  Particle tPart(tPType); 
  for (int i = 0 ; i < nsamples ; i++)  {
     double p = gRandom->Uniform(pmin, pmax) ;
     double cth = gRandom->Uniform(-1., 1.) ;
     double sth = sqrt(1. - cth*cth) ;
     double phi = gRandom->Uniform(0., 2.*M_PI) ;

     double m = tPType->GetMass() ;
     double ep = sqrt(p*p + m*m) ;
     double px = p*sth*cos(phi) ;
     double py = p*sth*sin(phi) ;
     double pz = p*cth;
     tPart.SetParticlePX(ep, px, py, pz, 0, 0, 0, 0) ;
     double temper = 0.160 ;

     const size_t ncases=5 ;
     double wFather[ncases]  = { } ;
     for (int icase = 0 ; icase < ncases ; icase++) {
       wFather[icase] = weight(&tPart, m, temper, icase)  ;
       double x = pz/p ;
       switch(icase) {
         case 0:
           h0input.Fill(p,wFather[icase]*legendre(0,x)) ;
           break;
         case 1:
           h1input.Fill(p,wFather[icase]*legendre(1,x)) ;
           break ;
         case 2:
           h2input.Fill(p, wFather[icase]*legendre(2,x)) ;
           break;
         case 3:
           h3input.Fill(p, wFather[icase]*legendre(3,x)) ;
           break ;
         case 4:
           h1input_phis.Fill(p, wFather[icase]*legendre(1,x)) ;
           break ;
         default:
           std::cout << "** Bad selection in switch icase = " << icase << std::endl;
           break ;
       }
     }

     tPDecay.DecayParticle(&tPart) ;

     for (auto daughter = tPList.begin() ;   daughter != tPList.end() ; ++daughter) {

       if (daughter->GetParticleType()->GetPDGCode()  != 211) {
         continue ;
       }

       double ep1, px1, py1, pz1 ;
       daughter->GetMomentum(&ep1, &px1, &py1, &pz1) ;
       
       double p1 = daughter->GetP() ;
       double cth1 = pz1/p1 ;

       double wD = 0.;
       for (int icase = 0; icase < ncases; icase++ ) {
         switch(icase) {
           case 0:
             wD = legendre(0,cth1) ;
             h0output.Fill(p1,wFather[icase] * wD) ;
             break;
           case 1:
             wD = legendre(1,cth1) ;
             h1output.Fill(p1,wFather[icase] * wD) ;
             break ;
           case 2:
             wD = legendre(2,cth1) ;
             h2output.Fill(p1, wFather[icase] * wD) ;
             break;
           case 3:
             wD = legendre(3,cth1) ;
             h3output.Fill(p1, wFather[icase] * wD) ;
             break ;
           case 4:
             wD = legendre(1,cth1) ;
             h1output_phis.Fill(p1, wFather[icase] * wD) ;
             break ;
           default:
             std::cout << "** Bad selection in switch icase = " << icase << std::endl;
             break ;
         }
       }
     }
     tPList.clear() ;
  }
  double sumw = h0input.GetSumOfWeights() ;
  h0output.Scale(1./sumw) ;
  h1output.Scale(1./sumw) ;
  h2output.Scale(1./sumw) ;
  h3output.Scale(1./sumw) ;
  h1output_phis.Scale(1./sumw) ;

  double volume = 4.0*M_PI*pmax ;
  double integral =sumw * volume/(1.0 * nsamples) ;
  std::cout << "The number of omegas's is " << integral << std::endl;
  
  fout.Write() ;
  fout.Close() ;
   
//  p.SetParticlePX(pt->GetMass(), 0, 0, 0, 0, 0, 0, 0) ;
  
     
  delete tParser ;
  delete sMainConfig ;
  delete tPartDB ;
  
  return 0;
}

