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

#include <sstream>
#include <TMath.h>
#include "THGlobal.h"
#include "Configurator.h"
#include "Parser.h"
#include "Model_HRG.h"

using namespace TMath;
using namespace std;


Model_HRG::Model_HRG()
{
  mName = "HRG";
//  ReadParameters();
  Description();
  mTemperature = 0.155 ; // Units GeV
  mSize = 30. / kHbarC ; // The size of the box, e.g. 30 fm,  in 1/GeV
  mHyperCube = 4.0 * M_PI * 1. * pow(mSize, 3.);

// calculate parameter hash
  ostringstream oss;
  oss << mName;
  oss << mTemperature << mSize ;
  CalculateHash(TString(oss.str()));

}

Model_HRG::~Model_HRG()
{
}

double Model_HRG::GetIntegrand(ParticleType* aPartType)
{
  
  double Spin   = aPartType->GetSpin();
  
// Generate momentum components. Here we 
// are sampling three uniform deviates
//
// Zet = [0..1 ] , a transformed value of P
// CosP = [-1..1] , the polar angle of P
// PhiP = [0...2 Pi], the azimuthal angle of P
// X,Y,Z = [0... Size]
//
// The Hyper volume is 2 * Pi * 2 * 1 * Size^3
//
  double P, dP;
  {
    double Zet = gRandom->Rndm();
    P  = Zet / (1.0 - Zet);
    dP = 1.0 / ( (1.0 - Zet) * (1.0 - Zet) );
  }
  double phiP  = 2.0 * Pi() * gRandom->Rndm();
  double cosP  = gRandom->Uniform(-1., 1.) ;
  double sinP  = sqrt(1. - cosP*cosP) ;


// Type of statistics Bose-Einstein or Fermi-Dirac
  double tStatistics = ( (Spin - static_cast<int>(Spin)) < 0.01 ? -1.0 : +1.0 );

  double PdotU = Hypot(aPartType->GetMass(), P) ;

  double Integrand = (2.0 * Spin + 1.0) * 1.0 / kTwoPi3 * pow(P,2) * dP * 1.0 / (Exp( ( PdotU) / mTemperature ) + tStatistics);
  

// Return values
  Xt = 0. ;
  Xx = gRandom->Uniform(0., mSize) ;
  Xy = gRandom->Uniform(0., mSize) ;
  Xz = gRandom->Uniform(0., mSize) ;
  Pe = PdotU;
  Px = P * sinP * cos(phiP);
  Py = P * sinP * sin(phiP) ;
  Pz = P * cosP;
  return Integrand;
}

void Model_HRG::Description()
{
  ostringstream oss;
  oss << "##################################################"<< endl;
  oss << MODEL_NAME(mName);
  oss << "# - parameter #1           : " <<MODEL_PAR_DESC(mTemperature,"[GeV]");
  oss << "# - parameter #2           : " <<MODEL_PAR_DESC(mSize,"[fm]");
  oss << "# Parameters hash (CRC32)  : " <<MODEL_PAR_DESC(mHash,  "");
  oss << "##################################################"<< endl;
  mDescription = oss.str();
}

void Model_HRG::AddParameterBranch(TTree* aTree)
{
  Model_t tPar;
  
  tPar.dummy = 0.0;
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_)->Fill();
}

