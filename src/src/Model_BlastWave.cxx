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

#include "Model_BlastWave.h"
#include "Configurator.h"
#include "Parser.h"
#include "THGlobal.h"
#include <TMath.h>
#include <sstream>

using namespace TMath;
using namespace std;

Model_BlastWave::Model_BlastWave()
    : Model(), mThermo(0), mRapPSRange(0.0), mRapSRange(0.0), mRhoMax(0.0), mTau(0.0), mVt(0.0) {}

Model_BlastWave::Model_BlastWave(const TString& modelini) {
  mName = "Blast-Wave";
  mThermo = new Thermodynamics();
  ReadParameters(modelini);
  Description();
  mHyperCube = mRapSRange * 2. * Pi() * mRhoMax * mRapPSRange * 1.0 * 2. * Pi();
}

Model_BlastWave::~Model_BlastWave() {
  delete mThermo;
}

double Model_BlastWave::GetHyperCubeVolume() { return mHyperCube; }

double Model_BlastWave::GetIntegrand(ParticleType* aPartType, ParticleCoor& coor) {
  double dSigmaP, PdotU;
  double Spin, Statistics;
  double Tau, Rho, PhiS, RapS;
  double Mt, Pt, PhiP, RapP;
  double dPt;

  // Type of statistics Bose-Einstein or Fermi-Dirac
  Spin = aPartType->GetSpin();
  Statistics = ((Spin - static_cast<int>(Spin)) < 0.01 ? -1.0 : +1.0);
  // Generate spatial components
  Rho = mRhoMax * gRandom->Rndm();
  PhiS = 2.0 * Pi() * gRandom->Rndm();
  RapS = mRapSRange * gRandom->Rndm() - 0.5 * mRapSRange;
  Tau = mTau;
  // Generate momentum components
  {
    double Zet = gRandom->Rndm();
    Pt = Zet / (1.0 - Zet);
    dPt = 1.0 / ((1.0 - Zet) * (1.0 - Zet));
  }
  PhiP = 2.0 * Pi() * gRandom->Rndm();
  {
    double RapPS = mRapPSRange * gRandom->Rndm() - 0.5 * mRapPSRange;
    RapP = RapS + RapPS;
  }
  Mt = Hypot(aPartType->GetMass(), Pt);
  // Invariants
  PdotU = 1.0 / Sqrt(1 - mVt * mVt) * (Mt * CosH(RapS - RapP) - mVt * Pt * Cos(PhiS - PhiP));
  dSigmaP = Tau * Rho * Mt * CosH(RapS - RapP);
  // integrand
  double integrand = (2.0 * Spin + 1.0) * 1.0 / kTwoPi3 * Pt * dPt * dSigmaP * 1.0 /
                     (Exp((PdotU - mThermo->GetChemicalPotential(aPartType)) / mThermo->GetTemperature()) +
                      Statistics);
  // particle X and P coordinates
  coor.mass = aPartType->GetMass();
  coor.t  = Tau * CosH(RapS);
  coor.x  = Rho * Cos(PhiS);
  coor.y  = Rho * Sin(PhiS);
  coor.z  = Tau * SinH(RapS);
  coor.e  = Mt * CosH(RapP);
  coor.px = Pt * Cos(PhiP);
  coor.py = Pt * Sin(PhiP);
  coor.pz = Mt * SinH(RapP);
  coor.w  = integrand * mHyperCube;
  return integrand;
}

void Model_BlastWave::Description() {
  ostringstream oss;
  oss << "##################################################" << endl;
  oss << MODEL_NAME(mName);
  oss << "# - rel. rapidity range    : " << MODEL_PAR_DESC(mRapPSRange, "[units]");
  oss << "# - spatial rapidity range : " << MODEL_PAR_DESC(mRapSRange, "[units]");
  oss << "# - max cylinder size      : " << MODEL_PAR_DESC(mRhoMax * kHbarC, "[fm]");
  oss << "# - Blast-Wave time        : " << MODEL_PAR_DESC(mTau * kHbarC, "[fm]");
  oss << "# - transverse velocity    : " << MODEL_PAR_DESC(mVt, "[c]");
  oss << "# - freeze-out temperature : "
      << MODEL_PAR_DESC(mThermo->GetTemperature() * 1000.0, "[MeV]");
  oss << "# - chem. potential Mu_B   : " << MODEL_PAR_DESC(mThermo->GetMuB() * 1000.0, "[MeV]");
  oss << "# - chem. potential Mu_I3  : " << MODEL_PAR_DESC(mThermo->GetMuI() * 1000.0, "[MeV]");
  oss << "# - chem. potential Mu_S   : " << MODEL_PAR_DESC(mThermo->GetMuS() * 1000.0, "[MeV]");
  oss << "# - chem. potential Mu_C   : " << MODEL_PAR_DESC(mThermo->GetMuC() * 1000.0, "[MeV]");
  oss << "# Parameters hash (CRC32)  : " << MODEL_PAR_DESC(mHash, "");
  oss << "##################################################" << endl;
  mDescription = oss.str();
}

void Model_BlastWave::AddParameterBranch(TTree* aTree) {
  Model_t_BlastWave tPar;

  tPar.RapPSRange = mRapPSRange;
  tPar.RapSRange = mRapSRange;
  tPar.RhoMax = mRhoMax * kHbarC;
  tPar.Tau = mTau * kHbarC;
  tPar.Vt = mVt;
  tPar.Temp = mThermo->GetTemperature() * 1000.0;
  tPar.MuB = mThermo->GetMuB() * 1000.0;
  tPar.MuI = mThermo->GetMuI() * 1000.0;
  tPar.MuS = mThermo->GetMuS() * 1000.0;
  tPar.MuC = mThermo->GetMuC() * 1000.0;
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_BLASTWAVE_)->Fill();
}

void Model_BlastWave::ReadParameters(const TString& modelini) {
  Configurator* tModelParam;
  Parser* tParser;

  tModelParam = new Configurator;
  tParser = new Parser(modelini.Data());
  tParser->ReadINI(tModelParam);
  delete tParser;

  try {
    mRapPSRange = tModelParam->GetParameter("RapPSRange").Atof();                     // [1]
    mRapSRange = tModelParam->GetParameter("RapSRange").Atof();                       // [1]
    mRhoMax = tModelParam->GetParameter("RhoMax").Atof() / kHbarC;                    // [GeV^-1]
    mTau = tModelParam->GetParameter("Tau").Atof() / kHbarC;                          // [GeV^-1]
    mVt = tModelParam->GetParameter("VelT").Atof();                                   // [c]
    mThermo->SetTemperature(tModelParam->GetParameter("Temperature").Atof() * 0.001); // [GeV]
    mThermo->SetChemistry(tModelParam->GetParameter("MuB").Atof() * 0.001,
                          tModelParam->GetParameter("MuI").Atof() * 0.001,
                          tModelParam->GetParameter("MuS").Atof() * 0.001,
                          tModelParam->GetParameter("MuC").Atof() * 0.001); // [GeV]
  } catch (TString tError) {
    PRINT_MESSAGE("<Model_BlastWave::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary model parameters.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }

  // calculate parameter hash
  ostringstream oss;
  oss << GetName();
  oss << mRapPSRange << mRapSRange << mRhoMax << mTau << mVt;
  oss << mThermo->GetTemperature() << mThermo->GetMuB() << mThermo->GetMuI() << mThermo->GetMuS()
      << mThermo->GetMuC();
  CalculateHash(TString(oss.str()));

  delete tModelParam;
}
