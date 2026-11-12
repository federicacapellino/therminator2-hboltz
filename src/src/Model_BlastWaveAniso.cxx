// clang-format off
/********************************************************************************
*                                                                              *
*             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
*                                                                              *
*  Model_BlastWaveAniso: modification of Model_BlastWave (Release 2.0.3)        *
*  implementing the anisotropic Blast-Wave of                                   *
*    F. Retiere and M. A. Lisa, Phys. Rev. C 70, 044907 (2004) [RL]             *
*                                                                              *
*  Changes with respect to Model_BlastWave:                                     *
*   - elliptic transverse shape: rt = sqrt(x^2/Rx^2 + y^2/Ry^2)      [RL 2]     *
*   - Woods-Saxon edge Omega(rt) = 1/(1+exp((rt-1)/a_s))             [RL 1]     *
*   - boost direction phi_b normal to the elliptic sub-shell:                   *
*       tan(phi_s) = (Ry/Rx)^2 tan(phi_b)                            [RL 3]     *
*   - flow rapidity rho = rt^n (rho0 + rho2 cos 2phi_b), n = 1 in RL [RL 4]     *
*   - Gaussian emission in proper time (Tau0, DeltaTau)              [RL 7]     *
*  Reaction plane = x-z plane; Ry > Rx is an out-of-plane extended source.      *
*                                                                              *
*  Limit Rx = Ry = RhoMax, Rho2 = 0, SurfaceDiff = 0, DeltaTau = 0,             *
*  FlowPower = 0, Rho0 = atanh(VelT) reproduces Model_BlastWave sample by       *
*  sample (same random numbers, same weights).                                  *
*                                                                              *
********************************************************************************/
// clang-format on

#include "Model_BlastWaveAniso.h"
#include "Configurator.h"
#include "Parser.h"
#include "THGlobal.h"
#include <TMath.h>
#include <sstream>

using namespace TMath;
using namespace std;

Model_BlastWaveAniso::Model_BlastWaveAniso()
    : Model(), mThermo(0), mRapPSRange(0.0), mRapSRange(0.0), mRx(0.0), mRy(0.0), mAs(0.0),
      mTau0(0.0), mDTau(0.0), mRho0(0.0), mRho2(0.0), mFlowPower(1.0), mRtMax(1.0),
      mHyperCube(0.0) {}

Model_BlastWaveAniso::Model_BlastWaveAniso(const TString& modelini) {
  mName = "Blast-Wave-Aniso";
  mThermo = new Thermodynamics();
  ReadParameters(modelini);
  Description();
  // CHANGED: the radial variable is rt in [0, mRtMax] instead of Rho in [0, RhoMax];
  // the Jacobian Rx*Ry*rt goes into the integrand.
  mHyperCube = mRapSRange * 2. * Pi() * mRtMax * mRapPSRange * 1.0 * 2. * Pi();
}

Model_BlastWaveAniso::~Model_BlastWaveAniso() {
  delete mThermo;
}

double Model_BlastWaveAniso::GetHyperCubeVolume() { return mHyperCube; }

double Model_BlastWaveAniso::GetIntegrand(ParticleType* aPartType, ParticleCoor& coor) {
  double dSigmaP, PdotU;
  double Spin, Statistics;
  double Tau, Rt, Theta, PhiB, RapS;
  double Xs, Ys, RhoF, Omega;
  double Mt, Pt, PhiP, RapP;
  double dPt;

  // Type of statistics Bose-Einstein or Fermi-Dirac
  Spin = aPartType->GetSpin();
  Statistics = ((Spin - static_cast<int>(Spin)) < 0.01 ? -1.0 : +1.0);
  // Generate spatial components
  // CHANGED: elliptic coordinates x = Rx rt cos(Theta), y = Ry rt sin(Theta),
  // dx dy = Rx Ry rt drt dTheta. Theta is the elliptic angle, not the spatial azimuth.
  Rt = mRtMax * gRandom->Rndm();
  Theta = 2.0 * Pi() * gRandom->Rndm();
  RapS = mRapSRange * gRandom->Rndm() - 0.5 * mRapSRange;
  // CHANGED: proper time from a normalized Gaussian truncated to tau > 0 (RL eq. 7).
  // The Gaussian is the sampling pdf, so only the factor Tau enters dSigmaP;
  // int Tau P(Tau) dTau = Tau0, i.e. the yield equals that of the sharp Tau0 surface.
  if (mDTau > 0.0) {
    do {
      Tau = gRandom->Gaus(mTau0, mDTau);
    } while (Tau <= 0.0);
  } else {
    Tau = mTau0;
  }
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
  // CHANGED: fluid element
  Xs = mRx * Rt * Cos(Theta);
  Ys = mRy * Rt * Sin(Theta);
  // boost along grad(rt) ~ (cos(Theta)/Rx, sin(Theta)/Ry): normal to the sub-shell (RL eq. 3)
  PhiB = ATan2(Sin(Theta) / mRy, Cos(Theta) / mRx);
  // flow rapidity (RL eq. 4, generalized to rt^n)
  RhoF = ((mFlowPower == 0.0) ? 1.0 : Power(Rt, mFlowPower)) * (mRho0 + mRho2 * Cos(2.0 * PhiB));
  // source weight (RL eq. 1)
  if (mAs > 0.0)
    Omega = 1.0 / (1.0 + Exp((Rt - 1.0) / mAs));
  else
    Omega = (Rt <= 1.0 ? 1.0 : 0.0);
  // Invariants
  // CHANGED: position-dependent flow; PdotU = -p.u (mostly plus) = energy in the local rest frame
  // (RL eq. 11)
  PdotU = Mt * CosH(RhoF) * CosH(RapS - RapP) - Pt * SinH(RhoF) * Cos(PhiB - PhiP);
  // CHANGED: transverse Jacobian Rx Ry rt replaces Rho
  dSigmaP = Tau * mRx * mRy * Rt * Mt * CosH(RapS - RapP);
  // integrand
  // CHANGED: extra factor Omega
  double integrand =
      (2.0 * Spin + 1.0) * 1.0 / kTwoPi3 * Pt * dPt * dSigmaP * Omega * 1.0 /
      (Exp((PdotU - mThermo->GetChemicalPotential(aPartType)) / mThermo->GetTemperature()) +
       Statistics);
  // particle X and P coordinates
  coor.mass = aPartType->GetMass();
  coor.t = Tau * CosH(RapS);
  coor.x = Xs;
  coor.y = Ys;
  coor.z = Tau * SinH(RapS);
  coor.e = Mt * CosH(RapP);
  coor.px = Pt * Cos(PhiP);
  coor.py = Pt * Sin(PhiP);
  coor.pz = Mt * SinH(RapP);
  coor.w = integrand * mHyperCube;
  return integrand;
}

void Model_BlastWaveAniso::Description() {
  ostringstream oss;
  oss << "##################################################" << endl;
  oss << MODEL_NAME(mName);
  oss << "# - rel. rapidity range    : " << MODEL_PAR_DESC(mRapPSRange, "[units]");
  oss << "# - spatial rapidity range : " << MODEL_PAR_DESC(mRapSRange, "[units]");
  oss << "# - in-plane radius Rx     : " << MODEL_PAR_DESC(mRx * kHbarC, "[fm]");
  oss << "# - out-of-plane radius Ry : " << MODEL_PAR_DESC(mRy * kHbarC, "[fm]");
  oss << "# - surface diffuseness    : " << MODEL_PAR_DESC(mAs, "[units]");
  oss << "# - mean proper time Tau0  : " << MODEL_PAR_DESC(mTau0 * kHbarC, "[fm]");
  oss << "# - emission duration      : " << MODEL_PAR_DESC(mDTau * kHbarC, "[fm]");
  oss << "# - flow rapidity Rho0     : " << MODEL_PAR_DESC(mRho0, "[units]");
  oss << "# - flow modulation Rho2   : " << MODEL_PAR_DESC(mRho2, "[units]");
  oss << "# - flow profile power     : " << MODEL_PAR_DESC(mFlowPower, "[units]");
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

void Model_BlastWaveAniso::AddParameterBranch(TTree* aTree) {
  Model_t_BlastWaveAniso tPar;

  tPar.RapPSRange = mRapPSRange;
  tPar.RapSRange = mRapSRange;
  tPar.Rx = mRx * kHbarC;
  tPar.Ry = mRy * kHbarC;
  tPar.As = mAs;
  tPar.Tau0 = mTau0 * kHbarC;
  tPar.DTau = mDTau * kHbarC;
  tPar.Rho0 = mRho0;
  tPar.Rho2 = mRho2;
  tPar.FlowPower = mFlowPower;
  tPar.Temp = mThermo->GetTemperature() * 1000.0;
  tPar.MuB = mThermo->GetMuB() * 1000.0;
  tPar.MuI = mThermo->GetMuI() * 1000.0;
  tPar.MuS = mThermo->GetMuS() * 1000.0;
  tPar.MuC = mThermo->GetMuC() * 1000.0;
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_BLASTWAVEANISO_)->Fill();
}

void Model_BlastWaveAniso::ReadParameters(const TString& modelini) {
  Configurator* tModelParam;
  Parser* tParser;

  tModelParam = new Configurator;
  tParser = new Parser(modelini.Data());
  tParser->ReadINI(tModelParam);
  delete tParser;

  try {
    mRapPSRange = tModelParam->GetParameter("RapPSRange").Atof();                     // [1]
    mRapSRange = tModelParam->GetParameter("RapSRange").Atof();                       // [1]
    mRx = tModelParam->GetParameter("Rx").Atof() / kHbarC;                            // [GeV^-1]
    mRy = tModelParam->GetParameter("Ry").Atof() / kHbarC;                            // [GeV^-1]
    mAs = tModelParam->GetParameter("SurfaceDiff").Atof();                            // [1]
    mTau0 = tModelParam->GetParameter("Tau0").Atof() / kHbarC;                        // [GeV^-1]
    mDTau = tModelParam->GetParameter("DeltaTau").Atof() / kHbarC;                    // [GeV^-1]
    mRho0 = tModelParam->GetParameter("Rho0").Atof();                                 // [1]
    mRho2 = tModelParam->GetParameter("Rho2").Atof();                                 // [1]
    mFlowPower = tModelParam->GetParameter("FlowPower").Atof();                       // [1]
    mThermo->SetTemperature(tModelParam->GetParameter("Temperature").Atof() * 0.001); // [GeV]
    mThermo->SetChemistry(tModelParam->GetParameter("MuB").Atof() * 0.001,
                          tModelParam->GetParameter("MuI").Atof() * 0.001,
                          tModelParam->GetParameter("MuS").Atof() * 0.001,
                          tModelParam->GetParameter("MuC").Atof() * 0.001); // [GeV]
  } catch (TString tError) {
    PRINT_MESSAGE("<Model_BlastWaveAniso::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary model parameters.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  if (mRx <= 0.0 || mRy <= 0.0 || mTau0 <= 0.0 || mDTau < 0.0 || mAs < 0.0) {
    PRINT_MESSAGE("<Model_BlastWaveAniso::ReadParameters>\tNeed Rx, Ry, Tau0 > 0 and DeltaTau, "
                  "SurfaceDiff >= 0.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  // sampled rt range: hard edge -> 1; Woods-Saxon -> 1 + 10 a_s (Omega < 5e-5 beyond)
  mRtMax = (mAs > 0.0) ? 1.0 + 10.0 * mAs : 1.0;

  // calculate parameter hash (all parameters: they change the multiplicities)
  ostringstream oss;
  oss << GetName();
  oss << mRapPSRange << mRapSRange << mRx << mRy << mAs << mTau0 << mDTau << mRho0 << mRho2
      << mFlowPower;
  oss << mThermo->GetTemperature() << mThermo->GetMuB() << mThermo->GetMuI() << mThermo->GetMuS()
      << mThermo->GetMuC();
  CalculateHash(TString(oss.str()));

  delete tModelParam;
}
