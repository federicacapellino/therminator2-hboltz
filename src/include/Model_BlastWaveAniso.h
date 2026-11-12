// clang-format off
/********************************************************************************
*  Model_BlastWaveAniso: anisotropic Blast-Wave of                              *
*    F. Retiere and M. A. Lisa, Phys. Rev. C 70, 044907 (2004)                  *
*  Derived from Model_BlastWave of THERMINATOR 2 (arXiv:1102.0273).             *
********************************************************************************/
// clang-format on

#ifndef _TH2_MODEL_BLASTWAVEANISO_H_
#define _TH2_MODEL_BLASTWAVEANISO_H_

#include "Model.h"
#include "Thermodynamics.h"

class Model_BlastWaveAniso : public Model {
public:
  Model_BlastWaveAniso();
  Model_BlastWaveAniso(const TString& modelini);
  ~Model_BlastWaveAniso();

  double GetIntegrand(ParticleType* aPartType, ParticleCoor& coor) override;
  double GetHyperCubeVolume() override;
  void AddParameterBranch(TTree* aTree) override;

private:
  void Description();
  void ReadParameters(const TString& modelini);

  Thermodynamics* mThermo;
  double mRapPSRange; // range of Y - Y_s                         [1]
  double mRapSRange;  // space-time rapidity range                [1]
  double mRx;         // in-plane radius                          [GeV^-1]
  double mRy;         // out-of-plane radius                      [GeV^-1]
  double mAs;         // surface diffuseness a_s (0 = box)        [1]
  double mTau0;       // mean freeze-out proper time              [GeV^-1]
  double mDTau;       // Gaussian emission duration               [GeV^-1]
  double mRho0;       // flow rapidity at rt = 1                  [1]
  double mRho2;       // cos(2 phi_b) modulation of flow rapidity [1]
  double mFlowPower;  // n in rho = rt^n (rho0 + rho2 cos 2phi_b) [1]
  double mRtMax;      // upper limit of sampled rt                [1]
  double mHyperCube;  // volume of the sampled parameter space
};

#endif
