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

#include "Model.h"
#include "Configurator.h"
#include "Crc32.h"
#include "THGlobal.h"
#include <sys/stat.h>

using namespace std;

Model::Model()
    : Xt(0.0), Xx(0.0), Xy(0.0), Xz(0.0), Pe(0.0), Px(0.0), Py(0.0), Pz(0.0), mHyperCube(0.0) {
  mName = "";
  mHash = "";
  mDescription = "";
}

Model::~Model() {}

void Model::AddParameterBranch(TTree* aTree) {
  Model_t tPar;

  tPar.dummy = -1.0;
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_);
}

void Model::SetParticlePX(Particle* aParticle) {
  aParticle->SetParticlePX(Pe, Px, Py, Pz, Xt, Xx, Xy, Xz);
}

double Model::GetHyperCubeVolume() {
  return mHyperCube;
}

const char* Model::GetHash() {
  return mHash.Data();
}

const char* Model::GetName() {
  return mName.Data();
}

const char* Model::GetDescription() {
  return mDescription.Data();
}

void Model::CalculateHash(TString aPreHash) {
  Crc32 tHash;
  tHash.Update(aPreHash.Data(), aPreHash.Length());
  tHash.Finish();
  mHash = tHash.GetValueHex();
  PRINT_DEBUG_1("<Model::Hash>\t" << aPreHash.Data() << " -> 0x" << mHash);
}
