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

#include "ParticleDB.h"
#include "DecayTable.h"
#include "Parser.h"
#include "THGlobal.h"
#include <TString.h>
#include <cstdlib>
#include <iostream>

void ReadSHARE(ParticleDB* aPartDB, const TString& aShareDir);
void CheckSHARE(ParticleDB* aPartDB) {
  ParticleType* tType;
  DecayTable* tDecTable;
  double SumBR;

  PRINT_DEBUG_2("<ParticleDB>\tRead " << aPartDB->GetParticleTypeCount() << " particle types.");
  for (int tPart = 0; tPart < aPartDB->GetParticleTypeCount(); tPart++) {
    tType = aPartDB->GetParticleType(tPart);
    PRINT_DEBUG_2("\tParticle " << tType->GetNumber() << ": " << tType->GetName() << ", Mass = "
                                << tType->GetMass() << ", Gamma = " << tType->GetGamma()
                                << ", Spin = " << tType->GetSpin() << ", I  = " << tType->GetI()
                                << ", I3 = " << tType->GetI3() << ", BarionN = "
                                << tType->GetBarionN() << ", StrangeN = " << tType->GetStrangeN()
                                << ", CharmN = " << tType->GetCharmN() << ", Charge = "
                                << tType->GetCharge() << ", MC# = " << tType->GetPDGCode());
    SumBR = 0.0;
    if ((tDecTable = tType->GetTable())) {
      for (int tChan = 0; tChan < tDecTable->GetChannelCount() + 1; tChan++) {
        if (tDecTable->GetDecayChannel(tChan)->Is3Particle()) {
          PRINT_DEBUG_2(
              "\t\tChannel "
              << tChan << ": "
              << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChan)->GetParticle1()))
                     ->GetName()
              << " + "
              << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChan)->GetParticle2()))
                     ->GetName()
              << " + "
              << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChan)->GetParticle3()))
                     ->GetName()
              << ", BR = " << tDecTable->GetDecayChannel(tChan)->GetBranchingRatio()
              << ", Step = " << tDecTable->GetDecayStep(tChan));
        } else {
          PRINT_DEBUG_2(
              "\t\tChannel "
              << tChan << ": "
              << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChan)->GetParticle1()))
                     ->GetName()
              << " + "
              << (aPartDB->GetParticleType(tDecTable->GetDecayChannel(tChan)->GetParticle2()))
                     ->GetName()
              << ", BR = " << tDecTable->GetDecayChannel(tChan)->GetBranchingRatio()
              << ", Step = " << tDecTable->GetDecayStep(tChan));
        }
        SumBR += tDecTable->GetDecayChannel(tChan)->GetBranchingRatio();
      }
      PRINT_DEBUG_2("\t\tSum BR = " << SumBR);
    }
  }
}

ParticleDB::ParticleDB(const TString& path) {
  mParticleTable.clear();
  mParticleNames.clear();

  TString tPath = path;
  if (tPath.IsNull()) {
    const char* tEnv = std::getenv("THERMINATOR_SHARE");
    if (tEnv)
      tPath = tEnv;
    else {
      PRINT_MESSAGE("<ParticleDB>\tNo share directory. Set THERMINATOR_SHARE or pass a path.");
      exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
    }
  }

  if (!tPath.EndsWith("/"))
    tPath += "/";

  ReadSHARE(this, tPath);
  CheckSHARE(this);
}

ParticleDB::~ParticleDB() {
  mParticleTable.clear();
  mParticleNames.clear();
}

int ParticleDB::AddParticleType(ParticleType* aPartType) {
  mParticleTable.push_back(*aPartType);
  mParticleNames[aPartType->GetName()] = mParticleTable.size() - 1;
  return mParticleTable.size() - 1;
}

ParticleType* ParticleDB::GetParticleType(int aIndex) {
  return &(mParticleTable[aIndex]);
}

ParticleType* ParticleDB::GetParticleType(TString aName) {
  auto it = mParticleNames.find(aName);
  if (it == mParticleNames.end()) {
    PRINT_MESSAGE("<ParticleDB::GetParticleType>\tParticle '" << aName << "' not found.");
    return nullptr;
  }
  return &(mParticleTable[it->second]);
}

int ParticleDB::GetParticleTypeIndex(TString aName) {
  return mParticleNames[aName];
}

int ParticleDB::GetParticleTypeCount() {
  return mParticleTable.size();
}

int ParticleDB::ExistsParticleType(TString aName) {
  return mParticleNames.count(aName);
}

/* @brief Reads in the particle database:

@param [in/out] aPartDB a previously allocated particle database. After
calling routine the database is filled up

@param [in] aShareDir path of the shared files
*/
void ReadSHARE(ParticleDB* aPartDB, const TString& aShareDir) {
  auto tParser = std::make_unique<Parser>((aShareDir + "particles.data").Data());
  tParser->ReadSHAREParticles(aPartDB);

  tParser = std::make_unique<Parser>((aShareDir + "decays.data").Data());
  tParser->ReadSHAREDecays(aPartDB);
}
