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

#include "CLI11.hpp"
#include "Configurator.h"
#include "EventGenerator.h"
#include "Model_BlastWave.h"
#include "Model_HRG.h"
#include "ParticleDB.h"
#include "THGlobal.h"
#include <TRandom.h>
#include <TString.h>
#include <memory>

static void MessageIntro();

// Returns the model code and model ini path from events.ini settings.
// Exits on unknown model name or missing parameters.
static std::pair<int, TString> ReadModelConfig(Configurator* config) {
  TString tModel;
  int tModelCode = -1;
  try {
    tModel = config->GetParameter("FreezeOutModel");
  } catch (TString& tError) {
    PRINT_MESSAGE("<therm2_events>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  TString tModelINI = config->GetParameter("FreezeOutModelINI", "");
  if      (tModel == "BlastWave") tModelCode = 1;
  else if (tModel == "HRG")       tModelCode = 99;
  else {
    PRINT_MESSAGE("<therm2_events>\tUnknown FreezeOutModel: " << tModel);
    exit(_ERROR_GENERAL_MODEL_UNKNOWN_);
  }
  return {tModelCode, tModelINI};
}

static std::unique_ptr<Model> CreateModel(int modelCode, const TString& modelINI) {
  switch (modelCode) {
  case 1:
    return std::make_unique<Model_BlastWave>(modelINI);
  case 99:
    return std::make_unique<Model_HRG>();
  default:
    PRINT_MESSAGE("<therm2_events>\tModel code " << modelCode << " is not implemented.");
    exit(_ERROR_GENERAL_MODEL_UNKNOWN_);
  }
}

int main(int argc, char** argv) {
  std::string iniFile = "events.ini";
  int seed = 0;

  CLI::App app{"Therminator 2 event generator"};
  app.add_option("ini", iniFile, "main settings file")->default_str("events.ini");
  app.add_option("--seed", seed, "random seed (0 = unique seed from clock, default)")->default_val(0);
  CLI11_PARSE(app, argc, argv);

  MessageIntro();

  auto tConfig = std::make_unique<Configurator>(iniFile);

  gRandom->SetSeed(seed);

  auto [modelCode, modelINI] = ReadModelConfig(tConfig.get());

  auto tPartDB = std::make_unique<ParticleDB>(tConfig->GetParameter("ShareDir", ""));
  auto tModel = CreateModel(modelCode, modelINI);
  auto tEventGen = std::make_unique<EventGenerator>(tConfig.get(), tPartDB.get(), tModel.get());

  tEventGen->AddLogEntry("[input]\t" + iniFile);
  tEventGen->AddLogEntry("[input]\t" + std::string(modelINI.Data()));

  tEventGen->GenerateEvents();

  return 0;
}

static void MessageIntro() {
  PRINT_MESSAGE("  ***********************************************************************");
  PRINT_MESSAGE("  *\t\tTHERMINATOR 2 EVENTS version " << _THERMINATOR2_VERSION_ << "\t\t\t*");
  PRINT_MESSAGE("  *\t\t\t\t\t\t\t\t\t*");
  PRINT_MESSAGE("  * authors: M.Chojnacki, A.Kisiel, W.Florkowski, W.Broniowski\t\t*");
  PRINT_MESSAGE("  * cite as: arXiv:1102.0273\t\t\t\t\t\t*");
  PRINT_MESSAGE("  * webpage: http://therminator2.ifj.edu.pl/\t\t\t\t*");
  PRINT_MESSAGE("  ***********************************************************************");
}
