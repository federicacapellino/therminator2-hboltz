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
#include "Model_BlastWaveAniso.h"
#include "Model_HRG.h"
#include "ParticleDB.h"
#include "StructModel.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TRandom.h>
#include <TString.h>
#include <TTree.h>
#include <filesystem>
#include <fstream>
#include <memory>

static void MessageIntro();

struct ExistingRunInfo {
  bool found = false;
  std::string hash;
  unsigned int storedEvents = 0;
};

// Scans aEventDir for an existing event_*.root file and, if found, reads the hash/event-count
// recorded in its "parameters" tree (written by EventGenerator::SaveAsRoot). Used to decide
// whether a fresh GenerateEvents() call would just duplicate an already-generated sample -- note
// EventGenerator::EventFileName() bakes a fresh UUID into every filename, so without this check
// every run silently adds a new batch of files rather than overwriting anything.
static ExistingRunInfo CheckExistingEvents(const TString& aEventDir) {
  ExistingRunInfo tInfo;

  if (!std::filesystem::exists(aEventDir.Data())) return tInfo;

  std::string tFoundFile;
  for (const auto& tEntry : std::filesystem::directory_iterator(aEventDir.Data())) {
    std::string tName = tEntry.path().filename().string();
    if (tName.rfind("event_", 0) == 0 && tEntry.path().extension() == ".root") {
      tFoundFile = tEntry.path().string();
      break;
    }
  }
  if (tFoundFile.empty()) return tInfo;

  std::unique_ptr<TFile> tFile(TFile::Open(tFoundFile.c_str()));
  if (!tFile || tFile->IsZombie()) return tInfo;

  TTree* tTree = (TTree*) tFile->Get(_PARAMETERS_TREE_);
  if (!tTree) return tInfo;

  char tHash[9] = {0};
  UInt_t tStored = 0;
  tTree->SetBranchAddress(_MODELHASH_BRANCH_, tHash);
  tTree->SetBranchAddress(_STOREDEVENTS_BRANCH_, &tStored);
  tTree->GetEntry(0);

  tInfo.found = true;
  tInfo.hash = std::string(tHash);
  tInfo.storedEvents = tStored;
  return tInfo;
}

// Writes a human-readable record of the resolved run configuration into aEventDir: a verbatim
// copy of the main ini file (and the model ini, if any), plus a small manifest built from the
// model's own description (which already includes the hash used by CheckExistingEvents above)
// and the event count.
static void WriteRunRecord(const TString& aEventDir, const std::string& aIniFile,
                           const TString& aModelINI, Model* aModel, int aNumberOfEvents) {
  namespace fs = std::filesystem;
  std::error_code tEc;

  fs::copy_file(aIniFile, fs::path(aEventDir.Data()) / fs::path(aIniFile).filename(),
               fs::copy_options::overwrite_existing, tEc);
  if (aModelINI.Length() > 0) {
    fs::copy_file(aModelINI.Data(), fs::path(aEventDir.Data()) / fs::path(aModelINI.Data()).filename(),
                 fs::copy_options::overwrite_existing, tEc);
  }

  std::ofstream tManifest((fs::path(aEventDir.Data()) / "run_manifest.txt").string());
  if (tManifest) {
    tManifest << aModel->GetDescription();
    tManifest << "NumberOfEvents\t= " << aNumberOfEvents << '\n';
  }
}

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
  if      (tModel == "BlastWave")      tModelCode = 1;
  else if (tModel == "BlastWaveAniso") tModelCode = 2;
  else if (tModel == "HRG")            tModelCode = 99;
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
  case 2:
    return std::make_unique<Model_BlastWaveAniso>(modelINI);
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
  bool force = false;

  CLI::App app{"Therminator 2 event generator"};
  app.add_option("ini", iniFile, "main settings file")->default_str("events.ini");
  app.add_option("--seed", seed, "random seed (0 = unique seed from clock, default)")->default_val(0);
  app.add_flag("--force", force, "regenerate even if EventDir already has a matching sample");
  CLI11_PARSE(app, argc, argv);

  MessageIntro();

  auto tConfig = std::make_unique<Configurator>(iniFile);

  gRandom->SetSeed(seed);

  auto [modelCode, modelINI] = ReadModelConfig(tConfig.get());

  auto tPartDB = std::make_unique<ParticleDB>(tConfig->GetParameter("ShareDir", ""));
  auto tModel = CreateModel(modelCode, modelINI);

  TString tEventDir = tConfig->GetParameter("EventDir");
  int tNumberOfEvents = tConfig->GetParameter("NumberOfEvents").Atoi();

  if (!force) {
    ExistingRunInfo tExisting = CheckExistingEvents(tEventDir);
    if (tExisting.found) {
      if (tExisting.hash == std::string(tModel->GetHash()) &&
          tExisting.storedEvents >= static_cast<unsigned int>(tNumberOfEvents)) {
        PRINT_MESSAGE("<therm2_events>\t\""
                      << tEventDir << "\" already has " << tExisting.storedEvents
                      << " events for this configuration (hash " << tExisting.hash
                      << "); skipping generation.");
        PRINT_MESSAGE("\tUse --force to regenerate anyway.");
        return 0;
      } else {
        PRINT_MESSAGE("<therm2_events>\t\"" << tEventDir << "\" already contains a different run:");
        PRINT_MESSAGE("\tcached:    hash " << tExisting.hash << ", " << tExisting.storedEvents
                                           << " events");
        PRINT_MESSAGE("\trequested: hash " << tModel->GetHash() << ", " << tNumberOfEvents
                                           << " events");
        PRINT_MESSAGE("\tRefusing to overwrite. Clear/rename \"" << tEventDir
                                                                 << "\" first, or pass --force.");
        return 1;
      }
    }
  }

  auto tEventGen = std::make_unique<EventGenerator>(tConfig.get(), tPartDB.get(), tModel.get());

  tEventGen->AddLogEntry("[input]\t" + iniFile);
  tEventGen->AddLogEntry("[input]\t" + std::string(modelINI.Data()));

  tEventGen->GenerateEvents();
  tEventGen.reset(); // flush/close all output files before copying the ini/manifest alongside them

  WriteRunRecord(tEventDir, iniFile, modelINI, tModel.get(), tNumberOfEvents);

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
