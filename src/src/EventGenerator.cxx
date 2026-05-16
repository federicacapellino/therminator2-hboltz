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

#include "EventGenerator.h"
#include "Configurator.h"
#include "Crc32.h"
#include "Event.h"
#include "StructEvent.h"
#include "THGlobal.h"
#include <TDatime.h>
#include <TUUID.h>
#include <fstream>
#include <sstream>

using namespace std;

EventGenerator::EventGenerator()
    : mDB(0), mInteg(0), mEvent(0), mParticleTree(0), mFile(0), mFileCounter(0), mEventCounter(0),
      mNumberOfEvents(0), mEventExportType(0), kEventsPerFile(_EVENTS_PER_FILE_),
      mIntegrateSample(0) {}

EventGenerator::EventGenerator(Configurator* config, ParticleDB* aDB, Model* model)
    : mDB(aDB), mParticleTree(0), mFile(0), mFileCounter(0), mEventCounter(0), mNumberOfEvents(0),
      mEventExportType(0), kEventsPerFile(_EVENTS_PER_FILE_) {
  TDatime tDate;
  tDate.Set();
  mTimeStamp = tDate.AsSQLString();

  // Generate an 8-char run ID — unique per invocation, shared across all
  // files from this run, safe for parallel jobs writing to the same directory.
  mRunID = TString(TUUID().AsString()).Remove(8);

  ReadParameters(config);

  mInteg = new Integrator(mIntegrateSample, model);
  mInteg->SetMultiplicities(mDB);
  AddLogEntry(mInteg->GetFileDescription());

  mEvent = new Event(mDB, mInteg);
}

TString EventGenerator::EventFileName(int index) const {
  return TString::Format("%sevent_%s_%05i.root", mEventDIR.Data(), mRunID.Data(), index);
}

EventGenerator::~EventGenerator() {
  if (mNumberOfEvents) {
    int tFirst = mFileCounter - (mNumberOfEvents / kEventsPerFile) -
                 (((mNumberOfEvents % kEventsPerFile) > 0) ? 1 : 0);
    TString tMsg;
    if (mFileCounter - tFirst == 1)
      tMsg = TString::Format("[output]\t\"%s\"\t[events]\t%i\t%i", EventFileName(tFirst).Data(),
                             mNumberOfEvents, kEventsPerFile);
    else
      tMsg = TString::Format("[output]\t\"%s\" - \"%s\"\t[events]\t%i\t%i",
                             EventFileName(tFirst).Data(), EventFileName(mFileCounter - 1).Data(),
                             mNumberOfEvents, kEventsPerFile);
    AddLogEntry(tMsg.Data());
  }

  delete mInteg;
  delete mEvent;
}

void EventGenerator::GenerateEvents() {
  TDatime tDate;

  tDate.Set();
  Print("[" + std::string(tDate.AsSQLString()) + "]\tGenerating events");
  for (int tIter = 0; tIter < mNumberOfEvents; tIter++) {
    mEventCounter = tIter + 1;
    SetEventID(tIter);

    mEvent->Reset();
    mEvent->GeneratePrimordials();
    mEvent->DecayParticles();

    SaveEvent();
  }
}

void EventGenerator::SaveEvent() {
  switch (mEventExportType) {
  case 0:
    SaveAsRoot();
    break;
  case 1:
    SaveAsText();
    break;
  case 2:
    SaveAsText();
    SaveAsRoot();
    break;
  }
}

void EventGenerator::SaveAsRoot() {
  TDatime tDate;

  char tEventFName[kFileNameMaxChar];
  char tTempFName[kFileNameMaxChar];
  static ParticleCoor tPartCoor;
  static StructEvent tStructEvent;
  list<Particle>::iterator tParticle_i;

  // open new file every _EVENTS_PER_FILE_ events
  if (!((mEventCounter - 1) % kEventsPerFile)) {
    sprintf(tTempFName, "%sevent_%s_%05i.tmp", mEventDIR.Data(), mRunID.Data(), mFileCounter);
    mFile = new TFile(tTempFName, "RECREATE");
    mFile->cd();

    mParticleTree = new TTree(_PARTICLES_TREE_, "particle tree");
    mEventTree = new TTree(_EVENTS_TREE_, "event tree");
    mParameterTree = new TTree(_PARAMETERS_TREE_, "parameters and model description tree");

    // (void*) cast to avoid some compilation errors on older ROOT versions
    char tTimeStamp[21];
    sprintf(tTimeStamp, "%s", mTimeStamp.Data());

    mParticleTree->Branch(_PARTICLE_BRANCH_, &tPartCoor, _PARTICLE_FORMAT_);

    mEventTree->Branch(_EVENTS_BRANCH_, &tStructEvent, _EVENTS_FORMAT_);

    // Model independent parameters
    mParameterTree->Branch(_INTEGRATESAMPLE_BRANCH_, (UInt_t*)&mIntegrateSample, "i");
    mParameterTree->Branch(_TIMESTAMP_BRANCH_, (Char_t*)tTimeStamp, _TIMESTAMP_FORMAT_);
    mParameterTree->Branch(_MODELNAME_BRANCH_,
                           (Char_t*)mEvent->GetIntegrator()->GetModel()->GetName(),
                           _MODELNAME_FORMAT_);
    mParameterTree->Branch(_MODELHASH_BRANCH_,
                           (Char_t*)mEvent->GetIntegrator()->GetModel()->GetHash(),
                           _MODELHASH_FORMAT_);
    mParameterTree->Branch(_MODELDESCRIPTION_BRANCH_,
                           (Char_t*)mEvent->GetIntegrator()->GetModel()->GetDescription(),
                           _MODELDESCRIPTION_FORMAT_);

    // Add Model specific parameters, The derived model class is responsible for this. And that the
    // parameter structure has a lifetime long enough for the tree to be filled.
    mEvent->GetIntegrator()->GetModel()->AddParameterBranch(mParameterTree);

    // Fill up the tree
    mParameterTree->Fill();
    PRINT_DEBUG_2("<EventGenerator::SaveAsRoot>\tCreated file " << tTempFName);
  }

  // add all Particle entries to file and Event information
  if (mFile) {
    for (tParticle_i = mEvent->GetParticleList().begin();
         tParticle_i != mEvent->GetParticleList().end(); tParticle_i++) {
      tPartCoor = (*tParticle_i);
      tPartCoor.t *= kHbarC;
      tPartCoor.x *= kHbarC;
      tPartCoor.y *= kHbarC;
      tPartCoor.z *= kHbarC;
      tPartCoor.eventid = GetEventID();
      mParticleTree->Fill();
    }
    tStructEvent.eventID = GetEventID();
    tStructEvent.entries = mEvent->GetParticleList().size();
    tStructEvent.entriesprev = 0; // not used here
    mEventTree->Fill();
    PRINT_DEBUG_2("<EventGenerator::SaveAsRoot>\tEvent " << mEventCounter << " saved.");
  }

  // close & rename file
  if (!(mEventCounter % kEventsPerFile) || (mEventCounter == mNumberOfEvents)) {
    mFile->Write();
    mFile->Close();
    mFile = 0;
    // rename temporary file to ROOT
    snprintf(tEventFName, kFileNameMaxChar, "%s", EventFileName(mFileCounter).Data());
    snprintf(tTempFName, kFileNameMaxChar, "%sevent_%s_%05i.tmp", mEventDIR.Data(), mRunID.Data(),
             mFileCounter);

    // Hopefully this operating system command works
    if (rename(tTempFName, tEventFName) != 0) {
      PRINT_MESSAGE("<EventGenerator::SaveAsRoot>\tUnable to rename temp file to " << tTempFName);
      exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
    } else {
      tDate.Set();
      Print("\n[" + std::string(tDate.AsSQLString()) + "]\tFile " + tEventFName + " written");
    }
    mFileCounter++;
  }
}

void EventGenerator::SaveAsText() {
  TDatime tDate;
  ofstream tFile;
  list<Particle>::iterator tParticle_i;

  tFile.open((mEventDIR + "event.txt").Data(), ios_base::app);
  if ((tFile) && (tFile.is_open())) {
    if (static_cast<long>(tFile.tellp()) == 0) {
      tDate.Set();
      PRINT_MESSAGE("\n[" << tDate.AsSQLString() << "]\tFile " << (mEventDIR + "event.txt").Data()
                          << " created.");
      tFile << "# THERMINATOR 2 text output" << endl;
      tFile << "#<EVENT_ENTRY>"
               "\tEID\tfatherEID\tPID\tfatherPID\trootPID\tdecayed\tmass\tE\tp_"
               "x\tp_y\tp_z\tt\tx\ty\tz\t</EVENT_ENTRY>"
            << endl;
    }
    tFile << "#<EVENT_ID>\t0x" << hex << uppercase << GetEventID() << nouppercase << dec
          << "\t</EVENT_ID>" << endl;
    tFile << "#<NO_OF_PARTICLES>\t" << mEvent->GetParticleList().size() << "\t</NO_OF_PARTICLES>"
          << endl;
    for (tParticle_i = mEvent->GetParticleList().begin();
         tParticle_i != mEvent->GetParticleList().end(); tParticle_i++) {
      tFile << tParticle_i->MakeTEXTEntry();
    }
    PRINT_DEBUG_2("<EventGenerator::SaveAsText>\tEventID " << GetEventID() << " saved.");
    tFile.close();
  } else {
    PRINT_MESSAGE("<EventGenerator::SaveAsText>\tUnable to create file "
                  << (mEventDIR + "event.txt").Data());
    exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
  }
}

void EventGenerator::ReadParameters(Configurator* config) {
  TString tExportType;
  try {
    mNumberOfEvents = (config->GetParameter("NumberOfEvents")).Atoi();
    mIntegrateSample = (config->GetParameter("IntegrateSamples")).Atoi();
    mEventDIR = config->GetParameter("EventDir");
    tExportType = config->GetParameter("EventFileType");
  } catch (TString tError) {
    PRINT_MESSAGE("<EventGenerator::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  mLogName = config->GetParameter("LogFile", "");
  if (tExportType == "root")
    mEventExportType = 0;
  else if (tExportType == "text")
    mEventExportType = 1;
  else if (tExportType == "root&text")
    mEventExportType = 2;
}

/*
   @brief Constructs Id number for the event and saves the result
*/
void EventGenerator::SetEventID(int aEventIter) {
  ostringstream oss;
  Crc32 tEventID;

  oss << mTimeStamp.Data() << "Event: " << aEventIter;
  tEventID.Update(oss.str().data(), oss.str().length());
  tEventID.Finish();

  mEventID = tEventID.GetValue();
}

/*
   @brief Returns the saved EventID
*/
unsigned int EventGenerator::GetEventID() const {
  return mEventID;
}

void EventGenerator::Print(const std::string& msg) const {
  std::cout << "[" << mRunID << "]\t" << msg << std::endl;
}

void EventGenerator::AddLogEntry(const std::string& aEntry) {
  TString tLogName;
  TDatime tDate;
  std::ofstream tFile;

  tFile.open(mLogName, std::ios_base::app);
  if (static_cast<long>(tFile.tellp()) == 0) {
    tFile << "# THERMINATOR 2 Log File" << std::endl;
  }
  tFile << '[' << tDate.AsSQLString() << "]\t" << mRunID << "\t" << aEntry << std::endl;
  tFile.close();
}
