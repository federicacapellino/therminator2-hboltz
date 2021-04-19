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

#include <stdio.h>
#include <fstream>
#include <TDatime.h>
#include <TString.h>
#include <TMath.h>
#include "Integrator.h"
#include "THGlobal.h"

#include "Model.h"

using namespace std;

Integrator::Integrator()
: mFOModel(0),  mNSamples(0)
{
}

Integrator::Integrator(int aNSamples, Model *FOModel)
: mNSamples(aNSamples)
{ 
  mFOModel = FOModel;
}

Integrator::~Integrator()
{
   // not owner of FOModel
};

Model* Integrator::GetModel()
{
  return mFOModel;
}

void Integrator::GenerateParticles(ParticleType* aPartType, int aPartCount, list<Particle>* aParticles)
{
  int      tIter = 0;
  double   tFMax;
  double   tVal;
  double   tValTest;
  Particle* tParticle; 
  
  PRINT_DEBUG_3("Integrator::GenerateParticles\t "<< aPartType->GetName() <<" B:"<< aPartType->GetBarionN() <<" I3:"<< aPartType->GetI3() <<" S:"<< aPartType->GetStrangeN()<<" C:"<< aPartType->GetCharmN());

  // Use accept/reject to sample the phase space.  This uses the previously
  // recorded maximum  of the integrand (see SetMultiplicities) which is
  // returned by get integrand.
  tFMax = aPartType->GetMaxIntegrand();  
  while (tIter < aPartCount) {
    tVal      = mFOModel->GetIntegrand(aPartType);
    tValTest  = gRandom->Rndm() * tFMax;
    if (tValTest<tVal) {
      tParticle = new Particle(aPartType);
      mFOModel->SetParticlePX(tParticle);
      aParticles->push_back(*tParticle);
      tIter++;
      delete tParticle;
    }
  }
}

bool Integrator::SetMultiplicities(ParticleDB *aDB, const TString &filename)
{
  // Make or read table with probabilities
  ifstream tFileIn;
  ofstream tFileOut;
  TDatime  tDate;
  char     tBuff[2*kFileNameMaxChar];
  char     tMultiName[kFileNameMaxChar];
  char     tPart[100];
  double   tMaxInt;
  double   tMulti;

  sprintf(tMultiName,"%s_%s.txt",filename.Data(),mFOModel->GetHash());

  tFileIn.open(tMultiName);
  if ((tFileIn) && (tFileIn.is_open())) {

    PRINT_DEBUG_1("<Integrator::SetMultiplicities>\tReading Max Integrand and Multiplicity values from " << tMultiName);
    
    // Read results from file
    while (!tFileIn.eof()) {
      *tPart = 0;
      tFileIn >> tPart;
      if(!(*tPart) || (*tPart == '#')) {
        tFileIn.getline(tPart,100);
        continue;
      }
      tFileIn >> tMaxInt >> tMulti;
      aDB->GetParticleType(tPart)->SetMaxIntegrand(tMaxInt);
      aDB->GetParticleType(tPart)->SetMultiplicity(tMulti);
      PRINT_DEBUG_2("\t"<<tPart << " " << tMaxInt << " " << tMulti);
    }

    tFileIn.close();
    sprintf(tBuff,"[input]\t%s",tMultiName);
    mFileDescription = tBuff;
    return true;

  } else {

    // Compute and write the data to a file
    char tTempName[kFileNameMaxChar];
    sprintf(tTempName,"%s_%s.tmp",filename.Data(),mFOModel->GetHash());
    tFileOut.open(tTempName);
    if ((tFileOut) && (tFileOut.is_open())) {
      PRINT_DEBUG_1("<Integrator::SetMultiplicities>\tMax Integrand and Multiplicity file " << tMultiName << " not found.");
      tDate.Set();
      PRINT_MESSAGE("["<<tDate.AsSQLString() << "]\tCalculating Max Integrand and Multiplicity");
      tFileOut << mFOModel->GetDescription();
      tFileOut << endl;
      tFileOut << "# Particle name\tMax integrand\tMultiplicity" << endl;

      // Loop over particles writing the FO function maximum and multiplicity
      for(int tIter = 0; tIter < aDB->GetParticleTypeCount(); tIter++) {
        // Main step of sampling the fo function many times for a given particle
        tMulti  = Integrate(aDB->GetParticleType(tIter));
        tMaxInt = aDB->GetParticleType(tIter)->GetMaxIntegrand();

        // Write the results to the multiplicity files
        tFileOut << aDB->GetParticleType(tIter)->GetName() << '\t' << tMaxInt << '\t' << tMulti << endl;
        cout << "\r\tparticle ("<< (tIter + 1) << "/" << aDB->GetParticleTypeCount() << "): " << aDB->GetParticleType(tIter)->GetName();
        cout.flush();
      }
      cout << endl;
      tFileOut.close();
      if(rename(tTempName,tMultiName) != 0) {
        PRINT_MESSAGE("<Integrator::SetMultiplicities>\tUnable to rename temp file to " << tMultiName);
        exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
      }
    } else {
      PRINT_MESSAGE("<Integrator::SetMultiplicities>\tUnable to create file " << tMultiName);
      exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
    }

    // The job is done. Record the result to the log file.
    sprintf(tBuff,"[output]\t%s_%s.txt",filename.Data(),mFOModel->GetHash());
    mFileDescription = tBuff;
    return false ;
  }
  return false ; // No path to this statement
}

double Integrator::Integrate(ParticleType* aPartType)
{
  double tMaxInt;
  double tMulti;
  double tVal;
  int    tIter;
  
  // For the given particle determine the Max of integrand and the compute the
  // multiplicity. This uses the member frunctions of the freezeout model.
  tMaxInt = 0.0;
  tMulti  = 0.0;

  // - Generate mNSamples over a given hypercube 
  // - Evaluate the integrand at each of these points
  // - The mean multiplicity is the HyperCubeVolume/NSamples
  for (tIter = 0; tIter < mNSamples; tIter++) {
    tVal = mFOModel->GetIntegrand(aPartType);
    if (tVal>tMaxInt)
      tMaxInt = tVal;
    tMulti += tVal;
  }
  tMulti *= mFOModel->GetHyperCubeVolume() / (1.0 * mNSamples); 
  aPartType->SetMaxIntegrand(tMaxInt);
  aPartType->SetMultiplicity(tMulti);

  return tMulti;
}
