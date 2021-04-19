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
#include <TDatime.h>
#include "Crc32.h"
#include "Configurator.h"
#include "ParticleDecayer.h"
#include "Event.h"
#include "THGlobal.h"


using namespace std;

Event::Event()
: mPartDB(0), mInteg(0),  mDistribution(0)
{
  mMultiplicities.clear();
  Reset();
}

Event::Event(ParticleDB* aDB, Integrator* aInteg)
: mPartDB(aDB), mInteg(aInteg), mDistribution(0)
{ 
  mMultiplicities.clear();
  mMultiplicities.resize(mPartDB->GetParticleTypeCount());
  Reset();
}

Event::~Event()
{
  mParticles.clear();
  mMultiplicities.clear();
}

void Event::Reset()
{
  mParticles.clear();
  Particle::ZeroEID();
}

list<Particle>* Event::GetParticleList()
{
  return &mParticles;
}

Integrator* Event::GetIntegrator() const
{
  return mInteg;
}

ParticleDB* Event::GetParticleDB() const
{
  return mPartDB;
}


void Event::GeneratePrimordials()
{ 
  GenerateMultiplicities();
  for (int tIter=0; tIter<mPartDB->GetParticleTypeCount(); tIter++) {
    if(! strstr(mPartDB->GetParticleType(tIter)->GetName(),"gam000zer")) { 
      // Disable primordial photon generation This is a not a photon, generate
      // the the particle,  and store them in mParticles
      mInteg->GenerateParticles(mPartDB->GetParticleType(tIter), mMultiplicities[tIter], &mParticles);

    } else {
      // disable primordial photons generation
      continue;
    }
  }
}

void Event::DecayParticles()
{
  list<Particle>::iterator tIter;
  ParticleType*    tFatherType;
  ParticleDecayer* tDecayer;
  
  tDecayer = new ParticleDecayer(mPartDB, &mParticles);

  tIter = mParticles.begin();
// as new particles are added from decays the end() of the list moves until all particles had decayed
  do {
    tFatherType = tIter->GetParticleType();
    // if not stable or stable but has a decay table with at least one decay channel
    if((tFatherType->GetGamma() >= 0.0) && (tFatherType->GetTable()) && (tFatherType->GetTable()->GetChannelCount() + 1 > 0))
      tDecayer->DecayParticle( &(*tIter) );
    tIter++;
  } while (tIter != mParticles.end());
  delete tDecayer;
}

void Event::GenerateMultiplicities()
{
  if(mDistribution == 0) { // Poisson
    for (int tIter=0; tIter<mPartDB->GetParticleTypeCount(); tIter++)
      mMultiplicities[tIter] = gRandom->Poisson(mPartDB->GetParticleType(tIter)->GetMultiplicity());
  } else if(mDistribution == 1) { // Negative Binomial
    for (int tIter=0; tIter<mPartDB->GetParticleTypeCount(); tIter++)
      mMultiplicities[tIter] = 0; // HOW?
  }
}

void Event::SetDistribution(const std::string &name) 
{
   if (name == "NegativeBinomial") {
      mDistribution=1;
   } else if (name == "Poisson") {
      mDistribution=0;
   } else {
      PRINT_MESSAGE("Event::SetDistribution: bad selection in switch. Setting the distribution to Poisson.") ;
      mDistribution=0;
   }
}

