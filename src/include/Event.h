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

#ifndef _TH2_EVENT_H_
#define _TH2_EVENT_H_

#include "Integrator.h"
#include "Particle.h"
#include "ParticleDB.h"
#include <TRandom2.h>
#include <list>
#include <vector>

class Event {

public:
  Event();
  Event(ParticleDB* aDB, Integrator* aInteg);
  ~Event();

  void Reset();
  std::list<Particle>& GetParticleList();
  Integrator* GetIntegrator() const;
  ParticleDB* GetParticleDB() const;

  void GeneratePrimordials();
  void DecayParticles();
  void SetDistribution(const std::string& name = "Poisson");

private:
  void GenerateMultiplicities();

  std::list<Particle> mParticles;
  std::vector<int> mMultiplicities;
  ParticleDB* mPartDB;
  Integrator* mInteg;
  int mDistribution; // type of multiplicity distribution: 0 = Poissonian, 1 - NegativeBinomial
};

#endif

/*! @file Event.h
@brief Definition of Event class. Generates primordial, decays resonances and stores particles.
 */

/*! @class Event
@brief Generates primordial particles and decays resonances. Created Particle is stored in a
list-type container class.

Primordial particles are generated with a Monte-Carlo method by a Model via Integrator class.

The type of distribution that is generated can be set to "Poisson" or "NegativeBinomial"

Event randomly generates the average multiplicity of primordial particles according to the given
distribution, by default with Poissonian.

@fn Event::Event()
@brief Default constructor.

@fn Event::Event(ParticleDB* aDB, Integrator* aInteg)
@brief Creates an Event with particles from database and generated with given Model via Integrator.
@param [in] aDB pointer to database with ParticleType to be generated in event.
@param [in] aInteg pointer to Integrator that uses Model to generate particles.

@fn Event::~Event()
@brief Destructor.

@fn void Event::Reset(int aEventIter=0)
@brief Resets the Event. Calculate eventID (CRC-32) from the current date and event iterator.
@param [in] aEventIter event iterator

@fn std::list<Particle>* Event::GetParticleList()
@brief Returns the pointer to Particle list created in Event.

@fn Integrator* Event::GetIntegrator() const
@brief Returns the pointer to Integrator.

@fn ParticleDB* Event::GetParticleDB() const
@brief Returns the pointer to ParticleDB.

@fn unsigned int Event::GetEventID() const;
@brief Returns the event ID.

@fn void Event::GeneratePrimordials()
@brief Generates primordial particles and stores them in the particle list.

See Integrator::GenerateParticles() for details.<br />

@fn void Event::DecayParticles(int aSeed=0)
@brief Uses a ParticleDecayer object to decay all unstable particles in the Particle list.

See ParticleDecayer::DecayParticle() for details. <br />
Created child particles are appended to the end of the list. Decay cascade continues until all
unstable  particles decay and ParticleDecayer reaches the end of the list.

@fn void Event::SetDistribution(const std::string &name)
@brief Sets the distribution for choosing the number of particles, given the mean number.

@param[in] name Selects the distribution, possible values are "Poisson" or "NegativeBinomial"
 */
