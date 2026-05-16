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

#ifndef _TH2_INTEGRATOR_H_
#define _TH2_INTEGRATOR_H_

#include "Model.h"
#include "Particle.h"
#include "ParticleDB.h"
#include "ParticleType.h"
#include <TRandom2.h>
#include <list>

class Integrator {
public:
  Integrator();
  Integrator(int aNSamples, Model* FOModel);
  ~Integrator();

  Model* GetModel();
  void GenerateParticles(ParticleType* aPartType, int aPartCount, std::list<Particle>* aParticles);
  bool SetMultiplicities(ParticleDB* aDB, const TString& filename = "fmultiplicity");
  std::string GetFileDescription() { return mFileDescription; }

private:
  double Integrate(ParticleType* aPartType);

  Model* mFOModel;
  int mNSamples;
  std::string mFileDescription;
};

#endif

/*! @file Integrator.h
 * @brief Definition of Integrator class. Integrates the Cooper-Frye formula and generates
 * primordial particles.
 */

/*! @class Integrator
 * @brief Integrates the Cooper-Frye formula using a Monte-Carlo method
 * and randomly generates primordial particles on a given freeze-out hypersurface.
 *
 * @fn Integrator::Integrator()
 * @brief Default constructor.
 *
 * @fn Integrator::Integrator(int aNSamples,Model *FOModel)
 * @brief Sets the number of Monte-Carlo samples to determine the maximal value of the integrand and
 average multiplicity  of the primordial particles.
 * @param [in] aNSamples number of Monte-Carlo samples
 * @param [in] FOModel, a previously allocated Model class
 *
 * @fn Integrator::~Integrator()
 * @brief Destructor
 *
 * @fn Model* Integrator::GetModel()
 * @brief Returns a pointer to the currently used Model object.
 *
 * @fn void Integrator::GenerateParticles(ParticleType* aPartType, int aPartCount,
 std::list<Particle>* aParticles)
 * @brief Generates a given number of primordial particles of a given type and stores them in the
 list-type container.
 * @param [in] aPartType type of particle
 * @param [in] aPartCount number of particles to be generated
 * @param [out] aParticles std::list container object that stores all generated particles
 *
 * @fn bool Integrator::SetMultiplicities(ParticleDB* aDB, const TString &filename)
 * @brief Read or write the particle database with appropriate average multiplicity and maximum
 value of the integrand. Returns true if the database is read, and false if the values are computed
 and the file is written.
 *
 Based on the model parameters a CRC32 hash (such as hash = "AAA1E9AE") is
 created. If a file with name,  filename_hash.txt, with such hash exists
 then the multiplicities and  max integrand data is taken from the  file.
 Otherwise these numbers are calculate with Monte-Carlo with the number
 of samples defined.

 * @param [in] aDB particle data-base. See ParticleDB class.
 * @param [in] a string containing the filename
 *
 */
