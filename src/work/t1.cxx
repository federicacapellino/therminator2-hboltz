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

#include "Configurator.h"
#include "Particle.h"
#include "ParticleDB.h"
#include "ParticleDecayer.h"
#include "ParticleType.h"
#include "THGlobal.h"
#include <TFile.h>
#include <TH1D.h>
#include <TRandom.h>
#include <cmath>
#include <iostream>
#include <list>
#include <memory>

double legendre(unsigned l, double x) {
  double p = 0.;
  switch (l) {
  case 0:
    p = 1.;
    break;
  case 1:
    p = x;
    break;
  case 2:
    p = (3. * x * x - 1.) / 2.;
    break;
  case 3:
    p = (5. * x * x * x - 3. * x) / 2.;
    break;
  default:
    std::cout << "Bad l in legendre: " << l << std::endl;
    break;
  }
  return p * std::sqrt(2. * l + 1.);
}

double weight(Particle* prt, double m, double temper, int icase) {
  double Spin = prt->GetParticleType()->GetSpin();
  double Statistics = ((Spin - static_cast<int>(Spin)) < 0.01) ? -1.0 : +1.0;
  double e = prt->e;
  double u = std::exp(-e / temper);
  double p = prt->GetP();
  double f = u / (1. - Statistics * u);
  double phis = 1. + Statistics * f;
  double val = (2. * Spin + 1.) * f * p * p / kTwoPi3;
  double x = prt->pz / p;
  switch (icase) {
  case 0:
    return val;
  case 1:
    return val * legendre(1, x) * p / e;
  case 2:
    return val * legendre(2, x) * p * p * phis;
  case 3:
    return val * legendre(3, x) * p * p * p / e * phis;
  case 4:
    return val * legendre(1, x) * 2. / 5. * p * p * p / e * phis;
  }
  return 0.;
}

int main(int argc, char** argv) {
  std::string iniFile = "./events.ini";
  auto tConfig = std::make_unique<Configurator>(iniFile);

  auto tPartDB = std::make_unique<ParticleDB>(tConfig->GetParameter("ShareDir", ""));

  std::list<Particle> tPList;
  ParticleDecayer tPDecay(tPartDB.get(), &tPList);

  auto tPType = tPartDB->GetParticleType("om0782zer");
  std::cout << "Processing decays of " << tPType->GetName() << std::endl;

  TFile fout("funcs.root", "RECREATE", "a test file");
  TH1::SetDefaultSumw2();

  TH1D h0input("om0782zer_dNdp0", "dN/dp l=0", 100, 0., 3.);
  TH1D h0output("pi0139plu_from_om0782zer_dNdp0", "dN/dp l=0", 100, 0., 3.);
  TH1D h1input("om0782zer_dNdp1", "dN/dp l=1", 100, 0., 3.);
  TH1D h1output("pi0139plu_from_om0782zer_dNdp1", "dN/dp l=1", 100, 0., 3.);
  TH1D h2input("om0782zer_dNdp2", "dN/dp l=2", 100, 0., 3.);
  TH1D h2output("pi0139plu_from_om0782zer_dNdp2", "dN/dp l=2", 100, 0., 3.);
  TH1D h3input("om0782zer_dNdp3", "dN/dp l=3", 100, 0., 3.);
  TH1D h3output("pi0139plu_from_om0782zer_dNdp3", "dN/dp l=3", 100, 0., 3.);
  TH1D h1input_phis("om0782zer_dNdp1_phis", "dN/dp l=1 phi3", 100, 0., 3.);
  TH1D h1output_phis("pi0139plu_from_om0782zer_dNdp1_phis", "dN/dp l=1 phis", 100, 0., 3.);

  constexpr double pmin = 0.;
  constexpr double pmax = 4.0;
  constexpr int nsamples = 40000000;
  constexpr double temper = 0.160;
  constexpr int ncases = 5;

  Particle tPart(tPType);
  for (int i = 0; i < nsamples; i++) {
    double p = gRandom->Uniform(pmin, pmax);
    double cth = gRandom->Uniform(-1., 1.);
    double sth = std::sqrt(1. - cth * cth);
    double phi = gRandom->Uniform(0., 2. * M_PI);
    double m = tPType->GetMass();
    double ep = std::sqrt(p * p + m * m);
    tPart.SetParticlePX(ep, p * sth * std::cos(phi), p * sth * std::sin(phi), p * cth, 0, 0, 0, 0);

    double wFather[ncases] = {};
    double x = p * cth / p; // = cth
    for (int icase = 0; icase < ncases; icase++) {
      wFather[icase] = weight(&tPart, m, temper, icase);
      switch (icase) {
      case 0:
        h0input.Fill(p, wFather[0] * legendre(0, x));
        break;
      case 1:
        h1input.Fill(p, wFather[1] * legendre(1, x));
        break;
      case 2:
        h2input.Fill(p, wFather[2] * legendre(2, x));
        break;
      case 3:
        h3input.Fill(p, wFather[3] * legendre(3, x));
        break;
      case 4:
        h1input_phis.Fill(p, wFather[4] * legendre(1, x));
        break;
      }
    }

    tPDecay.DecayParticle(&tPart);

    for (auto& daughter : tPList) {
      if (daughter.GetParticleType()->GetPDGCode() != 211)
        continue;

      double ep1, px1, py1, pz1;
      daughter.GetMomentum(&ep1, &px1, &py1, &pz1);
      double p1 = daughter.GetP();
      double cth1 = pz1 / p1;

      for (int icase = 0; icase < ncases; icase++) {
        double wD = legendre((icase == 4) ? 1 : icase, cth1);
        switch (icase) {
        case 0:
          h0output.Fill(p1, wFather[0] * wD);
          break;
        case 1:
          h1output.Fill(p1, wFather[1] * wD);
          break;
        case 2:
          h2output.Fill(p1, wFather[2] * wD);
          break;
        case 3:
          h3output.Fill(p1, wFather[3] * wD);
          break;
        case 4:
          h1output_phis.Fill(p1, wFather[4] * wD);
          break;
        }
      }
    }
    tPList.clear();
  }

  double sumw = h0input.GetSumOfWeights();
  double volume = 4. * M_PI * pmax;
  std::cout << "Number of omegas: " << sumw * volume / nsamples << std::endl;

  h0output.Scale(1. / sumw);
  h1output.Scale(1. / sumw);
  h2output.Scale(1. / sumw);
  h3output.Scale(1. / sumw);
  h1output_phis.Scale(1. / sumw);

  fout.Write();
  fout.Close();
  return 0;
}
