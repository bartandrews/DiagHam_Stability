////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//     class of hamiltonian with particles on the 4D manifold T2 x T2         //
//                     with 4D generic two-body interaction                   //
//                                                                            //
//                        last modification : 15/02/2017                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef PARTICLEONT2XT2GENERICTWOBODYHAMILTONIAN_H
#define PARTICLEONT2XT2GENERICTWOBODYHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereFullHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;
class Polynomial;


class ParticleOnT2xT2GenericTwoBodyHamiltonian : public AbstractQHEOnSphereFullHamiltonian
{

 protected:
  
  // number of flux quanta for the first torus
  int NbrFluxQuanta1;
  // number of flux quanta for the second torus
  int NbrFluxQuanta2;

  // ratio between the width in the x direction and the width in the y direction for the first torus
  double Ratio1;
  // inverse ratio between the width in the y direction and the width in the x direction for the first torus
  double InvRatio1;
  // ratio between the width in the x direction and the width in the y direction for the second torus
  double Ratio2;
  // inverse ratio between the width in the y direction and the width in the x direction for the second torus
  double InvRatio2;

  // number of pseudo-potentials
  int NbrPseudoPotentials;
  // pseudo-potential first sphere relative angular momenta 
  int* PseudoPotentialMomentum1;
  // pseudo-potential second sphere relative angular momenta
  int* PseudoPotentialMomentum2;
  //pseudo-potential coefficients
  double* PseudoPotentials;

  // Laguerre polynomial for the pseudopotentials
  Polynomial* LaguerrePolynomials;

 public:

  // default constructor
  //
  ParticleOnT2xT2GenericTwoBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuanta1 = number of flux quanta for the first torus
  // nbrFluxQuanta2 = number of flux quanta for the second torus  
  // ratio1 = ratio between the width in the x direction and the width in the y directiona for the first torus
  // ratio2 = ratio between the width in the x direction and the width in the y directiona for the second torus
  // nbrPseudoPotentials = number of pseudo-potentials
  // pseudoPotentialMomentum1= first torus pseudo-potential relative angular momenta
  // pseudoPotentialMomentum2 = second torus pseudo-potential relative angular momenta
  // pseudoPotentials = pseudo-potential coefficients
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnT2xT2GenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta1, int nbrFluxQuanta2, 
					   double ratio1, double ratio2, 
					   int nbrPseudoPotentials, int* pseudoPotentialMomentum1, int* pseudoPotentialMomentum2, 
					   double* pseudoPotentials, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnT2xT2GenericTwoBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // get a linearized index from the two angular momenta
  //
  // lz = z projection of the angular momentum for the first sphere
  // kz = z projection of the angular momentum for the second sphere
  // return value = linearized index 
  virtual int GetLinearizedIndex(int lz, int kz);

  // get the two angular momenta associated to a given linearized index
  //
  // index = linearized index 
  // lz = reference on the z projection of the angular momentum for the first sphere
  // kz = reference on the  z projection of the angular momentum for the second sphere
  virtual void GetLinearizedIndex(int index, int& lz, int& kz);

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrFluxQuanta = number of flux quanta
  // ratio = ratio between the width in the x direction and the width in the y direction
  // pseudopotentialMomentum = pseudo-potential relative angular momenta
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrFluxQuanta, 
						double ratio, int pseudopotentialMomentum);

  // evaluate the numerical coefficient  in front of the a+_(lz1,kz1) a+_(lz2,kz2) a_(lz3,kz3) a_(lz4,kz4) coupling term
  //
  // lz1 = first lz index
  // lz2 = second lz index
  // lz3 = third lz index
  // lz4 = fourth lz index
  // kz1 = first kz index
  // kz2 = second kz index
  // kz3 = third kz index
  // kz4 = fourth kz index
  // pseudopotentialIndex1 = pseudopotential index for the interaction on the first cylinder 
  // pseudopotentialIndex2 = pseudopotential index for the interaction on the second cylinder 
  // pseudoPotential = pseudopotential amplitude
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int lz1, int lz2, int lz3, int lz4, 
						int kz1, int kz2, int kz3, int kz4, 
						int pseudopotentialIndex1, int pseudopotentialIndex2, double pseudoPotential);
};

// get a linearized index from the two angular momenta
//
// lz = z projection of the angular momentum for the first sphere
// kz = z projection of the angular momentum for the second sphere
// return value = linearized index 

inline int ParticleOnT2xT2GenericTwoBodyHamiltonian::GetLinearizedIndex(int lz, int kz)
{
  return ((lz * this->NbrFluxQuanta2) + kz);
}

// get the two angular momenta associated to a given linearized index
//
// index = linearized index 
// lz = reference on the z projection of the angular momentum for the first sphere
// kz = reference on the  z projection of the angular momentum for the second sphere

inline void ParticleOnT2xT2GenericTwoBodyHamiltonian::GetLinearizedIndex(int index, int& lz, int& kz)
{
  lz = index / this->NbrFluxQuanta2;
  kz = index % this->NbrFluxQuanta2;
}

#endif
