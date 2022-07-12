////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//             class of hamiltonian with particles on the 4D manifold         //
//                   T2 x S2 with 4D two body generic interaction             //
//                                                                            //
//                        last modification : 09/03/2017                      //
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


#ifndef PARTICLEONT2XS2TWOBODYGENERICHAMILTONIAN_H
#define PARTICLEONT2XS2TWOBODYGENERICHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnS2xS2GenericTwoBodyHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Polynomial/Polynomial.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnT2xS2GenericTwoBodyHamiltonian : public ParticleOnS2xS2GenericTwoBodyHamiltonian
{

 protected:
  
  // ratio between the length in the x direction and the length in the y direction of the torus
  double Ratio; 

  // Laguerre polynomial for the pseudopotentials
  Polynomial* LaguerrePolynomials;

 public:

  // default constructor
  //
  ParticleOnT2xS2GenericTwoBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuantumTorus = number of flux quanta piercing the torus
  // ratio = ratio between the length in the x direction and the length in the y direction of the torus
  // nbrFluxQuantumSphere = number of flux quanta piercing the sphere
  // nbrPseudoPotentials = number of pseudo-potentials
  // pseudoPotentialMomentumTorus = torus pseudo-potential relative momenta
  // pseudoPotentialMomentumSphere = sphere pseudo-potential relative angular momenta
  // pseudoPotentials = pseudo-potential coefficients
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnT2xS2GenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuantumTorus, double ratio, int nbrFluxQuantumSphere, 
					   int nbrPseudoPotentials, int* pseudoPotentialMomentumTorus, int* pseudoPotentialMomentumSphere, 
					   double* pseudoPotentials, AbstractArchitecture* architecture, long memory = -1);
  
  // destructor
  //
  ~ParticleOnT2xS2GenericTwoBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // evaluate the numerical coefficient in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term for the torus part
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrFluxQuanta = number of flux quanta
  // ratio = ratio between the width in the x direction and the width in the y direction
  // pseudopotentialMomentum = pseudo-potential relative momentum
  // return value = numerical coefficient
  virtual double EvaluateTorusInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrFluxQuanta, 
						     double ratio, int pseudopotentialMomentum);

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term for the sphere [art
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // clebsch = reference on the Clebsch-Gordan coefficients for the sphere geometry
  // nbrFluxQuanta = number of flux quanta
  // pseudopotentialMomentum = pseudo-potential relative momentum
  // return value = numerical coefficient
  virtual double EvaluateSphereInteractionCoefficient(int m1, int m2, int m3, int m4, ClebschGordanCoefficients& clebsch,
						      int nbrFluxQuanta,int pseudopotentialMomentum);

  // evaluate the numerical coefficient  in front of the a+_(lz1,kz1) a+_(lz2,kz2) a_(lz3,kz3) a_(lz4,kz4) coupling term
  //
  // ky1 = first ky index
  // ky2 = second ky index
  // ky3 = third ky index
  // ky4 = fourth ky index
  // lz1 = first lz index
  // lz2 = second lz index
  // lz3 = third lz index
  // lz4 = fourth lz index
  // clebsch = reference on the Clebsch-Gordan coefficients for the sphere geometry
  // pseudopotentialIndex1 = pseudopotential index for the interaction on the torus 
  // pseudopotentialIndex2 = pseudopotential index for the interaction on the sphere
  // pseudoPotential = pseudopotential amplitude
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int ky1, int ky2, int ky3, int ky4, 
						int lz1, int lz2, int lz3, int lz4, 
						ClebschGordanCoefficients& clebsch,
						int pseudopotentialIndex1, int pseudopotentialIndex2, double pseudoPotential);

  // get all the indices that should appear in the annihilation/creation operators
  //  
  virtual void GetIndices();
 
};

#endif
