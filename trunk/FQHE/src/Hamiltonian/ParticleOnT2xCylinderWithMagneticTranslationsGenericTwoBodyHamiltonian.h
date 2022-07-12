////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hamiltonian associated to particles on a             //
//        4D space  torus x cylinder with a generic two-body interaction      //
//                          and magnetic translations                         //
//                                                                            //
//                        last modification : 07/03/2017                      //
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


#ifndef PARTICLEONT2XCYLINDERWITHMAGNETICTRANSLATIONSGENERICTWOBODYHAMILTONIAN_H
#define PARTICLEONT2XCYLINDERWITHMAGNETICTRANSLATIONSGENERICTWOBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian.h"
#include "Polynomial/Polynomial.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian : public ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian
{

 protected:

  // ratio between the length in the x direction and the length in the y direction of the cylinder
  double CylinderRatio; 

 public:

  // default constructor
  //
  ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuantumTorus = number of flux quanta piercing the torus
  // kxMomentum = momentum in the x direction for the torus
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // nbrFluxQuantumCylinder = number of flux quanta piercing the cylinder
  // cylinderRatio = ratio between the length in the x direction and the length in the y direction of the cylinder
  // nbrPseudoPotentials = number of pseudo-potentials
  // pseudoPotentialMomentumTorus = torus pseudo-potential relative momenta
  // pseudoPotentialMomentumCylinder = cylinder pseudo-potential relative angular momenta
  // pseudoPotentials = pseudo-potential coefficients
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int nbrFluxQuantumTorus, int xMomentum,
									 double ratio, int nbrFluxQuantumCylinder, double cylinderRatio, int nbrPseudoPotentials, 
									 int* pseudoPotentialMomentumTorus, int* pseudoPotentialMomentumCylinder, 
									 double* pseudoPotentials, AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

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
  // exponentialCoefficient = array that contains the precomputed exponential factors for the cylinder  
  // pseudopotentialIndex1 = pseudopotential index for the interaction on the torus 
  // pseudopotentialIndex2 = pseudopotential index for the interaction on the cylinder
  // pseudoPotential = pseudopotential amplitude
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int ky1, int ky2, int ky3, int ky4, 
						int lz1, int lz2, int lz3, int lz4, 
						double* exponentialCoefficient, 
						int pseudopotentialIndex1, int pseudopotentialIndex2, double pseudoPotential);

  // evaluate the numerical coefficient  in front of the a+_lz1 a+_lz2 a_lz3 a_lz4 coupling term for the cylinder part
  //
  // lz1 = first lz index
  // lz2 = second lz index
  // lz3 = third lz index
  // lz4 = fourth lz index
  // pseudopotentialIndex = pseudopotential index for the interaction on the cylinder 
  // exponentialCoefficient = array that contains the precomputed exponential factors for the cylinder  
  // return value = numerical coefficient
  virtual double EvaluateCylinderInteractionCoefficient(int lz1, int lz2, int lz3, int lz4, int pseudopotentialIndex, double* exponentialCoefficient);


};

#endif
