////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//             generic interaction defined by its pseudopotential             //
//                                                                            //
//                        last modification : 03/06/2004                      //
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


#ifndef PARTICLEONSPHEREGENERICANISOTROPICHAMILTONIAN_H
#define PARTICLEONSPHEREGENERICANISOTROPICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnSphereGenericAnisotropicHamiltonian : public AbstractQHEOnSphereHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  double* PseudoPotential;

  // array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
  double* OneBodyPotentials;

  // array with the anisotropic pseudo-potentials alpha=2 (last element corresponds to the delta interaction)
  double* AnisotropicPseudoPotentialAlpha2;

  // array with the anisotropic pseudo-potentials alpha=4 (last element corresponds to the delta interaction)
  double* AnisotropicPseudoPotentialAlpha4;

  //chirality of the anisotropic pseudopotential
  int Chirality;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // architecture = architecture to use for precalculation
  // pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
  // anisotropicPseudoPotentialAlpha2 = array with the anisotropic pseudo-potentials alpha=2 (ordered such that the first element corresponds to the delta interaction)
  // anisotropicPseudoPotentialAlpha4 = array with the anisotropic pseudo-potentials alpha=4 (ordered such that the first element corresponds to the delta interaction)
  // chirality = chirality of the anisotropic pseudopotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  // hermitianFlag = flag to indicate if hermitian symmetry of Hamiltonian shall be used
  ParticleOnSphereGenericAnisotropicHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, double* pseudoPotential, double* anisotropicPseudoPotential2, double* anisotropicPseudoPotentialAlpha4, int chirality, double l2Factor,
				     AbstractArchitecture* architecture, long memory = -1, 
				     bool onDiskCacheFlag = false, char* precalculationFileName = 0, bool hermitianFlag = false);

  // constructor with one body terms
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // architecture = architecture to use for precalculation
  // pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
  // oneBodyPotentials = array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
  // anisotropicPseudoPotentialAlpha2 = array with the anisotropic pseudo-potentials alpha=2 (ordered such that the first element corresponds to the delta interaction)
  // anisotropicPseudoPotentialAlpha4 = array with the anisotropic pseudo-potentials alpha=4 (ordered such that the first element corresponds to the delta interaction)
  // chirality = chirality of the anisotropic pseudopotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  // hermitianFlag = flag to indicate if hermitian symmetry of Hamiltonian shall be used
  ParticleOnSphereGenericAnisotropicHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
				     double* pseudoPotential, double* oneBodyPotentials, double* anisotropicPseudoPotentialAlpha2, double* anisotropicPseudoPotentialAlpha4, int chirality, double l2Factor,
				     AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,
				     char* precalculationFileName, bool hermitianFlag = false);
    
  // destructor
  //
  ~ParticleOnSphereGenericAnisotropicHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  List<Matrix*> RightInteractionOperators();  


 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();


};

#endif
