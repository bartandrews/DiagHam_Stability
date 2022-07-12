////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
// antiparallel magnetic field, coulomb interaction and magnetic translations //
//                                                                            //
//                        last modification : 30/05/2018                      //
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


#ifndef PARTICLEONTORUSCOULOMBWITHSPINANDMAGNETICTRANSLATIONSTIMEREVERSALSYMMETRICHAMILTONIAN_H
#define PARTICLEONTORUSCOULOMBWITHSPINANDMAGNETICTRANSLATIONSTIMEREVERSALSYMMETRICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Polynomial;


class ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian : public ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian
{

 protected:

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // ratio = ratio between the width in the x direction and the width in the y direction
  // layerSeparation = layer separation in units of magnetic length
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
							    double ratio, double layerSeparation, AbstractArchitecture* architecture, int memory = -1, char* precalculationFileName = 0);

  // constructor from pseudopotentials
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // ratio = ratio between the width in the x direction and the width in the y direction
  // nbrPseudopotentialsUpUp = number of pseudopotentials for up-up interaction
  // pseudopotentialsUpUp = pseudopotential coefficients for up-up interaction
  // nbrPseudopotentialsDownDown = number of pseudopotentials for down-down interaction
  // pseudopotentialsDownDown = pseudopotential coefficients for down-down interaction
  // nbrPseudopotentialsUpDown = number of pseudopotentials for up-down interaction
  // pseudopotentialsUpDown = pseudopotential coefficients for up-down interaction
  // spinFluxUp = additional inserted flux for spin up
  // spinFluxDown = additional inserted flux for spin down
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
								   double ratio, 
								   int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
								   int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
								   int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
								   double spinFluxUp, double spinFluxDown, 
								   AbstractArchitecture* architecture, long memory, char* precalculationFileName, 
								   double* oneBodyPotentielUpUp, double* oneBodyPotentielDownDown, double* oneBodyPotentielUpDown);
  
  // destructor
  //
  ~ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);
  
  
  // get all the indices that should appear in the annihilation/creation operators
  //
  virtual void GetIndices();

 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // layerSeparation = separation of layers
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficientDownDown(int m1, int m2, int m3, int m4, double layerSeparation=0.0);
  
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // layerSeparation = separation of layers
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficientUpDown(int m1, int m2, int m3, int m4, double layerSeparation=0.0);


};

#endif // PARTICLEONTORUSCOULOMBWITHSPINANDMAGNETICTRANSLATIONSHAMILTONIAN_H
