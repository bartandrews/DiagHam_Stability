////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                   Class author : Cecile Repellin                           //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//generic two body interaction and a time reversal symmetric degree of freedom//
//                                                                            //
//                        last modification : 10/02/2015                      //
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

#ifndef PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONSTIMEREVERSALSYMMETRICGENERICHAMILTONIAN_H
#define PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONSTIMEREVERSALSYMMETRICGENERICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Polynomial;


class ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian : public ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian
  {
  public:
      
  // default constructor    
  ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian();
      
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
  ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
								   double ratio, 
								   int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
								   int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
								   int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
								   double spinFluxUp, double spinFluxDown, 
								   AbstractArchitecture* architecture, long memory, char* precalculationFileName, 
								   double* oneBodyPotentielUpUp = 0, double* oneBodyPotentielDownDown = 0, double* oneBodyPotentielUpDown = 0);
  
  // destructor
  //
  ~ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian();

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
  
  protected:
  
  // get all the indices that should appear in the annihilation/creation operators
  //
  virtual void GetIndices();
  

 private:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotentials = pseudopotential coefficients
  // spinFluxM1 = additional inserted flux for m1
  // spinFluxM2 = additional inserted flux for m2
  // spinFluxM3 = additional inserted flux for m3
  // spinFluxM4 = additional inserted flux for m4
  // return value = numerical coefficient
  double EvaluateInteractionCoefficientUpUp(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
					double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4);
  
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotentials = pseudopotential coefficients
  // spinFluxM1 = additional inserted flux for m1
  // spinFluxM2 = additional inserted flux for m2
  // spinFluxM3 = additional inserted flux for m3
  // spinFluxM4 = additional inserted flux for m4
  // return value = numerical coefficient
  double EvaluateInteractionCoefficientDownDown(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
					double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4);
  
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotentials = pseudopotential coefficients
  // spinFluxM1 = additional inserted flux for m1
  // spinFluxM2 = additional inserted flux for m2
  // spinFluxM3 = additional inserted flux for m3
  // spinFluxM4 = additional inserted flux for m4
  // return value = numerical coefficient
  double EvaluateInteractionCoefficientUpDown(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
					double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4); 

};

#endif
