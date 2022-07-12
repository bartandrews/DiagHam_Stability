////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//      two body pseudopotential interaction and magnetic translations        //
//                                                                            //
//                        last modification : 07/03/2014                      //
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


#ifndef PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONSGENERICHAMILTONIAN_H
#define PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONSGENERICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Polynomial;


class ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian : public AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian
{

 protected:

  // Number of Pseudopotential for up-up interaction
  int NbrPseudopotentialsUpUp;
  // pseudopotential coefficients for up-up interaction
  double* PseudopotentialsUpUp;
  // Number of Pseudopotential for down-down interaction
  int NbrPseudopotentialsDownDown;
  // pseudopotential coefficients for down-down interaction
  double* PseudopotentialsDownDown;
  // Number of Pseudopotential for up-down interaction
  int NbrPseudopotentialsUpDown;
  // pseudopotential coefficients for up-down interaction
  double* PseudopotentialsUpDown;
  // maximum number of pseudopotentials
  int MaxNbrPseudopotentials;
  // Laguerre polynomial for the pseudopotentials
  Polynomial* LaguerrePolynomials;

  // additional inserted flux for spin up
  double SpinFluxUp;
  // additional inserted flux for spin down
  double SpinFluxDown;


 public:
   
   // default constructor
  // 
  ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian();


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
  ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
								   double ratio, 
								   int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
								   int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
								   int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
								   double spinFluxUp, double spinFluxDown, 
								   AbstractArchitecture* architecture, long memory, char* precalculationFileName, 
								   double* oneBodyPotentielUpUp = 0, double* oneBodyPotentielDownDown = 0, double* oneBodyPotentielUpDown = 0);
  
  // destructor
  //
  ~ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  virtual void ShiftHamiltonian (double shift);

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

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
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
						double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4);


};

#endif
