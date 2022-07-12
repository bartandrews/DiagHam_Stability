////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                  coulomb interaction and magnetic translations             //
//                                                                            //
//                        last modification : 02/10/2003                      //
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


#ifndef PARTICLEONTORUSCOULOMBWITHSPINANDMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTORUSCOULOMBWITHSPINANDMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Polynomial;


class ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian : public AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian
{

 protected:

  // parameter of layer separation if non-SU(2) invariant interaction
  double LayerSeparation;

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
  // 
  ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian();
     
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
  ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
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
  ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
								   double ratio, 
								   int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
								   int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
								   int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
								   double spinFluxUp, double spinFluxDown, 
								   AbstractArchitecture* architecture, long memory, char* precalculationFileName, 
								   double* oneBodyPotentielUpUp, double* oneBodyPotentielDownDown, double* oneBodyPotentielUpDown);
  
  // destructor
  //
  ~ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian();

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
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, double layerSeparation=0.0);

  // get fourier transform of interaction
  // Q2_half = one half of q² value
  // layerSeparation = layer separation
  virtual double GetVofQ(double Q2_half, double layerSeparation=0.0);

  
  // evaluate Wigner crystal energy per particle
  //
  // return value = Wigner crystal energy per particle
  virtual double EvaluateWignerCrystalEnergy ();

  // evaluate Misra function (integral of t^n exp (-xt) between 1 and +inf)
  //
  // n = index of the Misra function
  // x = point where the function has to be evaluated (> 0)
  // return value = value of the n-Misra function at x
  virtual double MisraFunction (double n, double x);

  // evaluate part of the integral needed in the Misra function (integral of t^n exp (-xt) between min and max)
  //
  // n = index of the Misra function
  // x = point where the function has to be evaluated (> 0)
  // min = lower bound of the integral
  // max = upper bound of the integral
  // nbrSubdivision = number of subdivision used for the integral
  // return value = value of the integral
  virtual double PartialMisraFunction (double n, double x, double min, double max, int nbrSubdivision);

};

#endif // PARTICLEONTORUSCOULOMBWITHSPINANDMAGNETICTRANSLATIONSHAMILTONIAN_H
