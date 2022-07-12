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
//                          at momentum points 0 or pi                        //
//                                                                            //
//                        last modification : 02/11/2019                      //
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


#ifndef PARTICLEONTORUSCOULOMBWITHMAGNETICTRANSLATIONSREALHAMILTONIAN_H
#define PARTICLEONTORUSCOULOMBWITHMAGNETICTRANSLATIONSREALHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnTorusWithMagneticTranslationsRealHamiltonian.h"
#include "Polynomial/Polynomial.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnTorusCoulombWithMagneticTranslationsRealHamiltonian : public AbstractQHEOnTorusWithMagneticTranslationsRealHamiltonian
{

 protected:

  // energy of wigner crystal reference
  double WignerEnergy;

  // landau Level index
  int LandauLevel;

  // Number of Pseudopotential
  int NbrPseudopotentials;

  // pseudopotential coefficients
  double *Pseudopotentials;

  // flag indicating whether Coulomb part is present
  bool HaveCoulomb;

  // form factor of the interaction (a single Laguerre polynomial for the Landau levels of GaAs)
  Polynomial FormFactor;
  
  // Laguerre polynomial for the pseudopotentials
  Polynomial *LaguerreM;

 public:

  // default constructor
  //
  ParticleOnTorusCoulombWithMagneticTranslationsRealHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // ratio = ratio between the width in the x direction and the width in the y direction
  // haveCoulomb = flag indicating whether a coulomb term is present
  // landauLevel = landauLevel to be simulated
  // nbrPseudopotentials = number of pseudopotentials indicated
  // pseudopotentials = pseudopotential coefficients
  // noWignerEnergy = do not consider the energy contribution from the Wigner crystal 
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusCoulombWithMagneticTranslationsRealHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum,
							    double ratio, bool haveCoulomb, int landauLevel, int nbrPseudopotentials, double* pseudopotentials, bool noWignerEnergy, AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusCoulombWithMagneticTranslationsRealHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // shift RealHamiltonian from a given energy
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
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

  // get fourier transform of interaction
  // Q2_half = one half of q² value
  virtual double GetVofQ(double Q2_half);

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

#endif
