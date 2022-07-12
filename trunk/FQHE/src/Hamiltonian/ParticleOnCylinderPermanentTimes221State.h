////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          laplacian delta interaction                       //
//                                                                            //
//                        last modification : 29/06/2010                      //
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


#ifndef PARTICLEONCYLINDERPERMANENTTIMES221STATE_H
#define PARTICLEONCYLINDERPERMANENTTIMES221STATE_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnCylinderThreeBodyHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnCylinderPermanentTimes221State : public AbstractQHEOnCylinderThreeBodyHamiltonian
{

 protected:

  //consider spinful Gaffnian (instead of permanent times 221) or not
  bool GaffnianFlag;
  //consider NASS state at nu=4/3 or not
  bool NASSFlag;

  // array containing all interaction factors 
  Complex* InteractionFactorsIntra;
  // number of interaction factors
  int NbrInteractionFactorsIntra;
  // arrays for indices attached to each interaction factor
  int* M1ValueIntra;
  int* M2ValueIntra;
  int* M3ValueIntra;
  int* M4ValueIntra;

  // array containing all interaction factors 
  Complex* InteractionFactorsInter;
  // number of interaction factors
  int NbrInteractionFactorsInter;
  // arrays for indices attached to each interaction factor
  int* M1ValueInter;
  int* M2ValueInter;
  int* M3ValueInter;
  int* M4ValueInter;

  // array containing all interaction factors 
  Complex* InteractionFactors12;
  // number of interaction factors
  int NbrInteractionFactors12;
  // arrays for indices attached to each interaction factor
  int* M1Value12;
  int* M2Value12;
  int* M3Value12;
  int* M4Value12;
  int* M5Value12;
  int* M6Value12;

  // array containing all interaction factors 
  Complex* InteractionFactors32;
  // number of interaction factors
  int NbrInteractionFactors32;
  // arrays for indices attached to each interaction factor
  int* M1Value32;
  int* M2Value32;
  int* M3Value32;
  int* M4Value32;
  int* M5Value32;
  int* M6Value32;

  //one body terms c_{mu}^+ c_{mu}, c_{md}^+ c_{md}, and c_{Xm}^+ c_{Xm}
  Complex* OneBodyUpUp;
  Complex* OneBodyUpDown;
  Complex* OneBodyDownDown;


 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // gaffnianFlag = consider spinful Gaffnian instead of permanent times 221
  // nassFlag = consider NASS state at nu=4/3
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnCylinderPermanentTimes221State(ParticleOnSphereWithSpin* particles, int nbrParticles, int maxMomentum, int totalSpin, double ratio, bool gaffnianFlag, bool nassFlag, AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnCylinderPermanentTimes221State();

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

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored

  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored

  ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								  int firstComponent, int nbrComponent);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element

  long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted

  void PartialEnableFastMultiplication(int firstComponent, int lastComponent);

 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a+_m3 a_m4 a_m5 a_m6 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // m5 = fifth index
  // m6 = sixth index
  // return value = numerical coefficient
  Complex EvaluateInteractionCoefficient12(int m1, int m2, int m3, int m4, int m5, int m6, bool gaffnianFlag, bool nassFlag);
  Complex EvaluateInteractionCoefficient32(int m1, int m2, int m3, int m4, int m5, int m6, bool gaffnianFlag, bool nassFlag);
  Complex EvaluateInteractionCoefficientIntra(int m1, int m2, int m3, int m4, bool gaffnianFlag, bool nassFlag);
  Complex EvaluateInteractionCoefficientInter(int m1, int m2, int m3, int m4, bool gaffnianFlag, bool nassFlag);

  int NumberOfPermutations32(int n1, int n2, int n3);
  int NumberOfPermutations12(int n1, int n2);
};

#endif
