////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                       coulombian interaction  and spin                     //
//                                                                            //
//                        last modification : 16/09/2002                      //
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


#ifndef PARTICLEONTORUSCOULOMBWITHSPINHAMILTONIAN_H
#define PARTICLEONTORUSCOULOMBWITHSPINHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpin.h"
#include "Hamiltonian/AbstractHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnTorusCoulombWithSpinHamiltonian : public AbstractHamiltonian
{

 protected:
  
  // Hilbert space associated to the system
  ParticleOnTorusWithSpin* Particles;

  // number of particles
  int NbrParticles;

  // energy of the Wigner crystal
  double WignerEnergy;

  // magnetic coupling constant times magnetic field 
  double MagneticG;

  // maximum momentum value reached by a particle in the state
  int MaxMomentum;
  // number of Lz values in a state
  int NbrLzValue;

  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;
  // layer separation
  double LayerSeparation;

  // array containing all interaction factors 
  double* InteractionFactors;
  // array containing all interaction factors 
  double* InterLayerInteractionFactors;  
  // number of interaction factors
  int NbrInteractionFactors;
  // arrays for indices attached to each interaction factor
  int* M1Value;
  int* M2Value;
  int* M3Value;
  int* M4Value;

  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // number of non-null term in the hamiltonian for each state
  int* NbrInteractionPerComponent;
  // index of the state obtained for each term of the hamiltonian when applying on a given state
  int** InteractionPerComponentIndex;
  // multiplicative coefficient obtained for each term of the hamiltonian when applying on a given state and with a given destination state
  double** InteractionPerComponentCoefficient;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // magneticG = magnetic coupling constant times magnetic field 
  ParticleOnTorusCoulombWithSpinHamiltonian(ParticleOnTorusWithSpin* particles, int nbrParticles, int maxMomentum, double ratio, double magneticG = 0.0, double layerSeparation=0.0, long memory=0);

  // destructor
  //
  ~ParticleOnTorusCoulombWithSpinHamiltonian();

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

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  bool SavePrecalculation (char* fileName);

  // set the value of the magnetic coupling constant times magnetic field
  //
  // magneticG = magnetic coupling constant times magnetic field 
  void SetMagneticG (double magneticG);

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

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
				  int firstComponent, int nbrComponent);

 
  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  List<Matrix*> RightInteractionOperators();  

  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, ParticleOnTorusCoulombWithSpinHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnTorusCoulombWithSpinHamiltonian& H);

 private:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // get fourier transform of interaction
  // Q2_half = one half of q² value
  // layerSeparation = layer separation
  double GetVofQ(double Q2_half, double layerSeparation=0.0);

  // test the amount of memory needed for fast multiplication algorithm
  //
  // return value = amount of memory needed
  long FastMultiplicationMemory();

  // enable fast multiplication algorithm
  //
  void EnableFastMultiplication();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, double layerSeparation = 0.0);

  // evaluate Wigner crystal energy per particle
  //
  // return value = Wigner crystal energy per particle
  double EvaluateWignerCrystalEnergy ();

  // evaluate Misra function (integral of t^n exp (-xt) between 1 and +inf)
  //
  // n = index of the Misra function
  // x = point where the function has to be evaluated (> 0)
  // return value = value of the n-Misra function at x
  double MisraFunction (double n, double x);

  // evaluate part of the integral needed in the Misra function (integral of t^n exp (-xt) between min and max)
  //
  // n = index of the Misra function
  // x = point where the function has to be evaluated (> 0)
  // min = lower bound of the integral
  // max = upper bound of the integral
  // nbrSubdivision = number of subdivision used for the integral
  // return value = value of the integral
  double PartialMisraFunction (double n, double x, double min, double max, int nbrSubdivision);

};

#endif
