////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//               class of Haldane model with interacting particles            //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 16/08/2011                      //
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


#ifndef PARTICLEONTWISTEDTORUSGENERICHAMILTONIAN_H
#define PARTICLEONTWISTEDTORUSGENERICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Polynomial/Polynomial.h"


#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnTwistedTorusGenericHamiltonian : public ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian
{

 protected:
  // Maximum angular momentum
  int MaxMomentum;
  // angle to use for twisted torus pseudopotentials
  double Angle;
  // aspect ratio of Torus
  double Ratio;
  // energy of wigner crystal reference
  double WignerEnergy;
  // cosine of the torus angle
  double CosTheta;
  // sine of the torus angle
  double SinTheta;
  // length of torus along x-direction
  double Lx;
  // length of torus along y-direction
  double Ly;
  // reciprocal lattice vector x-direction
  double Gx;
  // reciprocal lattice vector y-direction
  double Gy;
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
  ParticleOnTwistedTorusGenericHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the lengths of the two fundamental cycles of the torus L1 / L2
  // angle = angle between the two fundamental cycles of the torus, along (L1 sin, L1 cos) and (0, L2)
  // haveCoulomb = flag indicating whether a coulomb term is present
  // landauLevel = landauLevel to be simulated
  // nbrPseudopotentials = number of pseudopotentials indicated
  // pseudopotentials = pseudopotential coefficients
  // noWignerEnergy = do not consider the energy contribution from the Wigner crystal 
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnTwistedTorusGenericHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum,
										   
					   double ratio, double angle, bool haveCoulomb, int landauLevel, int nbrPseudopotentials, double* pseudopotentials, bool noWignerEnergy, AbstractArchitecture* architecture, long memory);

  // destructor
  //
  ~ParticleOnTwistedTorusGenericHamiltonian();
  

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
  //double RectangularEvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  //Complex TwistedEvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  Complex EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);
  
  // get fourier transform of interaction
  // Q2_half = one half of q² value
  double GetVofQ(double Q2_half);

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
