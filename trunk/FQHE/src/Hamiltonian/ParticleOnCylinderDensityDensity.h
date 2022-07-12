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


#ifndef PARTICLEONCYLINDERDENSITYDENSITY_H
#define PARTICLEONCYLINDERDENSITYDENSITY_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnCylinderHamiltonian.h"
#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"
#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnCylinderDensityDensity : public AbstractQHEOnCylinderHamiltonian
{

 protected:

 //Landau level index
 int LandauLevel;

 //(x0,y0) coordinate of the origin
 double X0Value;
 double Y0Value;

 //(x,y) coordinate where to evaluate the rho-rho correlation
 double XValue;
 double YValue;

 int HoppingCutoff;

 //one body basis
 ParticleOnCylinderFunctionBasis* Basis;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // landauLevel = LL index
  // x0,y0 = coordinates of the origin
  // x,y = coordinates of the point where to evaluate rho-rho
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnCylinderDensityDensity(ParticleOnSphere* particles, int nbrParticles, int maxMomentum, double ratio, int landauLevel, double x0, double y0, double x, double y, int hoppingCutoff, AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnCylinderDensityDensity();

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
  // return value = numerical coefficient
  Complex EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

};

#endif
