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


#ifndef PARTICLEONCYLINDERORBITALPROJECTION_H
#define PARTICLEONCYLINDERORBITALPROJECTION_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnCylinderHamiltonian.h"
#include "Polynomial/Polynomial.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnCylinderOrbitalProjection : public AbstractQHEOnCylinderHamiltonian
{

 protected:

  //orbital to be projected out
  int OrbitalIndex;

  //shape (anisotropy) parameter
  double Anisotropy;

  //position in real space
  double X0;
  double Y0;

  // array that contains all one-body interaction factors
  Complex* OneBodyInteractionFactors;

  //array the contains all one-body interaction indices
  int* OneBodyM1Values;
  int* OneBodyM2Values;

  //number of one-body terms
  int NbrOneBodyInteractionFactors;

  // Laguerre polynomial
  Polynomial *LaguerreM;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // orbitalIndex = index of the orbital to be projected out 
  // anisotropy = shape (anisotropy) parameter
  // x0, y0 = position in real space
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnCylinderOrbitalProjection(ParticleOnSphere* particles, int nbrParticles, int maxMomentum, double ratio, int orbitalIndex, double anisotropy, double x0, double y0,
					   AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnCylinderOrbitalProjection();

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

  double Integrand(double qx, void *p);

  double OrbitalProjectionMatrixElement(double xj1, double xj2, double x0, double y0, int orbitalIndex, double anisotropy, Polynomial* laguerreM, double kappa, int maxMomentum, double &error);

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

 protected:

  // evaluate the numerical coefficient  in front of the a+_m1 a_m2 coupling term
  //
  // m1 = first index
  // m2 = second index
  // return value = numerical coefficient
  Complex EvaluateInteractionCoefficient(int m1, int m2);

  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

};
#endif
