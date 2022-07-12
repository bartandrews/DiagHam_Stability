////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//       generic two body interaction without assuming ky conservation        //
//                                                                            //
//                        last modification : 11/09/2016                      //
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


#ifndef PARTICLEONTORUSGENERICHAMILTONIANALLMOMENTA_H
#define PARTICLEONTORUSGENERICHAMILTONIANALLMOMENTA_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnTorusHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Polynomial;


class ParticleOnTorusGenericHamiltonianAllMomenta : public AbstractQHEOnTorusHamiltonian
{

 protected:

  // Number of Pseudopotential
  int NbrPseudopotentials;
  // pseudopotential coefficients
  double* Pseudopotentials;
  // Laguerre polynomial for the pseudopotentials
  Polynomial* LaguerrePolynomials;

  // additional inserted flux along the x axis in flux quantum unit
  double XFlux;
  // additional inserted flux along the y axis in flux quantum unit
  double YFlux;

 public:

  // default constructor
  //
  ParticleOnTorusGenericHamiltonianAllMomenta();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotential = pseudopotential coefficients
  // oneBodyPotentials = array that contains the diagonal terms of the one-body potential
  // xFlux = additional inserted flux along the x axis in flux quantum unit
  // yFlux = additional inserted flux along the y axis in flux quantum unit
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusGenericHamiltonianAllMomenta(ParticleOnTorus* particles, int nbrParticles, int maxMomentum, double ratio, 
					      int nbrPseudopotentials, double* pseudopotentials, double* oneBodyPotentials,
					      double xFlux, double yFlux, 
					      AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);
    
  // destructor
  //
  ~ParticleOnTorusGenericHamiltonianAllMomenta();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

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
  // fluxM1 = additional inserted flux for m1
  // fluxM2 = additional inserted flux for m2
  // fluxM3 = additional inserted flux for m3
  // fluxM4 = additional inserted flux for m4
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4,
						double fluxM1, double fluxM2, double fluxM3, double fluxM4);

};

#endif
