////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//              class of bosons on the CP2 sphere with delta interaction      // 
//                                                                            //
//                        last modification : 10/01/2013                      //
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


#ifndef PARTICLEONCP2DELTAHAMILTONIAN_H
#define PARTICLEONCP2DELTAHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereFullHamiltonian.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/BosonOnCP2.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnCP2DeltaHamiltonian : public AbstractQHEOnSphereFullHamiltonian
{
  
 protected:
   
   //Number of flux quanta
  int NbrFluxQuanta;
  //Number of orbitals for a state with a given number of flux quanta
  int NbrLzValue;
  // array that gives the value of tz for one particle corresponding to the linearized index
  int* quantumNumberTz;
  // array that gives the value of y for one particle corresponding to the linearized index
  int* quantumNumberY;
  // array that gives the value of r for one particle corresponding to the linearized index
  int* quantumNumberR;
  // array that gives the value of s for one particle corresponding to the linearized index
  int* quantumNumberS;
  //Hilbert Space CP2
  BosonOnCP2* Particles2;
   // array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
  double* OneBodyPotentials;
  
 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuanta = number of flux quanta
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCP2DeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, 
						    AbstractArchitecture* architecture, long memory = -1);

  // constructor with one body terms
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuanta = number of flux quanta
  // oneBodyPotentials = array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCP2DeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, double* oneBodyPotentials, 
						    AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnCP2DeltaHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute //the matrix element for the two body delta interaction between two particles 
  
  // r1 = creation r
  // s1 = creation s
  // r2 = annihilation r
  // s2 = annihilation s
  // return value = corresponding matrix element
  double ComputeTwoBodyMatrixElement(int r1, int s1, int r2, int s2, int r3, int s3, int r4, int s4, FactorialCoefficient &Coef);
  
  

};



#endif
