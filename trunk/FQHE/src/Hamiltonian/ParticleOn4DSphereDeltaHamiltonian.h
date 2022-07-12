////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Cecile Repellin                   //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//              class of bosons on the 4d sphere with delta interaction       // 
//                                                                            //
//                        last modification : 24/10/2012                      //
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


#ifndef PARTICLEON4DSPHEREDELTAHAMILTONIAN_H
#define PARTICLEON4DSPHEREDELTAHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereFullHamiltonian.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOn4DSphereDeltaHamiltonian : public AbstractQHEOnSphereFullHamiltonian
{
  
 protected:
   
   //Number of flux quanta
  int NbrFluxQuanta;
  //Number of orbitals for a state with a given number of flux quanta
  int NbrLzValue;
  
  // array that gives the value of j for one particle corresponding to the linearized index
  int* quantumNumberJ;
  // array that gives the value of jz for one particle corresponding to the linearized index
  int* quantumNumberJz;
  // array that gives the value of kz for one particle corresponding to the linearized index
  int* quantumNumberKz;
 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuanta = number of flux quanta
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOn4DSphereDeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, 
						    AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOn4DSphereDeltaHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body delta interaction between two particles 
  //
  // j1 = creation j
  // jz1 = creation jz
  // kz1 = creation kz
  // j2 = annihilation j
  // jz2 = annihilation jz
  // kz2 = annihilation kz
  // return value = corresponding matrix element
  double ComputeTwoBodyMatrixElement(int j1, int jz1, int kz1, int j2, int jz2, int kz2, int j3, int jz3, int kz3, int j4, int jz4, int kz4, FactorialCoefficient &Coef);
  
  // get the quantum numbers j, jz, kz of a one particle state 
  //
  //quantumNumberJ = array that gives the quantum number j for a single particle state
  //quantumNumberJz = array that gives the quantum number j for a single particle stateDescription
  //quantumNumberKz = array that gives the quantum number j for a single particle state
  inline void GetQuantumNumbersFromLinearizedIndex(int* quantumNumberJ, int* quantumNumberJz, int* quantumNumberKz)
  {
    for (int j = 0; j <= this->NbrFluxQuanta; ++j)
      {
	for (int jz = 0; jz <= j; ++jz)
	  {
	    for (int kz = 0; kz <= this->NbrFluxQuanta - j; ++kz)
	      {
		int index = -((j-1)*j*(2*j-1))/6 + this->NbrFluxQuanta*j*(j-1)/2 + (this->NbrFluxQuanta+1)*j+ (this->NbrFluxQuanta + 1 - j)*jz +kz;
		quantumNumberJ[index] = j;
		quantumNumberJz[index] = jz;
		quantumNumberKz[index] = kz;
		
	      }
	  }
      }
  }


};



#endif
