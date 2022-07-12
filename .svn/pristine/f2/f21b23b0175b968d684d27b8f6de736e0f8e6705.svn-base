////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//   class of bosons on the CP2 sphere with generic two-body interaction      // 
//                                                                            //
//                        last modification : 14/02/2013                      //
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


#ifndef PARTICLEONCP2GENERICTWOBODYHAMILTONIAN_H
#define PARTICLEONCP2GENERICTWOBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereFullHamiltonian.h"
#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnCP2.h"
#include "HilbertSpace/FermionOnCP2.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnCP2GenericTwoBodyHamiltonian : public AbstractQHEOnSphereFullHamiltonian
{
  
 protected:
    //Number of flux quanta
  int NbrFluxQuanta;
  //Number of orbitals for a state with a given number of flux quanta
  int NbrLzValue;   
  // array containing the SU(3) pseudopotentials describing the interaction
  double* PseudoPotential;
  // index of the highest index non zero pseudopotential
  int PseudoPotentialIndexMax;
  // array that gives the value of tz for one particle corresponding to the linearized index
  int* quantumNumberTz;
  // array that gives the value of y for one particle corresponding to the linearized index
  int* quantumNumberY;
  // array that gives the value of r for one particle corresponding to the linearized index
  int* quantumNumberR;
  // array that gives the value of s for one particle corresponding to the linearized index
  int* quantumNumberS;
  //Hilbert Space CP2
  BosonOnCP2* ParticlesBosons;
  FermionOnCP2* ParticlesFermions;
  // array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
  double* OneBodyPotentials;
  int* RepresentationDimension;
  
 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuanta = number of flux quanta
  // pseudoPotential = pointer to an array containing the SU(3) pseudopotentials describing the interaction
  // nbrPseudoPotentials = number of pseudo-potentials
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCP2GenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, 
					 double* pseudoPotential, int nbrPseudoPotentials, 
					 AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnCP2GenericTwoBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
  
  

};



#endif
