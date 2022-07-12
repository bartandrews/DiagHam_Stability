////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//          class of a two body interaction projected onto two bands          //
//      from an ASCII file providing the two body matrix real elements        //
//                                                                            //
//                        last modification : 04/05/2020                      //
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


#ifndef PARTICLEONLATTICEFROMFILEINTERACTIONTWOBANDREALHAMILTONIAN_H
#define PARTICLEONLATTICEFROMFILEINTERACTIONTWOBANDREALHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullTwoBandRealHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian : public ParticleOnLatticeQuantumSpinHallFullTwoBandRealHamiltonian
{

 protected:
 
  // numerical factor for momentum along x
  double KxFactor;
  // numerical factor for momentum along y
  double KyFactor;
  
  // index of the filled first band
  int BandIndex1;
  // index of the filled second band
  int BandIndex2;
  
  // use flat band model
  bool FlatBand;
  
  // global rescaling factor for the two-body interaction term
  double InteractionRescalingFactor;

  // include an additional spin 1/2 degree of freedom, building an SU(2) invariant interaction
  bool AdditionalSpinFlag;
  
  // name of the ASCII file containing the matrix element for the generic two body interaction term
  char* MatrixElementsInteractionFile;
  
 public:

  // default constructor
  //
  ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // matrixElementsInteractionFile = name of the ASCII file containing the matrix element for the generic two body interaction term
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // interactionRescalingFactor = global rescaling factor for the two-body interaction term
  // spinFlag = include an additional spin 1/2 degree of freedom, building an SU(2) invariant interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,
							     char* matrixElementsInteractionFile,
							     Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag,
							     double interactionRescalingFactor, 
							     bool spinFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
  
  // evaluate all one-body factors
  //     
  virtual void EvaluateOneBodyFactors();

  // test which internal degrees of freedom are conserved in the matrix elements
  //   
  // nbrMatrixElements = number of matrix elements
  // sigmaIndices1 = array for internal degrees of freedom of the first creation operator
  // sigmaIndices2 = array for internal degrees of freedom of the second creation operator
  // sigmaIndices3 = array for internal degrees of freedom of the first annihilation operator
  // sigmaIndices4 = array for internal degrees of freedom of the second annihilation operator
  // return value = array that indicates which internal degrees of freedom are conserved  
  virtual bool**** TestMatrixElementsConservedDegreesOfFreedom (int nbrMatrixElements, int* sigmaIndices1, int* sigmaIndices2,
								int* sigmaIndices3, int* sigmaIndices4);
  
  // free the array tagging which internal degrees of freedom are conserved in the matrix elements
  //   
  // internalIndicesFlags = array that indicates which internal degrees of freedom are conserved
  void FreeMatrixElementsConservedDegreesOfFreedom (bool**** internalIndicesFlags);

};


#endif
