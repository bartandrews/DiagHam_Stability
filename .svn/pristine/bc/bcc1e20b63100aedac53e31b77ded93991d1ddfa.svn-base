////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//          class of a generic density-density two body interaction           //
//                     projected onto a single band and                       //
//              assuming a Bloch form for the tight binding model             //
//                                                                            //
//                        last modification : 25/09/2014                      //
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


#ifndef PARTICLEONLATTICEGENERICTWOBODYINTERACTIONSINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICEGENERICTWOBODYINTERACTIONSINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

 protected:
  
  
  // number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  int* NbrInteractingOrbitals;
  // orbital indices of the orbitals interacting with each orbital within the unit cell 
  // at the origin through a density-density term  
  int** InteractingOrbitalsOrbitalIndices;  
  // spatial indices (sorted as 2 consecutive integers) of the orbitals interacting 
  // with each orbital within the unit cell at the origin through a density-density term
  int** InteractingOrbitalsSpatialIndices;
  // intensity of each density-density term 
  double** InteractingOrbitalsPotentials; 
  
  // index of the band to be filled
  int BandIndex;
  
  // use flat band model
  bool FlatBand;
  
 public:

  // constructor
  //
  ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // bandIndex = index of the band in which the Hamiltonian is projected 
  // nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  // interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  // interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  // interactingOrbitalsPotentials = intensity of each density-density term 
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, 
									 int nbrSiteX, int nbrSiteY, int bandIndex, 
									 int* nbrInteractingOrbitals, int** interactingOrbitalsOrbitalIndices,
									 int** interactingOrbitalsSpatialIndices, double** interactingOrbitalsPotentials,
									 Abstract2DTightBindingModel* tightBindingModel,
									 bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

};



#endif
