////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//      class of checkerboard lattice model with interacting particles        //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 08/09/2011                      //
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


#ifndef PARTICLEONLATTICEHALFCONTINUOUSHOFSTADTERMODEL_H
#define PARTICLEONLATTICEHALFCONTINUOUSHOFSTADTERMODEL_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Tools/FTITightBinding/TightBindingModelHalfContinuousHofstadterModel.h"
#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeHalfContinuousHofstadterModelSingleBandHamiltonian : public  ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{
 protected:

  // Number of vectors of the reciprocal space in each direction to keep. It must be even!
  int NbrReciprocalVectors;
  
  // nearest neighbor density-density potential strength
  double UPotential;

  // index of the band to be filled
  int BandIndex;
  
  // use flat band model
  bool FlatBand;

  TightBindingModelHalfContinuousHofstadterModel * LocalTightBindingModel;  
  
public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrReciprocalVectors = number of vectors of the reciprocal space in each direction to keep. It must be even!
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeHalfContinuousHofstadterModelSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, int  nbrReciprocalVectors, double uPotential, Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, int bandIndex, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeHalfContinuousHofstadterModelSingleBandHamiltonian();
  
 protected:
  
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();  
};

#endif

