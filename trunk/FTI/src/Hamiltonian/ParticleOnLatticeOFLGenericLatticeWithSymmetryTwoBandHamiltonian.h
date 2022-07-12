////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                          class author: Gunnar Möller                       //
//                                                                            //
//      class of Hamiltonian of bosons (currently) in a generic optical       //
//                            flux lattice                                    //
//                                                                            //
//                        last modification : 01/07/2014                      //
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


#ifndef PARTICLEONLATTICEOFLGENERICLATTICEWITHSYMMETRYTWOBANDHAMILTONIAN_H
#define PARTICLEONLATTICEOFLGENERICLATTICEWITHSYMMETRYTWOBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnSquareLatticeTwoBandSimpleTIHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Tools/FTITightBinding/TightBindingModelOFLGenericLatticeWithSymmetry.h"
#include "Matrix/RealSymmetricMatrix.h"
#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian : public 
ParticleOnSquareLatticeTwoBandSimpleTIHamiltonian
{
 protected:
  
  // number of reciprocal lattice vectors along k1 / k2 direction in tight-binding model:
  int NbrReciprocalVectors1;
  int NbrReciprocalVectors2;

  // SymmetryMultipliers: enlargement of the unit cell along LatticeVectors 1,2
  int SymmetryMultiplier1;
  int SymmetryMultiplier2;

  // nearest neighbor density-density potential strength
  RealSymmetricMatrix UPotentialMatrix;
  // 
  double UPotential;
 
  // use flat band model
  bool FlatBand;
  bool NoDispersionFlag;

  TightBindingModelOFLGenericLatticeWithSymmetry * LocalTightBindingModel;  
  
public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // tightBindingModel = Generic Tight Binding Model With Symmetry
  // uPotential = spin-independent strength of the repulsive two body contact interactions
  // flatBandFlag = use flat band model
  // noDispersionFlag = remove the band dispersion and put a constant gap of 10 between the two lowest bands
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, TightBindingModelOFLGenericLatticeWithSymmetry* tightBindingModel, double uPotential, bool flatBandFlag, bool noDispersionFlag, AbstractArchitecture* architecture, long memory = -1);

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // tightBindingModel = Generic Tight Binding Model With Symmetry
  // uPotentialMatrix = spin-dependent strength of the repulsive two body contact interactions
  // flatBandFlag = use flat band model
  // noDispersionFlag = remove the band dispersion and put a constant gap of 10 between the two lowest bands
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, TightBindingModelOFLGenericLatticeWithSymmetry* tightBindingModel, RealSymmetricMatrix &uPotentialMatrix, bool flatBandFlag, bool noDispersionFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian();
  
 protected:
  
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();  

};

#endif

