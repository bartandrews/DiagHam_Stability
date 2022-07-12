////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of particle on lattice 1-body operator              //
//                                                                            //
//                        last modification : 09/04/2008                      //
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

#ifndef PARTICLEONLATTICEMOMENTUMOPERATOR_H
#define PARTICLEONLATTICEMOMENTUMOPERATOR_H

#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnLattice.h"

class ParticleOnLatticeMomentumOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  ParticleOnLattice* Particle;
  
  // X-component of the momentum
  int MomentumX;
  // Y-component of the momentum
  int MomentumY;

  // Lattice dimensions
  int Lx;
  int Ly;
  int NbrSubLattices;
  // number of single-particle terms
  int NbrTerms;

  // creation and annihilation indices of one-body terms
  int *CreationIndices;
  int *AnnihilationIndices;
  
  // table to store the phases associated with creation operators exp[i K.(Ri-Rf)]
  double *PhaseTableRe;
  double *PhaseTableIm;
  
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // lx = length in x-direction
  // ly = length in y-direction
  // subl = number of sublattices
  // momentumX = X-component of the momentum
  // momentumY = Y-component of the momentum
  // offsetX = absolute offset of momentum values along x-axis
  // offsetY = absolute offset of momentum values along y-axis
  ParticleOnLatticeMomentumOperator(ParticleOnLattice* particle, int lx, int ly, int subl = 1, int momentumX=0, int momentumY=0, double offsetX=0.0, double offsetY=0.0);

  // copy constructor
  //
  // oper = reference on the operator to copy
  ParticleOnLatticeMomentumOperator(const ParticleOnLatticeMomentumOperator& oper);

  // destructor
  //
  ~ParticleOnLatticeMomentumOperator();
  
  // clone operator without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractOperator* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which operator acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where operator acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();

  // change values of momentum represented
  // momentumX = X-component of the momentum
  // momentumY = Y-component of the momentum
  // offsetX = absolute offset of momentum values along x-axis
  // offsetY = absolute offset of momentum values along y-axis
  void SetMomentum (int momentumX, int momentumY, double offsetX=0.0, double offsetY=0.0);

  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent);

  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent);
   
  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
			       int firstComponent, int nbrComponent);
  
};

#endif 
