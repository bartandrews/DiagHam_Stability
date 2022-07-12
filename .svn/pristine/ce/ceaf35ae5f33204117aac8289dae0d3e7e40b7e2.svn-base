////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of particle on sphere square of total spin operator           //
//                                                                            //
//                        last modification : 01/01/2007                      //
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


#ifndef PARTICLEONSPHERESQUARETOTALSPINOPERATOR_H
#define PARTICLEONSPHERESQUARETOTALSPINOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"
#include "Matrix/RealMatrix.h"


class ParticleOnSphereSquareTotalSpinOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  ParticleOnSphereWithSU4Spin* Particle;

  //  twice the total spin projection (i.e Sz) value
  int TotalSz;
  // maximum Lz value reached by a particle
  int LzMax;
  //  number of particles
  int NbrParticles;

  // shift due to the Sz^2 and anticommutator terms
  double Shift;

 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // lzMax = maximum Lz value reached by a fermion
  // nbrParticles = number of particles
  // totalSz = twice the total spin projection (i.e Sz) value
  ParticleOnSphereSquareTotalSpinOperator(ParticleOnSphereWithSU4Spin* particle, int lzMax, int nbrParticles, int totalSz);

  // copy constructor
  //
  // oper = reference on the operator to copy
  ParticleOnSphereSquareTotalSpinOperator(const ParticleOnSphereSquareTotalSpinOperator& oper);

  // destructor
  //
  ~ParticleOnSphereSquareTotalSpinOperator();
  
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
  
  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent);

  // multiply a vector by the current operator for a given range of indices 
  // and add result to another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
				  int firstComponent, int nbrComponent);
  
};

#endif
