////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of particle on sphere lowering momentum L operator           //
//                                                                            //
//                        last modification : 06/06/2006                      //
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


#ifndef PARTICLEONSPHERELMINUSOPERATOR_H
#define PARTICLEONSPHERELMINUSOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Vector/RealVector.h"


class ParticleOnSphereLMinusOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles with totalLz  Lz momentum, target space has to be fixed to hilbert space totalLz - 1 Lz momentum
  ParticleOnSphere* Particle;

  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a particle
  int LzMax;

  // vector where all coefficents that come from the L- operator are stored
  RealVector Coefficients;

 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles with totalLz  Lz momentum, target space has to be fixed to hilbert space totalLz - 1 Lz momentum
  // totalLz = momentum total value of the source hilbert space associated
  // lzMax = maximum Lz value reached by a fermion
  ParticleOnSphereLMinusOperator(ParticleOnSphere* particle, int totalLz, int lzMax);

  // copy constructor
  //
  // oper = reference on the operator to copy
  ParticleOnSphereLMinusOperator(const ParticleOnSphereLMinusOperator& oper);

  // destructor
  //
  ~ParticleOnSphereLMinusOperator();
  
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
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);
   
  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
			       int firstComponent, int nbrComponent);
  
};

#endif
