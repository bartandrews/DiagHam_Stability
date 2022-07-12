////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of particle on sphere density operator               //
//                                                                            //
//                        last modification : 11/12/2002                      //
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


#ifndef PARTICLEONSPHEREDENSITYOPERATOR_H
#define PARTICLEONSPHEREDENSITYOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnSphere.h"

#include <iostream>


class MathematicaOutput;


class ParticleOnSphereDensityOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  ParticleOnSphere* Particle;
  
  // index of the creation operator that is part of the density operator (i.e. a+_i)
  int OperatorIndexDagger;
  // index of the annihilation operator that is part of the density operator (i.e. a_i)
  int OperatorIndex;
  
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // index = index of the density operator
  ParticleOnSphereDensityOperator(ParticleOnSphere* particle, int index);

  // constructor when dealing with two different Hilbert spaces
  //
  // particle = hilbert space associated to the right hand state (target space has to be fixed to the hilbert space associated to the left hand state)
  // indexDagger = index of the creation operator that is part of the density operator
  // index = index of the annihilation operator that is part of the density operator
  ParticleOnSphereDensityOperator(ParticleOnSphere* particle, int indexDagger, int index);

  // destructor
  //
  ~ParticleOnSphereDensityOperator();
  
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
