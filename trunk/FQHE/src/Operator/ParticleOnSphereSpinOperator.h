////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Thiebaut                  //
//                                                                            //
//                                                                            //
//                 class Sx operator for particle with spin                   //
//                                                                            //
//                        last modification : 29/03/2013                      //
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


#ifndef PARTICLEONSPHERESPINOPERATOR_H
#define PARTICLEONSPHERESPINOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

class ParticleOnSphereSpinOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  ParticleOnSphereWithSpin* Particle;

  // Maximum individual momentum
  int LMax;
   
  // Coordinate index : 0 for x, 1 for y and 2 for z
  int CoordinateIndex;
    
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // coordinate = Coordinate index : 0 for x, 1 for y and 2 for z
  // lzMax = maximum single particle momentum
  ParticleOnSphereSpinOperator(ParticleOnSphereWithSpin* particle, int coordinate, int lzMax);

  // destructor
  //
  ~ParticleOnSphereSpinOperator();
  
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
 
  // evaluate the SxSquare matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElementSquare (RealVector& V1, RealVector& V2);

  // evaluate part of the SxSquare matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElementSquare (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent);

  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent);

  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
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
