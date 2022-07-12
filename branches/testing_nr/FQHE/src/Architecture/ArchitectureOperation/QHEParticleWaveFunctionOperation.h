////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of QHE particle wave function evaluation operation         //
//                                                                            //
//                        last modification : 29/07/2004                      //
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


#ifndef QHEPARTICLEWAVEFUNCTIONOPERATION_H
#define QHEPARTICLEWAVEFUNCTIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractScalarSumOperation.h"
#include "HilbertSpace/AbstractQHEParticle.h"


class AbstractFunctionBasis;
class RealVector;


class QHEParticleWaveFunctionOperation: public AbstractScalarSumOperation
{

 protected:

  // pointer to the Hilbert space
  AbstractQHEParticle* HilbertSpace;

  // vector corresponding to the state in the Fock basis  
  RealVector* State;
  // array of vectors corresponding to the states in the Fock basis
  RealVector* States;

  // vector whose components give coordinates of the point where the wave function has to be evaluated
  RealVector* Position;
  // one body real space basis to use
  AbstractFunctionBasis* Basis;
  // indicate which coordinates will be change during next time step (-1 if no time coherence has to be used)
  int NextCoordinates;

 public:
  
  // constructor 
  //
  // space = pointer to the Hilbert space to use
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // nextCoordinates = indicate which coordinates will be change during next time step (-1 if no time coherence has to be used)
  QHEParticleWaveFunctionOperation(AbstractQHEParticle* space, RealVector* state, RealVector* position, 
				   AbstractFunctionBasis* basis, int nextCoordinates = -1);

  // constructor for multiple wave function evaluations
  //
  // space = pointer to the Hilbert space to use
  // states = array of vectors corresponding to the states in the Fock basis
  // nbrStates = number of states in the states array
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // nextCoordinates = indicate which coordinates will be change during next time step (-1 if no time coherence has to be used)
  QHEParticleWaveFunctionOperation(AbstractQHEParticle* space, RealVector* states, int nbrStates, RealVector* position, 
				   AbstractFunctionBasis* basis, int nextCoordinates = -1);

  // copy constructor 
  //
  // operation = reference on operation to copy
  QHEParticleWaveFunctionOperation(const QHEParticleWaveFunctionOperation& operation);
  
  // destructor
  //
  ~QHEParticleWaveFunctionOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // get dimension (i.e. Hilbert space dimension)
  //
  // return value = dimension
  int GetDimension ();

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
 protected:

  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

// get dimension (i.e. Hilbert space dimension)
//
// return value = dimension

inline int QHEParticleWaveFunctionOperation::GetDimension ()
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}

#endif
