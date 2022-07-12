////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of FQHE on disk quasihole propagator operation          //
//                                                                            //
//                        last modification : 08/03/2009                      //
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


#include "config.h"
#include "Architecture/ArchitectureOperation/FQHEDiskQuasiholePropagatorOperation.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"


// constructor 
//
// space = pointer to the Hilbert space to use
// jackPolynomial = vector corresponding to the Jack polynomial

FQHEDiskQuasiholePropagatorOperation::FQHEDiskQuasiholePropagatorOperation (ParticleOnSphere* space, RealVector* jackPolynomial)
{
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = space->GetLargeHilbertSpaceDimension();
  this->HilbertSpace = (ParticleOnSphere*) space->Clone();
  this->JackPolynomial = jackPolynomial;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Scalars = 0;
  this->NbrScalars = 1;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHEDiskQuasiholePropagatorOperation::FQHEDiskQuasiholePropagatorOperation(const FQHEDiskQuasiholePropagatorOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->JackPolynomial = operation.JackPolynomial;
  this->HilbertSpace = (ParticleOnSphere*) operation.HilbertSpace->Clone();
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Scalar = operation.Scalar;
  this->Scalars = 0; 
  this->NbrScalars = 1;
}
  
// destructor
//

FQHEDiskQuasiholePropagatorOperation::~FQHEDiskQuasiholePropagatorOperation()
{
  delete this->HilbertSpace;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHEDiskQuasiholePropagatorOperation::SetLargeIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->LargeFirstComponent = firstComponent;
  this->LargeNbrComponent = nbrComponent;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHEDiskQuasiholePropagatorOperation::Clone()
{
  return new FQHEDiskQuasiholePropagatorOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHEDiskQuasiholePropagatorOperation::RawApplyOperation()
{
  this->Scalar = this->HilbertSpace->JackSqrNormalization(*(this->JackPolynomial), this->LargeFirstComponent, this->LargeNbrComponent);
  return true;
}

