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
#include "Vector/LongRationalVector.h"
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
  this->RationalJackPolynomial = 0;
  this->JackPolynomial2 = 0;
  this->RationalJackPolynomial2 = 0;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Scalars = 0;
  this->LongRationalScalars = 0;
  this->NbrScalars = 1;
}

// constructor 
//
// space = pointer to the Hilbert space to use
// jackPolynomial = vector corresponding to the Jack polynomial

FQHEDiskQuasiholePropagatorOperation::FQHEDiskQuasiholePropagatorOperation (ParticleOnSphere* space, LongRationalVector* jackPolynomial)
{
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = space->GetLargeHilbertSpaceDimension();
  this->HilbertSpace = (ParticleOnSphere*) space->Clone();
  this->RationalJackPolynomial = jackPolynomial;
  this->JackPolynomial = 0;
  this->JackPolynomial2 = 0;
  this->RationalJackPolynomial2 = 0;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Scalars = 0;
  this->LongRationalScalars = 0;
  this->NbrScalars = 1;
}

// constructor for a scalar product
//
// space = pointer to the Hilbert space to use
// jackPolynomial1 = vector corresponding to the first Jack polynomial
// jackPolynomial2 = vector corresponding to the second Jack polynomial

FQHEDiskQuasiholePropagatorOperation::FQHEDiskQuasiholePropagatorOperation(ParticleOnSphere* space, RealVector* jackPolynomial1, RealVector* jackPolynomial2)
{
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = space->GetLargeHilbertSpaceDimension();
  this->HilbertSpace = (ParticleOnSphere*) space->Clone();
  this->JackPolynomial = jackPolynomial1;
  this->RationalJackPolynomial = 0;
  this->JackPolynomial2 = jackPolynomial2;
  this->RationalJackPolynomial2 = 0;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Scalars = 0;
  this->LongRationalScalars = 0;
  this->NbrScalars = 1;
}

// constructor  for a scalar product
//
// space = pointer to the Hilbert space to use
// jackPolynomial1 = vector corresponding to the first Jack polynomial
// jackPolynomial2 = vector corresponding to the second Jack polynomial

FQHEDiskQuasiholePropagatorOperation::FQHEDiskQuasiholePropagatorOperation(ParticleOnSphere* space, LongRationalVector* jackPolynomial1, LongRationalVector* jackPolynomial2)
{
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = space->GetLargeHilbertSpaceDimension();
  this->HilbertSpace = (ParticleOnSphere*) space->Clone();
  this->RationalJackPolynomial = jackPolynomial1;
  this->JackPolynomial = 0;
  this->JackPolynomial2 = 0;
  this->RationalJackPolynomial2 = jackPolynomial2;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Scalars = 0;
  this->LongRationalScalars = 0;
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
  this->RationalJackPolynomial = operation.RationalJackPolynomial;
  this->JackPolynomial2 = operation.JackPolynomial2;
  this->RationalJackPolynomial2 = operation.RationalJackPolynomial2;
  this->HilbertSpace = (ParticleOnSphere*) operation.HilbertSpace->Clone();
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Scalar = operation.Scalar;
  this->LongRationalScalar = operation.LongRationalScalar;
  this->Scalars = 0; 
  this->LongRationalScalars = 0;
  this->NbrScalars = 1;
}
  
// destructor
//

FQHEDiskQuasiholePropagatorOperation::~FQHEDiskQuasiholePropagatorOperation()
{
  delete this->HilbertSpace;
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHEDiskQuasiholePropagatorOperation::Clone()
{
  return new FQHEDiskQuasiholePropagatorOperation (*this);
}
  
// get dimension (i.e. Hilbert space dimension, nbr of subdivisions,...), return 0 if large number are required
// 
// return value = dimension  

int FQHEDiskQuasiholePropagatorOperation::GetDimension ()
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}


// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHEDiskQuasiholePropagatorOperation::RawApplyOperation()
{
  if ((this->RationalJackPolynomial2 == 0) && (this->JackPolynomial2 == 0))
    {
      if (this->RationalJackPolynomial == 0)
	{
	  this->Scalar = this->HilbertSpace->JackSqrNormalization(*(this->JackPolynomial), this->LargeFirstComponent, this->LargeNbrComponent);
	}
      else
	{
	  this->LongRationalScalar = this->HilbertSpace->JackSqrNormalization(*(this->RationalJackPolynomial), this->LargeFirstComponent, this->LargeNbrComponent);
	}
    }
  else
    {
      if (this->RationalJackPolynomial == 0)
	{
	  this->Scalar = this->HilbertSpace->JackScalarProduct(*(this->JackPolynomial), *(this->JackPolynomial2), this->LargeFirstComponent, this->LargeNbrComponent);
	}
      else
	{
	  this->LongRationalScalar = this->HilbertSpace->JackScalarProduct(*(this->RationalJackPolynomial), *(this->RationalJackPolynomial2), this->LargeFirstComponent, this->LargeNbrComponent);
	}
    }
  return true;
}

