////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of abstract operator                       //
//                                                                            //
//                        last modification : 22/03/2002                      //
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


#include "Operator/AbstractOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Complex.h"

#include <stdlib.h>


// destructor
//

AbstractOperator::~AbstractOperator() 
{
}

// store operator into an hermitian matrix
//
// M = reference on matrix where operator has to be stored
// return value = reference on  corresponding hermitian matrix

HermitianMatrix& AbstractOperator::GetOperator (HermitianMatrix& M)
{
  ComplexVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = Complex(1.0, 0.0);
      this->Multiply(TmpV1, TmpV2);
      for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	M.SetMatrixElement(i, j, TmpV2[j]);
      TmpV1[i] = Complex(0.0, 0.0);	
    }
  return M;
}
  
// store real part of operator into a real symmetric matrix
//
// M = reference on matrix where operator has to be stored
// return value = reference on  corresponding real symmetric matrix 

RealSymmetricMatrix& AbstractOperator::GetOperator (RealSymmetricMatrix& M)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->Multiply(TmpV1, TmpV2);
      for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	M.SetMatrixElement(i, j, TmpV2[j]);
      TmpV1[i] = 0.0;	
    }
  return M;
}
  
// store operator into a matrix
//
// M = reference on matrix where operator has to be stored
// return value = reference on  corresponding matrix 

Matrix& AbstractOperator::GetOperator (Matrix& M)
{
  switch (M.GetMatrixType())
    {
      case (Matrix::RealElements | Matrix::Symmetric):
	return this->GetOperator((RealSymmetricMatrix&) M);
      break;
      case (Matrix::ComplexElements | Matrix::Hermitian):
	return this->GetOperator((HermitianMatrix&) M);
      break;
    default:
      return M;
      break;
    }
  return M;
}
  
// return matrix representation of current operator
//
// return value = reference to representation

Matrix* AbstractOperator::GetOperator ()
{
  HermitianMatrix* TmpH = new HermitianMatrix(this->GetHilbertSpaceDimension());
  this->GetOperator(*TmpH);
  return TmpH;
}
  

