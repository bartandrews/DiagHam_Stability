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
#include "MathTools/Complex.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif

#include <stdlib.h>


using std::cout;
using std::endl;


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
  
// multiply a vector by the current operator and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractOperator::Multiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, vSource.GetVectorDimension());
}
   

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractOperator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  vDestination.ClearVector();
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current operator for a given range of indices 
// and add result to another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						  int firstComponent, int nbrComponent)
{
  return vDestination;
}


// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractOperator::Multiply(Vector& vSource, Vector& vDestination, 
				      int firstComponent, int nbrComponent)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    return vDestination;
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}

// multiply a vector by the current operator and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vector where result has been stored

ComplexVector& AbstractOperator::Multiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, vSource.GetVectorDimension());
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  vDestination.ClearVector();
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current operator for a given range of indices 
// and add result to another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  return vDestination;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, AbstractOperator& O)
{
  RealVector TmpV2 (O.GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [O.GetHilbertSpaceDimension()];
  for (int i = 0; i < O.GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(O.GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2[i - 1] = 0.0;
	}
      TmpV2[i] = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < O.GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < O.GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j][i];
	  Str << "   ";
	}
      Str << endl;
    }
  return Str;
}
 
#ifdef USE_OUTPUT
// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, AbstractOperator& O)
{
  RealVector TmpV2 (O.GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [O.GetHilbertSpaceDimension()];
  for (int i = 0; i < O.GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(O.GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2[i - 1] = 0.0;
	}
      TmpV2[i] = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (O.GetHilbertSpaceDimension() - 1); ++i)
    {
      Str << "{";
      for (int j = 0; j < (O.GetHilbertSpaceDimension() - 1); ++j)
	{
	  Str << TmpV[j][i];
	  Str << ",";
	}
      Str << TmpV[O.GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (O.GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][O.GetHilbertSpaceDimension() - 1];
      Str << ",";
    }
  Str << TmpV[O.GetHilbertSpaceDimension() - 1][O.GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}
#endif
