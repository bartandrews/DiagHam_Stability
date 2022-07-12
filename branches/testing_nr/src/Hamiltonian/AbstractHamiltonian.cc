////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract hamiltonian                      //
//                                                                            //
//                        last modification : 28/02/2001                      //
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


#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "MathTools/Complex.h"
#include "BitmapTools/BitmapPicture/AbstractBitmapPicture.h"
#include "BitmapTools/BitmapPicture/TgaFormat.h" 
#include "BitmapTools/Color/PicRGB.h"


#include <stdlib.h>


using std::cout;
using std::endl;


// destructor
//

AbstractHamiltonian::~AbstractHamiltonian() 
{
}

// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

HermitianMatrix& AbstractHamiltonian::GetHamiltonian (HermitianMatrix& M)
{
  ComplexVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1.Re(i) = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	{
	  M.SetMatrixElement(i, j, TmpV2[j]);
	}
      TmpV1.Re(i) = 0.0;
    }
  return M;  
}
  
// store real part of Hamiltonian into a real symmetric matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding real symmetric matrix 

RealSymmetricMatrix& AbstractHamiltonian::GetHamiltonian (RealSymmetricMatrix& M)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	M.SetMatrixElement(i, j, TmpV2[j]);
      TmpV1[i] = 0.0;	
    }
  return M;
}
  
// store Hamiltonian into a matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding matrix 

Matrix& AbstractHamiltonian::GetHamiltonian (Matrix& M)
{
  switch (M.GetMatrixType())
    {
      case (Matrix::RealElements | Matrix::Symmetric):
	return this->GetHamiltonian((RealSymmetricMatrix&) M);
      break;
      case (Matrix::ComplexElements | Matrix::Hermitian):
	return this->GetHamiltonian((HermitianMatrix&) M);
      break;
    default:
      return M;
      break;
    }
  return M;
}
  
// return matrix representation of current Hamiltonian
//
// return value = reference to representation

Matrix* AbstractHamiltonian::GetHamiltonian ()
{
  HermitianMatrix* TmpH = new HermitianMatrix(this->GetHilbertSpaceDimension());
  this->GetHamiltonian(*TmpH);
  return TmpH;
}
  
// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						  int firstComponent, int nbrComponent)
{
  return this->LowLevelMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							  int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vector by the current hamiltonian for a given range of indices 
// and store result in another set of vector, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			   destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vector sat which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							     int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}
 
// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result in another set of vectors, low level function (no architecture optimization)
//
// vSource = array of vectors to be multiplied
// vDestination = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelAddMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			      destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  return this->LowLevelMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int firstComponent, int nbrComponent)
{
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							     int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			   destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}
 

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
								int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelAddMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			      destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// store Hamiltonian into a picture (drawing non zero element with a color scale)
//
// error = absolute minimum value to be considered as non zero element
// return value = pointer to the picture associated to the matrix

AbstractBitmapPicture* AbstractHamiltonian::GetHamiltonianPicture (double error)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  Color BlackColor (0.0, 0.0, 0.0);
  Color WhiteColor (1.0, 1.0, 1.0);
  TgaFormat* TmpPicture = new TgaFormat (this->GetHilbertSpaceDimension(), this->GetHilbertSpaceDimension());
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < this->GetHilbertSpaceDimension(); j++)
	{
	  if (fabs(TmpV2[j]) < error)
	    {
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, WhiteColor);
	    }
	  else
	    {
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, BlackColor);
	    }
	}
      TmpV1[i] = 0.0;	
    }
  return TmpPicture;
}
  
// store Hamiltonian into a picture (drawing non zero element in black)
//
// error = absolute minimum value to be considered as non zero element
// return value = pointer to the picture associated to the matrix

AbstractBitmapPicture* AbstractHamiltonian::GetHamiltonianColorPicture (double error)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  Color RedColor (1.0, 0.0, 0.0);
  Color GreenColor (0.0, 1.0, 0.0);
  Color BlueColor (0.0, 0.0, 1.0);
  Color BlackColor (0.0, 0.0, 0.0);
  TgaFormat* TmpPicture = new TgaFormat (this->GetHilbertSpaceDimension(), this->GetHilbertSpaceDimension());
  double Max = 0.0;
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < this->GetHilbertSpaceDimension(); j++)
	{
	  if (fabs(TmpV2[j]) > error)
	    {
	      if (fabs(TmpV2[j]) > Max)
		{
		  Max = fabs(TmpV2[j]);
		}
	    }
	}
      TmpV1[i] = 0.0;	
    }
  Max = 4.0 / Max;
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < this->GetHilbertSpaceDimension(); j++)
	{
	  if (fabs(TmpV2[j]) < error)
	    {
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, BlackColor);
	    }
	  else
	    {
	      double Fac = (fabs(TmpV2[j]) * Max);
	      Color TmpColor;
	      if (Fac >= 4.0)
		{
		  TmpColor = RedColor;
		}
	      else
		if (Fac >= 3.0)
		  {
		    TmpColor = GreenColor * (4.0 - Fac) + RedColor;
		  }
		else
		  if (Fac >= 2.0)
		    {
		      TmpColor = GreenColor + RedColor * (Fac - 2.0);
		    }
		  else
		    if (Fac >= 1.0)
		      {
			TmpColor = BlueColor * (2.0 - Fac) + GreenColor;
		      }
		    else
		      {
			TmpColor = GreenColor * (Fac - 1.0) + BlueColor;
		      }
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, TmpColor);
	    }
	}
      TmpV1[i] = 0.0;	
    }
  return TmpPicture;
}
  
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractHamiltonian::LeftInteractionOperators()
{
  return List<Matrix*>();
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractHamiltonian::RightInteractionOperators()
{
  return List<Matrix*>();
}

// multiply a vector by the current hamiltonian and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::Multiply(Vector& vSource, Vector& vDestination)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    {
      return vDestination;
    }
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->LowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::Multiply(Vector& vSource, Vector& vDestination, 
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

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::AddMultiply(Vector& vSource, Vector& vDestination)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    {
      return vDestination;
    }
  if (vSource.GetVectorType() == Vector::RealDatas)
    {
      return this->LowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->LowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::AddMultiply(Vector& vSource, Vector& vDestination, 
					 int firstComponent, int nbrComponent)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    return vDestination;
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::MultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultipleMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors);
    }
  else
    {
      return this->LowLevelMultipleMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::MultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors, 
					      int firstComponent, int nbrComponent)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  cout << "check 0" << endl;
  for (int i = 1; i < nbrVectors; ++i)
    {
      if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
	return vDestinations;
    }
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      cout << "check real " << endl;
      return this->LowLevelMultipleMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelMultipleMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::MultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultipleAddMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors);
    }
  else
    {
      return this->LowLevelMultipleAddMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::MultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors,
						 int firstComponent, int nbrComponent)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultipleAddMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelMultipleAddMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
}

