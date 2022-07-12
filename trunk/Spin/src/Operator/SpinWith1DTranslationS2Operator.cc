////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of S^2 operator for aa 1D system with translation         //
//                                                                            //
//                        last modification : 23/03/2016                      //
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


#include "Operator/SpinWith1DTranslationS2Operator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


using std::cout;
using std::endl;


// constructor from default data
//
// chain = pointer to the Hilbert space
// nbrSpin = number of spins

SpinWith1DTranslationS2Operator::SpinWith1DTranslationS2Operator(AbstractSpinChainWithTranslations* chain, int nbrSpin)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
}

// destructor
//

SpinWith1DTranslationS2Operator::~SpinWith1DTranslationS2Operator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* SpinWith1DTranslationS2Operator::Clone ()
{
  return new SpinWith1DTranslationS2Operator(this->Chain, this->NbrSpin);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinWith1DTranslationS2Operator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (AbstractSpinChainWithTranslations*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinWith1DTranslationS2Operator::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int SpinWith1DTranslationS2Operator::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex SpinWith1DTranslationS2Operator::PartialMatrixElement (RealVector& V1, RealVector& V2, 
							       long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  double Tmp = 0.0;
  int TmpPos;
  int NbrTranslation;
  double TmpCoef;
  double TmpDiagonal = 0.0;
  double Coef = 2.0 * M_PI * ((double) this->Chain->GetMomentum()) / ((double) this->NbrSpin);
  double* ExponentialTable = new double[2 * this->NbrSpin];
  if (this->Chain->GetMomentum() == 0)
    {
      for (int i = 0; i < (2 * this->NbrSpin); ++i)
	{
	  ExponentialTable[i] = 1.0;
	}
    }
  else
    {
      ExponentialTable[0] = 1.0;
      for (int i = 1; i < (2 * this->NbrSpin); ++i)
	{
	  ExponentialTable[i] = -ExponentialTable[i - 1];
	}
    }
  for (int j = 0; j < this->NbrSpin; ++j)
    {      
      int TmpLocalSpin = this->Chain->GetLocalSpin(j);
      TmpDiagonal += 0.25 * ((double) this->Chain->GetLocalSpin(j)) * ((double) (this->Chain->GetLocalSpin(j) + 2));
    }

  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      for (int j = 0; j < this->NbrSpin; ++j)
	{
	  for (int k = 0; k < j; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef, NbrTranslation);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  Tmp += ExponentialTable[NbrTranslation] * V1[TmpPos] * V2[i] * TmpCoef;
		}
	      Tmp += 2.0 * V1[i] * V2[i] * this->Chain->SziSzj(j, k, i);
	    }
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef, NbrTranslation);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  Tmp += ExponentialTable[NbrTranslation] * V1[TmpPos] * V2[i] * TmpCoef;
		}
	    }
	}
      Tmp += V1[i] * V2[i] * TmpDiagonal;
    }
  delete[] ExponentialTable;
  return Complex(Tmp);
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex SpinWith1DTranslationS2Operator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, 
							       long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex Tmp = 0.0;
  int TmpPos;
  int NbrTranslation;
  double TmpCoef;
  double TmpDiagonal = 0.0;
  double Coef = 2.0 * M_PI * ((double) this->Chain->GetMomentum()) / ((double) this->NbrSpin);
  Complex* ExponentialTable = new Complex[2 * this->NbrSpin];
  for (int i = 0; i < (2 * this->NbrSpin); ++i)
    {
      ExponentialTable[i] = Phase(Coef * ((double) i));
    }
  for (int j = 0; j < this->NbrSpin; ++j)
    {      
      int TmpLocalSpin = this->Chain->GetLocalSpin(j);
      TmpDiagonal += 0.25 * ((double) this->Chain->GetLocalSpin(j)) * ((double) (this->Chain->GetLocalSpin(j) + 2));
    }
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      for (int j = 0; j < this->NbrSpin; ++j)
	{
	  for (int k = 0; k < j; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef, NbrTranslation);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  Tmp += ExponentialTable[NbrTranslation] * Conj(V1[TmpPos]) * V2[i] * TmpCoef;
		}
	      Tmp += 2.0 * Conj(V1[i]) * V2[i] * this->Chain->SziSzj(j, k, i);
	    }
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef, NbrTranslation);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  Tmp += ExponentialTable[NbrTranslation] * Conj(V1[TmpPos]) * V2[i] * TmpCoef;
		}
	    }
	}
      Tmp += Conj(V1[i]) * V2[i] * TmpDiagonal;
    }
  delete[] ExponentialTable;
  return Complex(Tmp);
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& SpinWith1DTranslationS2Operator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							      int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  int TmpPos;
  int NbrTranslation;
  double TmpCoef;
  double TmpDiagonal = 0.0;
  double* ExponentialTable = new double[2 * this->NbrSpin];
  if (this->Chain->GetMomentum() == 0)
    {
      for (int i = 0; i < (2 * this->NbrSpin); ++i)
	{
	  ExponentialTable[i] = 1.0;
	}
    }
  else
    {
      ExponentialTable[0] = 1.0;
      for (int i = 1; i < (2 * this->NbrSpin); ++i)
	{
	  ExponentialTable[i] = -ExponentialTable[i - 1];
	}
    }
  vDestination.ClearVector();
  for (int j = 0; j < this->NbrSpin; ++j)
    {      
      int TmpLocalSpin = this->Chain->GetLocalSpin(j);
      TmpDiagonal += 0.25 * ((double) this->Chain->GetLocalSpin(j)) * ((double) (this->Chain->GetLocalSpin(j) + 2));
    }
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      for (int j = 0; j < this->NbrSpin; ++j)
	{
	  for (int k = 0; k < j; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef, NbrTranslation);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  vDestination[TmpPos] += ExponentialTable[NbrTranslation] * vSource[i] * TmpCoef;
		}
	      vDestination[i] += 2.0 * vSource[i] * this->Chain->SziSzj(j, k, i);
	    }
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef, NbrTranslation);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  vDestination[TmpPos] += ExponentialTable[NbrTranslation] * vSource[i] * TmpCoef;
		}
	    }
	}
      vDestination[i] += vSource[i] * TmpDiagonal;
    }
  delete[] ExponentialTable;
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

ComplexVector& SpinWith1DTranslationS2Operator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								 int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  int TmpPos;
  int NbrTranslation;
  double TmpCoef;
  double TmpDiagonal = 0.0;
  Complex* ExponentialTable = new Complex[2 * this->NbrSpin];
  double Coef = 2.0 * M_PI * ((double) this->Chain->GetMomentum()) / ((double) this->NbrSpin);
  vDestination.ClearVector();
  for (int i = 0; i < (2 * this->NbrSpin); ++i)
    {
      ExponentialTable[i] = Phase(Coef * ((double) i));
    }
  for (int j = 0; j < this->NbrSpin; ++j)
    {      
      int TmpLocalSpin = this->Chain->GetLocalSpin(j);
      TmpDiagonal += 0.25 * ((double) this->Chain->GetLocalSpin(j)) * ((double) (this->Chain->GetLocalSpin(j) + 2));
    }
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      for (int j = 0; j < this->NbrSpin; ++j)
	{
	  for (int k = 0; k < j; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef, NbrTranslation);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  vDestination[TmpPos] += ExponentialTable[NbrTranslation] * vSource[i] * TmpCoef;
		}
	      vDestination[i] += 2.0 * vSource[i] * this->Chain->SziSzj(j, k, i);
	    }
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef, NbrTranslation);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  vDestination[TmpPos] += ExponentialTable[NbrTranslation] * vSource[i] * TmpCoef;
		}
	    }
	}
      vDestination[i] += vSource[i] * TmpDiagonal;
    }
  delete[] ExponentialTable;
  return vDestination;
}

