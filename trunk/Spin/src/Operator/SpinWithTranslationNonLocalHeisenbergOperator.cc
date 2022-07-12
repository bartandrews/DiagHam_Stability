////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of non-local symmetry operator for the Heisenberg model        //
//                            with translations                               //
//                                                                            //
//                        last modification : 26/04/2022                      //
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


#include "Operator/SpinWithTranslationNonLocalHeisenbergOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


using std::cout;
using std::endl;


// constructor from default datas
//
// chain = pointer to the Hilbert space
// nbrSpin = number of spins
// periodicBoundaryConditions = use periodic boundary conditions

SpinWithTranslationNonLocalHeisenbergOperator::SpinWithTranslationNonLocalHeisenbergOperator(AbstractSpinChainWithTranslations* chain, int nbrSpin)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
}

// destructor
//

SpinWithTranslationNonLocalHeisenbergOperator::~SpinWithTranslationNonLocalHeisenbergOperator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* SpinWithTranslationNonLocalHeisenbergOperator::Clone ()
{
  return new SpinWithTranslationNonLocalHeisenbergOperator(this->Chain, this->NbrSpin);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinWithTranslationNonLocalHeisenbergOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (AbstractSpinChainWithTranslations*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinWithTranslationNonLocalHeisenbergOperator::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int SpinWithTranslationNonLocalHeisenbergOperator::GetHilbertSpaceDimension ()
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

Complex SpinWithTranslationNonLocalHeisenbergOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, 
									     long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  double Tmp = 0.0;
  // int TmpPos;
  // double TmpCoef;
  // for (int i = (int) firstComponent; i < dim; ++i)
  //   {	
  //     for (int j = 0; j < this->NbrSpin; ++j)
  // 	{
  // 	  for (int k = j + 1; k < this->NbrSpin; ++k)
  // 	    {
  // 	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
  // 	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
  // 		{
  // 		  Tmp += V1[TmpPos] * V2[i] * TmpCoef;
  // 		}
  // 	      TmpPos = this->Chain->SmiSpj(k, j, i, TmpCoef);
  // 	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
  // 		{
  // 		  Tmp -= V1[TmpPos] * V2[i] * TmpCoef;
  // 		}
  // 	    }
  // 	}
  //   }
  return Complex(Tmp);
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex SpinWithTranslationNonLocalHeisenbergOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, 
									     long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex Tmp = 0.0;
  // int TmpPos;
  // double TmpCoef;
  // double TmpDiagonal = 0.0;

  // for (int i = (int) firstComponent; i < dim; ++i)
  //   {	 
  //     TmpDiagonal = 0.0;
  //     for (int j = 0; j < this->NbrSpin; ++j)
  // 	{
  // 	  for (int k = j + 1; k < this->NbrSpin; ++k)
  // 	    {
  // 	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
  // 	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
  // 		{
  // 		  Tmp += Conj(V1[TmpPos]) * V2[i] * TmpCoef;
  // 		}
  // 	      TmpPos = this->Chain->SmiSpj(k, j, i, TmpCoef);
  // 	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
  // 		{
  // 		  Tmp -= Conj(V1[TmpPos]) * V2[i] * TmpCoef;
  // 		}
  // 	    }
  // 	}
  //    }
  return Tmp;
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& SpinWithTranslationNonLocalHeisenbergOperator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
									    int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  int TmpPos;
  double TmpCoef;
  int TmpNbrTranslations;
  
  vDestination.ClearVector();

  for (int i = (int) firstComponent; i < dim; ++i)
    {
      for (int j = 0; j < (this->NbrSpin - 2); ++j)
	{
	  TmpPos = this->Chain->Pijk(j, j + 1, j + 2, i, TmpCoef, TmpNbrTranslations);
	  if (TmpPos < this->Chain->GetHilbertSpaceDimension())
	    {
	      vDestination[TmpPos] += vSource[i];
	    }
	  TmpPos = this->Chain->Pminusijk(j, j + 1, j + 2, i, TmpCoef, TmpNbrTranslations);
	  if (TmpPos < this->Chain->GetHilbertSpaceDimension())
	    {
	      vDestination[TmpPos] -= vSource[i];
	    }
	}
      TmpPos = this->Chain->Pijk(0, this->NbrSpin - 2, this->NbrSpin - 1, i, TmpCoef, TmpNbrTranslations);
      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
	{
	  vDestination[TmpPos] += vSource[i];
	}
      TmpPos = this->Chain->Pminusijk(0, this->NbrSpin - 2, this->NbrSpin - 1, i, TmpCoef, TmpNbrTranslations);
      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
	{
	  vDestination[TmpPos] -= vSource[i];
	}
      TmpPos = this->Chain->Pijk(0, 1, this->NbrSpin - 1, i, TmpCoef, TmpNbrTranslations);
      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
	{
	  vDestination[TmpPos] += vSource[i];
	}
      TmpPos = this->Chain->Pminusijk(0, 1, this->NbrSpin - 1, i, TmpCoef, TmpNbrTranslations);
      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
	{
	  vDestination[TmpPos] -= vSource[i];
	}
    }
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

ComplexVector& SpinWithTranslationNonLocalHeisenbergOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									       int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  int TmpPos;
  double TmpCoef;
  vDestination.ClearVector();
  return vDestination;
}

