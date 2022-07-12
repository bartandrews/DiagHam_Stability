////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                             class of S^2 operator                          //
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


#include "Operator/SpinS2Operator.h"
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

SpinS2Operator::SpinS2Operator(AbstractSpinChain* chain, int nbrSpin)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
}

// destructor
//

SpinS2Operator::~SpinS2Operator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* SpinS2Operator::Clone ()
{
  return new SpinS2Operator(this->Chain, this->NbrSpin);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinS2Operator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (AbstractSpinChain*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinS2Operator::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int SpinS2Operator::GetHilbertSpaceDimension ()
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

Complex SpinS2Operator::PartialMatrixElement (RealVector& V1, RealVector& V2, 
					      long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  double Tmp = 0.0;
  int TmpPos;
  double TmpCoef;
  double TmpDiagonal = 0.0;
  for (int i = (int) firstComponent; i < dim; ++i)
    {	
      TmpDiagonal = 0.0; 
      for (int j = 0; j < this->NbrSpin; ++j)
	{
	  int TmpLocalSpin = this->Chain->GetLocalSpin(j,i);
	  TmpDiagonal += 0.25 * ((double) TmpLocalSpin) * ((double) (TmpLocalSpin + 2));
//	  cout <<i<< " "<<j<< " "<< TmpLocalSpin<<endl;
	  for (int k = 0; k < j; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  Tmp += V1[TmpPos] * V2[i] * TmpCoef;
		}
//	      cout <<j<<" " <<k<<" "<<i<< " "<<this->Chain->SziSzj(j, k, i)<<endl;
	      Tmp += 2.0 * V1[i] * V2[i] * this->Chain->SziSzj(j, k, i);
	    }
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  Tmp += V1[TmpPos] * V2[i] * TmpCoef;
		}
	    }
	}
      Tmp += V1[i] * V2[i] * TmpDiagonal;
    }
  return Complex(Tmp);
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex SpinS2Operator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, 
					      long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex Tmp = 0.0;
  int TmpPos;
  double TmpCoef;
  double TmpDiagonal = 0.0;

  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      TmpDiagonal = 0.0;
      for (int j = 0; j < this->NbrSpin; ++j)
	{
	  int TmpLocalSpin = this->Chain->GetLocalSpin(j,i);
//	  cout << TmpLocalSpin <<endl;
	  TmpDiagonal += 0.25 * ((double) TmpLocalSpin) * ((double) (TmpLocalSpin + 2));
//	  cout <<i<< " "<<j<< " "<< TmpLocalSpin<<endl;
	  for (int k = 0; k < j; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  Tmp += Conj(V1[TmpPos]) * V2[i] * TmpCoef;
		}
	      Tmp += 2.0 * Conj(V1[i]) * V2[i] * this->Chain->SziSzj(j, k, i);
//	      cout <<j<<" " <<k<<" "<<i<< " "<<this->Chain->SziSzj(j, k, i)<<endl;
	    }
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  Tmp += Conj(V1[TmpPos]) * V2[i] * TmpCoef;
		}
	    }
	}
      Tmp += Conj(V1[i]) * V2[i] * TmpDiagonal;
    }
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

RealVector& SpinS2Operator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
					     int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  int TmpPos;
  double TmpCoef;
  double TmpDiagonal = 0.0;
  
  vDestination.ClearVector();
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      TmpDiagonal = 0.0;
      for (int j = 0; j < this->NbrSpin; ++j)
	{
	  int TmpLocalSpin = this->Chain->GetLocalSpin(j,i);
	  TmpDiagonal += 0.25 * ((double) TmpLocalSpin) * ((double) (TmpLocalSpin + 2));
	  
	  for (int k = 0; k < j; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  vDestination[TmpPos] += vSource[i] * TmpCoef;
		}
	      vDestination[i] += 2.0 * vSource[i] * this->Chain->SziSzj(j, k, i);
	    }
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  vDestination[TmpPos] += vSource[i] * TmpCoef;
		}
	    }
	}
      vDestination[i] += vSource[i] * TmpDiagonal;
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

ComplexVector& SpinS2Operator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  int TmpPos;
  double TmpCoef;
  double TmpDiagonal = 0.0;
  vDestination.ClearVector();
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      TmpDiagonal = 0.0;
      for (int j = 0; j < this->NbrSpin; ++j)
	{
	  int TmpLocalSpin = this->Chain->GetLocalSpin(j,i);
	  TmpDiagonal += 0.25 * ((double) TmpLocalSpin) * ((double) (TmpLocalSpin + 2));
	  for (int k = 0; k < j; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  vDestination[TmpPos] += vSource[i] * TmpCoef;
		}
	      vDestination[i] += 2.0 * vSource[i] * this->Chain->SziSzj(j, k, i);
	    }
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      TmpPos = this->Chain->SmiSpj(j, k, i, TmpCoef);
	      if (TmpPos < this->Chain->GetHilbertSpaceDimension())
		{
		  vDestination[TmpPos] += vSource[i] * TmpCoef;
		}
	    }
	}
      vDestination[i] += vSource[i] * TmpDiagonal;
    }
  return vDestination;
}

