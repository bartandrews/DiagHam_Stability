////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of H_plus operator used as building block of the pair       //
//                     hopping hamiltonian in the spin 1 language             //
//                                                                            //
//                        last modification : 27/03/2019                      //
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


#include "Operator/PairHoppingHPlusOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "HilbertSpace/PairHoppingP1AsSpin1Chain.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainLong.h"


using std::cout;
using std::endl;


// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// pValue = value that defines the filling factor p/(2p+1)
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

PairHoppingHPlusOperator::PairHoppingHPlusOperator(AbstractSpinChain* chain, int nbrSpin, int pValue, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PValue = pValue;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
}

// destructor
//

PairHoppingHPlusOperator::~PairHoppingHPlusOperator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* PairHoppingHPlusOperator::Clone ()
{
  return new PairHoppingHPlusOperator(this->Chain, this->NbrSpin, this->PValue, this->PeriodicBoundaryConditions);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PairHoppingHPlusOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (AbstractSpinChain*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* PairHoppingHPlusOperator::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int PairHoppingHPlusOperator::GetHilbertSpaceDimension ()
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

Complex PairHoppingHPlusOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, 
							long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  double Tmp = 0.0;
  int pos;
  int NbrUnitCells = this->NbrSpin / this->PValue;
  if (this->NbrSpin <= 32)
    {
      PairHoppingP1AsSpin1Chain* TmpSpace = (PairHoppingP1AsSpin1Chain*) this->Chain->Clone();
      for (int i = (int) firstComponent; i < dim; ++i)
	{
	  double TmpValue = V2[i];
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      pos = TmpSpace->PlusMinusOperatorPlus(j - 1, j, i);
	      if (pos != dim)
		{
		   Tmp += V1[pos] * TmpValue;
		}
	    }
	  if (this->PeriodicBoundaryConditions == true)
	    {
	      pos = TmpSpace->PlusMinusOperatorPlus(NbrUnitCells - 1, 0, i);
	      if (pos != dim)
		{
		  Tmp += V1[pos] * TmpValue;
		}
	    }      
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  pos = TmpSpace->SwapOperatorPlus(j, p - 1, i);
		  if (pos != dim)
		    {
		      Tmp += V1[pos] * TmpValue;
		    }
		}
	    }
	}
      delete TmpSpace;
   }
  else
    {
      cout << "PairHoppingHPlusOperator::PartialMatrixElement not defined for more than 32 spins" << endl;
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

RealVector& PairHoppingHPlusOperator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						       int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  double Tmp = 0.0;
  int pos;
  int NbrUnitCells = this->NbrSpin / this->PValue;
  if (this->NbrSpin <= 32)
    {
      PairHoppingP1AsSpin1Chain* TmpSpace = (PairHoppingP1AsSpin1Chain*) this->Chain->Clone();
      for (int i = (int) firstComponent; i < dim; ++i)
	{
	  double TmpValue = vSource[i];
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      pos = TmpSpace->PlusMinusOperatorPlus(j - 1, j, i);
	      if (pos != dim)
		{
		  vDestination[pos] += TmpValue;
		}
	    }
	  if (this->PeriodicBoundaryConditions == true)
	    {
	      pos = TmpSpace->PlusMinusOperatorPlus(NbrUnitCells - 1, 0, i);
	      if (pos != dim)
		{
		  vDestination[pos] += TmpValue;
		}
	    }      
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  pos = TmpSpace->SwapOperatorPlus(j, p - 1, i);
		  if (pos != dim)
		    {
		      vDestination[pos] += TmpValue;
		    }
		}
	    }
	}
      delete TmpSpace;
   }
  else
    {
      cout << "PairHoppingHPlusOperator::LowLevelMultiply not defined for more than 32 spins" << endl;
    }
  return vDestination;
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex PairHoppingHPlusOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, 
							long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex Tmp = 0.0;
  int pos;
  int NbrUnitCells = this->NbrSpin / this->PValue;
  if (this->NbrSpin <= 32)
    {
      PairHoppingP1AsSpin1Chain* TmpSpace = (PairHoppingP1AsSpin1Chain*) this->Chain->Clone();
      for (int i = (int) firstComponent; i < dim; ++i)
	{
	  Complex TmpValue = V2[i];
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      pos = TmpSpace->PlusMinusOperatorPlus(j - 1, j, i);
	      if (pos != dim)
		{
		   Tmp += Conj(V1[pos]) * TmpValue;
		}
	    }
	  if (this->PeriodicBoundaryConditions == true)
	    {
	      pos = TmpSpace->PlusMinusOperatorPlus(NbrUnitCells - 1, 0, i);
	      if (pos != dim)
		{
		  Tmp += Conj(V1[pos]) * TmpValue;
		}
	    }      
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  pos = TmpSpace->SwapOperatorPlus(j, p - 1, i);
		  if (pos != dim)
		    {
		      Tmp += Conj(V1[pos]) * TmpValue;
		    }
		}
	    }
	}
      delete TmpSpace;
   }
  else
    {
      cout << "PairHoppingHPlusOperator::PartialMatrixElement not defined for more than 32 spins" << endl;
    }
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

ComplexVector& PairHoppingHPlusOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							  int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  double Tmp = 0.0;
  int pos;
  int NbrUnitCells = this->NbrSpin / this->PValue;
  if (this->NbrSpin <= 32)
    {
      PairHoppingP1AsSpin1Chain* TmpSpace = (PairHoppingP1AsSpin1Chain*) this->Chain->Clone();
      for (int i = (int) firstComponent; i < dim; ++i)
	{
	  Complex TmpValue = vSource[i];
	  for (int j = 1; j < NbrUnitCells; ++j)
	    {
	      pos = TmpSpace->PlusMinusOperatorPlus(j - 1, j, i);
	      if (pos != dim)
		{
		  vDestination[pos] += TmpValue;
		}
	    }
	  if (this->PeriodicBoundaryConditions == true)
	    {
	      pos = TmpSpace->PlusMinusOperatorPlus(NbrUnitCells - 1, 0, i);
	      if (pos != dim)
		{
		  vDestination[pos] += TmpValue;
		}
	    }      
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      for (int p = 1; p < this->PValue; ++p)
		{
		  pos = TmpSpace->SwapOperatorPlus(j, p - 1, i);
		  if (pos != dim)
		    {
		      vDestination[pos] += TmpValue;
		    }
		}
	    }
	}
      delete TmpSpace;
   }
  else
    {
      cout << "PairHoppingHPlusOperator::LowLevelMultiply not defined for more than 32 spins" << endl;
    }
  return vDestination;
}

