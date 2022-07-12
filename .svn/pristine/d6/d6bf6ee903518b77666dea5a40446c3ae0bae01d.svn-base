////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of multiple spin chain hamiltonian                 //
//                                                                            //
//                        last modification : 12/09/2002                      //
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


#include "Hamiltonian/MultipleSpinChainHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::endl;
using std::ostream;


// constructor from default datas
//
// chains = pointer to Hilbert space of the associated system
// nbrSpinPerChain = number of spins per chain
// nbrChain = number of chains
// j1 = coupling constant between nearest neighbours in a chain
// j2 = coupling constant between nearest neighbours in two diffferent succesive chains
// j3 = first coupling constant between second nearest neighbours in two diffferent succesive chains 
// j4 = second coupling constant between second nearest neighbours in two diffferent succesive chains

MultipleSpinChainHamiltonian::MultipleSpinChainHamiltonian(AbstractSpinChain* chains, int nbrSpinPerChain, int nbrChain, double j1, double j2, double j3, double j4)
{
  this->Chains = chains;
  this->NbrSpinPerChain = nbrSpinPerChain;
  this->NbrChain = nbrChain;
  this->NbrSpin = this->NbrChain * this->NbrSpinPerChain;
  this->J1 = j1;
  this->J2 = j2;
  this->J3 = j3;
  this->J4 = j4;
  this->HalfJ1 = 0.5 * this->J1;
  this->HalfJ2 = 0.5 * this->J2;
  this->HalfJ3 = 0.5 * this->J3;
  this->HalfJ4 = 0.5 * this->J4;
  this->SzSzContributions = new double [this->Chains->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

MultipleSpinChainHamiltonian::~MultipleSpinChainHamiltonian() 
{
  delete[] this->SzSzContributions;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void MultipleSpinChainHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chains = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chains->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* MultipleSpinChainHamiltonian::GetHilbertSpace ()
{
  return this->Chains;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

MultipleSpinChainHamiltonian& MultipleSpinChainHamiltonian::SetChain(AbstractSpinChain* chain)
{  
  delete[] this->SzSzContributions;
  this->Chains = chain;
  this->SzSzContributions = new double [this->Chains->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int MultipleSpinChainHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chains->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void MultipleSpinChainHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chains->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs
bool MultipleSpinChainHamiltonian::SavePrecalculation (char* fileName)
{
  return false;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex MultipleSpinChainHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  double Coef;
  int Pos;
  int dim = this->Chains->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {

      // SzSz contributions
      x += V1[i] * this->SzSzContributions[i] * V2[i];
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex MultipleSpinChainHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& MultipleSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Chains->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& MultipleSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						   int firstComponent, int nbrComponent) 
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Chains->GetHilbertSpaceDimension();
  double Coef;
  int Pos;
  int ReducedNbrChain = this->NbrChain - 1;
  int ReducedNbrSpin = this->NbrSpin - 1;
  int ReducedNbrSpinPerChain = this->NbrSpinPerChain - 1;
  int Lim;
  double TmpHalfJ1;
  double TmpHalfJ2;
  double TmpHalfJ3;
  double TmpHalfJ4;
  for (int i = firstComponent; i < LastComponent; ++i)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      TmpHalfJ1 = this->HalfJ1 * vSource[i];
      TmpHalfJ2 = this->HalfJ2 * vSource[i];
      TmpHalfJ3 = this->HalfJ3 * vSource[i];
      TmpHalfJ4 = this->HalfJ4 * vSource[i];
      for (int j = 0; j < ReducedNbrSpinPerChain; ++j)
	{
	  Lim = j * this->NbrChain + ReducedNbrChain;
	  for (int k = j * this->NbrChain; k < Lim; ++k)
	    {
	      Pos = this->Chains->SmiSpj(k, k + 1, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ2 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + 1, k, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ2 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k, k + this->NbrChain, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ1 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + this->NbrChain, k, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ1 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k, k + this->NbrChain + 1, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ3 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + this->NbrChain + 1, k, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ3 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + 1, k + this->NbrChain, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ4 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + this->NbrChain, k + 1, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ4 * Coef;
		}
	    }
	  Pos = this->Chains->SmiSpj(Lim, j * this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ2 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(j * this->NbrChain, Lim, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ2 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(Lim, Lim + this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ1 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(Lim + this->NbrChain, Lim, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ1 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(Lim, (j + 1) * this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ3 * Coef;
	    }
	  Pos = this->Chains->SmiSpj((j + 1) * this->NbrChain, Lim, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ3 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(j * this->NbrChain, Lim + this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ4 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(Lim + this->NbrChain, j * this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ4 * Coef;
	    }
	}

      for (int j = this->NbrSpin - this->NbrChain; j < ReducedNbrSpin; ++j)
	{
	  Pos = this->Chains->SmiSpj(j, j + 1, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ2 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(j + 1, j, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ2 * Coef;
	    }
	}
      Pos = this->Chains->SmiSpj(ReducedNbrSpin, this->NbrSpin - this->NbrChain, i, Coef);     
      if (Pos != Dim)
	{
	  vDestination[Pos] += TmpHalfJ2 * Coef;
	}
      Pos = this->Chains->SmiSpj(this->NbrSpin - this->NbrChain, ReducedNbrSpin, i, Coef);     
      if (Pos != Dim)
	{
	  vDestination[Pos] += TmpHalfJ2 * Coef;
	}
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& MultipleSpinChainHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Chains->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& MultipleSpinChainHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							      int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  double Coef;
  int Pos;
  int Dim = this->Chains->GetHilbertSpaceDimension();
  int ReducedNbrChain = this->NbrChain - 1;
  int ReducedNbrSpin = this->NbrSpin - 1;
  int ReducedNbrSpinPerChain = this->NbrSpinPerChain - 1;
  int Lim;
  double TmpHalfJ1;
  double TmpHalfJ2;
  double TmpHalfJ3;
  double TmpHalfJ4;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      TmpHalfJ1 = this->HalfJ1 * vSource[i];
      TmpHalfJ2 = this->HalfJ2 * vSource[i];
      TmpHalfJ3 = this->HalfJ3 * vSource[i];
      TmpHalfJ4 = this->HalfJ4 * vSource[i];

      for (int j = 0; j < ReducedNbrSpinPerChain; ++j)
	{
	  Lim = j * this->NbrChain + ReducedNbrChain;
	  for (int k = j * this->NbrChain; k < Lim; ++k)
	    {
	      Pos = this->Chains->SmiSpj(k, k + 1, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ2 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + 1, k, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ2 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k, k + this->NbrChain, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ1 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + this->NbrChain, k, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ1 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k, k + this->NbrChain + 1, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ3 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + this->NbrChain + 1, k, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ3 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + 1, k + this->NbrChain, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ4 * Coef;
		}
	      Pos = this->Chains->SmiSpj(k + this->NbrChain, k + 1, i, Coef);
	      if (Pos != Dim)
		{
		  vDestination[Pos] += TmpHalfJ4 * Coef;
		}
	    }
	  Pos = this->Chains->SmiSpj(Lim, j * this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ2 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(j * this->NbrChain, Lim, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ2 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(Lim, Lim + this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ1 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(Lim + this->NbrChain, Lim, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ1 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(Lim, (j + 1) * this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ3 * Coef;
	    }
	  Pos = this->Chains->SmiSpj((j + 1) * this->NbrChain, Lim, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ3 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(j * this->NbrChain, Lim + this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ4 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(Lim + this->NbrChain, j * this->NbrChain, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ4 * Coef;
	    }
	}

      for (int j = this->NbrSpin - this->NbrChain; j < ReducedNbrSpin; ++j)
	{
	  Pos = this->Chains->SmiSpj(j, j + 1, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ2 * Coef;
	    }
	  Pos = this->Chains->SmiSpj(j + 1, j, i, Coef);
	  if (Pos != Dim)
	    {
	      vDestination[Pos] += TmpHalfJ2 * Coef;
	    }
	}
      Pos = this->Chains->SmiSpj(ReducedNbrSpin, this->NbrSpin - this->NbrChain, i, Coef);     
      if (Pos != Dim)
	{
	  vDestination[Pos] += TmpHalfJ2 * Coef;
	}
      Pos = this->Chains->SmiSpj(this->NbrSpin - this->NbrChain, ReducedNbrSpin, i, Coef);     
      if (Pos != Dim)
	{
	  vDestination[Pos] += TmpHalfJ2 * Coef;
	}
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& MultipleSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& MultipleSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						   int firstComponent, int nbrComponent)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& MultipleSpinChainHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
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
ComplexVector& MultipleSpinChainHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int firstComponent, int nbrComponent)
{
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> MultipleSpinChainHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  int Dim = this->Chains->GetHilbertSpaceDimension();
  for (int i = 0; i < this->NbrChain; ++i)
    {
      RealSymmetricMatrix* Sx = new RealSymmetricMatrix (Dim, true);
      RealAntisymmetricMatrix* Sy = new RealAntisymmetricMatrix (Dim, true);
      RealSymmetricMatrix* Sz = new RealSymmetricMatrix (Dim, true);
      this->Chains->Sxi(i, *Sx);
      this->Chains->Syi(i, *Sy);
      this->Chains->Szi(i, *Sz);
      TmpList += Sx;
      TmpList += Sy;
      TmpList += Sz;
    }
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> MultipleSpinChainHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  int Dim = this->Chains->GetHilbertSpaceDimension();
  for (int i = (this->NbrSpin - this->NbrChain); i < this->NbrSpin; ++i)
    {
      RealSymmetricMatrix* Sx = new RealSymmetricMatrix (Dim, true);
      RealAntisymmetricMatrix* Sy = new RealAntisymmetricMatrix (Dim, true);
      RealSymmetricMatrix* Sz = new RealSymmetricMatrix (Dim, true);
      this->Chains->Sxi(i, *Sx);
      this->Chains->Syi(i, *Sy);
      this->Chains->Szi(i, *Sz);
      TmpList += Sx;
      TmpList += Sy;
      TmpList += Sz;
    }
  return TmpList;
}

// evaluate diagonal matrix elements
// 

void MultipleSpinChainHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chains->GetHilbertSpaceDimension();
  int ReducedNbrChain = this->NbrChain - 1;
  int ReducedNbrSpin = this->NbrSpin - 1;
  int ReducedNbrSpinPerChain = this->NbrSpinPerChain - 1;
  int Lim;
  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < ReducedNbrSpinPerChain; j++)
	{
	  Lim = j * this->NbrChain + ReducedNbrChain;
	  for (int k = j * this->NbrChain; k < Lim; ++k)
	    {
	      this->SzSzContributions[i] += this->J2 * this->Chains->SziSzj(k, k + 1, i);
	      this->SzSzContributions[i] += this->J1 * this->Chains->SziSzj(k, k + this->NbrChain, i);
	      this->SzSzContributions[i] += this->J3 * this->Chains->SziSzj(k, k + this->NbrChain + 1, i);
	      this->SzSzContributions[i] += this->J4 * this->Chains->SziSzj(k + 1, k + this->NbrChain, i);
	    }
	  this->SzSzContributions[i] += this->J2 * this->Chains->SziSzj(Lim, j * this->NbrChain, i);
	  this->SzSzContributions[i] += this->J1 * this->Chains->SziSzj(Lim, Lim + this->NbrChain, i);
	  this->SzSzContributions[i] += this->J3 * this->Chains->SziSzj(Lim, (j + 1) * this->NbrChain, i);
	  this->SzSzContributions[i] += this->J4 * this->Chains->SziSzj(j * this->NbrChain, Lim + this->NbrChain, i);
	}

      for (int j = this->NbrSpin - this->NbrChain; j < ReducedNbrSpin; ++j)
	{
	  this->SzSzContributions[i] += this->J2 * this->Chains->SziSzj(j, j + 1, i);
	}
      this->SzSzContributions[i] += this->J2 * this->Chains->SziSzj(ReducedNbrSpin, this->NbrSpin - this->NbrChain, i);     
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, MultipleSpinChainHamiltonian& H) 
{
  RealVector TmpV2 (H.Chains->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Chains->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chains->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Chains->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Chains->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Chains->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j][i] << "    ";
	}
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, MultipleSpinChainHamiltonian& H) 
{
  RealVector TmpV2 (H.Chains->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Chains->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chains->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Chains->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Chains->GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Chains->GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Chains->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Chains->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Chains->GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Chains->GetHilbertSpaceDimension() - 1][H.Chains->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

