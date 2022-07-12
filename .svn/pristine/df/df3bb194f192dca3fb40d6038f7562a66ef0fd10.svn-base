////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of sand glass spin chain hamiltonian                //
//                                                                            //
//                        last modification : 02/06/2003                      //
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


#include "Hamiltonian/SandGlassSpinChainHamiltonian.h"
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
// chain = reference on Hilbert space of the associated system
// nbrSandGlass = number of sand glasses
// jInternal = coupling constant between spins into a sand glass
// jExternal = coupling constant between spins between sand glasses

SandGlassSpinChainHamiltonian::SandGlassSpinChainHamiltonian(AbstractSpinChain* chain, int nbrSandGlass, double jInternal, double jExternal)
{
  this->Chain = chain;
  this->NbrSandGlass = nbrSandGlass;
  this->NbrSpin = this->NbrSandGlass * 5;
  this->JInternal = jInternal;
  this->HalfJInternal = 0.5 * this->JInternal;
  this->JExternal = jExternal;
  this->HalfJExternal = 0.5 * this->JExternal;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

SandGlassSpinChainHamiltonian::~SandGlassSpinChainHamiltonian() 
{
  delete[] this->SzSzContributions;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SandGlassSpinChainHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SandGlassSpinChainHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

SandGlassSpinChainHamiltonian& SandGlassSpinChainHamiltonian::SetChain(AbstractSpinChain* chain)
{  
  delete[] this->SzSzContributions;
  this->Chain = chain;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SandGlassSpinChainHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SandGlassSpinChainHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex SandGlassSpinChainHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  double coef;
  int pos;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int TotalIndex = 0;
  int ReducedNbrSandGlass = this->NbrSandGlass - 1;
  for (int i = 0; i < dim; i++)
    {
     TotalIndex = 0;
      for (int j = 0; j < ReducedNbrSandGlass; j++)
	{
	  pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 2, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 2, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 3, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 4, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 3, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 4, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJInternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex + 5, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJExternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 5, TotalIndex + 2, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJExternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 8, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJExternal * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 8, TotalIndex + 4, i, coef);
	  if (pos != dim)
	    {
	      x += V1[pos] * this->HalfJExternal * coef * V2[i];
	    }
	  TotalIndex += 5;	  
	}
      pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 2, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 2, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 3, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 4, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 3, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 4, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  x += V1[pos] * this->HalfJInternal * coef * V2[i];
	}
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex SandGlassSpinChainHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& SandGlassSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Chain->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& SandGlassSpinChainHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Chain->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& SandGlassSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& SandGlassSpinChainHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of idinces 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& SandGlassSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
								    int firstComponent, int nbrComponent)
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  int lim = nbrComponent + firstComponent;
  double coef;
  int pos;
  int TotalIndex = 0;
  int ReducedNbrSandGlass = this->NbrSandGlass - 1;
  for (int i = firstComponent; i < lim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = firstComponent; i < lim; i++)
    {
      TotalIndex = 0;
      for (int j = 0; j < ReducedNbrSandGlass; j++)
	{
	  pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 3, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 4, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 3, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 4, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex + 5, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJExternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 5, TotalIndex + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJExternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 8, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJExternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 8, TotalIndex + 4, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJExternal * coef * vSource[i];
	    }
	  TotalIndex += 5;	  
	}
      pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 2, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 2, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 3, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 4, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 3, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 4, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
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

RealVector& SandGlassSpinChainHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								       int firstComponent, int nbrComponent)
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  int lim = nbrComponent + firstComponent;
  double coef;
  int pos;
  int TotalIndex = 0;
  int ReducedNbrSandGlass = this->NbrSandGlass - 1;

  for (int i = firstComponent; i < lim; i++)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      TotalIndex = 0;
      for (int j = 0; j < ReducedNbrSandGlass; j++)
	{
	  pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 3, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 4, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 3, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 4, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex + 5, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJExternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 5, TotalIndex + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJExternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 8, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJExternal * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(TotalIndex + 8, TotalIndex + 4, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJExternal * coef * vSource[i];
	    }
	  TotalIndex += 5;	  
	}
      pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 2, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 2, TotalIndex, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex, TotalIndex + 2, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 3, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 3, TotalIndex + 4, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 3, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 1, TotalIndex + 4, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(TotalIndex + 4, TotalIndex + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJInternal * coef * vSource[i];
	}
    }

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

ComplexVector& SandGlassSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								       int firstComponent, int nbrComponent)
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

ComplexVector& SandGlassSpinChainHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									  int firstComponent, int nbrComponent)
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> SandGlassSpinChainHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> SandGlassSpinChainHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate diagonal matrix elements
// 

void SandGlassSpinChainHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  int TotalIndex = 0;
  int ReducedNbrSandGlass = this->NbrSandGlass - 1;
  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] = 0.0;      
      TotalIndex = 0;
      for (int j = 0; j < ReducedNbrSandGlass; j++)
	{
	  this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex, TotalIndex + 1, i);
	  this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 1, TotalIndex + 2, i);
	  this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 2, TotalIndex, i);
	  this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 1, TotalIndex + 3, i);
	  this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 3, TotalIndex + 4, i);
	  this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 4, TotalIndex + 1, i);
	  this->SzSzContributions[i] += this->JExternal * this->Chain->SziSzj(TotalIndex + 2, TotalIndex + 5, i);
	  this->SzSzContributions[i] += this->JExternal * this->Chain->SziSzj(TotalIndex + 4, TotalIndex + 8, i);
	  TotalIndex += 5;	  
	}
      this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex, TotalIndex + 1, i);
      this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 1, TotalIndex + 2, i);
      this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 2, TotalIndex, i);
      this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 1, TotalIndex + 3, i);
      this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 3, TotalIndex + 4, i);
      this->SzSzContributions[i] += this->JInternal * this->Chain->SziSzj(TotalIndex + 4, TotalIndex + 1, i);
     }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, SandGlassSpinChainHamiltonian& H) 
{
  RealVector TmpV2 (H.Chain->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Chain->GetHilbertSpaceDimension(); j++)
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

MathematicaOutput& operator << (MathematicaOutput& Str, SandGlassSpinChainHamiltonian& H) 
{
  RealVector TmpV2 (H.Chain->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Chain->GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Chain->GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Chain->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Chain->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Chain->GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Chain->GetHilbertSpaceDimension() - 1][H.Chain->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

