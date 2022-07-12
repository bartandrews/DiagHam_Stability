////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of generic XXZ spin chain hamiltonian               //
//                                                                            //
//                        last modification : 10/05/2001                      //
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


#include "Hamiltonian/XXZSpinChainHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::endl;
using std::ostream;


// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin 1
// jx = coupling constants between Sx components
// jz = coupling constants between Sz components

XXZSpinChainHamiltonian::XXZSpinChainHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* jx, 
						 double* jz)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->Jz = new double [this->NbrSpin - 1];
  this->HalfJx = new double [this->NbrSpin - 1];
  for (int i = 0; i < (this->NbrSpin - 1); i++)
    {
      this->Jz[i] = jz[i];
      this->HalfJx[i] = jx[i] * 0.5;
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

XXZSpinChainHamiltonian::~XXZSpinChainHamiltonian() 
{
  delete[] this->Jz;
  delete[] this->HalfJx;
  delete[] this->SzSzContributions;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void XXZSpinChainHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* XXZSpinChainHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

XXZSpinChainHamiltonian& XXZSpinChainHamiltonian::SetChain(AbstractSpinChain* chain)
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

int XXZSpinChainHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void XXZSpinChainHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex XXZSpinChainHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  double coef;
  int MaxPos = this->NbrSpin - 1;
  int pos;
  int dim = this->Chain->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {

      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    x += V1[pos] * this->HalfJx[j] * coef * V2[i];
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    x += V1[pos] * this->HalfJx[j] * coef * V2[i];
	}

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

Complex XXZSpinChainHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& XXZSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  for (int i = 0; i < dim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = 0; i < dim; i++)
    {

      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJx[j] * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJx[j] * coef * vSource[i];
	    }
	}
    }
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

RealVector& XXZSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
					      int firstComponent, int nbrComponent)
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  int Lim = firstComponent + nbrComponent;
  for (int i = firstComponent; i < Lim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = firstComponent; i < Lim; i++)
    {

      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJx[j] * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJx[j] * coef * vSource[i];
	    }
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

ComplexVector& XXZSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
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

ComplexVector& XXZSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						   int firstComponent, int nbrComponent)
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> XXZSpinChainHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  int Dim = this->Chain->GetHilbertSpaceDimension();
  RealSymmetricMatrix* Sx = new RealSymmetricMatrix (Dim, true);
  RealAntisymmetricMatrix* Sy = new RealAntisymmetricMatrix (Dim, true);
  RealSymmetricMatrix* Sz = new RealSymmetricMatrix (Dim, true);
  this->Chain->Sxi(0, *Sx);
  this->Chain->Syi(0, *Sy);
  this->Chain->Szi(0, *Sz);
  TmpList += Sx;
  TmpList += Sy;
  TmpList += Sz;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> XXZSpinChainHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  int Dim = this->Chain->GetHilbertSpaceDimension();
  RealSymmetricMatrix* Sx = new RealSymmetricMatrix (Dim, true);
  RealAntisymmetricMatrix* Sy = new RealAntisymmetricMatrix (Dim, true);
  RealSymmetricMatrix* Sz = new RealSymmetricMatrix (Dim, true);
  this->Chain->Sxi(this->NbrSpin - 1, *Sx);
  this->Chain->Syi(this->NbrSpin - 1, *Sy);
  this->Chain->Szi(this->NbrSpin - 1, *Sz);
  TmpList += Sx;
  TmpList += Sy;
  TmpList += Sz;
  return TmpList;
}

// evaluate diagonal matrix elements
// 

void XXZSpinChainHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < (this->NbrSpin - 1); j++)
	{
	  this->SzSzContributions[i] += this->Jz[j] * this->Chain->SziSzj(j, j + 1, i);
	}
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, XXZSpinChainHamiltonian& H) 
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

MathematicaOutput& operator << (MathematicaOutput& Str, XXZSpinChainHamiltonian& H) 
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

