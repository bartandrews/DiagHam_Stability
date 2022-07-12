////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of periodic XXZ spin chain hamiltonian with second            //
//    nearest neighbor coupling in z direction and alternative coupling to    //
//                              the magnetic field                            //
//                                                                            //
//                        last modification : 04/02/2003                      //
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


#include "Hamiltonian/PeriodicXXZSpinChainWith2ndNeighborCouplingHamiltonian.h"
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
// nbrSpin = number of spin (rounded to the lower even number)
// jz = coupling constant between nearest neighbor Sz components (in jx coupling constant unit)
// jz2 = coupling constant between second nearest neighbor Sz components (in jx coupling constant unit)
// bField = magnetic field times the coupling constant (in jx coupling constant unit)

PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian(AbstractSpinChain* chain, 
																   int nbrSpin, double jz, 
																   double jz2, double bField)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->Jz = jz;
  this->Jz2 = jz2;
  this->BField = bField;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::~PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian() 
{
  delete[] this->SzSzContributions;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian& PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::SetChain(AbstractSpinChain* chain)
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

int PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
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
	    x += 0.5 * V1[pos] * coef * V2[i];
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    x += 0.5 * V1[pos] * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
      if (pos != dim)
	x += 0.5 * V1[pos] * coef * V2[i];
      pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
      if (pos != dim)
	x += 0.5 * V1[pos] * coef * V2[i];
      
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

Complex PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
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
	      vDestination[pos] += 0.5 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += 0.5 * coef * vSource[i];
	    }
	}
      pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += 0.5 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += 0.5 * coef * vSource[i];
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

RealVector& PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
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
	      vDestination[pos] += 0.5 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += 0.5 * coef * vSource[i];
	    }
	}
      pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += 0.5 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += 0.5 * coef * vSource[i];
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

ComplexVector& PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
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

ComplexVector& PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						   int firstComponent, int nbrComponent)
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate diagonal matrix elements
// 

void PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  int MaxIter = (this->NbrSpin / 2) - 1;
  int j = 0;
  int Twicej;
  double Coef;
  // SzSz part
  for (int i = 0; i < dim; ++i)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      for (j = 0; j < MaxIter; ++j)
	{
	  Twicej = 2 * j;
	  this->SzSzContributions[i] += this->Jz * this->Chain->SziSzj(Twicej, Twicej + 1, i);
	  this->SzSzContributions[i] += this->Jz * this->Chain->SziSzj(Twicej + 1, Twicej + 2, i);
	  this->SzSzContributions[i] += this->Jz2 * this->Chain->SziSzj(Twicej, Twicej + 2, i);
	  this->SzSzContributions[i] += this->Jz2 * this->Chain->SziSzj(Twicej + 1, Twicej + 3, i);
	  this->Chain->Szi(Twicej, i, Coef);
	  this->SzSzContributions[i] += this->BField * Coef;
	  this->Chain->Szi(Twicej + 1, i, Coef);
	  this->SzSzContributions[i] -= this->BField * Coef;
	}
      Twicej = 2 * j;
      this->SzSzContributions[i] += this->Jz * this->Chain->SziSzj(Twicej, Twicej + 1, i);
      this->SzSzContributions[i] += this->Jz * this->Chain->SziSzj(Twicej + 1, 0, i);
      this->SzSzContributions[i] += this->Jz2 * this->Chain->SziSzj(Twicej, 0, i);
      this->SzSzContributions[i] += this->Jz2 * this->Chain->SziSzj(Twicej + 1, 1, i);
      this->Chain->Szi(Twicej, i, Coef);
      this->SzSzContributions[i] += this->BField * Coef;
      this->Chain->Szi(Twicej + 1, i, Coef);
      this->SzSzContributions[i] -= this->BField * Coef;
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian& H) 
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

MathematicaOutput& operator << (MathematicaOutput& Str, PeriodicXXZSpinChainWithSecondNearestNeighborCouplingHamiltonian& H) 
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

