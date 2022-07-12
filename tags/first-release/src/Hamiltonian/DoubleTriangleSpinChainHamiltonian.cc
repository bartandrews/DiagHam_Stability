////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of double triangle spin chain  hamiltonian             //
//                                                                            //
//                        last modification : 13/12/2001                      //
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


#include "Hamiltonian/DoubleTriangleSpinChainHamiltonian.h"
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
// nbrTriangle = number of triangles
// j1 = nearest neighbour coupling constant
// j2 = second nearest neighbour coupling constant

DoubleTriangleSpinChainHamiltonian::DoubleTriangleSpinChainHamiltonian(AbstractSpinChain* chain, int nbrTriangle, double j1, double j2)
{
  this->Chain = chain;
  this->NbrTriangle = nbrTriangle;
  this->NbrSpin = 2 * this->NbrTriangle + 1;
  this->J1 = j1;
  this->HalfJ1 = this->J1 * 0.5;;
  this->J2 = j2;
  this->HalfJ2 = this->J2 * 0.5;;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

DoubleTriangleSpinChainHamiltonian::~DoubleTriangleSpinChainHamiltonian() 
{
  delete[] this->SzSzContributions;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void DoubleTriangleSpinChainHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* DoubleTriangleSpinChainHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

DoubleTriangleSpinChainHamiltonian& DoubleTriangleSpinChainHamiltonian::SetChain(AbstractSpinChain* chain)
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

int DoubleTriangleSpinChainHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void DoubleTriangleSpinChainHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex DoubleTriangleSpinChainHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  double coef;
  int ReducedNbrSpin = this->NbrSpin - 3;
  int SpinPos = 0;
  int pos;
  int dim = this->Chain->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
      SpinPos = 0;
      // J part of Hamiltonian      
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, coef);
      if (pos != dim)
	{
	  V1[pos] += (this->HalfJ2 - this->HalfJ1) * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, coef);
      if (pos != dim)
	{
	  V1[pos] += (this->HalfJ2 - this->HalfJ1) * coef * V2[i];
	}
      for (; SpinPos < ReducedNbrSpin; SpinPos += 2)
	{
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, coef);
	  if (pos != dim)
	    {
	      V1[pos] += this->HalfJ1 * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, coef);
	  if (pos != dim)
	    {
	      V1[pos] += this->HalfJ1 * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 2, i, coef);
	  if (pos != dim)
	    {
	      V1[pos] += this->HalfJ2 * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos, i, coef);
	  if (pos != dim)
	    {
	      V1[pos] += this->HalfJ2 * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 3, i, coef);
	  if (pos != dim)
	    {
	      V1[pos] += this->HalfJ2 * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 3, SpinPos + 1, i, coef);
	  if (pos != dim)
	    {
	      V1[pos] += this->HalfJ2 * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 2, i, coef);
	  if (pos != dim)
	    {
	      V1[pos] += this->HalfJ1 * coef * V2[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos + 1, i, coef);
	  if (pos != dim)
	    {
	      V1[pos] += this->HalfJ1 * coef * V2[i];
	    }
	}
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, coef);
      if (pos != dim)
	{
	  V1[pos] += this->HalfJ1 * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, coef);
      if (pos != dim)
	{
	  V1[pos] += this->HalfJ1 * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 2, i, coef);
      if (pos != dim)
	{
	  V1[pos] += this->HalfJ2 * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos, i, coef);
      if (pos != dim)
	{
	  V1[pos] += this->HalfJ2 * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 2, i, coef);
      if (pos != dim)
	{
	  V1[pos] += this->HalfJ2 * coef * V2[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos + 1, i, coef);
      if (pos != dim)
	{
	  V1[pos] += this->HalfJ2 * coef * V2[i];
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

Complex DoubleTriangleSpinChainHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& DoubleTriangleSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int SpinPos = 0;
  int ReducedNbrSpin = this->NbrSpin - 3;
  for (int i = 0; i < dim; i++)
    {
      vDestination[i] = this->SzSzContributions[i] * vSource[i];
    }
  for (int i = 0; i < dim; i++)
    {
      SpinPos = 0;
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += (this->HalfJ2 - this->HalfJ1) * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += (this->HalfJ2 - this->HalfJ1) * coef * vSource[i];
	}
      for (; SpinPos < ReducedNbrSpin; SpinPos += 2)
	{
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 3, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 3, SpinPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	    }
	}
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 2, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 2, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ2 * coef * vSource[i];
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

RealVector& DoubleTriangleSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						  int firstComponent, int nbrComponent) 
{
  int dim = firstComponent + nbrComponent;
  double coef;
  int pos;
  int SpinPos = 0;
  int ReducedNbrSpin = this->NbrSpin - 3;
  for (int i = firstComponent; i < dim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = firstComponent; i < dim; i++)
    {
      SpinPos = 0;
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += (this->HalfJ2 - this->HalfJ1) * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += (this->HalfJ2 - this->HalfJ1) * coef * vSource[i];
	}
      for (; SpinPos < ReducedNbrSpin; SpinPos += 2)
	{
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 3, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 3, SpinPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	    }
	}
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ1 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 2, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 2, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ2 * coef * vSource[i];
	}
      pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += this->HalfJ2 * coef * vSource[i];
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

ComplexVector& DoubleTriangleSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> DoubleTriangleSpinChainHamiltonian::LeftInteractionOperators()
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

List<Matrix*> DoubleTriangleSpinChainHamiltonian::RightInteractionOperators()
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

void DoubleTriangleSpinChainHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  // SzSz part
  int SpinPos = 0;
  int ReducedNbrSpin = this->NbrSpin - 3;
  for (int i = 0; i < dim; i++)
    {
       SpinPos = 0;
     // SzSz part
      this->SzSzContributions[i] = 0.0;
      this->SzSzContributions[i] += (this->J2 - this->J1) * this->Chain->SziSzj(SpinPos, SpinPos + 1, i);
      for (; SpinPos < ReducedNbrSpin; SpinPos += 2)
	{
	  this->SzSzContributions[i] += this->J1 * (this->Chain->SziSzj(SpinPos, SpinPos + 1, i) + 
						    this->Chain->SziSzj(SpinPos + 1, SpinPos + 2, i));
	  this->SzSzContributions[i] += this->J2 * (this->Chain->SziSzj(SpinPos, SpinPos + 2, i) + 
						    this->Chain->SziSzj(SpinPos + 1, SpinPos + 3, i));
	}
      this->SzSzContributions[i] += this->J1 * this->Chain->SziSzj(SpinPos, SpinPos + 1, i) + 
	this->J2 * this->Chain->SziSzj(SpinPos + 1, SpinPos + 2, i);
      this->SzSzContributions[i] += this->J2 * (this->Chain->SziSzj(SpinPos, SpinPos + 2, i));      
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, DoubleTriangleSpinChainHamiltonian& H) 
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

MathematicaOutput& operator << (MathematicaOutput& Str, DoubleTriangleSpinChainHamiltonian& H) 
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

