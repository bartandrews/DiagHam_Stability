////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of AKLT hamiltonian                         //
//                                                                            //
//                        last modification : 15/03/2001                      //
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


#include "Hamiltonian/AKLTHamiltonian.h"
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


// contructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin 1
// j = coupling constant between spin 1
// jp = coupling between ending spin 1/2 and spin 1

AKLTHamiltonian::AKLTHamiltonian(const Spin1AKLTChain& chain, int nbrSpin, double j, double jp)
{
  this->Chain = chain;
  this->J = j;
  this->HalfJ = 0.5 * this->J;
  this->Jp = jp;
  this->HalfJp = 0.5 * this->Jp;
  this->NbrSpin = nbrSpin;
  this->SzSzContributions = new double [this->Chain.GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

AKLTHamiltonian::~AKLTHamiltonian() 
{
  delete[] this->SzSzContributions;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

AKLTHamiltonian& AKLTHamiltonian::SetChain(const Spin1AKLTChain& chain)
{  
  delete[] this->SzSzContributions;
  this->Chain = chain;
  this->SzSzContributions = new double [this->Chain.GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  return *this;  
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AKLTHamiltonian::GetHilbertSpace ()
{
  return &(this->Chain);
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension
int AKLTHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain.GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AKLTHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain.GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AKLTHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  double coef;
  int pos;
  int dim = this->Chain.GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {

      // J part of Hamiltonian      
      for (int j = 1; j <= this->NbrSpin; j++)
	{
	  pos = this->Chain.SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    x+= V1[pos] * this->HalfJ * coef * V2[i];
	  pos = this->Chain.SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    x+= V1[pos] * this->HalfJ * coef * V2[i];
	}

      // Jp part of Hamiltonian      
      pos = this->Chain.SmiSpj(0, 1, i, coef);
      if (pos != dim)
	x+= V1[pos] * this->HalfJp * coef * V2[i];
      pos = this->Chain.SmiSpj(1, 0, i, coef);
      if (pos != dim)
	x+= V1[pos] * this->HalfJp * coef * V2[i];

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

Complex AKLTHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AKLTHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  int dim = this->Chain.GetHilbertSpaceDimension();
  double coef;
  int pos;
  for (int i = 0; i < dim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = 0; i < dim; i++)
    {

      // J part of Hamiltonian      
      for (int j = 1; j <= this->NbrSpin; j++)
	{
	  pos = this->Chain.SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    vDestination[pos] += this->HalfJ * coef * vSource[i];
	  pos = this->Chain.SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    vDestination[pos] += this->HalfJ * coef * vSource[i];
	}

      // Jp part of Hamiltonian      
      pos = this->Chain.SmiSpj(0, 1, i, coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJp * coef * vSource[i];
      pos = this->Chain.SmiSpj(1, 0, i, coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJp * coef * vSource[i];

    }
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AKLTHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AKLTHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  int Dim = this->Chain.GetHilbertSpaceDimension();
  RealSymmetricMatrix* Sx = new RealSymmetricMatrix (Dim, true);
  RealAntisymmetricMatrix* Sy = new RealAntisymmetricMatrix (Dim, true);
  RealSymmetricMatrix* Sz = new RealSymmetricMatrix (Dim, true);
  this->Chain.Sxi(this->NbrSpin, *Sx);
  this->Chain.Syi(this->NbrSpin, *Sy);
  this->Chain.Szi(this->NbrSpin, *Sz);
  TmpList += Sx;
  TmpList += Sy;
  TmpList += Sz;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AKLTHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  int Dim = this->Chain.GetHilbertSpaceDimension();
  RealSymmetricMatrix* Sx = new RealSymmetricMatrix (Dim, true);
  RealAntisymmetricMatrix* Sy = new RealAntisymmetricMatrix (Dim, true);
  RealSymmetricMatrix* Sz = new RealSymmetricMatrix (Dim, true);
  this->Chain.Sxi(this->NbrSpin, *Sx);
  this->Chain.Syi(this->NbrSpin, *Sy);
  this->Chain.Szi(this->NbrSpin, *Sz);
  TmpList += Sx;
  TmpList += Sy;
  TmpList += Sz;
  return TmpList;
}

// evaluate diagonal matrix elements
// 

void AKLTHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain.GetHilbertSpaceDimension();

  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      SzSzContributions[i] = this->Jp * this->Chain.SziSzj(0, 1, i);
      for (int j = 1; j < this->NbrSpin; j++)
	{
	  SzSzContributions[i] += this->J * this->Chain.SziSzj(j, j + 1, i);
	}
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, AKLTHamiltonian& H) 
{
  RealVector TmpV2 (H.Chain.GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Chain.GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain.GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Chain.GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Chain.GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Chain.GetHilbertSpaceDimension(); j++)
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

MathematicaOutput& operator << (MathematicaOutput& Str, AKLTHamiltonian& H) 
{
  RealVector TmpV2 (H.Chain.GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Chain.GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain.GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Chain.GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Chain.GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Chain.GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Chain.GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Chain.GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Chain.GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Chain.GetHilbertSpaceDimension() - 1][H.Chain.GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

