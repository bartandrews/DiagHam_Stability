////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin chain hamiltonian with translations             //
//                                                                            //
//                        last modification : 04/03/2002                      //
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


#include "Hamiltonian/SpinChainHamiltonianWithTranslations.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// j = coupling constant between spin

SpinChainHamiltonianWithTranslations::SpinChainHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->J = j;
  this->HalfJ = this->J * 0.5;
  this->Jz = this->J;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateOffDiagonalMatrixElements();
}

// destructor
//

SpinChainHamiltonianWithTranslations::~SpinChainHamiltonianWithTranslations() 
{
  delete[] this->SzSzContributions;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainHamiltonianWithTranslations::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateOffDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainHamiltonianWithTranslations::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

SpinChainHamiltonianWithTranslations& SpinChainHamiltonianWithTranslations::SetChain(AbstractSpinChainWithTranslations* chain)
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

int SpinChainHamiltonianWithTranslations::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainHamiltonianWithTranslations::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex SpinChainHamiltonianWithTranslations::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  double coef;
  int MaxPos = this->NbrSpin - 1;
  int pos;
  int dim = this->Chain->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
/*
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    x+= V1[pos] * this->HalfJ[j] * coef * V2[i];
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    x+= V1[pos] * this->HalfJ[j] * coef * V2[i];
	}
*/
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

Complex SpinChainHamiltonianWithTranslations::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& SpinChainHamiltonianWithTranslations::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  for (int i = 0; i < dim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = 0; i < dim; i++)
    {
/*
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ[j] * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ[j] * coef * vSource[i];
	    }
	}*/
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

RealVector& SpinChainHamiltonianWithTranslations::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							   int firstComponent, int nbrComponent) 
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  for (int i = 0; i < dim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = 0; i < dim; i++)
    {
/*
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ[j] * coef * vSource[i];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->HalfJ[j] * coef * vSource[i];
	    }
	}*/
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& SpinChainHamiltonianWithTranslations::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  for (int i = 0; i < dim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = 0; i < dim; i++)
    {
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

ComplexVector& SpinChainHamiltonianWithTranslations::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						   int firstComponent, int nbrComponent)
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> SpinChainHamiltonianWithTranslations::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  int Dim = this->Chain->GetHilbertSpaceDimension();
  RealSymmetricMatrix* Sx = new RealSymmetricMatrix (Dim, true);
  RealAntisymmetricMatrix* Sy = new RealAntisymmetricMatrix (Dim, true);
  RealSymmetricMatrix* Sz = new RealSymmetricMatrix (Dim, true);
  this->Chain->Sxi(this->NbrSpin - 1, *Sx);
  this->Chain->Syi(this->NbrSpin - 1, *Sy);
  this->Chain->Szi(this->NbrSpin - 1, *Sz);
//  this->Chain->Sxi(0, *Sx);
//  this->Chain->Syi(0, *Sy);
//  this->Chain->Szi(0, *Sz);
  TmpList += Sx;
  TmpList += Sy;
  TmpList += Sz;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> SpinChainHamiltonianWithTranslations::RightInteractionOperators()
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

void SpinChainHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < (this->NbrSpin - 1); j++)
	{
	  this->SzSzContributions[i] += this->Chain->SziSzj(j, j + 1, i);
	}
      this->SzSzContributions[i] += this->Chain->SziSzj(this->NbrSpin - 1, 0, i);
      this->SzSzContributions[i] *= this->Jz;
    }
}

// evaluate off diagonal matrix elements
// 

void SpinChainHamiltonianWithTranslations::EvaluateOffDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double* TmpMatrixElement = new double [dim];
  this->NbrNonZeroMatrixElement = new int [dim];
  int NbrMatrixElement = 0;
  double Coefficient = 0.0;
  int Translation = 0;
  int Index = 0;
  int Index2 =0;
  int Tmp;
  cout << dim << endl;

  // first pass : find all possible coefficients and the number of non zero matrix elements
  for (int i = 0; i < dim; i++)
    {
      this->NbrNonZeroMatrixElement[i] = 0;
      for (int j = 0; j < (this->NbrSpin - 1); j++)
	{
	  Index = this->Chain->SmiSpj(j, j + 1, i, Coefficient, Translation);
	  if (Coefficient != 0.0)
	    {
	      Tmp = 0;
	      while (Tmp < NbrMatrixElement)
		{
		  if (Coefficient == TmpMatrixElement[Tmp])
		    Tmp = NbrMatrixElement + 1;
		  ++Tmp;
		}
	      if (Tmp == NbrMatrixElement)
		{
		  TmpMatrixElement[NbrMatrixElement] = Coefficient;
		  ++NbrMatrixElement;
		}
	      ++this->NbrNonZeroMatrixElement[i];
	    }
	  Index = this->Chain->SmiSpj(j + 1, j, i, Coefficient, Translation);
	  if (Coefficient != 0.0)
	    {
	      Tmp = 0;
	      while (Tmp < NbrMatrixElement)
		{
		  if (Coefficient == TmpMatrixElement[Tmp])
		    Tmp = NbrMatrixElement + 1;
		  ++Tmp;
		}
	      if (Tmp == NbrMatrixElement)
		{
		  TmpMatrixElement[NbrMatrixElement] = Coefficient;
		  ++NbrMatrixElement;
		}
	      ++this->NbrNonZeroMatrixElement[i];
	    }
	}
      Index = this->Chain->SmiSpj(this->NbrSpin - 1, 0, i, Coefficient, Translation);
      if (Coefficient != 0.0)
	{
	  Tmp = 0;
	  while (Tmp < NbrMatrixElement)
	    {
	      if (Coefficient == TmpMatrixElement[Tmp])
		Tmp = NbrMatrixElement + 1;
	      ++Tmp;
	    }
	  if (Tmp == NbrMatrixElement)
	    {
	      TmpMatrixElement[NbrMatrixElement] = Coefficient;
	      ++NbrMatrixElement;
	    }
	  ++this->NbrNonZeroMatrixElement[i];
	}
      Index = this->Chain->SmiSpj(0, this->NbrSpin - 1, i, Coefficient, Translation);
      if (Coefficient != 0.0)
	{
	  Tmp = 0;
	  while (Tmp < NbrMatrixElement)
	    {
	      if (Coefficient == TmpMatrixElement[Tmp])
		Tmp = NbrMatrixElement + 1;
	      ++Tmp;
	    }
	  if (Tmp == NbrMatrixElement)
	    {
	      TmpMatrixElement[NbrMatrixElement] = Coefficient;
	      ++NbrMatrixElement;
	    }
	  ++this->NbrNonZeroMatrixElement[i];
	}
    }

  // create table containing all posible matrix element
  cout << "NbrMatrixElement= " << NbrMatrixElement << endl;
  this->MatrixElementArray = new double [NbrMatrixElement];
  for (int i = 0; i < NbrMatrixElement; ++i)
    {
      this->MatrixElementArray[i] = TmpMatrixElement[i];
      cout << TmpMatrixElement[i] << endl;
    }

  // evaluate all mask and shift needed to access informations about off-diagonal elements
  this->MatrixElementShift = 1;
  while (NbrMatrixElement > (1 << this->MatrixElementShift))
    {
      ++(this->MatrixElementShift);
    }
  this->MomentumShift = 1;
  while (this->NbrSpin > (1 << this->MomentumShift))
    ++(this->MomentumShift);
  Index = 1;
  while (dim > (1 << Index))
    ++(Index);
  this->IndexMask = (0xffffffff) >> (sizeof(int) * 8 - Index);
  this->MatrixElementMask = ((0xffffffff) >> (sizeof(int) * 8 - this->MatrixElementShift)) << Index;
  this->MomentumMask = ((0xffffffff) >> (sizeof(int) * 8 - this->MomentumShift)) << (Index + this->MatrixElementShift);  
  this->MomentumShift = this->MatrixElementShift + Index;
  this->MatrixElementShift = Index ;

  // second pass : store informations about non zero matrix element
  this->NonZeroMatrixElement = new int* [dim];
  for (int i = 0; i < dim; ++i)
    {
      if (this->NbrNonZeroMatrixElement[i] == 0)
	{
	  this->NonZeroMatrixElement[i] = 0;
	}
      else
	{
	  this->NonZeroMatrixElement[i] = new int [this->NbrNonZeroMatrixElement[i]];
	  Index = 0;
	  for (int j = 0; j < (this->NbrSpin - 1); j++)
	    {
	      Index2 = this->Chain->SmiSpj(j, j + 1, i, Coefficient, Translation);
	      if (Coefficient != 0.0)
		{
		  Tmp = 0;
		  while (Coefficient != this->MatrixElementArray[Tmp])
		    ++Tmp;
		  Tmp <<= this->MatrixElementShift;
		  Tmp |= Translation << this->MomentumShift;
		  Tmp |= Index2;
		  this->NonZeroMatrixElement[i][Index] = Tmp;
		  ++Index;
		}
	      Index2 = this->Chain->SmiSpj(j + 1, j, i, Coefficient, Translation);
	      if (Coefficient != 0.0)
		{
		  Tmp = 0;
		  while (Coefficient != this->MatrixElementArray[Tmp])
		    ++Tmp;
		  Tmp <<= this->MatrixElementShift;
		  Tmp |= Translation << this->MomentumShift;
		  Tmp |= Index2;
		  this->NonZeroMatrixElement[i][Index] = Tmp;
		  ++Index;
		}
	    }
	  Index2 = this->Chain->SmiSpj(this->NbrSpin - 1, 0, i, Coefficient, Translation);
	  if (Coefficient != 0.0)
	    {
	      Tmp = 0;
	      while (Coefficient != this->MatrixElementArray[Tmp])
		++Tmp;
	      Tmp <<= this->MatrixElementShift;
	      Tmp |= Translation << this->MomentumShift;
	      Tmp |= Index2;
	      this->NonZeroMatrixElement[i][Index] = Tmp;
	      ++Index;
	    }
	  Index2 = this->Chain->SmiSpj(0, this->NbrSpin - 1, i, Coefficient, Translation);
	  if (Coefficient != 0.0)
	    {
	      Tmp = 0;
	      while (Coefficient != this->MatrixElementArray[Tmp])
		++Tmp;
	      Tmp <<= this->MatrixElementShift;
	      Tmp |= Translation << this->MomentumShift;
	      Tmp |= Index2;
	      this->NonZeroMatrixElement[i][Index] = Tmp;
	      ++Index;
	    }
	}
    }

  for (int i = 0; i < NbrMatrixElement; ++i)
    {
      this->MatrixElementArray[i] *= this->HalfJ;
    }

  for (int i = 0; i < dim; ++i)
    {
      if (this->NbrNonZeroMatrixElement[i] == 0)
	{
	  cout << "state " << i << " not connected" << endl;
	  this->NonZeroMatrixElement[i] = 0;
	}
      else
	{
	  cout << "state " << i << " connected to :" << endl;
	  for (int j = 0; j < this->NbrNonZeroMatrixElement[i]; ++j)
	    {
	      Index = this->NonZeroMatrixElement[i][j];
	      cout << "    state " << (Index & this->IndexMask) << " with coefficient " 
		   << this->MatrixElementArray[((Index & this->MatrixElementMask) >> this->MatrixElementShift)]
		   << " and shift " << ((Index & this->MomentumMask) >> this->MomentumShift) << endl;
	    }
	}
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, SpinChainHamiltonianWithTranslations& H) 
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

MathematicaOutput& operator << (MathematicaOutput& Str, SpinChainHamiltonianWithTranslations& H) 
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

