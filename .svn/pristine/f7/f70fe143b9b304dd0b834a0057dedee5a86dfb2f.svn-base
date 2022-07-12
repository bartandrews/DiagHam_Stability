////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of  hamiltonian associated to trapped bosons            //
//                                                                            //
//                        last modification : 11/06/2002                      //
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


#include "Hamiltonian/TrappedBosonHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>


using std::endl;
using std::ostream;

#define FACTORIAL_MAX 100


// constructor from default datas
//
// bosons = Hilbert space associated to the system
// lzmax = maximum Lz value reached by a boson in the state

TrappedBosonHamiltonian::TrappedBosonHamiltonian(TrappedBosons* bosons, int lzmax)
{
  this->Bosons = bosons;
  this->LzMax = lzmax;
  this->EvaluateInteractionFactors();
}

// destructor
//

TrappedBosonHamiltonian::~TrappedBosonHamiltonian() 
{
  for (int i = 0; i <= this->LzMax; ++i)
    delete[] this->InteractionFactors[i];    
  delete[] this->InteractionFactors;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TrappedBosonHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  for (int i = 0; i <= this->LzMax; ++i)
    delete[] this->InteractionFactors[i];    
  delete[] this->InteractionFactors;
  this->Bosons = (TrappedBosons*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TrappedBosonHamiltonian::GetHilbertSpace ()
{
  return this->Bosons;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TrappedBosonHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Bosons->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TrappedBosonHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex TrappedBosonHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Bosons->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex TrappedBosonHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& TrappedBosonHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Bosons->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& TrappedBosonHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						      int firstComponent, int nbrComponent) 
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    vDestination[i] = 0.0;
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& TrappedBosonHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Bosons->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& TrappedBosonHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Bosons->GetHilbertSpaceDimension();
  int m1;
  int m2;
  int n1;
  int SumM;
  int Lim;
  int Index;
  double Coefficient;
  double* TmpInteractionFactors;
  double Factor;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (m1 = 0; m1 <= this->LzMax; ++m1)
	{
	  TmpInteractionFactors = this->InteractionFactors[m1];
	  for (m2 = 0; m2 < m1; ++m2)
	    {
	      Factor = 2 * TmpInteractionFactors[m2];
	      SumM = m1 + m2;
	      if (SumM >= this->LzMax)
		Lim = this->LzMax;
	      else
		Lim = SumM;
	      for (n1 = 0; n1 <= Lim; ++n1)
		{
		  Index = this->Bosons->AdAdAA(i, m1, m2, n1, SumM - n1, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += Coefficient * Factor * this->InteractionFactors[n1][SumM - n1] * vSource[i];
		    }
		}
	    }
	  Factor = TmpInteractionFactors[m1];
	  SumM = 2 * m1;
	  if (SumM >= this->LzMax)
	    Lim = this->LzMax;
	  else
	    Lim = SumM;
	  for (n1 = 0; n1 <= Lim; ++n1)
	    {
	      Index = this->Bosons->AdAdAA(i, m1, m1, n1, SumM - n1, Coefficient);
	      if (Index < Dim)
		{
		  vDestination[Index] += Coefficient * Factor * this->InteractionFactors[n1][SumM - n1] * vSource[i];
		}
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

ComplexVector& TrappedBosonHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
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

ComplexVector& TrappedBosonHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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

ComplexVector& TrappedBosonHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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
ComplexVector& TrappedBosonHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int firstComponent, int nbrComponent)
{
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> TrappedBosonHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> TrappedBosonHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void TrappedBosonHamiltonian::EvaluateInteractionFactors()
{
  this->InteractionFactors = new double* [this->LzMax + 1];
  FactorialCoefficient Coef;
  for (int i = 0; i <= this->LzMax; ++i)
    {	
      this->InteractionFactors[i] = new double [this->LzMax + 1];
      Coef.SetToOne();
      Coef.Power2Divide(i);
      this->InteractionFactors[i][0] = sqrt(Coef.GetNumericalValue());
      for (int j = 1; j <= this->LzMax; ++j)
	{      
	  Coef.SetToOne();
	  Coef.PartialFactorialMultiply(i + 1, i + j);
	  Coef.FactorialDivide(j);
	  Coef.Power2Divide(i+j);
	  this->InteractionFactors[i][j] = sqrt(Coef.GetNumericalValue());
	}
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, TrappedBosonHamiltonian& H) 
{
  RealVector TmpV2 (H.Bosons->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Bosons->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Bosons->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Bosons->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Bosons->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Bosons->GetHilbertSpaceDimension(); j++)
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

MathematicaOutput& operator << (MathematicaOutput& Str, TrappedBosonHamiltonian& H) 
{
  RealVector TmpV2 (H.Bosons->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Bosons->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Bosons->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Bosons->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Bosons->GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Bosons->GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Bosons->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Bosons->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Bosons->GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Bosons->GetHilbertSpaceDimension() - 1][H.Bosons->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

