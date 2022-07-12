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
#include "MathTools/Complex.h"
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
  this->EvaluateCosinusTable();
}

// destructor
//

SpinChainHamiltonianWithTranslations::~SpinChainHamiltonianWithTranslations() 
{
  delete[] this->SzSzContributions;
  delete[] this->CosinusTable;
  delete[] this->SinusTable;
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
  return Complex(0);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex SpinChainHamiltonianWithTranslations::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  Complex Z (0.0, 0.0);
  Complex TmpZ;
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i++)
    {
      Z += this->SzSzContributions[i] * (Conj(V1[i]) * V2[i]);
    }
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
				(V2.Im(i) * this->SinusTable[NbrTranslation]));
	      TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
				(V2.Im(i) * this->CosinusTable[NbrTranslation]));
	      Z += Conj(V1[pos]) * TmpZ;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
				(V2.Im(i) * this->SinusTable[NbrTranslation]));
	      TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
				(V2.Im(i) * this->CosinusTable[NbrTranslation]));
	      Z += Conj(V1[pos]) * TmpZ;
	    }
	}    
      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
			    (V2.Im(i) * this->SinusTable[NbrTranslation]));
	  TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
			    (V2.Im(i) * this->CosinusTable[NbrTranslation]));
	  Z += Conj(V1[pos]) * TmpZ;
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
			    (V2.Im(i) * this->SinusTable[NbrTranslation]));
	  TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
			    (V2.Im(i) * this->CosinusTable[NbrTranslation]));
	  Z += Conj(V1[pos]) * TmpZ;
	}      
    }
  return Z;
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
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  for (int i = firstComponent; i < Last; i++)
    {
      vDestination.Re(i) = this->SzSzContributions[i] * vSource.Re(i);
      vDestination.Im(i) = this->SzSzContributions[i] * vSource.Im(i);
    }
  for (int i = firstComponent; i < Last; i++)
    {
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	      vDestination.Im(pos) +=Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					     (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	    }
	}    
      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
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

ComplexVector& SpinChainHamiltonianWithTranslations::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									 int firstComponent, int nbrComponent)
{
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  for (int i = firstComponent; i < Last; i++)
    {
       vDestination.Re(i) += this->SzSzContributions[i] * vSource.Re(i);
       vDestination.Im(i) += this->SzSzContributions[i] * vSource.Im(i);
     for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	    }
	}    
      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	}      
    }
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> SpinChainHamiltonianWithTranslations::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> SpinChainHamiltonianWithTranslations::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void SpinChainHamiltonianWithTranslations::EvaluateCosinusTable()
{
  this->CosinusTable = new double [this->NbrSpin];
  this->SinusTable = new double [this->NbrSpin];
  double Coef = 2.0 * M_PI / ((double) this->NbrSpin) * ((double) this->Chain->GetMomentum());
  for (int i = 0; i < this->NbrSpin ; ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = sin(Coef * ((double) i));
    }
}

// evaluate diagonal matrix elements
// 

void SpinChainHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < (this->NbrSpin - 1); j++)
	{
	  this->SzSzContributions[i] += this->Chain->SziSzj(j, j + 1, i);
	}
      this->SzSzContributions[i] += this->Chain->SziSzj(this->NbrSpin - 1, 0, i);
      this->SzSzContributions[i] *= this->Jz;
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, SpinChainHamiltonianWithTranslations& H) 
{
  ComplexVector TmpV2 (H.Chain->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [H.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(H.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2.Re(i - 1) = 0.0;
	  TmpV2.Im(i - 1) = 0.0;
	}
      TmpV2.Re(i) = 1.0;
      TmpV2.Im(i) = 0.0;
      ((AbstractHamiltonian*) &H)->LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Chain->GetHilbertSpaceDimension(); j++)
	{
	  Str << "(" << TmpV[j].Re(i) << ", " << TmpV[j].Im(i) << ")    ";
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
  ComplexVector TmpV2 (H.Chain->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [H.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(H.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      ((AbstractHamiltonian*) &H)->LowLevelMultiply (TmpV2, TmpV[i]);
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

