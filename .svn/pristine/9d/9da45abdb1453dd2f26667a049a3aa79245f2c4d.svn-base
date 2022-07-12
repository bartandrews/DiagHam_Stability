////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of Haldane-Shasry hamiltonian with translations         //
//                                                                            //
//                        last modification : 26/05/2011                      //
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


#include "Hamiltonian/HaldaneShastryHamiltonianWithTranslations.h"
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


// default constructor
//

HaldaneShastryHamiltonianWithTranslations::HaldaneShastryHamiltonianWithTranslations()
{
}

// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin

HaldaneShastryHamiltonianWithTranslations::HaldaneShastryHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->JCoupling = new double [this->NbrSpin];
  double PiOnN = M_PI / ((double) this->NbrSpin);
  for (int i = 1; i < this->NbrSpin; ++i)
    this->JCoupling[i] = PiOnN * PiOnN / (sin (PiOnN * ((double) i)) * sin (PiOnN * ((double) i)));
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// destructor
//

HaldaneShastryHamiltonianWithTranslations::~HaldaneShastryHamiltonianWithTranslations() 
{
  delete[] this->SzSzContributions;
  delete[] this->CosinusTable;
  delete[] this->SinusTable;
  delete[] this->ExponentialTable;
  delete[] this->JCoupling;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void HaldaneShastryHamiltonianWithTranslations::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  delete[] this->JCoupling;
  this->Chain = (AbstractSpinChainWithTranslations*) hilbertSpace;
  this->JCoupling = new double [this->NbrSpin];
  double PiOnN = M_PI / ((double) this->NbrSpin);
  for (int i = 1; i < this->NbrSpin; ++i)
    this->JCoupling[i] = PiOnN * PiOnN / (sin (PiOnN * ((double) i)) * sin (PiOnN * ((double) i)));
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* HaldaneShastryHamiltonianWithTranslations::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

HaldaneShastryHamiltonianWithTranslations& HaldaneShastryHamiltonianWithTranslations::SetChain(AbstractSpinChainWithTranslations* chain)
{  
  delete[] this->SzSzContributions;
  this->Chain = chain;
  delete[] this->JCoupling;
  this->JCoupling = new double [this->NbrSpin];
  double PiOnN = M_PI / ((double) this->NbrSpin);
  for (int i = 1; i < this->NbrSpin; ++i)
    this->JCoupling[i] = PiOnN * PiOnN / (sin (PiOnN * ((double) i)) * sin (PiOnN * ((double) i)));
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int HaldaneShastryHamiltonianWithTranslations::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void HaldaneShastryHamiltonianWithTranslations::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs
bool HaldaneShastryHamiltonianWithTranslations::SavePrecalculation (char* fileName)
{
  return false;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex HaldaneShastryHamiltonianWithTranslations::MatrixElement (RealVector& V1, RealVector& V2) 
{
  return Complex(0);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex HaldaneShastryHamiltonianWithTranslations::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
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
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      pos = this->Chain->SmiSpj(j, k, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= 0.5 * this->JCoupling[k - j];
		  TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
				    (V2.Im(i) * this->SinusTable[NbrTranslation]));
		  TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
				    (V2.Im(i) * this->CosinusTable[NbrTranslation]));
		  Z += Conj(V1[pos]) * TmpZ;
		}
	      pos = this->Chain->SmiSpj(k, j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= 0.5 * this->JCoupling[k - j];
		  TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
				    (V2.Im(i) * this->SinusTable[NbrTranslation]));
		  TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
				    (V2.Im(i) * this->CosinusTable[NbrTranslation]));
		  Z += Conj(V1[pos]) * TmpZ;
		}
	    }    
	}
    }
  return Z;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& HaldaneShastryHamiltonianWithTranslations::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      pos = this->Chain->SmiSpj(j, k, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= 0.5 * this->JCoupling[k - j];
		  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		}
	      pos = this->Chain->SmiSpj(k, j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= 0.5 * this->JCoupling[k - j];
		  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		}
	    }    
	}
    }
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> HaldaneShastryHamiltonianWithTranslations::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> HaldaneShastryHamiltonianWithTranslations::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void HaldaneShastryHamiltonianWithTranslations::EvaluateCosinusTable()
{
  this->CosinusTable = new double [this->NbrSpin];
  this->SinusTable = new double [this->NbrSpin];
  this->ExponentialTable = new Complex [this->NbrSpin];
  double Coef = 2.0 * M_PI / ((double) this->NbrSpin) * ((double) this->Chain->GetMomentum());
  for (int i = 0; i < this->NbrSpin ; ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = sin(Coef * ((double) i));
      this->ExponentialTable[i].Re = this->CosinusTable[i];
      this->ExponentialTable[i].Im = this->SinusTable[i];
    }
}

// evaluate diagonal matrix elements
// 

void HaldaneShastryHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < (this->NbrSpin - 1); j++)
	{
	  for (int k = j + 1; k < this->NbrSpin; ++k)
	    {
	      this->SzSzContributions[i] += this->JCoupling[k - j] * this->Chain->SziSzj(j, k, i);
	    }
	}
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, HaldaneShastryHamiltonianWithTranslations& H) 
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

MathematicaOutput& operator << (MathematicaOutput& Str, HaldaneShastryHamiltonianWithTranslations& H) 
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

