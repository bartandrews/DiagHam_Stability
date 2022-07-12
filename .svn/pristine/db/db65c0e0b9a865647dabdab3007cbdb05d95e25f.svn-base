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


#include "Hamiltonian/SpinChainAKLTHamiltonianWithTranslations.h"
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

SpinChainAKLTHamiltonianWithTranslations::SpinChainAKLTHamiltonianWithTranslations()
{
}

// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// squareFactor = numerical factor in front of the 1/3 (S_i S_i+1)^2 term

SpinChainAKLTHamiltonianWithTranslations::SpinChainAKLTHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double squareFactor)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->J = 1.0;
  this->HalfJ = this->J * 0.5;
  this->Jz = this->J;
  this->LinearFactor = 1.0;
  this->SquareFactor = squareFactor / 3.0;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// constructor from default datas
//
// chain = pointer to Hilbert space of the associated system
// nbrSpin = number of spin
// linearFactor = numerical factor in front of the (S_i S_i+1) term
// quadraticFactor = numerical factor in front of the (S_i S_i+1)^2 term

SpinChainAKLTHamiltonianWithTranslations::SpinChainAKLTHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double linearFactor, double quadraticFactor)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->J = 1.0;
  this->HalfJ = this->J * 0.5;
  this->Jz = this->J;
  this->LinearFactor = linearFactor;
  this->SquareFactor = quadraticFactor;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// destructor
//

SpinChainAKLTHamiltonianWithTranslations::~SpinChainAKLTHamiltonianWithTranslations() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainAKLTHamiltonianWithTranslations::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainAKLTHamiltonianWithTranslations::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

SpinChainAKLTHamiltonianWithTranslations& SpinChainAKLTHamiltonianWithTranslations::SetChain(AbstractSpinChainWithTranslations* chain)
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

int SpinChainAKLTHamiltonianWithTranslations::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainAKLTHamiltonianWithTranslations::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs
bool SpinChainAKLTHamiltonianWithTranslations::SavePrecalculation (char* fileName)
{
  return false;
}
  
// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void SpinChainAKLTHamiltonianWithTranslations::EvaluateCosinusTable()
{
  this->CosinusTable = new double [2 * this->NbrSpin];
  this->SinusTable = new double [2 * this->NbrSpin];
  this->ExponentialTable = new Complex [2 * this->NbrSpin];
  double Coef = 2.0 * M_PI * ((double) this->Chain->GetMomentum()) / ((double) this->NbrSpin);
  for (int i = 0; i < (2 * this->NbrSpin); ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = sin(Coef * ((double) i));
      this->ExponentialTable[i] = Phase(Coef * ((double) i));
    }
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& SpinChainAKLTHamiltonianWithTranslations::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									     int firstComponent, int nbrComponent)
{
  double Coef;
  Complex Coef2;
  double Coef3;
  int NbrTranslation;
  int NbrTranslation2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  for (int i = firstComponent; i < Last; i++)
    {
      Complex TmpValue = vSource[i];
      vDestination[i] += this->SzSzContributions[i] * TmpValue;
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      Coef2 = 0.5 * Coef * TmpValue;
	      vDestination[pos] += Coef2 * this->LinearFactor * this->ExponentialTable[NbrTranslation];
	      Coef2 *= this->SquareFactor;
	      vDestination[pos] += Coef2 * this->ExponentialTable[NbrTranslation] * this->Chain->SziSzj(j, j + 1, i);
	    }
	  pos = this->Chain->SmiSpjSmiSpj(j, j + 1, j, j+1, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	    }	  
	  pos = this->Chain->SmiSpjSmiSpj(j + 1, j, j, j+1, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	    }	  
	  pos = this->Chain->SziSzjSmiSpj(j, j + 1, j, j+1, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.5 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	    }	  

	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      Coef2 = 0.5 * Coef * TmpValue;
	      vDestination[pos] += Coef2 * this->LinearFactor * this->ExponentialTable[NbrTranslation];
	      Coef2 *= this->SquareFactor;
	      vDestination[pos] += Coef2 * this->Chain->SziSzj(j, j + 1, i) * this->ExponentialTable[NbrTranslation];
	    }
	  pos = this->Chain->SmiSpjSmiSpj(j, j + 1, j + 1, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	    }	  
	  pos = this->Chain->SmiSpjSmiSpj(j + 1, j, j + 1, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	    }	  
	  pos = this->Chain->SziSzjSmiSpj(j, j + 1, j + 1, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      vDestination[pos] += (0.5 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	    }	  
	}    

      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  Coef2 = 0.5 * Coef * TmpValue;
	  vDestination[pos] += Coef2 * this->LinearFactor * this->ExponentialTable[NbrTranslation];
	  Coef2 *= this->SquareFactor;
	  vDestination[pos] += Coef2 * this->Chain->SziSzj(MaxPos, 0, i) * this->ExponentialTable[NbrTranslation];
	}
      pos = this->Chain->SmiSpjSmiSpj(MaxPos, 0, MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  vDestination[pos] += (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	}	  
      pos = this->Chain->SmiSpjSmiSpj(0, MaxPos, MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  vDestination[pos] += (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	}	  
      pos = this->Chain->SziSzjSmiSpj(0, MaxPos, MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  vDestination[pos] += (0.5 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	}	  

      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  Coef2 = 0.5 * Coef * TmpValue;
	  vDestination[pos] += Coef2 * this->LinearFactor * this->ExponentialTable[NbrTranslation];
	  Coef2 *= this->SquareFactor;
	  vDestination[pos] += Coef2 * this->Chain->SziSzj(MaxPos, 0, i) * this->ExponentialTable[NbrTranslation];
	}
      pos = this->Chain->SmiSpjSmiSpj(MaxPos, 0, 0, MaxPos, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  vDestination[pos] += (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	}	  
      pos = this->Chain->SmiSpjSmiSpj(0, MaxPos, 0, MaxPos, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  vDestination[pos] += (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	}	  
      pos = this->Chain->SziSzjSmiSpj(MaxPos, 0, 0, MaxPos, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  vDestination[pos] += (0.5 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation] * TmpValue;	      
	}	  
    }
  return vDestination;
}
 
// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* SpinChainAKLTHamiltonianWithTranslations::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										  int firstComponent, int nbrComponent)
{
  double Coef;
  double Coef3;
  int NbrTranslation;
  int NbrTranslation2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  Complex* TmpValues = new Complex[nbrVectors];
  Complex TmpCoef;
  for (int i = firstComponent; i < Last; i++)
    {
      for (int k = 0; k < nbrVectors; ++k)
	{
	  TmpValues[k] = vSources[k][i];
	  vDestinations[k][i] += this->SzSzContributions[i] * TmpValues[k];
	}
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      TmpCoef = 0.5 * Coef * this->ExponentialTable[NbrTranslation] * (this->LinearFactor + (this->SquareFactor * this->Chain->SziSzj(j, j + 1, i)));
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];
		}
	    }
	  pos = this->Chain->SmiSpjSmiSpj(j, j + 1, j, j+1, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      TmpCoef = (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
		}
	    }	  
	  pos = this->Chain->SmiSpjSmiSpj(j + 1, j, j, j+1, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      TmpCoef = (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
		}
	    }	  
	  pos = this->Chain->SziSzjSmiSpj(j, j + 1, j, j+1, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      TmpCoef = (0.5 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
		}
	    }	  

	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    { 
	      TmpCoef = 0.5 * Coef * this->ExponentialTable[NbrTranslation] * (this->LinearFactor  + (this->SquareFactor * this->Chain->SziSzj(j, j + 1, i)));
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];
		}
	    }
	  pos = this->Chain->SmiSpjSmiSpj(j, j + 1, j + 1, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      TmpCoef = (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
		}
	    }	  
	  pos = this->Chain->SmiSpjSmiSpj(j + 1, j, j + 1, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      TmpCoef = (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
		}
	    }	  
	  pos = this->Chain->SziSzjSmiSpj(j, j + 1, j + 1, j, i, Coef, NbrTranslation);
	  if (pos != dim)
	    {
	      TmpCoef = (0.5 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
		}
	    }	  
	}    

      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = 0.5 * Coef * this->ExponentialTable[NbrTranslation] * (this->LinearFactor + (this->SquareFactor * this->Chain->SziSzj(MaxPos, 0, i)));
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}
      pos = this->Chain->SmiSpjSmiSpj(MaxPos, 0, MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
	    }
	}	  
      pos = this->Chain->SmiSpjSmiSpj(0, MaxPos, MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
	    }
	}	  
      pos = this->Chain->SziSzjSmiSpj(0, MaxPos, MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = (0.5 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
	    }
	}	  

      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = 0.5 * Coef * this->ExponentialTable[NbrTranslation] * (this->LinearFactor + this->SquareFactor * this->Chain->SziSzj(MaxPos, 0, i));
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];
	    }
	}
      pos = this->Chain->SmiSpjSmiSpj(MaxPos, 0, 0, MaxPos, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
	    }
	}	  
      pos = this->Chain->SmiSpjSmiSpj(0, MaxPos, 0, MaxPos, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = (0.25 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
	    }
	}	  
      pos = this->Chain->SziSzjSmiSpj(MaxPos, 0, 0, MaxPos, i, Coef, NbrTranslation);
      if (pos != dim)
	{
	  TmpCoef = (0.5 * this->SquareFactor * Coef) * this->ExponentialTable[NbrTranslation];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][pos] += TmpCoef * TmpValues[k];	      
	    }
	}	  
    }
  delete[] TmpValues;
  return vDestinations;
}
 
// evaluate diagonal matrix elements
// 

void SpinChainAKLTHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double Tmp;
  int MaxSite = this->NbrSpin - 1;
  for (int i = 0; i < dim; ++i)
    {
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < MaxSite; j++)
	{
	  Tmp = this->Chain->SziSzj(j, j + 1, i);
	  this->SzSzContributions[i] += (this->LinearFactor + (this->SquareFactor * Tmp)) * Tmp;
	}
      Tmp = this->Chain->SziSzj(MaxSite, 0, i);
      this->SzSzContributions[i] += (this->LinearFactor + (this->SquareFactor * Tmp)) * Tmp;
    }
}

