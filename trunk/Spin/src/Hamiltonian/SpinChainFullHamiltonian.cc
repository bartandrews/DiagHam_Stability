////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of spin chain hamiltonian with fully anisotropic          //
//                     coupling constants and magnetic field                  //
//                                                                            //
//                        last modification : 15/08/2015                      //
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


#include "Hamiltonian/SpinChainFullHamiltonian.h"
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


// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// nbrSpin = number of spin
// j = array containing coupling constants between spins
// hxFactor = array containing the amplitude of the Zeeman term along x
// hyFactor = array containing the amplitude of the Zeeman term along y
// hzFactor = array containing the amplitude of the Zeeman term along z
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChainFullHamiltonian::SpinChainFullHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* jx, double* jy, double* jz, 
						   double* hx, double* hy, double* hz, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  if (this->PeriodicBoundaryConditions == false)
    {
      this->JPlusFactor = new double [this->NbrSpin - 1];
      this->JMinusFactor = new double [this->NbrSpin - 1];
      this->JzFactor = new double [this->NbrSpin - 1];
      for (int i = 0; i < (this->NbrSpin - 1); i++)
	{
	  this->JPlusFactor[i] = 0.25 * (jx[i] + jy[i]);
	  this->JMinusFactor[i] = 0.25 * (jx[i] - jy[i]);
	  this->JzFactor[i] = jz[i];
	}
    }
  else
    {
      this->JPlusFactor = new double [this->NbrSpin];
      this->JMinusFactor = new double [this->NbrSpin];
      this->JzFactor = new double [this->NbrSpin];
      for (int i = 0; i < this->NbrSpin; i++)
	{
	  this->JPlusFactor[i] = 0.25 * (jx[i] + jy[i]);
	  this->JMinusFactor[i] = 0.25 * (jx[i] - jy[i]);
	  this->JzFactor[i] = jz[i];
	}
   }
  this->HPlusFactor = new Complex [this->NbrSpin];
  this->HMinusFactor = new Complex [this->NbrSpin];
  this->HzFactor = new double [this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; i++)
    {
      this->HPlusFactor[i] = 0.5 * Complex(hx[i], hy[i]);
      this->HMinusFactor[i] = 0.5 * Complex(hx[i], -hy[i]);
      this->HzFactor[i] = hz[i];
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

SpinChainFullHamiltonian::~SpinChainFullHamiltonian() 
{
  delete[] this->JPlusFactor;
  delete[] this->JMinusFactor;
  delete[] this->JzFactor;
  delete[] this->SzSzContributions;
  if (this->HzFactor != 0)
    {
      delete[] this->HPlusFactor;  
      delete[] this->HMinusFactor;  
      delete[] this->HzFactor;  
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainFullHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainFullHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SpinChainFullHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainFullHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& SpinChainFullHamiltonian::SpinChainFullHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      Complex& TmpValue = vSource[i];
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->JPlusFactor[j] * TmpValue;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->JPlusFactor[j] * TmpValue;
	    }
	  pos = this->Chain->SpiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->JMinusFactor[j] * TmpValue;
	    }
	  pos = this->Chain->SmiSmj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->JMinusFactor[j] * TmpValue;
	    }
	  pos = this->Chain->Spi(j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->HPlusFactor[j] * TmpValue;
	    }
	  pos = this->Chain->Smi(j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->HMinusFactor[j] * TmpValue;
	    }
	}
      if (this->PeriodicBoundaryConditions == true)
	{
	  pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->JPlusFactor[MaxPos] * TmpValue;
	    }
	  pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->JPlusFactor[MaxPos] * TmpValue;
	    }
	  pos = this->Chain->SpiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->JMinusFactor[MaxPos] * TmpValue;
	    }
	  pos = this->Chain->SmiSmj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->JMinusFactor[MaxPos] * TmpValue;
	    }
	}
      pos = this->Chain->Spi(MaxPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += coef * this->HPlusFactor[MaxPos] * TmpValue;
	}
      pos = this->Chain->Smi(MaxPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += coef * this->HMinusFactor[MaxPos] * TmpValue;
	}
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* SpinChainFullHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  Complex* TmpValues = new Complex[nbrVectors];
  for (int k = 0; k < nbrVectors; ++k)
    {
      ComplexVector& TmpSource = vSources[k];
      ComplexVector& TmpDestination = vDestinations[k];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpDestination[i] += this->SzSzContributions[i] * TmpSource[i];
	}
    }
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int k = 0; k < nbrVectors; ++k)
	{
	  TmpValues[k] = vSources[k][i];
	}
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->JPlusFactor[j] * TmpValues[k];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->JPlusFactor[j] * TmpValues[k];
	    }
	  pos = this->Chain->SpiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->JMinusFactor[j] * TmpValues[k];
	    }
	  pos = this->Chain->SmiSmj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->JMinusFactor[j] * TmpValues[k];
	    }
	  pos = this->Chain->Spi(j, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->HPlusFactor[j] * TmpValues[k];
	    }
	  pos = this->Chain->Smi(j, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->HMinusFactor[j] * TmpValues[k];
	    }
	}
      if (this->PeriodicBoundaryConditions == true)
	{
	  pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->JPlusFactor[MaxPos] * TmpValues[k];
	    }
	  pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->JPlusFactor[MaxPos] * TmpValues[k];
	    }
	  pos = this->Chain->SpiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->JMinusFactor[MaxPos] * TmpValues[k];
	    }
	  pos = this->Chain->SmiSmj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->JMinusFactor[MaxPos] * TmpValues[k];
	    }
	}
      pos = this->Chain->Spi(MaxPos, i, coef);
      if (pos != dim)
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][pos] += coef * this->JPlusFactor[MaxPos] * TmpValues[k];
	}
      pos = this->Chain->Smi(MaxPos, i, coef);
      if (pos != dim)
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][pos] += coef * this->JMinusFactor[MaxPos] * TmpValues[k];
	}
    }
  delete[] TmpValues;
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void SpinChainFullHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  int MaxSite = this->NbrSpin - 1;
  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < MaxSite; j++)
	{
	  this->SzSzContributions[i] += this->JzFactor[j] * this->Chain->SziSzj(j, j + 1, i);
	}
      if (this->PeriodicBoundaryConditions == true)
	{
	  this->SzSzContributions[i] += this->JzFactor[MaxSite] * this->Chain->SziSzj(MaxSite, 0, i);
	}
    }
  if (this->HzFactor != 0)
    {
      double Coefficient;
      // Sz part
      for (int i = 0; i < dim; i++)
	{
	  for (int j = 0; j < this->NbrSpin; j++)
	    {
	      this->Chain->Szi(j, i, Coefficient);
	      this->SzSzContributions[i] += this->HzFactor[j] * Coefficient;
	    }
	}      
    }
}

