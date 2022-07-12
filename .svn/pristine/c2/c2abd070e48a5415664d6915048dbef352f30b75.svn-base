////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of two dimension Ising hamiltonian with transverse field     //
//                                                                            //
//                        last modification : 17/12/2015                      //
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


#include "Hamiltonian/TwoDimensionalTransverseFieldIsingHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/RealVector.h"
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

TwoDimensionalTransverseFieldIsingHamiltonian::TwoDimensionalTransverseFieldIsingHamiltonian()
{
  this->HxFactors = 0;
  this->HzFactors = 0;
}


// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// nbrSpinX = number of spin along the x direction
// nbrSpinY = number of spin along the y direction
// jFactor = amplitude of the Ising term
// hxFactor = amplitudes of the Zeeman term along x
// hzFactor = amplitudes of the Zeeman term along z
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

TwoDimensionalTransverseFieldIsingHamiltonian::TwoDimensionalTransverseFieldIsingHamiltonian(AbstractSpinChain* chain, int nbrSpinX, int nbrSpinY, double jFactor,
											     double** hxFactor, double** hzFactor, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->JFactor = jFactor;
  this->HxFactors = new double* [this->NbrSpinX];
  this->HzFactors = new double* [this->NbrSpinX];
  for (int i = 0; i < this->NbrSpinX; i++)
    {
      this->HxFactors[i] = new double [this->NbrSpinY];
      this->HzFactors[i] = new double [this->NbrSpinY];     
      for (int j = 0; j < this->NbrSpinY; j++)
	{
	  this->HxFactors[i][j] = 0.5 * hxFactor[i][j];
	  this->HzFactors[i][j] = hzFactor[i][j];
	}
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

TwoDimensionalTransverseFieldIsingHamiltonian::~TwoDimensionalTransverseFieldIsingHamiltonian() 
{
  if (this->HxFactors != 0)
    delete[] this->HxFactors; 
  if (this->HzFactors != 0)
    delete[] this->HzFactors;  
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalTransverseFieldIsingHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalTransverseFieldIsingHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalTransverseFieldIsingHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalTransverseFieldIsingHamiltonian::ShiftHamiltonian (double shift)
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

RealVector& TwoDimensionalTransverseFieldIsingHamiltonian::TwoDimensionalTransverseFieldIsingHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
      double& TmpValue = vSource[i];
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      pos = this->Chain->Spi(this->GetLinearizedIndex(j, k), i, coef);
	      if (pos != dim)
		{
		  vDestination[pos] += coef * this->HxFactors[j][k] * TmpValue;
		}
	      pos = this->Chain->Smi(this->GetLinearizedIndex(j, k), i, coef);
	      if (pos != dim)
		{
		  vDestination[pos] += coef * this->HxFactors[j][k] * TmpValue;
		}
	    }
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

RealVector* TwoDimensionalTransverseFieldIsingHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  double* TmpValues = new double[nbrVectors];
  for (int k = 0; k < nbrVectors; ++k)
    {
      RealVector& TmpSource = vSources[k];
      RealVector& TmpDestination = vDestinations[k];
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
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      pos = this->Chain->Spi(this->GetLinearizedIndex(j, k), i, coef);
	      if (pos != dim)
		{
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += coef * this->HxFactors[j][k] * TmpValues[l];
		}
	      pos = this->Chain->Smi(this->GetLinearizedIndex(j, k), i, coef);
	      if (pos != dim)
		{
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += coef * this->HxFactors[j][k] * TmpValues[l];
		}
	    }
	}
    }
  delete[] TmpValues;
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void TwoDimensionalTransverseFieldIsingHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 1; k < this->NbrSpinY; k++)
	    {
	      this->SzSzContributions[i] += this->Chain->SziSzj(this->GetLinearizedIndex(j, k - 1), 
								this->GetLinearizedIndex(j, k), i);
	    }
	}
      for (int k = 0; k < this->NbrSpinY; k++)
	{
	  for (int j = 1; j < this->NbrSpinX; j++)
	    {
	      this->SzSzContributions[i] += this->Chain->SziSzj(this->GetLinearizedIndex(j - 1, k), 
								this->GetLinearizedIndex(j, k), i);
	    }
	}
      for (int k = 0; k < this->NbrSpinY; k++)
	{
	  this->SzSzContributions[i] += this->Chain->SziSzj(this->GetLinearizedIndex(0, k), 
							    this->GetLinearizedIndex(this->NbrSpinX -1, k), i);
	}
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  this->SzSzContributions[i] += this->Chain->SziSzj(this->GetLinearizedIndex(j, 0), 
							    this->GetLinearizedIndex(j, this->NbrSpinY - 1), i);
	}
    }
  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] *= this->JFactor;
    }
  if (this->HzFactors != 0)
    {
      double Coefficient;
      for (int i = 0; i < dim; i++)
	{
	  for (int j = 0; j < this->NbrSpinX; j++)
	    {
	      for (int k = 0; k < this->NbrSpinY; k++)
		{
		  this->Chain->Szi(this->GetLinearizedIndex(j, k), i, Coefficient);
		  this->SzSzContributions[i] += this->HzFactors[j][k] * Coefficient;
		}
	    }
	}      
    }
}

