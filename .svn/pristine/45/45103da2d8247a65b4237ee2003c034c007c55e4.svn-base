////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of two dimension Heisenberg model                  //
//                                                                            //
//                        last modification : 13/05/2018                      //
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


#include "Hamiltonian/TwoDimensionalHeisenbergHamiltonian.h"
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

TwoDimensionalHeisenbergHamiltonian::TwoDimensionalHeisenbergHamiltonian()
{
}


// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// nbrSpinX = number of spin along the x direction
// nbrSpinY = number of spin along the y direction
// jFactor = Heisenberg XX coupling constant between nearest neighbors
// jzFactor = Heisenberg Z coupling constant between nearest neighbors
// openBoundaryConditions = flag indicating that open boundary conditions should be used 
TwoDimensionalHeisenbergHamiltonian::TwoDimensionalHeisenbergHamiltonian(AbstractSpinChain* chain, int nbrSpinX, int nbrSpinY, double jFactor, double jzFactor, bool openBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  if (openBoundaryConditions)
    {
      this->NbrLinksX = nbrSpinX - 1;
      this->NbrLinksY = nbrSpinY - 1;
    }
  else
    {
      this->NbrLinksX = nbrSpinX;
      this->NbrLinksY = nbrSpinY;
    }

  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->JFactor = jFactor;
  this->JzFactor = jzFactor;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

TwoDimensionalHeisenbergHamiltonian::~TwoDimensionalHeisenbergHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalHeisenbergHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalHeisenbergHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalHeisenbergHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalHeisenbergHamiltonian::ShiftHamiltonian (double shift)
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

RealVector& TwoDimensionalHeisenbergHamiltonian::TwoDimensionalHeisenbergHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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

      for (int j = 0; j < this->NbrLinksX; ++j)
	{
	  for (int k = 0; k < this->NbrLinksY; ++k)
	    {
	      int TmpIndex1 = this->GetLinearizedIndex(j, k);
	      int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
	      int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef);
 	      if (pos != dim)
		{
		  vDestination[pos] += 0.5 * coef * this->JFactor * TmpValue;
		}
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, coef);
 	      if (pos != dim)
		{
		  vDestination[pos] += 0.5 * coef * this->JFactor * TmpValue;
		}
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, coef);
 	      if (pos != dim)
		{
		  vDestination[pos] += 0.5 * coef * this->JFactor * TmpValue;
		}
	      pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex1, i, coef);
 	      if (pos != dim)
		{
		  vDestination[pos] += 0.5 * coef * this->JFactor * TmpValue;
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

RealVector* TwoDimensionalHeisenbergHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
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
//       for (int j = 0; j < this->NbrSpinX; j++)
// 	{
// 	  for (int k = 0; k < this->NbrSpinY; k++)
// 	    {
// 	      pos = this->Chain->Spi(this->GetLinearizedIndex(j, k), i, coef);
// 	      if (pos != dim)
// 		{
// 		  for (int l = 0; l < nbrVectors; ++l)
// 		    vDestinations[l][pos] += 0.5 * coef * this->HxFactors[j][k] * TmpValues[l];
// 		}
// 	      pos = this->Chain->Smi(this->GetLinearizedIndex(j, k), i, coef);
// 	      if (pos != dim)
// 		{
// 		  for (int l = 0; l < nbrVectors; ++l)
// 		    vDestinations[l][pos] += 0.5 * coef * this->HxFactors[j][k] * TmpValues[l];
// 		}
// 	    }
// 	}
    }
  delete[] TmpValues;
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void TwoDimensionalHeisenbergHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      double Tmp = 0.0;
      for (int j = 0; j < this->NbrLinksX; j++)
	{
	  for (int k = 0; k < this->NbrLinksY; k++)
	    {
	      Tmp += this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					 this->GetSafeLinearizedIndex(j, k + 1), i);
	      Tmp += this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					 this->GetSafeLinearizedIndex(j + 1, k), i);
	    }
	}
      Tmp *= this->JzFactor; 
      this->SzSzContributions[i] += Tmp;
   }
}

