////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of two dimension spin model that could host             //
//                            a Read-Rezayi Z3 phase                          //
//                                                                            //
//                        last modification : 24/05/2018                      //
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


#include "Hamiltonian/TwoDimensionalRRHamiltonian.h"
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

TwoDimensionalRRHamiltonian::TwoDimensionalRRHamiltonian()
{
}


// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// nbrSpinX = number of spin along the x direction
// nbrSpinY = number of spin along the y direction
// j1Factor = amplitude of the Heisenberg coupling between nearest neighbors
// j2Factor = amplitude of the (S_i S_j)^2 nearest neighbor coupling
// j3Factor = amplitude of the (S_i S_j)^3 nearest neighbor coupling
// jcFactor = amplitude of the chiral term
// periodicBoundaryConditions = assume periodic boundary conditions

TwoDimensionalRRHamiltonian::TwoDimensionalRRHamiltonian(AbstractSpinChain* chain, int nbrSpinX, int nbrSpinY, double j1Factor, 
							 double j2Factor, double j3Factor, double jcFactor, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->J1Factor = j1Factor;
  this->JFactor = j1Factor;
  this->JzFactor = j1Factor;
  this->J2Factor = j2Factor;
  this->J3Factor = j3Factor;
  this->JcFactor = jcFactor;
  this->HalfJcFactor = 0.5 * this->JcFactor;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

TwoDimensionalRRHamiltonian::~TwoDimensionalRRHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalRRHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalRRHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalRRHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalRRHamiltonian::ShiftHamiltonian (double shift)
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

RealVector& TwoDimensionalRRHamiltonian::TwoDimensionalRRHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
											  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  if (this->JcFactor != 0.0)
    {
      cout << "error in TwoDimensionalRRHamiltonian::LowLevelMultipleAddMultiply : cannot use real vectors in presence of the chiral term" << endl;
      return vDestination;
    }
  if (this->PeriodicBoundaryConditions == true)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  double& TmpValue = vSource[i];
	  vDestination[i] += this->SzSzContributions[i] * TmpValue;
	  
	  for (int j = 0; j < this->NbrSpinX; ++j)
	    {
	      for (int k = 0; k < this->NbrSpinY; ++k)
		{
		  int TmpIndex1 = this->GetLinearizedIndex(j, k);
		  int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
		  int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
		  
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
		  
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex2, TmpIndex1, i, dim, vDestination, TmpValue);
		  
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex3, i, dim, vDestination, TmpValue);
		  
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);
		}
	    }
	}
    }
  else
    {      
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  double& TmpValue = vSource[i];
	  vDestination[i] += this->SzSzContributions[i] * TmpValue;
	  
	  for (int j = 0; j < (this->NbrSpinX - 1); ++j)
	    {
	      for (int k = 0; k < (this->NbrSpinY - 1); ++k)
		{
		  int TmpIndex1 = this->GetLinearizedIndex(j, k);
		  int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
		  int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
		  
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex2, TmpIndex1, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex3, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);
		}
	    }
	  for (int j = 0; j < (this->NbrSpinX - 1); ++j)
	    {
	      int TmpIndex1 = this->GetLinearizedIndex(j, this->NbrSpinY - 1);
	      int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, this->NbrSpinY - 1);
	      
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex2, TmpIndex1, i, dim, vDestination, TmpValue);
	    }

	  for (int k = 0; k < (this->NbrSpinY - 1); ++k)
	    {
	      int TmpIndex1 = this->GetLinearizedIndex(this->NbrSpinX - 1, k);
	      int TmpIndex3 = this->GetSafeLinearizedIndex(this->NbrSpinX - 1, k + 1);
		  
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex3, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);
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

RealVector* TwoDimensionalRRHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
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
  if (this->JcFactor != 0.0)
    {
      cout << "error in TwoDimensionalRRHamiltonian::LowLevelMultipleAddMultiply : cannot use real vectors in presence of the chiral term" << endl;
      return vDestinations;
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

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& TwoDimensionalRRHamiltonian::TwoDimensionalRRHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int MaxNbrSpinX = this->NbrSpinX;
  if (this->PeriodicBoundaryConditions == true)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Complex& TmpValue = vSource[i];
	  vDestination[i] += this->SzSzContributions[i] * TmpValue;
	  
	  for (int j = 0; j < this->NbrSpinX; ++j)
	    {
	      for (int k = 0; k < this->NbrSpinY; ++k)
		{
		  int TmpIndex1 = this->GetLinearizedIndex(j, k);
		  int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
		  int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
		  int TmpIndex4 = this->GetSafeLinearizedIndex(j + 1, k + 1);
		  
 		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex2, TmpIndex1, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex3, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);

		  this->EvaluateOffDiagonalChiralContribution(TmpIndex1, TmpIndex2, TmpIndex4, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalChiralContribution(TmpIndex2, TmpIndex4, TmpIndex3, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalChiralContribution(TmpIndex4, TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalChiralContribution(TmpIndex3, TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
		}
	    }
	}
    }
  else
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Complex& TmpValue = vSource[i];
	  vDestination[i] += this->SzSzContributions[i] * TmpValue;
	  
	  for (int j = 0; j < (this->NbrSpinX - 1); ++j)
	    {
	      for (int k = 0; k < (this->NbrSpinY - 1); ++k)
		{
		  int TmpIndex1 = this->GetLinearizedIndex(j, k);
		  int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
		  int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
		  int TmpIndex4 = this->GetSafeLinearizedIndex(j + 1, k + 1);
		  
 		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
 		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex2, TmpIndex1, i, dim, vDestination, TmpValue);
 		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex3, i, dim, vDestination, TmpValue);
 		  this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);

		  this->EvaluateOffDiagonalChiralContribution(TmpIndex1, TmpIndex2, TmpIndex4, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalChiralContribution(TmpIndex2, TmpIndex4, TmpIndex3, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalChiralContribution(TmpIndex4, TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);
		  this->EvaluateOffDiagonalChiralContribution(TmpIndex3, TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
		}
	    }

	  for (int j = 0; j < (this->NbrSpinX - 1); ++j)
	    {
	      int TmpIndex1 = this->GetLinearizedIndex(j, this->NbrSpinY - 1);
	      int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, this->NbrSpinY - 1);
	      
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex2, TmpIndex1, i, dim, vDestination, TmpValue);
	    }

	  for (int k = 0; k < (this->NbrSpinY - 1); ++k)
	    {
	      int TmpIndex1 = this->GetLinearizedIndex(this->NbrSpinX - 1, k);
	      int TmpIndex3 = this->GetSafeLinearizedIndex(this->NbrSpinX - 1, k + 1);
		  
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex3, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);
	    }
	}
    }
  return vDestination;
}

// evaluate diagonal matrix elements
// 

void TwoDimensionalRRHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  // SzSz part
  if (this->PeriodicBoundaryConditions == true)
    {
      for (int i = 0; i < dim; i++)
	{
	  // SzSz part
	  this->SzSzContributions[i] = 0.0;
	  double Tmp = 0.0;
	  for (int j = 0; j < this->NbrSpinX; j++)
	    {
	      for (int k = 0; k < this->NbrSpinY; k++)
		{
		  double Tmp2;
		  Tmp2 = this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					     this->GetSafeLinearizedIndex(j, k + 1), i);
		  Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
		  Tmp2 = this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					     this->GetSafeLinearizedIndex(j + 1, k), i);
		  Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
		}
	    }
	  this->SzSzContributions[i] += Tmp;
	}
    }
  else
    {
      for (int i = 0; i < dim; i++)
	{
	  // SzSz part
	  this->SzSzContributions[i] = 0.0;
	  double Tmp = 0.0;
	  for (int j = 0; j < (this->NbrSpinX - 1); j++)
	    {
	      for (int k = 0; k < (this->NbrSpinY - 1); k++)
		{
		  double Tmp2;
		  Tmp2 = this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					     this->GetSafeLinearizedIndex(j, k + 1), i);
		  Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
 		  Tmp2 = this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
 					     this->GetSafeLinearizedIndex(j + 1, k), i);
 		  Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
		}
	    }
	  for (int j = 0; j < (this->NbrSpinX - 1); j++)
	    {
	      double Tmp2;
	      Tmp2 = this->Chain->SziSzj(this->GetLinearizedIndex(j, this->NbrSpinY - 1), 
					 this->GetSafeLinearizedIndex(j + 1, this->NbrSpinY - 1), i);
	      Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
	    }
	  for (int k = 0; k < (this->NbrSpinY - 1); k++)
	    {
	      double Tmp2;
	      Tmp2 = this->Chain->SziSzj(this->GetLinearizedIndex(this->NbrSpinX - 1, k), 
					 this->GetSafeLinearizedIndex(this->NbrSpinX - 1, k + 1), i);
	      Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
	    }
	  this->SzSzContributions[i] += Tmp;
	}
    }
}

