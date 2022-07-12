////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                         Class author: Cecile Repellin                      //
//                                                                            //
//                                                                            //
//       class of two dimensional SU(3) spin model on the square lattice      //
//                                                                            //
//                        last modification : 05/02/2018                      //
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


#include "Hamiltonian/TwoDimensionalSquareLatticeSU3Hamiltonian.h"
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
using std::min;
using std::max;



// default constructor
//

TwoDimensionalSquareLatticeSU3Hamiltonian::TwoDimensionalSquareLatticeSU3Hamiltonian()
{
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

TwoDimensionalSquareLatticeSU3Hamiltonian::TwoDimensionalSquareLatticeSU3Hamiltonian (AbstractSpinChain* chain, int nbrSpinX, int nbrSpinY, double jFactor, double jSquareExchangeFactor, bool periodicBoundaryConditionsX, bool periodicBoundaryConditionsY, int offset)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->PeriodicBoundaryConditionsX = periodicBoundaryConditionsX;
  this->PeriodicBoundaryConditionsY = periodicBoundaryConditionsY;
  this->JFactor = jFactor;
  this->JSquareExchangeFactor = jSquareExchangeFactor;
  
  this->Offset = offset;
  
  this->HermitianSymmetryFlag = false;
  
  // diagonal terms
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); ++i)
    this->SzSzContributions[i] = 0.0;
}

// destructor
//

TwoDimensionalSquareLatticeSU3Hamiltonian::~TwoDimensionalSquareLatticeSU3Hamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalSquareLatticeSU3Hamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalSquareLatticeSU3Hamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalSquareLatticeSU3Hamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalSquareLatticeSU3Hamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] = shift;
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& TwoDimensionalSquareLatticeSU3Hamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int pos;
  int TmpIndex;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  
  int MaxPosX = this->NbrSpinX;
  if (!this->PeriodicBoundaryConditionsX)
    MaxPosX = this->NbrSpinX - 1;
  int MaxPosY = this->NbrSpinY;
  if (!this->PeriodicBoundaryConditionsY)
    MaxPosY = this->NbrSpinY - 1;
      
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      
      for (int j = 0; j < MaxPosX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe (j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe (j + 1, k);
	      
	      pos = this->Chain->Pij(TmpIndex, TmpIndex1, i);
	      vDestination[pos] += vSource[i] * this->JFactor;
	    }
	 }
	 
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < MaxPosY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe (j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe (j - this->Offset, k + 1);
	      
	      pos = this->Chain->Pij(TmpIndex, TmpIndex1, i);
	      vDestination[pos] += vSource[i] * this->JFactor;
	    }
	 }
	 
      if (this->JSquareExchangeFactor != 0.0)
      {
// 	cout << MaxPosX << " " << MaxPosY << endl;
	for (int j = 0; j < MaxPosX; j++)
	{
	  for (int k = 0; k < MaxPosY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe (j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe (j + 1, k);
	      TmpIndex2 = this->GetLinearizedIndexSafe (j + 1 - this->Offset, k + 1);
	      TmpIndex3 = this->GetLinearizedIndexSafe (j - this->Offset, k + 1);
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex1, TmpIndex2, TmpIndex3, i);
	      vDestination[pos] += vSource[i] * this->JSquareExchangeFactor;
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex3, TmpIndex2, TmpIndex1, i);
	      vDestination[pos] += vSource[i] * this->JSquareExchangeFactor;
	    }
	 }
	}
      }

//   for (int i = firstComponent; i <LastComponent; ++i)
//     cout << vDestination[i] << " ";
//   cout << endl;
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

RealVector& TwoDimensionalSquareLatticeSU3Hamiltonian::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int pos;
  int TmpIndex;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  double TmpSum;
  
  int MaxPosX = this->NbrSpinX;
  if (!this->PeriodicBoundaryConditionsX)
    MaxPosX = this->NbrSpinX - 1;
  int MaxPosY = this->NbrSpinY;
  if (!this->PeriodicBoundaryConditionsY)
    MaxPosY = this->NbrSpinY - 1;
    
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      TmpSum = 0.0;
      for (int j = 0; j < MaxPosX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j + 1, k);
	      
	      pos = this->Chain->Pij(TmpIndex, TmpIndex1, i);
	      if (pos < i)
	      {
		vDestination[pos] += vSource[i] * this->JFactor;
		TmpSum += vSource[pos] * this->JFactor;
	      }
	      if (pos == i)
		TmpSum += vSource[pos] * this->JFactor;
	    }     
	 }
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < MaxPosY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j - this->Offset, k + 1);
	      
	      pos = this->Chain->Pij(TmpIndex, TmpIndex1, i);
	      if (pos < i)
	      {
		vDestination[pos] += vSource[i] * this->JFactor;
		TmpSum += vSource[pos] * this->JFactor;
	      }
	      if (pos == i)
		TmpSum += vSource[pos] * this->JFactor;
	    }     
	 }
      if (this->JSquareExchangeFactor != 0.0)
      {
	for (int j = 0; j < MaxPosX; j++)
	{
	  for (int k = 0; k < MaxPosY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe (j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe (j + 1, k);
	      TmpIndex2 = this->GetLinearizedIndexSafe (j + 1 - this->Offset, k + 1);
	      TmpIndex3 = this->GetLinearizedIndexSafe (j - this->Offset, k + 1);
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex1, TmpIndex2, TmpIndex3, i);
	      if (pos < dim)
	      	vDestination[pos] += vSource[i] * this->JSquareExchangeFactor;
	      
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex3, TmpIndex2, TmpIndex1, i);
	      if (pos < dim)
		vDestination[pos] += vSource[i] * this->JSquareExchangeFactor;
	    }
	  }
	}
      vDestination[i] += (this->SzSzContributions[i] * vSource[i] + TmpSum);
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

RealVector* TwoDimensionalSquareLatticeSU3Hamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int pos;
  
  int TmpIndex;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  
  int MaxPosX = this->NbrSpinX;
  if (!this->PeriodicBoundaryConditionsX)
    MaxPosX = this->NbrSpinX - 1;
  int MaxPosY = this->NbrSpinY;
  if (!this->PeriodicBoundaryConditionsY)
    MaxPosY = this->NbrSpinY - 1;
  
  double* TmpValues = new double[nbrVectors];
  for (int l = 0; l < nbrVectors; ++l)
    {
      RealVector& TmpSource = vSources[l];
      RealVector& TmpDestination = vDestinations[l];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpDestination[i] += this->SzSzContributions[i] * TmpSource[i];
	}
    }
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int l = 0; l < nbrVectors; ++l)
	{
	  TmpValues[l] = vSources[l][i];
	}
      for (int j = 0; j < MaxPosX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j + 1, k);
	      
	      pos = this->Chain->Pij(TmpIndex, TmpIndex1, i);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][pos] += TmpValues[l] * this->JFactor;

	  }
	}
	  
	for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < MaxPosY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j - this->Offset, k + 1);
	      
	      pos = this->Chain->Pij(TmpIndex, TmpIndex1, i);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][pos] += TmpValues[l] * this->JFactor;

	  }
	}
	
      if (this->JSquareExchangeFactor != 0.0)
      {
	for (int j = 0; j < MaxPosX; j++)
	{
	  for (int k = 0; k < MaxPosY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe (j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe (j + 1, k);
	      TmpIndex2 = this->GetLinearizedIndexSafe (j + 1 - this->Offset, k + 1);
	      TmpIndex3 = this->GetLinearizedIndexSafe (j - this->Offset, k + 1);
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex1, TmpIndex2, TmpIndex3, i);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][pos] += TmpValues[l] * this->JSquareExchangeFactor;
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex3, TmpIndex2, TmpIndex1, i);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][pos] += TmpValues[l] * this->JSquareExchangeFactor;
	    }
	 }	
	}
      }
  delete[] TmpValues;
  return vDestinations;
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

RealVector* TwoDimensionalSquareLatticeSU3Hamiltonian::HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int pos;
  
  int TmpIndex;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  
  int MaxPosX = this->NbrSpinX;
  if (!this->PeriodicBoundaryConditionsX)
    MaxPosX = this->NbrSpinX - 1;
  int MaxPosY = this->NbrSpinY;
  if (!this->PeriodicBoundaryConditionsY)
    MaxPosY = this->NbrSpinY - 1;
  
  double* TmpValues = new double[nbrVectors];
  double* TmpSums = new double[nbrVectors];
  
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int l = 0; l < nbrVectors; ++l)
	{
	  TmpValues[l] = vSources[l][i];
	  TmpSums[l] = 0.0;
	}
      for (int j = 0; j < MaxPosX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j + 1, k);
	      
	      
	      pos = this->Chain->Pij(TmpIndex, TmpIndex1, i);
	      if (pos < i)
	      {
		for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][pos] += TmpValues[l] * this->JFactor;
		  TmpSums[l] += vSources[l][pos] * this->JFactor;
		}
	      }
	      if (pos == i)
		for (int l = 0; l < nbrVectors; ++l)
		  TmpSums[l] += vSources[l][pos] * this->JFactor;
	    }
	}
	
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < MaxPosY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j - this->Offset, k + 1);
	      
	      
	      pos = this->Chain->Pij(TmpIndex, TmpIndex1, i);
	      if (pos < i)
	      {
		for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][pos] += TmpValues[l] * this->JFactor;
		  TmpSums[l] += vSources[l][pos] * this->JFactor;
		}
	      }
	      if (pos == i)
		for (int l = 0; l < nbrVectors; ++l)
		  TmpSums[l] += vSources[l][pos] * this->JFactor;
	    }
	}
	      
      if (this->JSquareExchangeFactor != 0.0)
      {
	for (int j = 0; j < MaxPosX; j++)
	{
	  for (int k = 0; k < MaxPosY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe (j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe (j + 1, k);
	      TmpIndex2 = this->GetLinearizedIndexSafe (j + 1 - this->Offset, k + 1);
	      TmpIndex3 = this->GetLinearizedIndexSafe (j - this->Offset, k + 1);
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex1, TmpIndex2, TmpIndex3, i);
	      if (pos < dim)
	      {
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += TmpValues[l] * this->JSquareExchangeFactor;
	      }
		
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex3, TmpIndex2, TmpIndex1, i);
	      if (pos < dim)
	      {
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += TmpValues[l] * this->JSquareExchangeFactor;
	      }
	    }
	 }
	}
	
	for (int l = 0; l < nbrVectors; ++l)
	  vDestinations[l][i] += (TmpSums[l]);
    }
  delete[] TmpValues;
  delete[] TmpSums;
  return vDestinations;
}


// ask if Hamiltonian implements hermitian symmetry operations
//

bool TwoDimensionalSquareLatticeSU3Hamiltonian::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}
