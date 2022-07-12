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
//                            with 2D translations                            //
//                                                                            //
//                        last modification : 07/02/2018                      //
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

#include "Hamiltonian/TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian.h"
#include "Vector/ComplexVector.h"
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
using std::min;
using std::max;


// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// nbrSpinX = number of spin along the x direction
// nbrSpinY = number of spin along the y direction
// jFactor = amplitude of the Ising term
// hxFactor = amplitudes of the Zeeman term along x
// hzFactor = amplitudes of the Zeeman term along z
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian (AbstractSpinChain* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, double jFactor, double jSquareExchangeFactor, int offset, long memory)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->JFactor = jFactor;
  this->JSquareExchangeFactor = jSquareExchangeFactor;
  
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  
  this->Offset = offset;
//   this->HermitianSymmetryFlag = true;
  this->HermitianSymmetryFlag = false;
  
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); ++i)
    this->SzSzContributions[i] = 0.0;
  this->EvaluateExponentialFactors();
  /*
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024l)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1l << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1l << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
      cout << endl;
      this->EnableFastMultiplication();
    }*/
}


// destructor
//

TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::~TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian() 
{
  for (int i = 0; i < this->NbrSpinX; ++i)
    delete[] this->ExponentialFactors[i];
  delete[] this->ExponentialFactors;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::ShiftHamiltonian (double shift)
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

ComplexVector& TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int pos;
  int TmpIndex;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j + 1, k);
	      
	      pos = this->Chain->Pij(TmpIndex1, TmpIndex, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += vSource[i] * TmpCoefficient * this->JFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      
	      TmpIndex1 = this->GetLinearizedIndexSafe(j - this->Offset, k + 1);      
	      pos = this->Chain->Pij(TmpIndex1, TmpIndex, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += vSource[i] * TmpCoefficient * this->JFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      
	      
	      TmpIndex1 = this->GetLinearizedIndexSafe (j + 1, k);
	      TmpIndex2 = this->GetLinearizedIndexSafe (j + 1 - this->Offset, k + 1);
	      TmpIndex3 = this->GetLinearizedIndexSafe (j - this->Offset, k + 1);
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex1, TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += vSource[i]* TmpCoefficient * this->JSquareExchangeFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex3, TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += vSource[i]* TmpCoefficient * this->JSquareExchangeFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
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

ComplexVector& TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int pos;
  int TmpIndex;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  double TmpCoefficient;
  Complex TmpSum;
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      TmpSum = 0.0;
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j + 1, k);
	      
	      pos = this->Chain->Pij(TmpIndex1, TmpIndex, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < i)
	      {
		vDestination[pos] += vSource[i] * TmpCoefficient * this->JFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		TmpSum += vSource[pos] * TmpCoefficient * this->JFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	      }
	      if (pos == i)
		TmpSum += vSource[pos] * TmpCoefficient * this->JFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	      
	      TmpIndex1 = this->GetLinearizedIndexSafe(j - this->Offset, k + 1);      
	      pos = this->Chain->Pij(TmpIndex1, TmpIndex, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < i)
	      {
		vDestination[pos] += vSource[i] * TmpCoefficient * this->JFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		TmpSum += vSource[pos] * TmpCoefficient * this->JFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	      }
	      if (pos == i)
		TmpSum += vSource[pos] * TmpCoefficient * this->JFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	      
	      
	      TmpIndex1 = this->GetLinearizedIndexSafe (j + 1, k);
	      TmpIndex2 = this->GetLinearizedIndexSafe (j + 1 - this->Offset, k + 1);
	      TmpIndex3 = this->GetLinearizedIndexSafe (j - this->Offset, k + 1);
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex1, TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < dim)
		vDestination[pos] += vSource[i]* TmpCoefficient * this->JSquareExchangeFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex3, TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < dim)
		vDestination[pos] += vSource[i]* TmpCoefficient * this->JSquareExchangeFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
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

ComplexVector* TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int pos;
  int TmpIndex;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  Complex* TmpValues = new Complex[nbrVectors];
  for (int l = 0; l < nbrVectors; ++l)
    {
      ComplexVector& TmpSource = vSources[l];
      ComplexVector& TmpDestination = vDestinations[l];
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
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {      
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j + 1, k);
	      
	      pos = this->Chain->Pij(TmpIndex1, TmpIndex, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * this->JFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      
	      TmpIndex1 = this->GetLinearizedIndexSafe(j - this->Offset, k + 1);  
	      pos = this->Chain->Pij(TmpIndex1, TmpIndex, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * this->JFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      
	      
	      TmpIndex1 = this->GetLinearizedIndexSafe (j + 1, k);
	      TmpIndex2 = this->GetLinearizedIndexSafe (j + 1 - this->Offset, k + 1);
	      TmpIndex3 = this->GetLinearizedIndexSafe (j - this->Offset, k + 1);
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex1, TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * this->JSquareExchangeFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex3, TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * this->JSquareExchangeFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
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

ComplexVector* TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int pos;
  int TmpIndex;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  Complex* TmpValues = new Complex[nbrVectors];
  Complex* TmpSums = new Complex[nbrVectors];
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int l = 0; l < nbrVectors; ++l)
	{
	  TmpValues[l] = vSources[l][i];
	  TmpSums[l] = 0.0;
	}
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      TmpIndex = this->GetLinearizedIndexSafe(j, k);
	      TmpIndex1 = this->GetLinearizedIndexSafe(j + 1, k);
	      
	      pos = this->Chain->Pij(TmpIndex1, TmpIndex, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < i)
		for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * this->JFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSums[l] += vSources[l][pos] * TmpCoefficient * this->JFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
	      
	      TmpIndex1 = this->GetLinearizedIndexSafe(j - this->Offset, k + 1);  
	      pos = this->Chain->Pij(TmpIndex1, TmpIndex, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < i)
		for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * this->JFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSums[l] += vSources[l][pos] * TmpCoefficient * this->JFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
	      if (pos == i)
		for (int l = 0; l < nbrVectors; ++l)
		  TmpSums[l] += vSources[l][pos] * TmpCoefficient * this->JFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		
		
	      TmpIndex1 = this->GetLinearizedIndexSafe (j + 1, k);
	      TmpIndex2 = this->GetLinearizedIndexSafe (j + 1 - this->Offset, k + 1);
	      TmpIndex3 = this->GetLinearizedIndexSafe (j - this->Offset, k + 1);
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex1, TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * this->JSquareExchangeFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	
	      pos = this->Chain->Pijkl(TmpIndex, TmpIndex3, TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * this->JSquareExchangeFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	    }
	}
	      
	for (int l = 0; l < nbrVectors; ++l)
	  vDestinations[l][i] += (this->SzSzContributions[i] * vSources[l][i] + TmpSums[l]);
    }
  delete[] TmpValues;
  delete[] TmpSums;
  return vDestinations;
}


// ask if Hamiltonian implements hermitian symmetry operations
//

bool TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// evaluate all exponential factors
//   

void TwoDimensionalSquareLatticeSU3And2DTranslationHamiltonian::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex*[this->NbrSpinX];
  for (int i = 0; i < this->NbrSpinX; ++i)
    { 
      this->ExponentialFactors[i] = new Complex[this->NbrSpinY];
      for (int j = 0; j < this->NbrSpinY; ++j)
	{ 
	  this->ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->NbrSpinX))
							       + (this->YMomentum * ((double) j) / ((double) this->NbrSpinY))));
	}
    }
}
