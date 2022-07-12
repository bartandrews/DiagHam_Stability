////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                         Class author: Cecile Repellin                      //
//                                                                            //
//                                                                            //
//         class of two dimensional spin model on the kagome lattice          //
//                                                                            //
//                        last modification : 18/06/2016                      //
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


#include "Hamiltonian/TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian.h"
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

TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian (AbstractSpinChain* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, double jFactor, double jDownFactor, double jEasyPlaneFactor, double jDownEasyPlaneFactor, int offset, long memory)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->JFactor = jFactor;
  this->JDownFactor = jDownFactor;
  this->JEasyPlaneFactor = jEasyPlaneFactor;
  this->JDownEasyPlaneFactor = jDownEasyPlaneFactor;
  this->JBreakD3Factor = 1.0;
  
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  
  this->Offset = offset;
  this->HermitianSymmetryFlag = true;
//   this->HermitianSymmetryFlag = false;
  
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateExponentialFactors();
  this->EvaluateDiagonalMatrixElements();
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



// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// nbrSpinX = number of spin along the x direction
// nbrSpinY = number of spin along the y direction
// jFactor = amplitude of the Ising term
// hxFactor = amplitudes of the Zeeman term along x
// hzFactor = amplitudes of the Zeeman term along z
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian (AbstractSpinChain* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, double jFactor, double jDownFactor, double jBreakD3Factor, int offset, long memory)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->JFactor = jFactor;
  this->JDownFactor = jDownFactor;
  this->JEasyPlaneFactor = jFactor;
  this->JDownEasyPlaneFactor = jDownFactor;
  this->JBreakD3Factor = jBreakD3Factor;
  
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  
  this->Offset = offset;
  this->HermitianSymmetryFlag = true;
//   this->HermitianSymmetryFlag = false;
  
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateExponentialFactors();
  this->EvaluateDiagonalMatrixElements();
}


// destructor
//

TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::~TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian() 
{
  for (int i = 0; i < this->NbrSpinX; ++i)
    delete[] this->ExponentialFactors[i];
  delete[] this->ExponentialFactors;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::ShiftHamiltonian (double shift)
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

ComplexVector& TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int pos3;
  int MaxPos = this->NbrSpin - 1;
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
	      //upward triangles
	      TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
	      TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
	      TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
	      //AB
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      //AC
	      pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      //BC
	      pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      
	      //downward triangles
	      if (this->NbrSpinX > 1)
	      {
		//BA
		pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j - 1, k, 1), this->GetLinearizedIndex(j, k, 0), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      }
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//CA
		pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - this->Offset, k - 1, 2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		if (this->NbrSpinX > 1)
		{
		  //BC
		  pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		    vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j - 1, k, 1), this->GetLinearizedIndex(j - this->Offset, k - 1, 2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		    vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		}
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

ComplexVector& TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int pos3;
  int MaxPos = this->NbrSpin - 1;
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
	      //upward triangles
	      TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
	      TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
	      TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
	      //AB
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < i)
	      {
		vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	      }
	      //AC
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	       if (pos < i)
	      {
		vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	      }
	      //BC
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	       if (pos < i)
	      {
		vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	      }
	      
	      //downward triangles
	      if (this->NbrSpinX > 1)
	      {
		//BA
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k, 1);
		TmpIndex2 = this->GetLinearizedIndex(j, k, 0);
		pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos < i)
		{
		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
	      }
	      
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//CA
		TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
		TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1, 2);
		pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos < i)
		{
		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JDownEasyPlaneFactor * this->JBreakD3Factor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
		if (this->NbrSpinX > 1)
		{
		  //BC
		  TmpIndex1 = this->GetLinearizedIndex(j - 1, k, 1);
		  TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1, 2);
		  pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos < i)
		  {
		    vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		    TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		  }
		}
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

ComplexVector* TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
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
	      //upward triangles
	      TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
	      TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
	      TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
	      //AB
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      //AC
	      pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      //BC
	      pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor  * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor  * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      
	      //downward triangles
	      if (this->NbrSpinX > 1)
	      {
		//BA
		pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j - 1, k, 1), this->GetLinearizedIndex(j, k, 0), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	      }
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//CA
		pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JDownEasyPlaneFactor * this->JBreakD3Factor  * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - this->Offset, k - 1, 2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JDownEasyPlaneFactor  * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		if (this->NbrSpinX > 1)
		{
		  //BC
		  pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  pos = this->Chain->SmiSpj(this->GetLinearizedIndex(j - 1, k, 1), this->GetLinearizedIndex(j - this->Offset, k - 1, 2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 0.5 * TmpValues[l] *TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		}
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

ComplexVector* TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
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
	      //upward triangles
	      TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
	      TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
	      TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
	      //AB
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos < i)
	      {
		for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JEasyPlaneFactor  * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSums[l] += 0.5 * vSources[l][pos] * TmpCoefficient * this->JEasyPlaneFactor  * this->JBreakD3Factor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
	      }
	      //AC
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	       if (pos < i)
	      {
		for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSums[l] += 0.5 * vSources[l][pos] * TmpCoefficient * this->JEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
	      }
	      //BC
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	       if (pos < i)
	      {
		for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JEasyPlaneFactor  * this->JBreakD3Factor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSums[l] += 0.5 * vSources[l][pos] * TmpCoefficient * this->JEasyPlaneFactor * this->JBreakD3Factor  * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
	      }
	      
	      //downward triangles
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//BA
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k, 1);
		TmpIndex2 = this->GetLinearizedIndex(j, k, 0);
		pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos < i)
		{
		  for (int l = 0; l < nbrVectors; ++l)
		  {
		    vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		    TmpSums[l] += 0.5 * vSources[l][pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		  }
		}
	      }
	      
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//CA
		TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
		TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1, 2);
		pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos < i)
		{
		  for (int l = 0; l < nbrVectors; ++l)
		  {
		    vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JDownEasyPlaneFactor * this->JBreakD3Factor  * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		    TmpSums[l] += 0.5 * vSources[l][pos] * TmpCoefficient * this->JDownEasyPlaneFactor  * this->JBreakD3Factor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		  }
		}
		if (this->NbrSpinX > 1)
		{
		  //BC
		  TmpIndex1 = this->GetLinearizedIndex(j - 1, k, 1);
		  TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1, 2);
		  pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos < i)
		  {
		    for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][pos] += 0.5 * TmpValues[l] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		      TmpSums[l] += 0.5 * vSources[l][pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		    }
		  }
		}
	      }
	    }
	}
	for (int l = 0; l < nbrVectors; ++l)
	  vDestinations[l][i] += (this->SzSzContributions[i] * vSources[l][i] + TmpSums[l]);
    }
  delete[] TmpValues;
  delete[] TmpSums;
  return vDestinations;
}


// evaluate diagonal matrix elements
// 

void TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoefficient;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] = 0.0;
      
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      //upward triangles
	      TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
	      TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
	      TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
	      //AB
	      TmpCoefficient = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
	      this->SzSzContributions[i] += this->JFactor  * this->JBreakD3Factor * TmpCoefficient;
	      //AC
	      TmpCoefficient = this->Chain->SziSzj(TmpIndex1, TmpIndex3, i);
	      this->SzSzContributions[i] += this->JFactor * TmpCoefficient;
	      //BC
	      TmpCoefficient = this->Chain->SziSzj(TmpIndex2, TmpIndex3, i);
	      this->SzSzContributions[i] += this->JFactor  * this->JBreakD3Factor * TmpCoefficient;
	      
	      //downward triangles
	      if (this->NbrSpinX > 1)
	      {
		//BA
		TmpCoefficient = this->Chain->SziSzj(this->GetLinearizedIndex(j - 1, k, 1), this->GetLinearizedIndex(j, k, 0), i);
		this->SzSzContributions[i] += this->JDownFactor * TmpCoefficient;
	      }
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//CA
		TmpCoefficient = this->Chain->SziSzj(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - this->Offset, k - 1, 2), i);
		this->SzSzContributions[i] += this->JDownFactor  * this->JBreakD3Factor * TmpCoefficient;
		if (this->NbrSpinX > 1)
		{
		  //BC
		  TmpCoefficient = this->Chain->SziSzj(this->GetLinearizedIndex(j - 1, k, 1), this->GetLinearizedIndex(j - this->Offset, k - 1, 2), i);
		  this->SzSzContributions[i] += this->JDownFactor * TmpCoefficient;
		}
	      }
	    }
	}
    }
}
/*

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::HermitianLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  double TmpCoefficient;
  double TmpCoef;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  AbstractSpinChain* TmpParticles = (AbstractSpinChain*) this->Particles->Clone();
  
  int* TmpIndexArray;
  Complex* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      Coefficient = vSource[l];
      for (j = 0; j < TmpNbrInteraction; ++j)
	vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
      vDestination[l] += this->SzSzContributions[l] * Coefficient;
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  int Index;  
  Complex TmpSum;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	{
	  TmpSum = 0.0;
	  for (int j = 0; j < this->NbrSpinX; j++)
	  {
	    for (int k = 0; k < this->NbrSpinY; k++)
	      {
		//upward triangles
		TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
		TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
		TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
		//AB
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos < i)
		{
		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
		//AC
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos < i)
		{
		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
		//BC
		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos < i)
		{
		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		  TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		}
	      
		//downward triangles
		if (this->NbrSpinX > 1)
		{
		  //BA
		  TmpIndex1 = this->GetLinearizedIndex(j - 1, k, 1);
		  TmpIndex2 = this->GetLinearizedIndex(j, k, 0);
		  pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos < i)
		  {
		    vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		    TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		  }
		}
	      
		if ((this->NbrSpinY > 1) || (this->Offset != 0))
		{
		  //CA
		  TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
		  TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1, 2);
		  pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos < i)
		  {
		    vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		    TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		  }
		  if (this->NbrSpinX > 1)
		  {
		    //BC
		    TmpIndex1 = this->GetLinearizedIndex(j - 1, k, 1);
		    TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1, 2);
		    pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		    if (pos < i)
		    {
		      vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
		      TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
		    }
		  }
		}
	      
	      }
	  }
	  vDestination[i] += (this->SzSzContributions[i] * vSource[i] + TmpSum);
	}
      }
      
}*/
/*
// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted
// return value = number of non-zero matrix element

long TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  double TmpCoefficient;
  long Memory = 0l;
  int NbrTranslationsX;
  int NbrTranslationsY;
  
  AbstractSpinChain* TmpParticles = (AbstractSpinChain*) this->Particles->Clone();
  int LastComponent = nbrComponent + firstComponent;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->NbrInteractionPerComponent[i - this->PrecalculationShift] = 0;
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      //upward triangles
	      TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
	      TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
	      TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
	      //AB
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      //AC
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      //BC
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      
	      //downward triangles
	      if (this->NbrSpinX > 1)
	      {
		//BA
		TmpIndex1 = min(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1));
		TmpIndex2 = max(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1));
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  if ((this->HermitianSymmetryFlag == false) || (pos < i))
		  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
		}
		if (this->HermitianSymmetryFlag == false)
		{
		  pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
		}
	      }
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//CA
		TmpIndex1 = min(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0));
		TmpIndex2 = max(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0));
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  if ((this->HermitianSymmetryFlag == false) || (pos < i))
		  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
		}
		if (this->HermitianSymmetryFlag == false)
		{
		  pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
		}
		if (this->NbrSpinX > 1)
		{
		  //BC
		  TmpIndex1 = min(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1));
		  TmpIndex2 = max(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1));
		  pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    if ((this->HermitianSymmetryFlag == false) || (pos < i))
		    {
		      ++Memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  }
		  if (this->HermitianSymmetryFlag == false)
		  {
		    pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		    if (pos != dim)
		    {
		      ++Memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  }
		}
	      }
	      
	    }
	}
    }    
  delete TmpParticles;
  return Memory;
}


// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  double TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int LastComponent = nbrComponent + firstComponent;
  AbstractSpinChain* TmpParticles = (AbstractSpinChain*) this->Particles->Clone();
  
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
    
  int CurrentNbrCounting = 0;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      //upward triangles
	      TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
	      TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
	      TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
	      //AB
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      //AC
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      //BC
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      
	      //downward triangles
	      if (this->NbrSpinX > 1)
	      {
		//BA
		TmpIndex1 = min(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1));
		TmpIndex2 = max(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1));
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  if ((this->HermitianSymmetryFlag == false) || (pos < i))
		  {
		    this->InteractionPerComponentIndex[Pos][position] = pos;
		    this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		    ++position;
		  }
		}
		if (this->HermitianSymmetryFlag == false)
		{
		  pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    this->InteractionPerComponentIndex[Pos][position] = pos;
		    this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		    ++position;
		  }
		}
	      }
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//CA
		TmpIndex1 = min(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0));
		TmpIndex2 = max(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0));
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  if ((this->HermitianSymmetryFlag == false) || (pos < i))
		  {
		    this->InteractionPerComponentIndex[Pos][position] = pos;
		    this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		    ++position;
		  }
		}
		if (this->HermitianSymmetryFlag == false)
		{
		  pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    this->InteractionPerComponentIndex[Pos][position] = pos;
		    this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		    ++position;
		  }
		}
		if (this->NbrSpinX > 1)
		{
		  //BC
		  TmpIndex1 = min(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1));
		  TmpIndex2 = max(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1));
		  pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    if ((this->HermitianSymmetryFlag == false) || (pos < i))
		    {
		      this->InteractionPerComponentIndex[Pos][position] = pos;
		      this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		      ++position;
		    }
		  }
		  if (this->HermitianSymmetryFlag == false)
		  {
		    pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		    if (pos != dim)
		    {
		      this->InteractionPerComponentIndex[Pos][position] = pos;
		      this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		      ++position;
		    }
		  }
		}
	      }
	      
	    }
	}
      ++Pos;
    }
}*/

// ask if Hamiltonian implements hermitian symmetry operations
//

bool TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// evaluate all exponential factors
//   

void TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::EvaluateExponentialFactors()
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
