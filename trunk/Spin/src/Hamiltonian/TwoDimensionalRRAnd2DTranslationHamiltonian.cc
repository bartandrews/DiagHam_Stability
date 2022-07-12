////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of two dimension spin model that could host             //
//                a Read-Rezayi Z3 phase with 2d translations                 //
//                                                                            //
//                        last modification : 27/07/2018                      //
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


#include "Hamiltonian/TwoDimensionalRRAnd2DTranslationHamiltonian.h"
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

TwoDimensionalRRAnd2DTranslationHamiltonian::TwoDimensionalRRAnd2DTranslationHamiltonian()
{
  this->J1Factor = 0.0;
  this->J2Factor = 0.0;
  this->J3Factor = 0.0;
  this->JcFactor = 0.0;
}

// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// xMomentum = momentum along the x direction
// nbrSpinX = number of spin along the x direction
// yMomentum = momentum along the y direction
// nbrSpinY = number of spin along the y direction
// j1Factor = amplitude of the Heisenberg coupling between nearest neighbors
// j2Factor = amplitude of the (S_i S_j)^2 nearest neighbor coupling
// j3Factor = amplitude of the (S_i S_j)^3 nearest neighbor coupling
// jcFactor = amplitude of the chiral term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

TwoDimensionalRRAnd2DTranslationHamiltonian::TwoDimensionalRRAnd2DTranslationHamiltonian(AbstractSpinChain* chain, int xMomentum, int nbrSpinX, 
											 int yMomentum, int nbrSpinY, double j1Factor, double j2Factor, double j3Factor, double jcFactor,
											 AbstractArchitecture* architecture, long memory)
{
  this->Chain = chain;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
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
  this->HalfJFactor = 0.5 * this->JFactor;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = true;
  this->EvaluateExponentialFactors();
  this->EvaluateDiagonalMatrixElements();
  if (memory == 0l)
    {
      this->FastMultiplicationFlag = false;
    }
  else
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
    }
}

// destructor
//

TwoDimensionalRRAnd2DTranslationHamiltonian::~TwoDimensionalRRAnd2DTranslationHamiltonian() 
{
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& TwoDimensionalRRAnd2DTranslationHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      Complex& TmpValue = vSource[i];
      for (int j = 0; j < this->NbrSpinX; j++)
 	{
 	  for (int k = 0; k < this->NbrSpinY; k++)
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

ComplexVector* TwoDimensionalRRAnd2DTranslationHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
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
  int NbrTranslationsX;
  int NbrTranslationsY;
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
      for (int j = 0; j < this->NbrSpinX; j++)
 	{
 	  for (int k = 0; k < this->NbrSpinY; k++)
 	    {
	      int TmpIndex1 = this->GetLinearizedIndex(j, k);
	      int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
	      int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
	      int TmpIndex4 = this->GetSafeLinearizedIndex(j + 1, k + 1);
	      
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex2, i, dim, vDestinations, nbrVectors, TmpValues);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex2, TmpIndex1, i, dim, vDestinations, nbrVectors, TmpValues);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex3, i, dim, vDestinations, nbrVectors, TmpValues);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex3, TmpIndex1, i, dim, vDestinations, nbrVectors, TmpValues);
	      
	      this->EvaluateOffDiagonalChiralContribution(TmpIndex1, TmpIndex2, TmpIndex4, i, dim, vDestinations, nbrVectors, TmpValues);
	      this->EvaluateOffDiagonalChiralContribution(TmpIndex2, TmpIndex4, TmpIndex3, i, dim, vDestinations, nbrVectors, TmpValues);
	      this->EvaluateOffDiagonalChiralContribution(TmpIndex4, TmpIndex3, TmpIndex1, i, dim, vDestinations, nbrVectors, TmpValues);
	      this->EvaluateOffDiagonalChiralContribution(TmpIndex3, TmpIndex1, TmpIndex2, i, dim, vDestinations, nbrVectors, TmpValues);
	    }
	}
    }
  delete[] TmpValues;
  return vDestinations;
}

// core part of the AddMultiply method
// 
// chain = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added  

void TwoDimensionalRRAnd2DTranslationHamiltonian::HermitianEvaluateAddMultiplyComponent(AbstractSpinChain* chain, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int NbrTranslationsX;
  int NbrTranslationsY;
  vDestination[index] += this->SzSzContributions[index - this->PrecalculationShift] * vSource[index];
  for (int j = 0; j < this->NbrSpinX; j++)
    {
      for (int k = 0; k < this->NbrSpinY; k++)
	{
	  int TmpIndex1 = this->GetLinearizedIndex(j, k);
	  int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
	  int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
	  int TmpIndex4 = this->GetSafeLinearizedIndex(j + 1, k + 1);
	  
	  this->HermitianEvaluateOffDiagonalPowerHeisenbergContribution(chain, TmpIndex1, TmpIndex2, index, vSource, vDestination);
	  this->HermitianEvaluateOffDiagonalPowerHeisenbergContribution(chain, TmpIndex2, TmpIndex1, index, vSource, vDestination);
	  this->HermitianEvaluateOffDiagonalPowerHeisenbergContribution(chain, TmpIndex1, TmpIndex3, index, vSource, vDestination);
	  this->HermitianEvaluateOffDiagonalPowerHeisenbergContribution(chain, TmpIndex3, TmpIndex1, index, vSource, vDestination);
	  
	  this->HermitianEvaluateOffDiagonalChiralContribution(chain, TmpIndex1, TmpIndex2, TmpIndex4, index, vSource, vDestination);
	  this->HermitianEvaluateOffDiagonalChiralContribution(chain, TmpIndex2, TmpIndex4, TmpIndex3, index, vSource, vDestination);
	  this->HermitianEvaluateOffDiagonalChiralContribution(chain, TmpIndex4, TmpIndex3, TmpIndex1, index, vSource, vDestination);
	  this->HermitianEvaluateOffDiagonalChiralContribution(chain, TmpIndex3, TmpIndex1, TmpIndex2, index, vSource, vDestination);
	}
    }
}

// core part of the AddMultiply method for a set of vectors
// 
// chain = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

void TwoDimensionalRRAnd2DTranslationHamiltonian::HermitianEvaluateAddMultiplyComponent(AbstractSpinChain* chain, int index, ComplexVector* vSources, 
											ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
}

// evaluate diagonal matrix elements
// 

void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateDiagonalMatrixElements()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int dim = this->Chain->GetHilbertSpaceDimension();
  int TmpNbrComponents = ((int) (MaxIndex - MinIndex)) + 1;
  // SzSz part
  for (int i = 0; i < TmpNbrComponents; i++)
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
					 this->GetSafeLinearizedIndex(j, k + 1), i + this->PrecalculationShift);
	      Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
	      Tmp2 = this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					 this->GetSafeLinearizedIndex(j + 1, k), i + this->PrecalculationShift);
	      Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
	    }
	}
      this->SzSzContributions[i] += Tmp;
   }
}

// evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
//
// i = linearized position of the first spin
// j = linearized position of the second spin
// index = index of the many-body state to act on
// dimension = total Hilbert space dimension
// vDestination = vector at which result has to be added
// coefficient = global multiplicative coefficient

void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, ComplexVector& vDestination, Complex& coefficient)
{
  double TmpCoefficient;
  double TmpCoefficient2;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int pos = this->Chain->SmiSpj(i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J1Factor * TmpValue2;
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      TmpValue2 *= TmpValue3;
      vDestination[pos] += (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
    }

  pos = this->Chain->SmiSpjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
    }
  pos = this->Chain->SmiSpjSmkSpl(j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
    }
  pos = this->Chain->SziSzjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
    }

  pos = this->Chain->SmiSpjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 =  (0.125 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SmiSpjSmkSplSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 =  (0.125 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SmiSpjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 =  (0.125 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SmiSpjSmkSplSmmSpn(j, i, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 =  (0.125 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J3Factor * TmpValue2;	      
    }

  pos = this->Chain->SziSzjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SziSzjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }

  pos = this->Chain->SmiSpjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SmiSpjSzkSzlSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }

  pos = this->Chain->SziSzjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }
}

// evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
//
// i = linearized position of the first spin
// j = linearized position of the second spin
// index = index of the many-body state to act on
// dimension = total Hilbert space dimension
// vDestinations = vectors to which results have to be added
// nbrVectors = number of vectors that have to be evaluated together
// coefficients = global multiplicative coefficients

void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, ComplexVector* vDestinations, int nbrVectors, Complex* coefficients)
{
  double TmpCoefficient;
  double TmpCoefficient2;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int pos = this->Chain->SmiSpj(i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];
	  vDestinations[k][pos] += (this->J1Factor + ((this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue3)) * TmpValue2;
	}
    }

  pos = this->Chain->SmiSpjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];
	  vDestinations[k][pos] +=  (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
	}
    }
  pos = this->Chain->SmiSpjSmkSpl(j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];
	  vDestinations[k][pos] +=  (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
	}
    }
  pos = this->Chain->SziSzjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];
	  vDestinations[k][pos] +=  (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
	}
    }

  pos = this->Chain->SmiSpjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.125;
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  vDestinations[k][pos] +=  TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];	      
	}
    }
  pos = this->Chain->SmiSpjSmkSplSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.125;
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  vDestinations[k][pos] +=  TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];	      
	}
    }
  pos = this->Chain->SmiSpjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.125;
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  vDestinations[k][pos] +=  TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];	      
	}
    }
  pos = this->Chain->SmiSpjSmkSplSmmSpn(j, i, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.125;
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  vDestinations[k][pos] +=  TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];	      
	}
    }

  pos = this->Chain->SziSzjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.25;
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  vDestinations[k][pos] +=  TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];	      
	}
    }
  pos = this->Chain->SziSzjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.25;
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  vDestinations[k][pos] +=  TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];	      
	}
    }

  pos = this->Chain->SmiSpjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.25;
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  vDestinations[k][pos] +=  TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];	      
	}
    }
  pos = this->Chain->SmiSpjSzkSzlSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.25;
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  vDestinations[k][pos] +=  TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];	      
	}
    }

  pos = this->Chain->SziSzjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.5;
      TmpCoefficient *= this->J3Factor;
      for (int k = 0; k < nbrVectors; ++k)
	{	  
	  vDestinations[k][pos] += TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficients[k];	      
	}
    }  
}

// evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
//
// i = linearized position of the first spin
// j = linearized position of the second spin
// k = linearized position of the second spin
// index = index of the many-body state to act on
// dimension = total Hilbert space dimension
// vDestination = vector at which result has to be added
// coefficient = global multiplicative coefficient

void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateOffDiagonalChiralContribution(int i, int j, int k, int index, int dimension, ComplexVector& vDestination, Complex& coefficient)
{
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int pos = this->Chain->SpiSmjSzk(i, j, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[pos].Re -= TmpCoefficient * Tmp.Im;
      vDestination[pos].Im += TmpCoefficient * Tmp.Re;
    }
  pos = this->Chain->SpiSmjSzk(i, k, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[pos].Re += TmpCoefficient * Tmp.Im;
      vDestination[pos].Im -= TmpCoefficient * Tmp.Re;
    }
  pos = this->Chain->SpiSmjSzk(j, k, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[pos].Re -= TmpCoefficient * Tmp.Im;
      vDestination[pos].Im += TmpCoefficient * Tmp.Re;
    }
  pos = this->Chain->SpiSmjSzk(k, j, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[pos].Re += TmpCoefficient * Tmp.Im;
      vDestination[pos].Im -= TmpCoefficient * Tmp.Re;
    }
  pos = this->Chain->SpiSmjSzk(k, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[pos].Re -= TmpCoefficient * Tmp.Im;
      vDestination[pos].Im += TmpCoefficient * Tmp.Re;
    }
  pos = this->Chain->SpiSmjSzk(j, i, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[pos].Re += TmpCoefficient * Tmp.Im;
      vDestination[pos].Im -= TmpCoefficient * Tmp.Re;
    }
}

// evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
//
// i = linearized position of the first spin
// j = linearized position of the second spin
// k = linearized position of the second spin
// index = index of the many-body state to act on
// dimension = total Hilbert space dimension
// vDestinations = vectors to which results have to be added
// nbrVectors = number of vectors that have to be evaluated together
// coefficients = global multiplicative coefficients

void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateOffDiagonalChiralContribution(int i, int j, int k, int index, int dimension, 
											       ComplexVector* vDestinations, int nbrVectors, Complex* coefficients)
{
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int pos = this->Chain->SpiSmjSzk(i, j, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      for (int k = 0; k < nbrVectors; ++k)
	{
	  Complex Tmp = coefficients[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  vDestinations[k][pos].Re -= TmpCoefficient * Tmp.Im;
	  vDestinations[k][pos].Im += TmpCoefficient * Tmp.Re;
	}
    }
  pos = this->Chain->SpiSmjSzk(i, k, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      for (int k = 0; k < nbrVectors; ++k)
	{
	  Complex Tmp = coefficients[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  vDestinations[k][pos].Re += TmpCoefficient * Tmp.Im;
	  vDestinations[k][pos].Im -= TmpCoefficient * Tmp.Re;
	}
    }
  pos = this->Chain->SpiSmjSzk(j, k, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      for (int k = 0; k < nbrVectors; ++k)
	{
	  Complex Tmp = coefficients[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  vDestinations[k][pos].Re -= TmpCoefficient * Tmp.Im;
	  vDestinations[k][pos].Im += TmpCoefficient * Tmp.Re;
	}
    }
  pos = this->Chain->SpiSmjSzk(k, j, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      for (int k = 0; k < nbrVectors; ++k)
	{
	  Complex Tmp = coefficients[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  vDestinations[k][pos].Re += TmpCoefficient * Tmp.Im;
	  vDestinations[k][pos].Im -= TmpCoefficient * Tmp.Re;
	}
    }
  pos = this->Chain->SpiSmjSzk(k, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      for (int k = 0; k < nbrVectors; ++k)
	{
	  Complex Tmp = coefficients[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  vDestinations[k][pos].Re -= TmpCoefficient * Tmp.Im;
	  vDestinations[k][pos].Im += TmpCoefficient * Tmp.Re;
	}
    }
  pos = this->Chain->SpiSmjSzk(j, i, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      for (int k = 0; k < nbrVectors; ++k)
	{
	  Complex Tmp = coefficients[k] * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  vDestinations[k][pos].Re += TmpCoefficient * Tmp.Im;
	  vDestinations[k][pos].Im -= TmpCoefficient * Tmp.Re;
	}
    }
}

// evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
//
// chain = pointer to the Hilbert space
// i = linearized position of the first spin
// j = linearized position of the second spin
// dimension = total Hilbert space dimension
// vSource = vector  to be multiplied
// vDestination = vector to which result has to be added

void TwoDimensionalRRAnd2DTranslationHamiltonian::HermitianEvaluateOffDiagonalPowerHeisenbergContribution(AbstractSpinChain* chain, int i, int j, int index,  ComplexVector& vSource, ComplexVector& vDestination)
{
  double TmpCoefficient;
  double TmpCoefficient2;
  int NbrTranslationsX;
  int NbrTranslationsY;
  Complex Coefficient = vSource[index];
  Complex TmpSum = 0.0;
  int TmpIndex = chain->SmiSpj(i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      double TmpValue3 = chain->SziSzj(i, j, index);
      vDestination[TmpIndex] += (this->J1Factor + ((this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue3)) * TmpValue2 * Coefficient;
      if (TmpIndex < index)
	{
	  TmpSum += (this->J1Factor + ((this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue3)) * Conj(TmpValue2) * vSource[TmpIndex];
	}
    }
  TmpIndex = chain->SmiSpjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      double TmpValue3 = chain->SziSzj(i, j, index);
      vDestination[TmpIndex] += (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += (this->J2Factor + (this->J3Factor * TmpValue3)) * Conj(TmpValue2) * vSource[TmpIndex];
	}
    }
  TmpIndex = chain->SmiSpjSmkSpl(j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      double TmpValue3 = chain->SziSzj(i, j, index);
      vDestination[TmpIndex] += (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += (this->J2Factor + (this->J3Factor * TmpValue3)) * Conj(TmpValue2) * vSource[TmpIndex];
	}
    }
  TmpIndex = chain->SziSzjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      double TmpValue3 = chain->SziSzj(i, j, index);
      vDestination[TmpIndex] += (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += (this->J2Factor + (this->J3Factor * TmpValue3)) * Conj(TmpValue2) * vSource[TmpIndex];
	}
    }

  TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 =  (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex] += TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += Conj(TmpValue2) * vSource[TmpIndex];
	}
    }
  TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 =  (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex] += TmpValue2 * Coefficient;
      if (TmpIndex < index)
	{
	  TmpSum += Conj(TmpValue2) * vSource[TmpIndex];
	}
    }
  TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 =  (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex] += TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += Conj(TmpValue2) * vSource[TmpIndex];
	}
    }
  TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 =  (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex] += TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += Conj(TmpValue2) * vSource[TmpIndex];
	}
    }

  TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex] +=  TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += Conj(TmpValue2) * vSource[TmpIndex];
	}
    }
  TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex] +=  TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += Conj(TmpValue2) * vSource[TmpIndex];
	}
    }

  TmpIndex = chain->SmiSpjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex] +=  TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += Conj(TmpValue2) * vSource[TmpIndex];
	}
    }
  TmpIndex = chain->SmiSpjSzkSzlSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex] +=  TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += Conj(TmpValue2) * vSource[TmpIndex];
	}
    }

  TmpIndex = chain->SziSzjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      Complex TmpValue2 = (0.5 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex] +=  TmpValue2 * Coefficient;	      
      if (TmpIndex < index)
	{
	  TmpSum += Conj(TmpValue2) * vSource[TmpIndex];
	}
    }

  vDestination[index] += TmpSum;
}


// evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
//
// chain = pointer to the Hilbert space
// i = linearized position of the first spin
// j = linearized position of the second spin
// k = linearized position of the second spin
// index = index of the many-body state to act on
// vDestination = vector at which result has to be added
// coefficient = global multiplicative coefficient

void TwoDimensionalRRAnd2DTranslationHamiltonian::HermitianEvaluateOffDiagonalChiralContribution(AbstractSpinChain* chain, int i, int j, int k, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  Complex Coefficient = vSource[index];
  Complex TmpSum = 0.0;
  int TmpIndex = chain->SpiSmjSzk(i, j, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex].Re -= TmpCoefficient * Tmp.Im;
      vDestination[TmpIndex].Im += TmpCoefficient * Tmp.Re;
      if (TmpIndex < index)
	{
	  Tmp = vSource[TmpIndex] * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	  TmpSum.Re += TmpCoefficient * Tmp.Im;
	  TmpSum.Im -= TmpCoefficient * Tmp.Re;
	}
    }
  TmpIndex = chain->SpiSmjSzk(i, k, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex].Re += TmpCoefficient * Tmp.Im;
      vDestination[TmpIndex].Im -= TmpCoefficient * Tmp.Re;
      if (TmpIndex < index)
	{
	  Tmp = vSource[TmpIndex] * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	  TmpSum.Re -= TmpCoefficient * Tmp.Im;
	  TmpSum.Im += TmpCoefficient * Tmp.Re;
	}
    }
  TmpIndex = chain->SpiSmjSzk(j, k, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex].Re -= TmpCoefficient * Tmp.Im;
      vDestination[TmpIndex].Im += TmpCoefficient * Tmp.Re;
      if (TmpIndex < index)
	{
	  Tmp = vSource[TmpIndex] * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	  TmpSum.Re += TmpCoefficient * Tmp.Im;
	  TmpSum.Im -= TmpCoefficient * Tmp.Re;
	}
    }
  TmpIndex = chain->SpiSmjSzk(k, j, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex].Re += TmpCoefficient * Tmp.Im;
      vDestination[TmpIndex].Im -= TmpCoefficient * Tmp.Re;
      if (TmpIndex < index)
	{
	  Tmp = vSource[TmpIndex] * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	  TmpSum.Re -= TmpCoefficient * Tmp.Im;
	  TmpSum.Im += TmpCoefficient * Tmp.Re;
	}
    }
  TmpIndex = chain->SpiSmjSzk(k, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex].Re -= TmpCoefficient * Tmp.Im;
      vDestination[TmpIndex].Im += TmpCoefficient * Tmp.Re;
      if (TmpIndex < index)
	{
	  Tmp = vSource[TmpIndex] * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	  TmpSum.Re += TmpCoefficient * Tmp.Im;
	  TmpSum.Im -= TmpCoefficient * Tmp.Re;
	}
    }
  TmpIndex = chain->SpiSmjSzk(j, i, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (TmpIndex <= index)
    {
      TmpCoefficient *= this->HalfJcFactor;
      Complex Tmp = Coefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
      vDestination[TmpIndex].Re += TmpCoefficient * Tmp.Im;
      vDestination[TmpIndex].Im -= TmpCoefficient * Tmp.Re;
      if (TmpIndex < index)
	{
	  Tmp = vSource[TmpIndex] * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
	  TmpSum.Re -= TmpCoefficient * Tmp.Im;
	  TmpSum.Im += TmpCoefficient * Tmp.Re;
	}
    }
  vDestination[index] += TmpSum;
}

// core part of the FastMultiplication method
// 
// chain = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray  

void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateFastMultiplicationComponent(AbstractSpinChain* chain, int index, 
										      int* indexArray, Complex* coefficientArray, long& position)
{
  for (int j = 0; j < this->NbrSpinX; j++)
    {
      for (int k = 0; k < this->NbrSpinY; k++)
	{
	  int TmpIndex1 = this->GetLinearizedIndex(j, k);
	  int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
	  int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
	  int TmpIndex4 = this->GetSafeLinearizedIndex(j + 1, k + 1);
	  
	  this->CoreEvaluateFastMultiplicationComponentPowerHeisenberg(chain, TmpIndex1, TmpIndex2, index, indexArray, coefficientArray, position);
	  this->CoreEvaluateFastMultiplicationComponentPowerHeisenberg(chain, TmpIndex2, TmpIndex1, index, indexArray, coefficientArray, position);
	  this->CoreEvaluateFastMultiplicationComponentPowerHeisenberg(chain, TmpIndex1, TmpIndex3, index, indexArray, coefficientArray, position);
	  this->CoreEvaluateFastMultiplicationComponentPowerHeisenberg(chain, TmpIndex3, TmpIndex1, index, indexArray, coefficientArray, position);
	  
	  this->CoreEvaluateFastMultiplicationComponentChiral(chain, TmpIndex1, TmpIndex2, TmpIndex4, index, indexArray, coefficientArray, position);
	  this->CoreEvaluateFastMultiplicationComponentChiral(chain, TmpIndex2, TmpIndex4, TmpIndex3, index, indexArray, coefficientArray, position);
	  this->CoreEvaluateFastMultiplicationComponentChiral(chain, TmpIndex4, TmpIndex3, TmpIndex1, index, indexArray, coefficientArray, position);
	  this->CoreEvaluateFastMultiplicationComponentChiral(chain, TmpIndex3, TmpIndex1, TmpIndex2, index, indexArray, coefficientArray, position);
	}
    }
}


// core part of the PartialFastMultiplication for all Hamiltonian terms (S_i S_j)^n 
//
// chain = pointer to the Hilbert space
// i = linearized position of the first spin
// j = linearized position of the second spin
// index = index of the many-body state to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray  

void TwoDimensionalRRAnd2DTranslationHamiltonian::CoreEvaluateFastMultiplicationComponentPowerHeisenberg(AbstractSpinChain* chain, int i, int j, int index, 
													 int* indexArray, Complex* coefficientArray, long& position)
{
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int AbsoluteIndex = index + this->PrecalculationShift;
  if (this->HermitianSymmetryFlag == false)
    {
      int Dim = chain->GetHilbertSpaceDimension();
      int TmpIndex = chain->SmiSpj(i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  double TmpValue3 = chain->SziSzj(i, j, AbsoluteIndex);
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (this->J1Factor + ((this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue3)) * TmpValue2;
	  ++position;
	}
      TmpIndex = chain->SmiSpjSmkSpl(i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  double TmpValue3 = chain->SziSzj(i, j, AbsoluteIndex);
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
	  ++position;
	}
      TmpIndex = chain->SmiSpjSmkSpl(j, i, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  double TmpValue3 = chain->SziSzj(i, j, AbsoluteIndex);
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
	  ++position;
	}
      TmpIndex = chain->SziSzjSmkSpl(i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  double TmpValue3 = chain->SziSzj(i, j, AbsoluteIndex);
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
	  ++position;
	}
      
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  ++position;
	}
      
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  ++position;
	}
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, j, i, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  ++position;
	}
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, j, i, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  ++position;
	}
      
      TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  ++position;
	}
      TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, j, i, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  ++position;
	}
      
      TmpIndex = chain->SmiSpjSzkSzlSmmSpn(i, j, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  ++position;
	}
      TmpIndex = chain->SmiSpjSzkSzlSmmSpn(j, i, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  ++position;
	}
      
      TmpIndex = chain->SziSzjSzkSzlSmmSpn(i, j, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.5 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  ++position;
	}
    }
  else
    {
      int TmpIndex = chain->SmiSpj(i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  double TmpValue3 = chain->SziSzj(i, j, AbsoluteIndex);
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (this->J1Factor + ((this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue3)) * TmpValue2;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SmiSpjSmkSpl(i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  double TmpValue3 = chain->SziSzj(i, j, AbsoluteIndex);
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SmiSpjSmkSpl(j, i, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  double TmpValue3 = chain->SziSzj(i, j, AbsoluteIndex);
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SziSzjSmkSpl(i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  double TmpValue3 = chain->SziSzj(i, j, AbsoluteIndex);
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, j, i, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, j, i, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.125 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      
      TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, j, i, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      
      TmpIndex = chain->SmiSpjSzkSzlSmmSpn(i, j, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SmiSpjSzkSzlSmmSpn(j, i, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.25 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      
      TmpIndex = chain->SziSzjSzkSzlSmmSpn(i, j, i, j, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  indexArray[position] = TmpIndex;
	  coefficientArray[position] = (0.5 * this->J3Factor * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
    }
}

//  core part of the PartialFastMultiplication for a single term ( S_i (S_j ^ S_k) )
//
// chain = pointer to the Hilbert space
// i = linearized position of the first spin
// j = linearized position of the second spin
// k = linearized position of the second spin
// index = index of the many-body state to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray  

void TwoDimensionalRRAnd2DTranslationHamiltonian::CoreEvaluateFastMultiplicationComponentChiral(AbstractSpinChain* chain, int i, int j, int k, int index, 
												int* indexArray, Complex* coefficientArray, long& position)
{
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int AbsoluteIndex = index + this->PrecalculationShift;

  if (this->HermitianSymmetryFlag == false)
    {
      int Dim = chain->GetHilbertSpaceDimension();
      int TmpIndex = chain->SpiSmjSzk(i, j, k, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(i, k, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(j, k, i, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(k, j, i, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(k, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(j, i, k, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  ++position;
	}
    }
  else
    {
      int TmpIndex = chain->SpiSmjSzk(i, j, k, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(i, k, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(j, k, i, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(k, j, i, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(k, i, j, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
      TmpIndex = chain->SpiSmjSzk(j, i, k, AbsoluteIndex, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= AbsoluteIndex)
	{
	  TmpCoefficient *= this->HalfJcFactor;
	  indexArray[position] = TmpIndex;
	  coefficientArray[position].Re = TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Im;
	  coefficientArray[position].Im = -TmpCoefficient * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY].Re;
	  if (TmpIndex == AbsoluteIndex)
	    {
	      coefficientArray[position] *= 0.5;
	    }
	  ++position;
	}
    }
}


// core part of the PartialFastMultiplicationMemory
// 
// chain = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations  

void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateFastMultiplicationMemoryComponent(AbstractSpinChain* chain, int firstComponent, int lastComponent, long& memory)
{
  for (int i = firstComponent; i < lastComponent; ++i)
    {
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      int TmpIndex1 = this->GetLinearizedIndex(j, k);
	      int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
	      int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
	      int TmpIndex4 = this->GetSafeLinearizedIndex(j + 1, k + 1);
	      
	      this->CoreEvaluateFastMultiplicationMemoryComponentPowerHeisenberg(chain, TmpIndex1, TmpIndex2, i, this->NbrInteractionPerComponent, memory);
	      this->CoreEvaluateFastMultiplicationMemoryComponentPowerHeisenberg(chain, TmpIndex2, TmpIndex1, i, this->NbrInteractionPerComponent, memory);
	      this->CoreEvaluateFastMultiplicationMemoryComponentPowerHeisenberg(chain, TmpIndex1, TmpIndex3, i, this->NbrInteractionPerComponent, memory);
	      this->CoreEvaluateFastMultiplicationMemoryComponentPowerHeisenberg(chain, TmpIndex3, TmpIndex1, i, this->NbrInteractionPerComponent, memory);
	            
	      this->CoreEvaluateFastMultiplicationMemoryComponentChiral(chain, TmpIndex1, TmpIndex2, TmpIndex4, i, this->NbrInteractionPerComponent, memory);
	      this->CoreEvaluateFastMultiplicationMemoryComponentChiral(chain, TmpIndex2, TmpIndex4, TmpIndex3, i, this->NbrInteractionPerComponent, memory);
	      this->CoreEvaluateFastMultiplicationMemoryComponentChiral(chain, TmpIndex4, TmpIndex3, TmpIndex1, i, this->NbrInteractionPerComponent, memory);
	      this->CoreEvaluateFastMultiplicationMemoryComponentChiral(chain, TmpIndex3, TmpIndex1, TmpIndex2, i, this->NbrInteractionPerComponent, memory);
	    }
	}
    }
}

// core part of the PartialFastMultiplicationMemory for all Hamiltonian terms (S_i S_j)^n 
//
// chain = pointer to the Hilbert space
// i = linearized position of the first spin
// j = linearized position of the second spin
// index = index of the many-body state to act on
// nbrInteractionPerComponent = array that contains the number of interaction per component
// memory = reference on the amount of memory required for precalculations  

void TwoDimensionalRRAnd2DTranslationHamiltonian::CoreEvaluateFastMultiplicationMemoryComponentPowerHeisenberg(AbstractSpinChain* chain, int i, int j, int index, int* nbrInteractionPerComponent, long& memory)
{
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  if (this->HermitianSymmetryFlag == false)
    {
      int Dim = chain->GetHilbertSpaceDimension();
      int TmpIndex = chain->SmiSpj(i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSmkSpl(j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SziSzjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SmiSpjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSzkSzlSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SziSzjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
    }
  else
    {
      int TmpIndex = chain->SmiSpj(i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSmkSpl(j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SziSzjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSmkSplSmmSpn(j, i, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SziSzjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SmiSpjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SmiSpjSzkSzlSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      
      TmpIndex = chain->SziSzjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
    }
}

//  core part of the PartialFastMultiplicationMemory for a single term ( S_i (S_j ^ S_k) )
//
// chain = pointer to the Hilbert space
// i = linearized position of the first spin
// j = linearized position of the second spin
// k = linearized position of the second spin
// index = index of the many-body state to act on
// nbrInteractionPerComponent = array that contains the number of interaction per component
// memory = reference on the amount of memory required for precalculations  

void TwoDimensionalRRAnd2DTranslationHamiltonian::CoreEvaluateFastMultiplicationMemoryComponentChiral(AbstractSpinChain* chain, int i, int j, int k, int index, int* nbrInteractionPerComponent, long& memory)
{
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;

  if (this->HermitianSymmetryFlag == false)
    {
      int Dim = chain->GetHilbertSpaceDimension();
      int TmpIndex = chain->SpiSmjSzk(i, j, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(i, k, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(j, k, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(k, j, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(k, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(j, i, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex < Dim)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
    }
  else
    {
      int TmpIndex = chain->SpiSmjSzk(i, j, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(i, k, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(j, k, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(k, j, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(k, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
      TmpIndex = chain->SpiSmjSzk(j, i, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
      if (TmpIndex <= index)
	{
	  ++memory;
	  ++nbrInteractionPerComponent[index - this->PrecalculationShift];
	}
    }
}
