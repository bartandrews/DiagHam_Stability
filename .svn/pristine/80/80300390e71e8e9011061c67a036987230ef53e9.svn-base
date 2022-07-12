////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                              class of XYZ chain                            //
//                                                                            //
//                        last modification : 16/12/2013                      //
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


#include "Hamiltonian/SpinChainXYZHamiltonian.h"
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

SpinChainXYZHamiltonian::SpinChainXYZHamiltonian()
{
}

// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// jxFactor = coupling along the x direction
// jyFactor = coupling along the y direction
// jzFactor = coupling along the z direction
// hFactor = Zeeman term 
// boundaryCondition = boundary condition to apply (0 for open chain, 1 for periodic, -1 for antiperiodic)

SpinChainXYZHamiltonian::SpinChainXYZHamiltonian(Spin1_2Chain* chain, int nbrSpin, 
						 double jxFactor, double jyFactor, double jzFactor, double hFactor, 
						 double boundaryCondition)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->FixedParityFlag = false;
  this->Parities = new double [this->Chain->GetHilbertSpaceDimension()];
  this->JxFactor = -jxFactor;
  this->JyFactor = -jyFactor;
  this->JzFactor = -jzFactor;
  this->FFactors = new double[this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; ++i)
    this->FFactors[i] = hFactor;
  this->BoundaryCondition = boundaryCondition;
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

SpinChainXYZHamiltonian::~SpinChainXYZHamiltonian() 
{
  if (this->FixedParityFlag == false)
    delete[] this->Parities;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainXYZHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (Spin1_2Chain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->Parities = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainXYZHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SpinChainXYZHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainXYZHamiltonian::ShiftHamiltonian (double shift)
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

RealVector& SpinChainXYZHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
      double TmpValue1 = (this->JxFactor - this->JyFactor) * vSource[i];
      double TmpValue2 = (this->JxFactor + this->JyFactor) * vSource[i];
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * TmpValue2;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * TmpValue2;
	    }
	  pos = this->Chain->SpiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * TmpValue1;
	    }
	  pos = this->Chain->SmiSmj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * TmpValue1;
	    }
	}
      if (this->BoundaryCondition != 0.0)
	{
	  pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValue2;
	    }
	  pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValue2;
	    }
	  pos = this->Chain->SpiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValue1;
	    }
	  pos = this->Chain->SmiSmj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValue1;
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

RealVector* SpinChainXYZHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  double* TmpValues1 = new double[nbrVectors];
  double* TmpValues2 = new double[nbrVectors];
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
	  TmpValues1[k] = (this->JxFactor - this->JyFactor) * vSources[k][i];
	  TmpValues2[k] = (this->JxFactor + this->JyFactor) * vSources[k][i];
	}
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues2[k];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues2[k];
	    }
	  pos = this->Chain->SpiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues1[k];
	    }
	  pos = this->Chain->SmiSmj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues1[k];
	    }
	}
      if (this->BoundaryCondition != 0.0)
	{
	  pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValues2[k];
	    }
	  pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValues2[k];
	    }
	  pos = this->Chain->SpiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValues1[k];
	    }
	  pos = this->Chain->SmiSmj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValues1[k];
	    }
	}
    }
  delete[] TmpValues1;
  delete[] TmpValues2;
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void SpinChainXYZHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      double Tmp = 0.0;
      double Tmp2;
      double Tmp3;
      this->Chain->Szi(0, i, Tmp2);
      Tmp2 *= this->FFactors[0]; 
      for (int j = 1; j < this->NbrSpin; ++j)
	{
	  Tmp += this->Chain->SziSzj(j - 1, j, i);
	  this->Chain->Szi(j, i, Tmp3);
	  Tmp2 += this->FFactors[j] * Tmp3;
	}
      this->SzSzContributions[i] = this->JzFactor * 4.0 * Tmp + Tmp2;
      this->Parities[i] = 1.0 - 2.0 * ((double) this->Chain->Parity(i));  
    }
  if (this->BoundaryCondition != 0.0)
    {
      for (int i = 0; i < dim; i++)
	{
	  this->SzSzContributions[i] += 4.0 * this->JzFactor * this->Chain->SziSzj(this->NbrSpin - 1, 0, i);
	}
    }
}

