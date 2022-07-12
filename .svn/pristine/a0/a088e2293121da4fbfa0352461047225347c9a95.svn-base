////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of Z2 interacting chain                       //
//                                                                            //
//                        last modification : 11/12/2013                      //
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


#include "Hamiltonian/SpinChainZ2InteractingHamiltonian.h"
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

SpinChainZ2InteractingHamiltonian::SpinChainZ2InteractingHamiltonian()
{
  this->Chain = 0;
  this->NbrSpin = 0;
  this->SzSzContributions = 0;
  this->Parities = 0;
  this->JFactor = 0.0;
  this->FFactors = 0;
  this->InteractionStrength = 0.0;
  this->BoundaryCondition = 0.0;
}

// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// jFactor = coupling along the z direction
// fFactor = Zeeman term
// interactionStrength = coupling along the x direction
// boundaryCondition = boundary condition to apply (0 for open chain, 1 for periodic, -1 for antiperiodic)

SpinChainZ2InteractingHamiltonian::SpinChainZ2InteractingHamiltonian(Spin1_2Chain* chain, int nbrSpin, 
								     double jFactor, double fFactor, double interactionStrength, 
								     double boundaryCondition)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->Parities = new double [this->Chain->GetHilbertSpaceDimension()];
  this->JFactor = -jFactor;
  this->FFactors = new double[this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; ++i)
    this->FFactors[i] = - 2.0 * fFactor;
  this->InteractionStrength = -interactionStrength;
  this->BoundaryCondition = boundaryCondition;
  this->EvaluateDiagonalMatrixElements();
}

// constructor
//
// chain = pointer to the Hilbert space of the system
// nbrSpin = number of spins
// jFactor = coupling along the z direction
// fFactors = Zeeman term for each site
// interactionStrength = coupling along the x direction
// boundaryCondition = boundary condition to apply (0 for open chain, 1 for periodic, -1 for antiperiodic)
SpinChainZ2InteractingHamiltonian::SpinChainZ2InteractingHamiltonian(Spin1_2Chain* chain, int nbrSpin, 
								     double jFactor, double* fFactors, double interactionStrength, 
								     double boundaryCondition)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->Parities = new double [this->Chain->GetHilbertSpaceDimension()];
  this->JFactor = -jFactor;
  this->FFactors = new double[this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; ++i)
    this->FFactors[i] = - 2.0 * fFactors[i];
  this->InteractionStrength = -interactionStrength;
  this->BoundaryCondition = boundaryCondition;
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

SpinChainZ2InteractingHamiltonian::~SpinChainZ2InteractingHamiltonian() 
{
  delete[] this->FFactors;
  delete[] this->SzSzContributions;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainZ2InteractingHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
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

AbstractHilbertSpace* SpinChainZ2InteractingHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

SpinChainZ2InteractingHamiltonian& SpinChainZ2InteractingHamiltonian::SetChain(AbstractSpinChain* chain)
{  
  delete[] this->SzSzContributions;
  delete[] this->Parities;
  this->Chain =  (Spin1_2Chain*) chain;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->Parities = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SpinChainZ2InteractingHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainZ2InteractingHamiltonian::ShiftHamiltonian (double shift)
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

RealVector& SpinChainZ2InteractingHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
      double TmpValue = this->JFactor * vSource[i];
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * TmpValue;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * TmpValue;
	    }
	  pos = this->Chain->SpiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * TmpValue;
	    }
	  pos = this->Chain->SmiSmj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * TmpValue;
	    }
	}
      if (this->BoundaryCondition != 0.0)
	{
	  pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValue;
	    }
	  pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValue;
	    }
	  pos = this->Chain->SpiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValue;
	    }
	  pos = this->Chain->SmiSmj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValue;
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

RealVector* SpinChainZ2InteractingHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
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
	TmpValues[k] = this->JFactor * vSources[k][i];
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues[k];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues[k];
	    }
	  pos = this->Chain->SpiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues[k];
	    }
	  pos = this->Chain->SmiSmj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues[k];
	    }
	}
      if (this->BoundaryCondition != 0.0)
	{
	  pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValues[k];
	    }
	  pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValues[k];
	    }
	  pos = this->Chain->SpiSpj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValues[k];
	    }
	  pos = this->Chain->SmiSmj(MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValues[k];
	    }
	}
    }
  delete[] TmpValues;
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void SpinChainZ2InteractingHamiltonian::EvaluateDiagonalMatrixElements()
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
      this->SzSzContributions[i] = this->InteractionStrength * 4.0 * Tmp + Tmp2;
      this->Parities[i] = 1.0 - 2.0 * ((double) this->Chain->Parity(i));  
    }
  if (this->BoundaryCondition != 0.0)
    {
      for (int i = 0; i < dim; i++)
	{
	  this->SzSzContributions[i] += 4.0 * this->InteractionStrength * this->Chain->SziSzj(this->NbrSpin - 1, 0, i);
	}
    }
}

