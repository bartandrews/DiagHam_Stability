////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of Potts 3 chain hamiltonian with translations            //
//                                                                            //
//                        last modification : 01/01/2014                      //
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


#include "Hamiltonian/Potts3ChainHamiltonianWithTranslations.h"
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


// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// momentum = momentum sector
// jFactor = magnitude of the coupling term
// phiJ = phase of the coupling term (in PI units)
// fFactor = magnitude of the Zeeman term
// phiF = phase of the Zeeman term (in PI units)
// boundaryCondition = type of boundary conditions if the chain is periodic (0 for 1, 1 for exp(i 2 \pi / 3), -1 1 for exp(i 2 \pi / 3)) 

Potts3ChainHamiltonianWithTranslations::Potts3ChainHamiltonianWithTranslations(Potts3ChainWithTranslations* chain, int nbrSpin, int momentum, double jFactor, double phiJ, double fFactor, double phiF, double boundaryCondition)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->Momentum = momentum;
  this->ReducedNbrSpin = this->NbrSpin -1;
  this->BoundaryCondition = boundaryCondition;
  this->PhiJ = phiJ;
  this->JFactor = -jFactor;
  this->JFullFactor = Polar(this->JFactor, 2.0 * M_PI * this->PhiJ);
  this->JFullFactors =  new Complex[this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; ++i)
    {
//       this->JFullFactors[i] = Polar(this->JFactor, 2.0 * M_PI * (this->PhiJ + 
// 								 ((((double) this->Momentum) + this->BoundaryCondition / 3.0) * ((double) i)) / ((double) this->NbrSpin)));
      this->JFullFactors[i] = Polar(this->JFactor, 2.0 * M_PI * this->PhiJ);
    }
  this->PhiF = phiF;
  this->FFactor = -fFactor;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->BoundaryFactors = new Complex [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

Potts3ChainHamiltonianWithTranslations::~Potts3ChainHamiltonianWithTranslations() 
{
  delete[] this->SzSzContributions;
  delete[] this->JFullFactors;
  delete[] this->BoundaryFactors;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void Potts3ChainHamiltonianWithTranslations::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  delete[] this->BoundaryFactors;
  this->Chain = (Potts3ChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->BoundaryFactors = new Complex [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* Potts3ChainHamiltonianWithTranslations::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int Potts3ChainHamiltonianWithTranslations::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void Potts3ChainHamiltonianWithTranslations::ShiftHamiltonian (double shift)
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

ComplexVector& Potts3ChainHamiltonianWithTranslations::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int NbrTranslations;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      Complex TmpValue = vSource[i];
      vDestination[i] += this->SzSzContributions[i] * TmpValue;
      for (int j = 0; j < this->ReducedNbrSpin; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      cout << i << " " << j << " 1 : " << Phase(2.0 * M_PI * (((double) this->Momentum) + (this->BoundaryCondition - 1.0) / 3.0) * ((double) NbrTranslations) / ((double) this->NbrSpin)) * this->JFullFactors[NbrTranslations] << " " << coef << endl;
	      vDestination[pos] += coef * Phase(2.0 * M_PI * (((double) this->Momentum) + (this->BoundaryCondition - 1.0) / 3.0) * ((double) NbrTranslations) / ((double) this->NbrSpin)) * this->JFullFactors[NbrTranslations] * TmpValue;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      cout << i << " " << j << " 2 : " << Phase(2.0 * M_PI * (((double) this->Momentum) + (this->BoundaryCondition - 1.0) / 3.0) * ((double) NbrTranslations) / ((double) this->NbrSpin)) * this->JFullFactors[NbrTranslations] << " " << coef << endl;
	      vDestination[pos] += coef * Phase(2.0 * M_PI * (((double) this->Momentum) + (this->BoundaryCondition - 1.0) / 3.0) * ((double) NbrTranslations) / ((double) this->NbrSpin)) * ConjugateProduct(this->JFullFactors[NbrTranslations], TmpValue);
	    }
	}
      pos = this->Chain->SmiSpj(0, this->ReducedNbrSpin, i, coef, NbrTranslations);
      if (pos != dim)
	{
	      cout << i << " " << this->ReducedNbrSpin << " 3 : " << Phase(2.0 * M_PI * (((double) this->Momentum) + (this->BoundaryCondition - 1.0) / 3.0) * ((double) NbrTranslations) / ((double) this->NbrSpin)) * this->JFullFactors[NbrTranslations] << " " << coef << endl;
	  vDestination[pos] += coef * Phase(2.0 * M_PI * (((double) this->Momentum) + (this->BoundaryCondition - 1.0) / 3.0) * ((double) NbrTranslations) / ((double) this->NbrSpin)) * this->JFullFactors[NbrTranslations] * this->BoundaryFactors[pos] * TmpValue;
	}
      pos = this->Chain->SmiSpj(this->ReducedNbrSpin, 0, i, coef, NbrTranslations);
      if (pos != dim)
	{
	      cout << i << " " << this->ReducedNbrSpin << " 4 : " << Phase(2.0 * M_PI * (((double) this->Momentum) + (this->BoundaryCondition - 1.0) / 3.0) * ((double) NbrTranslations) / ((double) this->NbrSpin)) * this->JFullFactors[NbrTranslations] << " " << coef << endl;
	  vDestination[pos] += coef * Phase(2.0 * M_PI * (((double) this->Momentum) + (this->BoundaryCondition - 1.0) / 3.0) * ((double) NbrTranslations) / ((double) this->NbrSpin)) * ConjugateProduct(this->JFullFactors[NbrTranslations] * this->BoundaryFactors[pos], TmpValue);
	}
    }
  return vDestination;
}

// evaluate diagonal matrix elements
// 

void Potts3ChainHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int Dimension = this->Chain->GetHilbertSpaceDimension();
  double Coefficient;
  double CoefficientCosine1 = 2.0 * this->FFactor * cos (2.0 * M_PI * this->PhiF);
  double CoefficientCosine2 = (cos (2.0 * M_PI / 3.0) - 1.0);
  double CoefficientSine1 = 2.0 * this->FFactor * sin (2.0 * M_PI * this->PhiF);
  double CoefficientSine2 = sin (2.0 * M_PI / 3.0);
  for (int i = 0; i < Dimension; i++)
    {
      double Tmp = 0.0;
      for (int j = 0; j < this->NbrSpin; j++)
	{
	  this->Chain->Szi(j, i, Coefficient);

	  Tmp += CoefficientCosine1 * (CoefficientCosine2 * Coefficient * Coefficient  + 1.0) - (Coefficient * CoefficientSine1 * CoefficientSine2);
	}
      this->SzSzContributions[i] = Tmp;
//       this->BoundaryFactors[i] = Polar(this->JFactor,
// 				       2.0 * M_PI * (this->PhiJ + ((this->BoundaryCondition - 1.0 + this->Chain->QValue(i)) / 3.0)));
      this->BoundaryFactors[i] = Polar(2.0 * M_PI * (((this->BoundaryCondition - 1.0 + this->Chain->QValue(i)) / 3.0)));
    }
}

