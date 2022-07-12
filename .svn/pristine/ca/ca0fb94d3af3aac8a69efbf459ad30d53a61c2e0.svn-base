////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of spin chain generalized AKLT hamiltonian with translations     //
//             at inversion symmetric points using local projectors           //
//                      on the spin 3 and spin 4 sectors                      //
//                                                                            //
//                        last modification : 12/11/2016                      //
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


#include "Hamiltonian/SpinChainAKLTP3P4RealHamiltonianWithTranslations.h"
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

SpinChainAKLTP3P4RealHamiltonianWithTranslations::SpinChainAKLTP3P4RealHamiltonianWithTranslations()
{
}

// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// p3Factor = numerical factor in front of the spin 3 projector
// p4Factor = numerical factor in front of the spin 4 projector

SpinChainAKLTP3P4RealHamiltonianWithTranslations::SpinChainAKLTP3P4RealHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double p3Factor, double p4Factor)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->J = 1.0;
  this->HalfJ = this->J * 0.5;
  this->Jz = this->J;
  this->SquareFactor = 0.0;
  this->P3Factor = p3Factor;
  this->P4Factor = p4Factor;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// destructor
//

SpinChainAKLTP3P4RealHamiltonianWithTranslations::~SpinChainAKLTP3P4RealHamiltonianWithTranslations() 
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

RealVector& SpinChainAKLTP3P4RealHamiltonianWithTranslations::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									  int firstComponent, int nbrComponent)
{
  double Coef;
  double Coef2;
  double Coef3;
  int NbrTranslation;
  int NbrTranslation2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoefficients[25];
  int TmpNbrTranslations[25];
  int TmpIndices[25];
  for (int i = firstComponent; i < Last; i++)
    {
      double TmpValue = vSource[i];
      vDestination[i] += this->SzSzContributions[i] * TmpValue;
      for (int j = 0; j < MaxPos; j++)
	{
	  int Tmp = this->Chain->Spin3Projector(j, j + 1, i, TmpIndices, TmpCoefficients, TmpNbrTranslations);
	  for (int k = 0; k < Tmp; ++k)
	    {
	      vDestination[TmpIndices[k]] +=  (this->P3Factor * TmpCoefficients[k]) * this->ExponentialTable[TmpNbrTranslations[k]] * TmpValue;
	    }
	  Tmp = this->Chain->Spin4Projector(j, j + 1, i, TmpIndices, TmpCoefficients, TmpNbrTranslations);
	  for (int k = 0; k < Tmp; ++k)
	    {
	      vDestination[TmpIndices[k]] +=  (this->P4Factor * TmpCoefficients[k]) * this->ExponentialTable[TmpNbrTranslations[k]] * TmpValue;
	    }
	}    
      int Tmp = this->Chain->Spin3Projector(MaxPos, 0, i, TmpIndices, TmpCoefficients, TmpNbrTranslations);
      for (int k = 0; k < Tmp; ++k)
	{
	  vDestination[TmpIndices[k]] +=  (this->P3Factor * TmpCoefficients[k]) * this->ExponentialTable[TmpNbrTranslations[k]] * TmpValue;
	}
      Tmp = this->Chain->Spin4Projector(MaxPos, 0, i, TmpIndices, TmpCoefficients, TmpNbrTranslations);
      for (int k = 0; k < Tmp; ++k)
	{
	  vDestination[TmpIndices[k]] +=  (this->P4Factor * TmpCoefficients[k]) * this->ExponentialTable[TmpNbrTranslations[k]] * TmpValue;
	}
    }
  return vDestination;
}
 
// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* SpinChainAKLTP3P4RealHamiltonianWithTranslations::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									       int firstComponent, int nbrComponent)
{
  double Coef;
  double Coef3;
  int NbrTranslation;
  int NbrTranslation2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double* TmpValues = new double[nbrVectors];
  double TmpCoef;
  double TmpCoefficients[25];
  int TmpNbrTranslations[25];
  int TmpIndices[25];
  for (int i = firstComponent; i < Last; i++)
    {
      for (int k = 0; k < nbrVectors; ++k)
	{
	  TmpValues[k] = vSources[k][i];
	  vDestinations[k][i] += this->SzSzContributions[i] * TmpValues[k];
	}
      for (int j = 0; j < MaxPos; j++)
	{
	  int Tmp = this->Chain->Spin3Projector(j, j + 1, i, TmpIndices, TmpCoefficients, TmpNbrTranslations);
	  for (int k = 0; k < Tmp; ++k)
	    {
	      TmpCoef =(this->P3Factor * TmpCoefficients[k]) * this->ExponentialTable[TmpNbrTranslations[k]];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][TmpIndices[k]] += TmpCoef * TmpValues[l];
		}
	    }
	  Tmp = this->Chain->Spin4Projector(j, j + 1, i, TmpIndices, TmpCoefficients, TmpNbrTranslations);
	  for (int k = 0; k < Tmp; ++k)
	    {
	      TmpCoef =(this->P4Factor * TmpCoefficients[k]) * this->ExponentialTable[TmpNbrTranslations[k]];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][TmpIndices[k]] += TmpCoef * TmpValues[l];
		}
	    }
	}    
      int Tmp = this->Chain->Spin3Projector(MaxPos, 0, i, TmpIndices, TmpCoefficients, TmpNbrTranslations);
      for (int k = 0; k < Tmp; ++k)
	{
	  TmpCoef =(this->P3Factor * TmpCoefficients[k]) * this->ExponentialTable[TmpNbrTranslations[k]];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      vDestinations[l][TmpIndices[k]] += TmpCoef * TmpValues[l];
	    }
	}
      Tmp = this->Chain->Spin4Projector(MaxPos, 0, i, TmpIndices, TmpCoefficients, TmpNbrTranslations);
      for (int k = 0; k < Tmp; ++k)
	{
	  TmpCoef =(this->P4Factor * TmpCoefficients[k]) * this->ExponentialTable[TmpNbrTranslations[k]];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      vDestinations[l][TmpIndices[k]] += TmpCoef * TmpValues[l];
	    }
	}
    }
  delete[] TmpValues;
  return vDestinations;
}
 
// evaluate diagonal matrix elements
// 

void SpinChainAKLTP3P4RealHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double Tmp;
  int MaxSite = this->NbrSpin - 1;
  for (int i = 0; i < dim; ++i)
    {
      this->SzSzContributions[i] = 0.0;
    }
}

