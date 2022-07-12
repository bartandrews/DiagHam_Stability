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


#include "Hamiltonian/SpinChainZ2InteractingHamiltonianWithTranslations.h"
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

SpinChainZ2InteractingHamiltonianWithTranslations::SpinChainZ2InteractingHamiltonianWithTranslations()
{
  this->Chain = 0;
  this->NbrSpin = 0;
  this->SzSzContributions = 0;
  this->Parities = 0;
  this->JFactor = 0.0;
  this->FFactor = 0;
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
// boundaryCondition = boundary condition to apply (1 for periodic, -1 for antiperiodic)

SpinChainZ2InteractingHamiltonianWithTranslations::SpinChainZ2InteractingHamiltonianWithTranslations(Spin1_2ChainWithTranslations* chain, int nbrSpin, 
												     double jFactor, double fFactor, double interactionStrength, 
												     double boundaryCondition)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->Parities = new double [this->Chain->GetHilbertSpaceDimension()];
  this->JFactor = -jFactor;
  this->FFactor = -  fFactor;
  this->InteractionStrength = -interactionStrength;
  this->BoundaryCondition = boundaryCondition;
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosineTable();
}

// destructor
//

SpinChainZ2InteractingHamiltonianWithTranslations::~SpinChainZ2InteractingHamiltonianWithTranslations() 
{
  delete[] this->SzSzContributions;
  delete[] this->ExponentialTable;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainZ2InteractingHamiltonianWithTranslations::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (Spin1_2ChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->Parities = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainZ2InteractingHamiltonianWithTranslations::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SpinChainZ2InteractingHamiltonianWithTranslations::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainZ2InteractingHamiltonianWithTranslations::ShiftHamiltonian (double shift)
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

ComplexVector& SpinChainZ2InteractingHamiltonianWithTranslations::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int NbrTranslations;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      Complex TmpValue = this->JFactor * vSource[i];
      cout << "TmpValue=" << TmpValue << endl;
      cout << "diag part = " << this->SzSzContributions[i] << endl;
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * (TmpValue * this->ExponentialTable[NbrTranslations]);
	      cout << "o1 : " <<  coef << " t=" << NbrTranslations << " " << pos << endl;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * (TmpValue * this->ExponentialTable[NbrTranslations]);
	      cout << "o2 : " <<  coef << " t=" << NbrTranslations << " " << pos << endl;
	    }
	  pos = this->Chain->SpiSpj(j, j + 1, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * (TmpValue * this->ExponentialTable[NbrTranslations]);
	      cout << "o3 : " <<  coef << " t=" << NbrTranslations << " " << pos << endl;
	    }
	  pos = this->Chain->SmiSmj(j, j + 1, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      vDestination[pos] += coef * (TmpValue * this->ExponentialTable[NbrTranslations]);
	      cout << "o4 : " <<  coef << " t=" << NbrTranslations << " " << pos << endl;
	    }
	}
      pos = this->Chain->SmiSpj(MaxPos, 0, i, coef, NbrTranslations);
      if (pos != dim)
	{
	  vDestination[pos] -= coef * this->Parities[pos] * this->BoundaryCondition * this->ExponentialTable[NbrTranslations] * TmpValue;
	  cout << "p1 : " <<  coef << " " <<  this->Parities[pos] << " " <<  this->BoundaryCondition << " t=" << NbrTranslations << " " << pos << endl;
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, coef, NbrTranslations);
      if (pos != dim)
	{
	  vDestination[pos] -= coef * this->Parities[pos] * this->BoundaryCondition * this->ExponentialTable[NbrTranslations] * TmpValue;
	  cout << "p2 : " <<  coef << " " <<  this->Parities[pos] << " " <<  this->BoundaryCondition << " t=" << NbrTranslations << " " << pos << endl;
	}
      pos = this->Chain->SpiSpj(MaxPos, 0, i, coef, NbrTranslations);
      if (pos != dim)
	{
	  vDestination[pos] -= coef * this->Parities[pos] * this->BoundaryCondition * this->ExponentialTable[NbrTranslations] * TmpValue;
	  cout << "p3 : " <<  coef << " " <<  this->Parities[pos] << " " <<  this->BoundaryCondition << " t=" << NbrTranslations << " " << pos << endl;
	}
      pos = this->Chain->SmiSmj(MaxPos, 0, i, coef, NbrTranslations);
      if (pos != dim)
	{
	  vDestination[pos] -= coef * this->Parities[pos] * this->BoundaryCondition  * this->ExponentialTable[NbrTranslations] * TmpValue;
	  cout << "p4 : " <<  coef << " " <<  this->Parities[pos] << " " <<  this->BoundaryCondition << " t=" << NbrTranslations << " " << pos << endl;
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

ComplexVector* SpinChainZ2InteractingHamiltonianWithTranslations::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
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
  int NbrTranslations;
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
	TmpValues[k] = this->JFactor * vSources[k][i];
      // J part of Hamiltonian      
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues[k] * this->ExponentialTable[NbrTranslations];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues[k] * this->ExponentialTable[NbrTranslations];
	    }
	  pos = this->Chain->SpiSpj(j, j + 1, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues[k] * this->ExponentialTable[NbrTranslations];
	    }
	  pos = this->Chain->SmiSmj(j, j + 1, i, coef, NbrTranslations);
	  if (pos != dim)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][pos] += coef * TmpValues[k] * this->ExponentialTable[NbrTranslations];
	    }
	}
      pos = this->Chain->SmiSpj(MaxPos, 0, i, coef, NbrTranslations);
      if (pos != dim)
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValues[k] * this->ExponentialTable[NbrTranslations];
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, coef, NbrTranslations);
      if (pos != dim)
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValues[k] * this->ExponentialTable[NbrTranslations];
	}
      pos = this->Chain->SpiSpj(MaxPos, 0, i, coef, NbrTranslations);
      if (pos != dim)
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValues[k] * this->ExponentialTable[NbrTranslations];
	}
      pos = this->Chain->SmiSmj(MaxPos, 0, i, coef, NbrTranslations);
      if (pos != dim)
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][pos] -= coef * this->Parities[pos] * this->BoundaryCondition * TmpValues[k] * this->ExponentialTable[NbrTranslations];
	}
    }
  delete[] TmpValues;
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void SpinChainZ2InteractingHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      double Tmp = 0.0;
      for (int j = 1; j < this->NbrSpin; ++j)
	{
	  Tmp += this->Chain->SziSzj(j - 1, j, i);
	}
      this->SzSzContributions[i] = this->InteractionStrength * 4.0 * Tmp + this->FFactor * this->Chain->TotalSz(i);
      this->Parities[i] = 1.0 - 2.0 * ((double) this->Chain->Parity(i));  
    }
  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] += 4.0 * this->InteractionStrength * this->Chain->SziSzj(this->NbrSpin - 1, 0, i);
    }
}

// evaluate all cosine/sine that are needed when computing matrix elements
//

void SpinChainZ2InteractingHamiltonianWithTranslations::EvaluateCosineTable()
{
  this->ExponentialTable = new Complex [this->NbrSpin];
  double Coef = 2.0 * M_PI / ((double) this->NbrSpin) * ((double) this->Chain->GetMomentum());
  for (int i = 0; i < this->NbrSpin ; ++i)
    {
      this->ExponentialTable[i] = Phase((Coef * ((double) i)));
    }
}
