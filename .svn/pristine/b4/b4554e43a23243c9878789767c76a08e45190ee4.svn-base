////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of Potts 3 chain hamiltonian                     //
//                     for the dual O'Brien-Fendley model                     //
//                                                                            //
//                        last modification : 02/12/2019                      //
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


#include "Hamiltonian/Potts3ChainDualOBrienFendleyHamiltonian.h"
#include "Vector/RealVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "GeneralTools/StringTools.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

Potts3ChainDualOBrienFendleyHamiltonian::Potts3ChainDualOBrienFendleyHamiltonian()
{
}

// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// periodicFlag = true if the chain is periodic
// memory = amount of memory that can be used from precalculations (in bytes)

Potts3ChainDualOBrienFendleyHamiltonian::Potts3ChainDualOBrienFendleyHamiltonian(Potts3Chain* chain, int nbrSpin, bool periodicFlag, long memory)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->ReducedNbrSpin = this->NbrSpin -1;
  this->PeriodicFlag = periodicFlag;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->FastMultiplicationFlag = false;
  this->NbrInteractionPerComponent = this->ReducedNbrSpin;
  if (this->PeriodicFlag == true)
    {
      this->NbrInteractionPerComponent++;
    }
  long RequireMemory = 2l * (((long) this->Chain->GetHilbertSpaceDimension()) * 
			     (this->NbrInteractionPerComponent * sizeof(int) + sizeof(double)));
  cout << "Precalculations requires ";
  // PrintMemorySize (cout, RequireMemory) << endl;
  // if (memory > RequireMemory)
  //   {
  //     cout << "use fast multiplication" << endl;
  //     this->FastMultiplicationFlag = true;
  //   }
  // else
  //   {
  //     cout << "cannot use fast multiplication" << endl;
  //     this->FastMultiplicationFlag = false;
  //   }
  this->FastMultiplicationFlag = false;
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

Potts3ChainDualOBrienFendleyHamiltonian::~Potts3ChainDualOBrienFendleyHamiltonian() 
{
  delete[] this->SzSzContributions;
  if (this->FastMultiplicationFlag == true)
    delete[] this->InteractionPerComponentIndex;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void Potts3ChainDualOBrienFendleyHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  if (this->FastMultiplicationFlag == true)
    delete[] this->InteractionPerComponentIndex;
  this->Chain = (Potts3Chain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* Potts3ChainDualOBrienFendleyHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int Potts3ChainDualOBrienFendleyHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void Potts3ChainDualOBrienFendleyHamiltonian::ShiftHamiltonian (double shift)
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

RealVector& Potts3ChainDualOBrienFendleyHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									    int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == true)
    {
      long Position = (2l * this->NbrInteractionPerComponent) * firstComponent;
      if (this->PeriodicFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      double TmpValue = vSource[i];
	      vDestination[i] += this->SzSzContributions[i] * TmpValue;
	      for (int j = 0; j < this->ReducedNbrSpin; ++j)
		{
		  // vDestination[this->InteractionPerComponentIndex[Position++]] += this->JFullFactor * TmpValue;
		  // vDestination[this->InteractionPerComponentIndex[Position++]] += ConjugateProduct(this->JFullFactor, TmpValue);
		}
	      // vDestination[this->InteractionPerComponentIndex[Position++]] += this->BoundaryFactor * TmpValue;
	      // vDestination[this->InteractionPerComponentIndex[Position++]] += ConjugateProduct(this->BoundaryFactor, TmpValue);
	    }
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      double TmpValue = vSource[i];
	      vDestination[i] += this->SzSzContributions[i] * TmpValue;
	      for (int j = 0; j < this->ReducedNbrSpin; ++j)
		{
		  // vDestination[this->InteractionPerComponentIndex[Position++]] += this->JFullFactor * TmpValue;
		  // vDestination[this->InteractionPerComponentIndex[Position++]] += ConjugateProduct(this->JFullFactor, TmpValue);
		}
	    }
	}
    }
  else
    {
      int dim = this->Chain->GetHilbertSpaceDimension();
      double Coefficient;
      int pos;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  double TmpValue = vSource[i];
	  vDestination[i] += this->SzSzContributions[i] * TmpValue;
	  for (int j = 0; j < this->ReducedNbrSpin; ++j)
	    {
	      pos = this->Chain->SmiSpj(j, j + 1, i, Coefficient);
	      if (pos != dim)
		{
		  this->Chain->Szi(j + 1, i, Coefficient);
		  if (Coefficient == 1.0)
		    {
		      vDestination[pos] += 2.0 * TmpValue;
		    }
		  else
		    {
		      vDestination[pos] -= TmpValue;
		    }
		  this->Chain->Szi(j, pos, Coefficient);
		  if (Coefficient == 1.0)
		    {
		      vDestination[pos] += 2.0 * TmpValue;
		    }
		  else
		    {
		      vDestination[pos] -= TmpValue;
		    }
		}
	      pos = this->Chain->SmiSpj(j + 1, j, i, Coefficient);
	      if (pos != dim)
		{
		  this->Chain->Szi(j + 1, i, Coefficient);
		  if (Coefficient == -1.0)
		    {
		      vDestination[pos] += 2.0 * TmpValue;
		    }
		  else
		    {
		      vDestination[pos] -= TmpValue;
		    }
		  this->Chain->Szi(j, pos, Coefficient);
		  if (Coefficient == -1.0)
		    {
		      vDestination[pos] += 2.0 * TmpValue;
		    }
		  else
		    {
		      vDestination[pos] -= TmpValue;
		    }
		}
	    }
	}
      if (this->PeriodicFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      double TmpValue = vSource[i];
	      pos = this->Chain->SmiSpj(this->ReducedNbrSpin, 0, i, Coefficient);
	      if (pos != dim)
		{
		  this->Chain->Szi(0, i, Coefficient);
		  if (Coefficient == 1.0)
		    {
		      vDestination[pos] += 2.0 * TmpValue;
		    }
		  else
		    {
		      vDestination[pos] -= TmpValue;
		    }
		  this->Chain->Szi(this->ReducedNbrSpin, pos, Coefficient);
		  if (Coefficient == 1.0)
		    {
		      vDestination[pos] += 2.0 * TmpValue;
		    }
		  else
		    {
		      vDestination[pos] -= TmpValue;
		    }
		}
	      pos = this->Chain->SmiSpj(0, this->ReducedNbrSpin, i, Coefficient);
	      if (pos != dim)
		{
		  this->Chain->Szi(0, i, Coefficient);
		  if (Coefficient == -1.0)
		    {
		      vDestination[pos] += 2.0 * TmpValue;
		    }
		  else
		    {
		      vDestination[pos] -= TmpValue;
		    }
		  this->Chain->Szi(this->ReducedNbrSpin, pos, Coefficient);
		  if (Coefficient == -1.0)
		    {
		      vDestination[pos] += 2.0 * TmpValue;
		    }
		  else
		    {
		      vDestination[pos] -= TmpValue;
		    }
		}
	    }
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

RealVector*  Potts3ChainDualOBrienFendleyHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  double* TmpValues = new double[nbrVectors];
  if (this->FastMultiplicationFlag == true)
    {
      long Position = (2l * this->NbrInteractionPerComponent) * firstComponent;
      int pos1;
      int pos2;
      if (this->PeriodicFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  TmpValues[k] = vSources[k][i];
		  vDestinations[k][i] += this->SzSzContributions[i] * TmpValues[k];
		}
	      for (int j = 0; j < this->ReducedNbrSpin; ++j)
		{
		  pos1 = this->InteractionPerComponentIndex[Position++];
		  pos2 = this->InteractionPerComponentIndex[Position++];
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      // vDestinations[k][pos1] += this->JFullFactor * TmpValues[k];
		      // vDestinations[k][pos2] += ConjugateProduct(this->JFullFactor, TmpValues[k]);
		    }
		}
	      pos1 = this->InteractionPerComponentIndex[Position++];
	      pos2 = this->InteractionPerComponentIndex[Position++];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  // vDestinations[k][pos1] += this->BoundaryFactor * TmpValues[k];
		  // vDestinations[k][pos2] += ConjugateProduct(this->BoundaryFactor, TmpValues[k]);
		}
	    }
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  TmpValues[k] = vSources[k][i];
		  vDestinations[k][i] += this->SzSzContributions[i] * TmpValues[k];
		}
	      for (int j = 0; j < this->ReducedNbrSpin; ++j)
		{
		  pos1 = this->InteractionPerComponentIndex[Position++];
		  pos2 = this->InteractionPerComponentIndex[Position++];
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      // vDestinations[k][pos1] += this->JFullFactor * TmpValues[k];
		      // vDestinations[k][pos2] += ConjugateProduct(this->JFullFactor, TmpValues[k]);
		    }
		}
	    }
	}
    }
  else
    {
      int dim = this->Chain->GetHilbertSpaceDimension();
      double Coefficient;
      int pos;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      TmpValues[k] = vSources[k][i];
	      vDestinations[k][i] += this->SzSzContributions[i] * TmpValues[k];
	    }
	  for (int j = 0; j < this->ReducedNbrSpin; ++j)
	    {
	      pos = this->Chain->SmiSpj(j, j + 1, i, Coefficient);
	      if (pos != dim)
		{
		  this->Chain->Szi(j + 1, i, Coefficient);
		  if (Coefficient == 1.0)
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += 2.0 * TmpValues[k];
			}
		    }
		  else
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] -= TmpValues[k];
			}
		    }
		  this->Chain->Szi(j, pos, Coefficient);
		  if (Coefficient == 1.0)
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += 2.0 * TmpValues[k];
			}
		    }
		  else
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] -= TmpValues[k];
			}
		    }
		}
	      pos = this->Chain->SmiSpj(j + 1, j, i, Coefficient);
	      if (pos != dim)
		{
		  this->Chain->Szi(j + 1, i, Coefficient);
		  if (Coefficient == -1.0)
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += 2.0 * TmpValues[k];
			}
		    }
		  else
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] -= TmpValues[k];
			}
		    }
		  this->Chain->Szi(j, pos, Coefficient);
		  if (Coefficient == -1.0)
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += 2.0 * TmpValues[k];
			}
		    }
		  else
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] -= TmpValues[k];
			}
		    }
		}
	    }
	}
      if (this->PeriodicFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  TmpValues[k] = vSources[k][i];
		}
	      pos = this->Chain->SmiSpj(this->ReducedNbrSpin, 0, i, Coefficient);
	      if (pos != dim)
		{
		  this->Chain->Szi(0, i, Coefficient);
		  if (Coefficient == 1.0)
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += 2.0 * TmpValues[k];
			}
		    }
		  else
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] -= TmpValues[k];
			}
		    }
		  this->Chain->Szi(this->ReducedNbrSpin, pos, Coefficient);
		  if (Coefficient == 1.0)
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += 2.0 * TmpValues[k];
			}
		    }
		  else
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] -= TmpValues[k];
			}
		    }
		}
	      pos = this->Chain->SmiSpj(0, this->ReducedNbrSpin, i, Coefficient);
	      if (pos != dim)
		{
		  this->Chain->Szi(0, i, Coefficient);
		  if (Coefficient == -1.0)
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += 2.0 * TmpValues[k];
			}
		    }
		  else
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] -= TmpValues[k];
			}
		    }
		  this->Chain->Szi(this->ReducedNbrSpin, pos, Coefficient);
		  if (Coefficient == -1.0)
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] += 2.0 * TmpValues[k];
			}
		    }
		  else
		    {
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][pos] -= TmpValues[k];
			}
		    }
		}
	    }
	}
    }
  delete[] TmpValues;
  return vDestinations;
}


 // evaluate diagonal matrix elements
// 

void Potts3ChainDualOBrienFendleyHamiltonian::EvaluateDiagonalMatrixElements()
{
  int Dimension = this->Chain->GetHilbertSpaceDimension();
  double Coefficient;
  double CoefficientCosine1 = 1.0;
  double CoefficientCosine2 = (cos (2.0 * M_PI / 3.0) - 1.0);
  double CoefficientSine1 = 0.0;
  double CoefficientSine2 = sin (2.0 * M_PI / 3.0);
  for (int i = 0; i < Dimension; i++)
    {
      double Tmp = 0.0;
      for (int j = 0; j < this->NbrSpin; j++)
	{
	  this->Chain->Szi(j, i, Coefficient);
	  Tmp += -2.0 + 3.0 * (Coefficient * Coefficient);
	}
      this->SzSzContributions[i] = Tmp;
    }
  if (this->FastMultiplicationFlag == true)
    {
      long Position = 0l;
      double Coef = 0.0;
      this->InteractionPerComponentIndex = new int[2 * this->NbrInteractionPerComponent * Dimension];
      for (int i = 0; i < Dimension; ++i)
	{
	  for (int j = 0; j < this->ReducedNbrSpin; ++j)
	    {
	      this->InteractionPerComponentIndex[Position++] = this->Chain->SmiSpj(j, j + 1, i, Coef);
	      this->InteractionPerComponentIndex[Position++] = this->Chain->SmiSpj(j + 1, j, i, Coef);
	    }
	  if (this->PeriodicFlag == true)
	    {
	      this->InteractionPerComponentIndex[Position++] = this->Chain->SmiSpj(0, this->ReducedNbrSpin, i, Coef);
	      this->InteractionPerComponentIndex[Position++] = this->Chain->SmiSpj(this->ReducedNbrSpin, 0, i, Coef);
	    }
	}      
    }
}

