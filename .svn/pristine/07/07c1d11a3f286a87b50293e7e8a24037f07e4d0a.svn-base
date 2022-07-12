////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of spin chain hamiltonian with a triple product term       //
//                                                                            //
//                        last modification : 21/06/2015                      //
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


#include "Hamiltonian/SpinChainTripleProductHamiltonian.h"
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
// j = array containing coupling constants between spins
// chi = array containing the coupling constants of the triple product
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChainTripleProductHamiltonian::SpinChainTripleProductHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* chi, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
      this->Hz = 0;
  if (this->PeriodicBoundaryConditions == false)
    {
      this->J = new double [this->NbrSpin - 1];
      this->Jz = new double [this->NbrSpin - 1];
      this->HalfJ = new double [this->NbrSpin - 1];
      this->Chi = new double [this->NbrSpin - 1];
      this->HalfChi = new double [this->NbrSpin - 1];
      for (int i = 0; i < (this->NbrSpin - 1); i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = j[i];
	  this->HalfJ[i] = j[i] * 0.5;
	  this->Chi[i] = chi[i]; 
	  this->HalfChi[i] = chi[i] * 0.5; 
	}
    }
  else
    {
      this->J = new double [this->NbrSpin];
      this->Jz = new double [this->NbrSpin];
      this->HalfJ = new double [this->NbrSpin];
      this->Chi = new double [this->NbrSpin];
      this->HalfChi = new double [this->NbrSpin];
      for (int i = 0; i < this->NbrSpin; i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = j[i];
	  this->HalfJ[i] = j[i] * 0.5;
	  this->Chi[i] = chi[i]; 
	  this->HalfChi[i] = chi[i] * 0.5; 
	}
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// j = array containing coupling constants between spins along x and z
// jz = array containing coupling constants between spins along z
// chi = array containing the coupling constants of the triple product
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChainTripleProductHamiltonian::SpinChainTripleProductHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz, double* chi, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->Hz = 0;
  if (this->PeriodicBoundaryConditions == false)
    {
      this->J = new double [this->NbrSpin - 1];
      this->Jz = new double [this->NbrSpin - 1];
      this->HalfJ = new double [this->NbrSpin - 1];
      this->Chi = new double [this->NbrSpin - 1];
      this->HalfChi = new double [this->NbrSpin - 1];
      for (int i = 0; i < (this->NbrSpin - 1); i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = jz[i];
	  this->HalfJ[i] = j[i] * 0.5;
	  this->Chi[i] = chi[i]; 
	  this->HalfChi[i] = chi[i] * 0.5; 
	}
    }
  else
    {
      this->J = new double [this->NbrSpin];
      this->Jz = new double [this->NbrSpin];
      this->HalfJ = new double [this->NbrSpin];
      this->Chi = new double [this->NbrSpin];
      this->HalfChi = new double [this->NbrSpin];
      for (int i = 0; i < this->NbrSpin; i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = jz[i];
	  this->HalfJ[i] = j[i] * 0.5;
	  this->Chi[i] = chi[i]; 
	  this->HalfChi[i] = chi[i] * 0.5; 
	}
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// j = array containing the coupling constants between spins along x and z
// jz = array containing the coupling constants between spins along z
// chi = array containing the coupling constants of the triple product
// hz = array containing the amplitude of the Zeeman term along z
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

SpinChainTripleProductHamiltonian::SpinChainTripleProductHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* j, double* jz, double* chi, double* hz, bool periodicBoundaryConditions)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  if (this->PeriodicBoundaryConditions == false)
    {
      this->J = new double [this->NbrSpin - 1];
      this->Jz = new double [this->NbrSpin - 1];
      this->HalfJ = new double [this->NbrSpin - 1];
      this->Chi = new double [this->NbrSpin - 1];
      this->HalfChi = new double [this->NbrSpin - 1];
      for (int i = 0; i < (this->NbrSpin - 1); i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = jz[i];
	  this->HalfJ[i] = j[i] * 0.5;
	  this->Chi[i] = chi[i]; 
	  this->HalfChi[i] = chi[i] * 0.5; 
	}
    }
  else
    {
      this->J = new double [this->NbrSpin];
      this->Jz = new double [this->NbrSpin];
      this->HalfJ = new double [this->NbrSpin];
      this->Chi = new double [this->NbrSpin];
      this->HalfChi = new double [this->NbrSpin];
      for (int i = 0; i < this->NbrSpin; i++)
	{
	  this->J[i] = j[i];
	  this->Jz[i] = jz[i];
	  this->HalfJ[i] = j[i] * 0.5;
	  this->Chi[i] = chi[i]; 
	  this->HalfChi[i] = chi[i] * 0.5; 
	}
    }
  this->Hz = new double [this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; i++)
    {
      this->Hz[i] = hz[i];
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

SpinChainTripleProductHamiltonian::~SpinChainTripleProductHamiltonian() 
{
  delete[] this->J;
  delete[] this->Jz;
  delete[] this->HalfJ;
  delete[] this->SzSzContributions;
  delete[] this->Chi; 
  delete[] this->HalfChi; 
  if (this->Hz != 0)
    {
      delete[] this->Hz;  
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainTripleProductHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainTripleProductHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

SpinChainTripleProductHamiltonian& SpinChainTripleProductHamiltonian::SetChain(AbstractSpinChain* chain)
{  
  delete[] this->SzSzContributions;
  this->Chain = chain;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SpinChainTripleProductHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainTripleProductHamiltonian::ShiftHamiltonian (double shift)
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

ComplexVector& SpinChainTripleProductHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int MaxPos = this->NbrSpin - 2;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      Complex TmpValue = vSource[i];
      vDestination[i] += this->SzSzContributions[i] * TmpValue;
      // J part of Hamiltonian     
      for (int j = 0; j < MaxPos; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += (this->HalfJ[j] * coef) * TmpValue;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += (this->HalfJ[j] * coef) * TmpValue;
	    }
	  pos = this->Chain->SpiSmjSzk(j, j + 1, j + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im += (this->HalfChi[j] * coef) * TmpValue.Re;
	      vDestination[pos].Re -= (this->HalfChi[j] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(j + 2, j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im += (this->HalfChi[j] * coef) * TmpValue.Re;
	      vDestination[pos].Re -= (this->HalfChi[j] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(j + 1, j + 2, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im += (this->HalfChi[j] * coef) * TmpValue.Re;
	      vDestination[pos].Re -= (this->HalfChi[j] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(j + 1, j, j + 2, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im -= (this->HalfChi[j] * coef) * TmpValue.Re;
	      vDestination[pos].Re += (this->HalfChi[j] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(j, j + 2, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im -= (this->HalfChi[j] * coef) * TmpValue.Re;
	      vDestination[pos].Re += (this->HalfChi[j] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(j + 2, j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im -= (this->HalfChi[j] * coef) * TmpValue.Re;
	      vDestination[pos].Re += (this->HalfChi[j] * coef) * TmpValue.Im;
	    }
	}
      pos = this->Chain->SmiSpj(MaxPos, MaxPos + 1, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += (this->HalfJ[MaxPos] * coef) * TmpValue;
	}
      pos = this->Chain->SmiSpj(MaxPos + 1, MaxPos, i, coef);
      if (pos != dim)
	{
	  vDestination[pos] += (this->HalfJ[MaxPos] * coef) * TmpValue;
	}
      if (this->PeriodicBoundaryConditions == true)
	{
	  pos = this->Chain->SmiSpj(MaxPos + 1, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += (this->HalfJ[MaxPos + 1] * coef) * TmpValue;
	    }
	  pos = this->Chain->SmiSpj(0, MaxPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += (this->HalfJ[MaxPos + 1] * coef) * TmpValue;
	    }

	  pos = this->Chain->SpiSmjSzk(MaxPos, MaxPos + 1, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im += (this->HalfChi[MaxPos] * coef) * TmpValue.Re;
	      vDestination[pos].Re -= (this->HalfChi[MaxPos] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(0, MaxPos, MaxPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im += (this->HalfChi[MaxPos] * coef) * TmpValue.Re;
	      vDestination[pos].Re -= (this->HalfChi[MaxPos] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(MaxPos + 1, 0, MaxPos, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im += (this->HalfChi[MaxPos] * coef) * TmpValue.Re;
	      vDestination[pos].Re -= (this->HalfChi[MaxPos] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(MaxPos + 1, MaxPos, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im -= (this->HalfChi[MaxPos] * coef) * TmpValue.Re;
	      vDestination[pos].Re += (this->HalfChi[MaxPos] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(MaxPos, 0, MaxPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im -= (this->HalfChi[MaxPos] * coef) * TmpValue.Re;
	      vDestination[pos].Re += (this->HalfChi[MaxPos] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(0, MaxPos + 1, MaxPos, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im -= (this->HalfChi[MaxPos] * coef) * TmpValue.Re;
	      vDestination[pos].Re += (this->HalfChi[MaxPos] * coef) * TmpValue.Im;
	    }

	  pos = this->Chain->SpiSmjSzk(MaxPos + 1, 0, 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im += (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Re;
	      vDestination[pos].Re -= (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(1, MaxPos + 1, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im += (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Re;
	      vDestination[pos].Re -= (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(0, 1, MaxPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im += (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Re;
	      vDestination[pos].Re -= (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(0, MaxPos + 1, 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im -= (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Re;
	      vDestination[pos].Re += (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(MaxPos + 1, 1, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im -= (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Re;
	      vDestination[pos].Re += (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Im;
	    }
	  pos = this->Chain->SpiSmjSzk(1, 0, MaxPos + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos].Im -= (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Re;
	      vDestination[pos].Re += (this->HalfChi[MaxPos + 1] * coef) * TmpValue.Im;
	    }
	}
   }
  return vDestination;
}
 
// evaluate diagonal matrix elements
// 

void SpinChainTripleProductHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  int MaxSite = this->NbrSpin - 1;

  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < MaxSite; j++)
	{
	  this->SzSzContributions[i] += this->Jz[j] * this->Chain->SziSzj(j, j + 1, i);
	}
      if (this->PeriodicBoundaryConditions == true)
	{
	  this->SzSzContributions[i] += this->Jz[MaxSite] * this->Chain->SziSzj(MaxSite, 0, i);
	}
    }
  if (this->Hz != 0)
    {
      double Coefficient;
      // Sz part
      for (int i = 0; i < dim; i++)
	{
	  for (int j = 0; j < this->NbrSpin; j++)
	    {
	      this->Chain->Szi(j, i, Coefficient);
	      this->SzSzContributions[i] += this->Hz[j] * Coefficient;
	    }
	}      
    }
}

