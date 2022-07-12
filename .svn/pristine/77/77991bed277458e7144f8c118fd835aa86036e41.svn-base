////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin chain hamiltonian with translations             //
//                                                                            //
//                        last modification : 04/03/2002                      //
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


#include "Hamiltonian/SpinChainRealHamiltonianWithTranslations.h"
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

SpinChainRealHamiltonianWithTranslations::SpinChainRealHamiltonianWithTranslations()
{
}

// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// j = coupling constant between spin
// nnCoupling = term to add to ZZ nearest-neighbour interaction
// nnnCoupling = nearest-neighbour interaction in Z direction

SpinChainRealHamiltonianWithTranslations::SpinChainRealHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j, double nnCoupling, double nnnCoupling)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->J = j;
  this->HalfJ = this->J * 0.5;
  this->Jz = this->J;
  if (nnCoupling != 0)
    {
      this->Jz += nnCoupling;
      cout << "Adjusting ZZ interaction: " << this->Jz << endl;
    }

  this->NNNCoupling = nnnCoupling;
  if (this->NNNCoupling != 0)
    cout << "Adding NNN interaction: " << this->NNNCoupling << endl;

  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// destructor
//

SpinChainRealHamiltonianWithTranslations::~SpinChainRealHamiltonianWithTranslations() 
{
  delete[] this->SzSzContributions;
  delete[] this->CosinusTable;
  delete[] this->SinusTable;
  delete[] this->ExponentialTable;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainRealHamiltonianWithTranslations::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainRealHamiltonianWithTranslations::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

SpinChainRealHamiltonianWithTranslations& SpinChainRealHamiltonianWithTranslations::SetChain(AbstractSpinChainWithTranslations* chain)
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

int SpinChainRealHamiltonianWithTranslations::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainRealHamiltonianWithTranslations::ShiftHamiltonian (double shift)
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

RealVector& SpinChainRealHamiltonianWithTranslations::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									  int firstComponent, int nbrComponent)
{
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  for (int i = firstComponent; i < Last; i++)
    {
       vDestination[i] += this->SzSzContributions[i] * vSource[i];
       for (int j = 0; j < MaxPos; j++)
	 {
	   pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	   if (pos != this->Chain->GetHilbertSpaceDimension())
	     {
	       Coef *= this->HalfJ;
	       vDestination[pos] += Coef * vSource[i] * this->CosinusTable[NbrTranslation];
	     }
	   pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	   if (pos != this->Chain->GetHilbertSpaceDimension())
	     {
	       Coef *= this->HalfJ;
	       vDestination[pos] += Coef * vSource[i] * this->CosinusTable[NbrTranslation];
	     }
	 }    
       pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
       if (pos != this->Chain->GetHilbertSpaceDimension())
	 {
	   Coef *= this->HalfJ;
	   vDestination[pos] += Coef * vSource[i] * this->CosinusTable[NbrTranslation];
	 }
       pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
       if (pos != this->Chain->GetHilbertSpaceDimension())
	 {
	   Coef *= this->HalfJ;
	   vDestination[pos] += Coef * vSource[i] * this->CosinusTable[NbrTranslation];
	 }      
    }
  return vDestination;
}
 
// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void SpinChainRealHamiltonianWithTranslations::EvaluateCosinusTable()
{
  this->CosinusTable = new double [this->NbrSpin];
  this->SinusTable = new double [this->NbrSpin];
  this->ExponentialTable = new double [this->NbrSpin];
  double Coef = 2.0 * M_PI / ((double) this->NbrSpin) * ((double) this->Chain->GetMomentum());
  for (int i = 0; i < this->NbrSpin ; ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = 0.0;
      this->ExponentialTable[i] = this->CosinusTable[i];
    }
}

// evaluate diagonal matrix elements
// 

void SpinChainRealHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < (this->NbrSpin - 1); j++)
	{
	  this->SzSzContributions[i] += this->Chain->SziSzj(j, j + 1, i);
	}
      this->SzSzContributions[i] += this->Chain->SziSzj(this->NbrSpin - 1, 0, i);
      this->SzSzContributions[i] *= this->Jz;
      if (this->NNNCoupling != 0)
        {
          for (int j = 0; j < (this->NbrSpin - 2); j++)
   	     this->SzSzContributions[i] += this->NNNCoupling * this->Chain->SziSzj(j, j + 2, i);
   	  this->SzSzContributions[i] += this->NNNCoupling * this->Chain->SziSzj(this->NbrSpin - 2, 0, i);
   	  this->SzSzContributions[i] += this->NNNCoupling * this->Chain->SziSzj(this->NbrSpin - 1, 1, i);
        }
    }
}

