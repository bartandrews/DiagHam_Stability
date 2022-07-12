////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of spin chain with only an external magnetic field         //
//                                                                            //
//                        last modification : 28/04/2015                      //
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


#include "Hamiltonian/SpinChainPureHFieldHamiltonian.h"
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

SpinChainPureHFieldHamiltonian::SpinChainPureHFieldHamiltonian()
{
}

// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// hxFactor = magnetic field along the x direction
// hyFactor = magnetic field along the y direction
// hzFactor = magnetic field along the z direction

SpinChainPureHFieldHamiltonian::SpinChainPureHFieldHamiltonian(Spin1_2Chain* chain, int nbrSpin, 
							       double hxFactor, double hyFactor, double hzFactor)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->HxFactor = hxFactor;
  this->HyFactor = hyFactor;
  this->HzFactor = hzFactor;
  this->HamiltonianShift = 0.0;
}

// destructor
//

SpinChainPureHFieldHamiltonian::~SpinChainPureHFieldHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainPureHFieldHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (Spin1_2Chain*) hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainPureHFieldHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SpinChainPureHFieldHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainPureHFieldHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& SpinChainPureHFieldHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Chain->GetHilbertSpaceDimension();
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  if ((this->HzFactor != 0.0) && (this->HxFactor == 0.0))
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  vDestination[i] += (this->HamiltonianShift + this->HzFactor * this->Chain->TotalSz(i)) * vSource[i];
	}
    }
  else
    {
      if ((this->HzFactor == 0.0) && (this->HxFactor != 0.0))
	{
 	  double Coefficient;
	  int TmpPos;
	  int Dim = this->Chain->GetHilbertSpaceDimension();
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      double Coefficient2 = this->HxFactor * vSource[i];
	      for (int j = 0; j < this->NbrSpin; ++j)
		{
		  TmpPos = this->Chain->Smi(j, i, Coefficient);
		  if (TmpPos != Dim)
		    {
		      vDestination[TmpPos] += Coefficient * Coefficient2;
		    }
		  TmpPos = this->Chain->Spi(j, i, Coefficient);
		  if (TmpPos != Dim)
		    {
		      vDestination[TmpPos] += Coefficient * Coefficient2;
		    }
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }	  
	}
      else
	{
 	  double Coefficient;
	  int TmpPos;
	  int Dim = this->Chain->GetHilbertSpaceDimension();
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      double Coefficient2 = this->HxFactor * vSource[i];
	      for (int j = 0; j < this->NbrSpin; ++j)
		{
		  TmpPos = this->Chain->Smi(j, i, Coefficient);
		  if (TmpPos != Dim)
		    {
		      vDestination[TmpPos] += Coefficient * Coefficient2;
		    }
		  TmpPos = this->Chain->Spi(j, i, Coefficient);
		  if (TmpPos != Dim)
		    {
		      vDestination[TmpPos] += Coefficient * Coefficient2;
		    }
		}
	      vDestination[i] += (this->HamiltonianShift + this->HzFactor * this->Chain->TotalSz(i)) * vSource[i];
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

RealVector* SpinChainPureHFieldHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
//   int dim = this->Chain->GetHilbertSpaceDimension();
//   double coef;
//   double coef2;
//   int pos;
//   int pos2;
//   int MaxPos = this->NbrSpin - 1;
//   double* TmpValues1 = new double[nbrVectors];
//   double* TmpValues2 = new double[nbrVectors];
//   for (int k = 0; k < nbrVectors; ++k)
//     {
//       RealVector& TmpSource = vSources[k];
//       RealVector& TmpDestination = vDestinations[k];
//       for (int i = firstComponent; i < LastComponent; ++i)
// 	{
// 	  TmpDestination[i] += this->SzSzContributions[i] * TmpSource[i];
// 	}
//     }
//   for (int i = firstComponent; i < LastComponent; ++i)
//     {
//       for (int k = 0; k < nbrVectors; ++k)
// 	{
// 	  TmpValues1[k] = (this->JxFactor - this->JyFactor) * vSources[k][i];
// 	  TmpValues2[k] = (this->JxFactor + this->JyFactor) * vSources[k][i];
// 	}
//       // J part of Hamiltonian      
//       for (int j = 0; j < MaxPos; ++j)
// 	{
// 	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
// 	  if (pos != dim)
// 	    {
// 	      for (int k = 0; k < nbrVectors; ++k)
// 		vDestinations[k][pos] += coef * TmpValues2[k];
// 	    }
// 	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
// 	  if (pos != dim)
// 	    {
// 	      for (int k = 0; k < nbrVectors; ++k)
// 		vDestinations[k][pos] += coef * TmpValues2[k];
// 	    }
// 	  pos = this->Chain->SpiSpj(j, j + 1, i, coef);
// 	  if (pos != dim)
// 	    {
// 	      for (int k = 0; k < nbrVectors; ++k)
// 		vDestinations[k][pos] += coef * TmpValues1[k];
// 	    }
// 	  pos = this->Chain->SmiSmj(j, j + 1, i, coef);
// 	  if (pos != dim)
// 	    {
// 	      for (int k = 0; k < nbrVectors; ++k)
// 		vDestinations[k][pos] += coef * TmpValues1[k];
// 	    }
// 	}
//       if (this->BoundaryCondition != 0.0)
// 	{
// 	  pos = this->Chain->SmiSpj(MaxPos, 0, i, coef);
// 	  if (pos != dim)
// 	    {
// 	      for (int k = 0; k < nbrVectors; ++k)
// 		vDestinations[k][pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValues2[k];
// 	    }
// 	  pos = this->Chain->SmiSpj(0, MaxPos, i, coef);
// 	  if (pos != dim)
// 	    {
// 	      for (int k = 0; k < nbrVectors; ++k)
// 		vDestinations[k][pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValues2[k];
// 	    }
// 	  pos = this->Chain->SpiSpj(MaxPos, 0, i, coef);
// 	  if (pos != dim)
// 	    {
// 	      for (int k = 0; k < nbrVectors; ++k)
// 		vDestinations[k][pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValues1[k];
// 	    }
// 	  pos = this->Chain->SmiSmj(MaxPos, 0, i, coef);
// 	  if (pos != dim)
// 	    {
// 	      for (int k = 0; k < nbrVectors; ++k)
// 		vDestinations[k][pos] += coef * this->Parities[pos] * this->BoundaryCondition * TmpValues1[k];
// 	    }
// 	}
//     }
//   delete[] TmpValues1;
//   delete[] TmpValues2;
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& SpinChainPureHFieldHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Chain->GetHilbertSpaceDimension();
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  if (this->HyFactor == 0.0)
    {
      if ((this->HzFactor != 0.0) && (this->HxFactor == 0.0))
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      vDestination[i] += (this->HamiltonianShift + this->HzFactor * this->Chain->TotalSz(i)) * vSource[i];
	    }
	}
      else
	{
	  if ((this->HzFactor == 0.0) && (this->HxFactor != 0.0))
	    {
	      double Coefficient;
	      int TmpPos;
	      int Dim = this->Chain->GetHilbertSpaceDimension();
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Complex Coefficient2 = this->HxFactor * vSource[i];
		  for (int j = 0; j < this->NbrSpin; ++j)
		    {
		      TmpPos = this->Chain->Smi(j, i, Coefficient);
		      if (TmpPos != Dim)
			{
			  vDestination[TmpPos] += Coefficient * Coefficient2;
			}
		      TmpPos = this->Chain->Spi(j, i, Coefficient);
		      if (TmpPos != Dim)
			{
			  vDestination[TmpPos] += Coefficient * Coefficient2;
			}
		    }
		  vDestination[i] += this->HamiltonianShift * vSource[i];
		}	  
	    }
	  else
	    {
	      double Coefficient;
	      int TmpPos;
	      int Dim = this->Chain->GetHilbertSpaceDimension();
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Complex Coefficient2 = this->HxFactor * vSource[i];
		  for (int j = 0; j < this->NbrSpin; ++j)
		    {
		      TmpPos = this->Chain->Smi(j, i, Coefficient);
		      if (TmpPos != Dim)
			{
			  vDestination[TmpPos] += Coefficient * Coefficient2;
			}
		      TmpPos = this->Chain->Spi(j, i, Coefficient);
		      if (TmpPos != Dim)
			{
			  vDestination[TmpPos] += Coefficient * Coefficient2;
			}
		    }
		  vDestination[i] += (this->HamiltonianShift + this->HzFactor * this->Chain->TotalSz(i)) * vSource[i];
		}
	    }
	}
    }
  else
    {
      if (this->HzFactor == 0.0)
	{
	  double Coefficient;
	  int TmpPos;
	  int Dim = this->Chain->GetHilbertSpaceDimension();
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Complex Coefficient2 = Complex(this->HxFactor, this->HyFactor) * vSource[i];
	      Complex Coefficient3 = Complex(this->HxFactor, -this->HyFactor) * vSource[i];
	      for (int j = 0; j < this->NbrSpin; ++j)
		{
		  TmpPos = this->Chain->Smi(j, i, Coefficient);
		  if (TmpPos != Dim)
		    {
		      vDestination[TmpPos] += Coefficient * Coefficient3;
		    }
		  TmpPos = this->Chain->Spi(j, i, Coefficient);
		  if (TmpPos != Dim)
		    {
		      vDestination[TmpPos] += Coefficient * Coefficient2;
		    }
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }	  
	}
      else
	{
	  double Coefficient;
	  int TmpPos;
	  int Dim = this->Chain->GetHilbertSpaceDimension();
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Complex Coefficient2 = Complex(this->HxFactor, this->HyFactor) * vSource[i];
	      Complex Coefficient3 = Complex(this->HxFactor, -this->HyFactor) * vSource[i];
	      for (int j = 0; j < this->NbrSpin; ++j)
		{
		  TmpPos = this->Chain->Smi(j, i, Coefficient);
		  if (TmpPos != Dim)
		    {
		      vDestination[TmpPos] += Coefficient * Coefficient3;
		    }
		  TmpPos = this->Chain->Spi(j, i, Coefficient);
		  if (TmpPos != Dim)
		    {
		      vDestination[TmpPos] += Coefficient * Coefficient2;
		    }
		}
	      vDestination[i] += (this->HamiltonianShift + this->HzFactor * this->Chain->TotalSz(i)) * vSource[i];
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

ComplexVector* SpinChainPureHFieldHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  return vDestinations;
}
