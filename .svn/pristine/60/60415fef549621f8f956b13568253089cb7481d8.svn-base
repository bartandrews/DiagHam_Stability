////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of unconventional 3D toric code                  //
//                                                                            //
//                        last modification : 23/05/2017                      //
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


#include "Hamiltonian/Unconventional3DToricCodeHamiltonian.h"
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

Unconventional3DToricCodeHamiltonian::Unconventional3DToricCodeHamiltonian()
{
}


// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// nbrSpinX = number of spin along the x direction
// nbrSpinY = number of spin along the y direction
// nbrSpinZ = number of spin along the z direction
// periodicBoundaryConditionsX = true if periodic boundary conditions have to be usedalong the x direction
// periodicBoundaryConditionsY = true if periodic boundary conditions have to be usedalong the y direction
// periodicBoundaryConditionsZ = true if periodic boundary conditions have to be usedalong the z direction

Unconventional3DToricCodeHamiltonian::Unconventional3DToricCodeHamiltonian(AbstractSpinChain* chain, int nbrSpinX, int nbrSpinY, int nbrSpinZ,
									   bool periodicBoundaryConditionsX, bool periodicBoundaryConditionsY, bool periodicBoundaryConditionsZ)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpinZ = nbrSpinZ;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY * this->NbrSpinZ;
  this->PeriodicBoundaryConditionsX = periodicBoundaryConditionsX;
  this->PeriodicBoundaryConditionsY = periodicBoundaryConditionsY;
  this->PeriodicBoundaryConditionsZ = periodicBoundaryConditionsZ;
  this->HamiltonianShift = 0.0;
}

// destructor
//

Unconventional3DToricCodeHamiltonian::~Unconventional3DToricCodeHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void Unconventional3DToricCodeHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (AbstractSpinChain*) hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* Unconventional3DToricCodeHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int Unconventional3DToricCodeHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void Unconventional3DToricCodeHamiltonian::ShiftHamiltonian (double shift)
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

RealVector& Unconventional3DToricCodeHamiltonian::Unconventional3DToricCodeHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
													    int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int MaxXIndex = this->NbrSpinX;
  if (this->PeriodicBoundaryConditionsX == false)
    {
      --MaxXIndex;
    }
  int MaxYIndex = this->NbrSpinY;
  if (this->PeriodicBoundaryConditionsY == false)
    {
      --MaxYIndex;
    }
  int MaxZIndex = this->NbrSpinZ;
  if (this->PeriodicBoundaryConditionsZ == false)
    {
      --MaxZIndex;
    }
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->HamiltonianShift * vSource[i];
      for (int j = 0; j < MaxXIndex; ++j)
	{
	  for (int k = 0; k < MaxYIndex; ++k)
	    {
	      for (int l = 0; l < MaxZIndex; ++l)
		{
		  this->ApplyCubeOperator(j, k, l, vDestination, i, vSource[i]);
//		  this->ApplyPlaquetteOperator(j, k, l, vDestination, i, vSource[i]);
		}
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

RealVector* Unconventional3DToricCodeHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
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
  int MaxXIndex = this->NbrSpinX;
  if (this->PeriodicBoundaryConditionsX == false)
    {
      --MaxXIndex;
    }
  int MaxYIndex = this->NbrSpinY;
  if (this->PeriodicBoundaryConditionsY == false)
    {
      --MaxYIndex;
    }
  int MaxZIndex = this->NbrSpinZ;
  if (this->PeriodicBoundaryConditionsZ == false)
    {
      --MaxZIndex;
    }
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int k = 0; k < nbrVectors; ++k)
	{
	  vDestinations[k][i] += this->HamiltonianShift * vSources[k][i];
	  TmpValues[k] = vSources[k][i];
	}
      for (int j = 0; j < MaxXIndex; ++j)
	{
	  for (int k = 0; k < MaxYIndex; ++k)
	    {
	      for (int l = 0; l < MaxZIndex; ++l)
		{		  
		  for (int m = 0; m < nbrVectors; ++m)	      
		    {
		      this->ApplyCubeOperator(j, k, l, vDestinations[m], i, vSources[m][i]);
		    }
		}
	    }
	}
    }
  delete[] TmpValues;
  return vDestinations;
}


// apply the cube operator to one basis state
//
// cornerX = lower leftmost frontmost cube corner x coordinate
// cornerY = lower leftmost frontmost cube corner y coordinate
// cornerZ = lower leftmost frontmost cube corner z coordinate
// vDestination = reference ont the vector to which result has to be added
// index = basis state index
// coefficient = basis state coefficient

void Unconventional3DToricCodeHamiltonian::ApplyCubeOperator(int cornerX, int cornerY, int cornerZ, RealVector& vDestination, int index, double coefficient)
{
  int IncCornerX = (cornerX + 1) % this->NbrSpinX;
  int IncCornerY = (cornerY + 1) % this->NbrSpinY;
  int IncCornerZ = (cornerZ + 1) % this->NbrSpinZ;
  int TmpDim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoefficient1 = (this->Chain->SziSzj(this->GetLinearizedIndex(IncCornerX, cornerY, cornerZ), this->GetLinearizedIndex(cornerX, IncCornerY, cornerZ), index) 
			   * this->Chain->SziSzj(this->GetLinearizedIndex(IncCornerX, cornerY, IncCornerZ), this->GetLinearizedIndex(cornerX, IncCornerY, IncCornerZ), index)) * 16.0 * coefficient;
  double TmpCoefficient2 = 0.0;
  double TmpCoefficient3 = 0.0;  
  int TmpPos = this->Chain->SpiSpj(this->GetLinearizedIndex(cornerX, cornerY, cornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, cornerZ)), index, TmpCoefficient2);
  if (TmpPos < TmpDim)
    {
      TmpCoefficient2 *= TmpCoefficient1;
      int TmpPos2 = this->Chain->SpiSpj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SmiSmj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SpiSmj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SmiSpj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
    }
  TmpPos = this->Chain->SmiSmj(this->GetLinearizedIndex(cornerX, cornerY, cornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, cornerZ)), index, TmpCoefficient2);
  if (TmpPos < TmpDim)
    {
      TmpCoefficient2 *= TmpCoefficient1;
      int TmpPos2 = this->Chain->SpiSpj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SmiSmj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SpiSmj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SmiSpj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
    }
  TmpPos = this->Chain->SpiSmj(this->GetLinearizedIndex(cornerX, cornerY, cornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, cornerZ)), index, TmpCoefficient2);
  if (TmpPos < TmpDim)
    {
      TmpCoefficient2 *= TmpCoefficient1;
      int TmpPos2 = this->Chain->SpiSpj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SmiSmj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SpiSmj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SmiSpj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
    }
  TmpPos = this->Chain->SmiSpj(this->GetLinearizedIndex(cornerX, cornerY, cornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, cornerZ)), index, TmpCoefficient2);
  if (TmpPos < TmpDim)
    {
      TmpCoefficient2 *= TmpCoefficient1;
      int TmpPos2 = this->Chain->SpiSpj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SmiSmj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SpiSmj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
      TmpPos2 = this->Chain->SmiSpj(this->GetLinearizedIndex(cornerX, cornerY, IncCornerZ), (this->GetLinearizedIndex(IncCornerX, IncCornerY, IncCornerZ)), TmpPos, TmpCoefficient3);
      if (TmpPos2 < TmpDim)
	{
	  vDestination[TmpPos2] -= TmpCoefficient2 * TmpCoefficient3;
	}
    }
}

// apply the plaquette operator in the xy plane to one basis state
//
// cornerX = leftmost frontmost plaquette corner x coordinate
// cornerY = leftmost frontmost plaquette corner y coordinate
// cornerZ = plaquette z coordinate
// vDestination = reference on the vector to which result has to be added
// index = basis state index
// coefficient = basis state coefficient

void Unconventional3DToricCodeHamiltonian::ApplyPlaquetteOperator(int cornerX, int cornerY, int cornerZ, RealVector& vDestination, int index, double coefficient)
{
  int IncCornerX = (cornerX + 1) % this->NbrSpinX;
  int IncCornerY = (cornerY + 1) % this->NbrSpinY;
  int TmpDim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoefficient1 = this->Chain->SziSzj(this->GetLinearizedIndex(IncCornerX, cornerY, cornerZ), this->GetLinearizedIndex(cornerX, IncCornerY, cornerZ), index) * 4.0 * coefficient;			   
  double TmpCoefficient2 = 0.0;
  int TmpPos = this->Chain->SpiSpj(this->GetLinearizedIndex(cornerX, cornerY, cornerZ), this->GetLinearizedIndex(IncCornerX, IncCornerY, cornerZ), index, TmpCoefficient2);
  if (TmpPos < TmpDim)
    {
      vDestination[TmpPos] -= TmpCoefficient1 * TmpCoefficient2;
    }
  TmpPos = this->Chain->SmiSmj(this->GetLinearizedIndex(cornerX, cornerY, cornerZ), this->GetLinearizedIndex(IncCornerX, IncCornerY, cornerZ), index, TmpCoefficient2);
  if (TmpPos < TmpDim)
    {
      vDestination[TmpPos] -= TmpCoefficient1 * TmpCoefficient2;
    }
  TmpPos = this->Chain->SpiSmj(this->GetLinearizedIndex(cornerX, cornerY, cornerZ), this->GetLinearizedIndex(IncCornerX, IncCornerY, cornerZ), index, TmpCoefficient2);
  if (TmpPos < TmpDim)
    {
      vDestination[TmpPos] -= TmpCoefficient1 * TmpCoefficient2;
    }
  TmpPos = this->Chain->SmiSpj(this->GetLinearizedIndex(cornerX, cornerY, cornerZ), this->GetLinearizedIndex(IncCornerX, IncCornerY, cornerZ), index, TmpCoefficient2);
  if (TmpPos < TmpDim)
    {
      vDestination[TmpPos] -= TmpCoefficient1 * TmpCoefficient2;
    }
}

