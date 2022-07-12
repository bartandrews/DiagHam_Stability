////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of periodic anisotropic spin chain hamiltonian           //
//                                                                            //
//                        last modification : 18/03/2002                      //
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


#include "Hamiltonian/PeriodicAnisotropicSpinChainHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::endl;
using std::ostream;


// constructor from default datas
//
// chain = pointer to Hilbert space of the associated system
// nbrSpin = number of spin
// jxx = array containing Jxx coupling constants between spins
// jyy = array containing Jyy coupling constants between spins
// jzz = array containing Jzz coupling constants between spins
// jxy = array containing Jxy coupling constants between spins
// jxz = array containing Jxz coupling constants between spins
// jyz = array containing Jyz coupling constants between spins
// sxFactor = array containing factors in front of Sx term (g mu B)
// syFactor = array containing factors in front of Sy term (g mu B)
// szFactor = array containing factors in front of Sz term (g mu B)

PeriodicAnisotropicSpinChainHamiltonian::PeriodicAnisotropicSpinChainHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* jxx,
										 double* jyy, double* jzz, double* jxy, double* jxz,
										 double* jyz, double* sxFactor, double* syFactor,
										 double* szFactor)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->Jzz = new double [this->NbrSpin];
  this->Jpm = new double [this->NbrSpin];
  this->JppRe = new double [this->NbrSpin];
  this->JppIm = new double [this->NbrSpin];
  this->JpzRe = new double [this->NbrSpin];
  this->JpzIm = new double [this->NbrSpin];
  this->SpFactorRe = new double [this->NbrSpin];
  this->SpFactorIm = new double [this->NbrSpin];
  this->SzFactor = new double [this->NbrSpin];
  for (int i = 0; i < (this->NbrSpin); i++)
    {
      this->Jzz[i] = jzz[i];
      this->Jpm[i] = 0.25 * (jxx[i] + jyy[i]);
      this->JppRe[i] = 0.25 * (jxx[i] - jyy[i]);
      this->JppIm[i] = - 0.5 * jxy[i];
      this->JpzRe[i] = 0.5 * jxz[i];
      this->JpzIm[i] = - 0.5 * jyz[i];
      this->SpFactorRe[i] = 0.5 * sxFactor[i];
      this->SpFactorIm[i] = - 0.5 * syFactor[i];
      this->SzFactor[i] = szFactor[i];
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

PeriodicAnisotropicSpinChainHamiltonian::~PeriodicAnisotropicSpinChainHamiltonian() 
{
  delete[] this->Jzz;
  delete[] this->Jpm;
  delete[] this->JppRe;
  delete[] this->JppIm;
  delete[] this->JpzRe;
  delete[] this->JpzIm;
  delete[] this->SpFactorRe;
  delete[] this->SpFactorIm;
  delete[] this->SzSzContributions;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PeriodicAnisotropicSpinChainHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* PeriodicAnisotropicSpinChainHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

PeriodicAnisotropicSpinChainHamiltonian& PeriodicAnisotropicSpinChainHamiltonian::SetChain(AbstractSpinChain* chain)
{  
  delete[] this->SzSzContributions;
  this->Chain = chain;                                                                                                                        this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int PeriodicAnisotropicSpinChainHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void PeriodicAnisotropicSpinChainHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicAnisotropicSpinChainHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  return Complex();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicAnisotropicSpinChainHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  double x = 0.0;
  double y = 0.0;
  double coef;
  int MaxPos = this->NbrSpin - 1;  
  int pos;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoef1;
  double TmpCoef2;
  double Tmpx;
  double Tmpy;
  int jinc;
  int j = 0;
  for (; j < MaxPos; ++j)
    {
      if (this->Jpm[j] != 0)
	{
	  TmpCoef1 = this->Jpm[j];
	  jinc = j + 1;
	  for (int i = 0; i < dim; ++i)
	    {	 	    
	      pos = this->Chain->SmiSpj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  coef *= TmpCoef1;
		  x += coef * (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  y += coef * (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		}
	      pos = this->Chain->SmiSpj(jinc, j, i, coef);
	      if (pos != dim)
		{
		  coef *= TmpCoef1;
		  x += coef * (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  y += coef * (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		}
	    }
	}
      if ((this->JppRe[j] != 0) || (this->JppIm[j] != 0))
	{
	  TmpCoef1 = this->JppRe[j];
	  TmpCoef2 = this->JppIm[j];
	  jinc = j + 1;
	  for (int i = 0; i < dim; ++i)
	    {	 	    
	      pos = this->Chain->SpiSpj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 - Tmpy * TmpCoef2);
		  y += coef * (Tmpx * TmpCoef2 + Tmpy * TmpCoef1);
		}
	      pos = this->Chain->SmiSmj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 + Tmpy * TmpCoef2);
		  y += coef * (Tmpy * TmpCoef1 - Tmpx * TmpCoef2);
		}
	    }
	}
      if ((this->JpzRe[j] != 0) || (this->JpzIm[j] != 0))
	{
	  TmpCoef1 = this->JpzRe[j];
	  TmpCoef2 = this->JpzIm[j];
	  jinc = j + 1;
	  for (int i = 0; i < dim; ++i)
	    {	 	    
	      pos = this->Chain->SpiSzj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 - Tmpy * TmpCoef2);
		  y += coef * (Tmpx * TmpCoef2 + Tmpy * TmpCoef1);
		}
	      pos = this->Chain->SpiSzj(jinc, j, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 - Tmpy * TmpCoef2);
		  y += coef * (Tmpx * TmpCoef2 + Tmpy * TmpCoef1);
		}
	      pos = this->Chain->SmiSzj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 + Tmpy * TmpCoef2);
		  y += coef * (Tmpy * TmpCoef1 - Tmpx * TmpCoef2);
		}
	      pos = this->Chain->SmiSzj(jinc, j, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 + Tmpy * TmpCoef2);
		  y += coef * (Tmpy * TmpCoef1 - Tmpx * TmpCoef2);
		}
	    }
	}
      if ((this->SpFactorRe[j] != 0) || (this->SpFactorIm[j] != 0))
	{
	  TmpCoef1 = this->SpFactorRe[j];
	  TmpCoef2 = this->SpFactorIm[j];
	  for (int i = 0; i < dim; ++i)
	    {	 	    
	      pos = this->Chain->Spi(j, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 - Tmpy * TmpCoef2);
		  y += coef * (Tmpx * TmpCoef2 + Tmpy * TmpCoef1);
		}
	      pos = this->Chain->Smi(j, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 + Tmpy * TmpCoef2);
		  y += coef * (Tmpy * TmpCoef1 - Tmpx * TmpCoef2);
		}
	    }
	}
    }
  if (this->Jpm[j] != 0)
    {
      TmpCoef1 = this->Jpm[j];
      for (int i = 0; i < dim; ++i)
	{	 	    
	  pos = this->Chain->SmiSpj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      coef *= TmpCoef1;
	      x += coef * (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      y += coef * (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	    }
	  pos = this->Chain->SmiSpj(0, j, i, coef);
	  if (pos != dim)
	    {
	      coef *= TmpCoef1;
	      x += coef * (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      y += coef * (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	    }
	}
    }
  if ((this->JppRe[j] != 0) || (this->JppIm[j] != 0))
    {
      TmpCoef1 = this->JppRe[j];
      TmpCoef2 = this->JppIm[j];
      for (int i = 0; i < dim; ++i)
	{	 	    
	  pos = this->Chain->SpiSpj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	      x += coef * (Tmpx * TmpCoef1 - Tmpy * TmpCoef2);
	      y += coef * (Tmpx * TmpCoef2 + Tmpy * TmpCoef1);
	    }
	  pos = this->Chain->SmiSmj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	      x += coef * (Tmpx * TmpCoef1 + Tmpy * TmpCoef2);
	      y += coef * (Tmpy * TmpCoef1 - Tmpx * TmpCoef2);
	    }
	}
    }
  if ((this->JpzRe[j] != 0) || (this->JpzIm[j] != 0))
    {
      TmpCoef1 = this->JpzRe[j];
      TmpCoef2 = this->JpzIm[j];
      for (int i = 0; i < dim; ++i)
	{	 	    
	  pos = this->Chain->SpiSzj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	      x += coef * (Tmpx * TmpCoef1 - Tmpy * TmpCoef2);
	      y += coef * (Tmpx * TmpCoef2 + Tmpy * TmpCoef1);
	    }
	  pos = this->Chain->SpiSzj(0, j, i, coef);
	  if (pos != dim)
	    {
	      Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	      x += coef * (Tmpx * TmpCoef1 - Tmpy * TmpCoef2);
	      y += coef * (Tmpx * TmpCoef2 + Tmpy * TmpCoef1);
	    }
	  pos = this->Chain->SmiSzj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	      x += coef * (Tmpx * TmpCoef1 + Tmpy * TmpCoef2);
	      y += coef * (Tmpy * TmpCoef1 - Tmpx * TmpCoef2);
	    }
	  pos = this->Chain->SmiSzj(0, j, i, coef);
	  if (pos != dim)
	    {
	      Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	      x += coef * (Tmpx * TmpCoef1 + Tmpy * TmpCoef2);
	      y += coef * (Tmpy * TmpCoef1 - Tmpx * TmpCoef2);
	    }
	}
    }
  if ((this->SpFactorRe[j] != 0) || (this->SpFactorIm[j] != 0))
    {
      TmpCoef1 = this->SpFactorRe[j];
      TmpCoef2 = this->SpFactorIm[j];
      for (int i = 0; i < dim; ++i)
	{	 	    
	  pos = this->Chain->Spi(j, i, coef);
	  if (pos != dim)
	    {
	      Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	      x += coef * (Tmpx * TmpCoef1 - Tmpy * TmpCoef2);
	      y += coef * (Tmpx * TmpCoef2 + Tmpy * TmpCoef1);
	    }
	  pos = this->Chain->Smi(j, i, coef);
	  if (pos != dim)
	    {
	      Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	      Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
	      x += coef * (Tmpx * TmpCoef1 + Tmpy * TmpCoef2);
	      y += coef * (Tmpy * TmpCoef1 - Tmpx * TmpCoef2);
	    }
	}
    }
  // SzSz contributions
  for (int i = 0; i < dim; i++)
    {	 	    
      x += this->SzSzContributions[i] * (V1.Re(i) * V2.Re(i) + V1.Im(i) * V2.Im(i));
      y += this->SzSzContributions[i] * (V1.Re(i) * V2.Im(i) - V1.Im(i) * V2.Re(i));
    }
  return Complex(x, y);
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& PeriodicAnisotropicSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& PeriodicAnisotropicSpinChainHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							      int firstComponent, int nbrComponent) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& PeriodicAnisotropicSpinChainHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& PeriodicAnisotropicSpinChainHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									 int firstComponent, int nbrComponent) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicAnisotropicSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Chain->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PeriodicAnisotropicSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									 int firstComponent, int nbrComponent) 
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  double TmpCoef1;
  double TmpCoef2;
  int jinc;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination.Re(i) = this->SzSzContributions[i] * vSource.Re(i);
      vDestination.Im(i) = this->SzSzContributions[i] * vSource.Im(i);
    }
  int j = 0;
  for (; j < MaxPos; ++j)
    {
      if (this->Jpm[j] != 0)
	{
	  TmpCoef1 = this->Jpm[j];
	  jinc = j + 1;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {	 	    
	      pos = this->Chain->SmiSpj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  coef *= TmpCoef1;
		  vDestination.Re(pos) += coef * vSource.Re(i);
		  vDestination.Im(pos) += coef * vSource.Im(i);
		}
	      pos = this->Chain->SmiSpj(jinc, j, i, coef);
	      if (pos != dim)
		{
		  coef *= TmpCoef1;
		  vDestination.Re(pos) += coef * vSource.Re(i);
		  vDestination.Im(pos) += coef * vSource.Im(i);
		}
	    }
	}
      if ((this->JppRe[j] != 0) || (this->JppIm[j] != 0))
	{
	  TmpCoef1 = this->JppRe[j];
	  TmpCoef2 = this->JppIm[j];
	  jinc = j + 1;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {	 	    
	      pos = this->Chain->SpiSpj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->SmiSmj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
		}
	    }
	}
      if ((this->JpzRe[j] != 0) || (this->JpzIm[j] != 0))
	{
	  TmpCoef1 = this->JpzRe[j];
	  TmpCoef2 = this->JpzIm[j];
	  jinc = j + 1;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {	 	    
	      pos = this->Chain->SpiSzj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->SpiSzj(jinc, j, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->SmiSzj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->SmiSzj(jinc, j, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
		}
	    }
	}
      if ((this->SpFactorRe[j] != 0.0) || (this->SpFactorIm[j] != 0.0))
	{
	  TmpCoef1 = this->SpFactorRe[j];
	  TmpCoef2 = this->SpFactorIm[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {	 	    
	      pos = this->Chain->Spi(j, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->Smi(j, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
		}
	    }
	}
    }
  if (this->Jpm[j] != 0)
    {
      TmpCoef1 = this->Jpm[j];
      for (int i = firstComponent; i < LastComponent; ++i)
	{	 	    
	  pos = this->Chain->SmiSpj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      coef *= TmpCoef1;
	      vDestination.Re(pos) += coef * vSource.Re(i);
	      vDestination.Im(pos) += coef * vSource.Im(i);
	    }
	  pos = this->Chain->SmiSpj(0, j, i, coef);
	  if (pos != dim)
	    {
	      coef *= TmpCoef1;
	      vDestination.Re(pos) += coef * vSource.Re(i);
	      vDestination.Im(pos) += coef * vSource.Im(i);
	    }
	}
    }
  if ((this->JppRe[j] != 0) || (this->JppIm[j] != 0))
    {
      TmpCoef1 = this->JppRe[j];
      TmpCoef2 = this->JppIm[j];
      jinc = 0;
      for (int i = firstComponent; i < LastComponent; ++i)
	{	 	    
	  pos = this->Chain->SpiSpj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->SmiSmj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
	    }
	}
    }
  if ((this->JpzRe[j] != 0) || (this->JpzIm[j] != 0))
    {
      TmpCoef1 = this->JpzRe[j];
      TmpCoef2 = this->JpzIm[j];
      jinc = 0;
      for (int i = firstComponent; i < LastComponent; ++i)
	{	 	    
	  pos = this->Chain->SpiSzj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->SpiSzj(0, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->SmiSzj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->SmiSzj(0, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
	    }
	}
    }
  if ((this->SpFactorRe[j] != 0.0) || (this->SpFactorIm[j] != 0.0))
    {
      TmpCoef1 = this->SpFactorRe[j];
      TmpCoef2 = this->SpFactorIm[j];
      for (int i = firstComponent; i < LastComponent; ++i)
	{	 	    
	  pos = this->Chain->Spi(j, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->Smi(j, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
	    }
	}
    }
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> PeriodicAnisotropicSpinChainHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  int Dim = this->Chain->GetHilbertSpaceDimension();
  RealSymmetricMatrix* Sx = new RealSymmetricMatrix (Dim, true);
  RealAntisymmetricMatrix* Sy = new RealAntisymmetricMatrix (Dim, true);
  RealSymmetricMatrix* Sz = new RealSymmetricMatrix (Dim, true);
  this->Chain->Sxi(0, *Sx);
  this->Chain->Syi(0, *Sy);
  this->Chain->Szi(0, *Sz);
  TmpList += Sx;
  TmpList += Sy;
  TmpList += Sz;
  return TmpList;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicAnisotropicSpinChainHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Chain->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PeriodicAnisotropicSpinChainHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									    int firstComponent, int nbrComponent) 
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  double TmpCoef1;
  double TmpCoef2;
  int jinc;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination.Re(i) += this->SzSzContributions[i] * vSource.Re(i);
      vDestination.Im(i) += this->SzSzContributions[i] * vSource.Im(i);
    }
  int j = 0;
  for (; j < MaxPos; ++j)
    {
      if (this->Jpm[j] != 0)
	{
	  TmpCoef1 = this->Jpm[j];
	  jinc = j + 1;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {	 	    
	      pos = this->Chain->SmiSpj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  coef *= TmpCoef1;
		  vDestination.Re(pos) += coef * vSource.Re(i);
		  vDestination.Im(pos) += coef * vSource.Im(i);
		}
	      pos = this->Chain->SmiSpj(jinc, j, i, coef);
	      if (pos != dim)
		{
		  coef *= TmpCoef1;
		  vDestination.Re(pos) += coef * vSource.Re(i);
		  vDestination.Im(pos) += coef * vSource.Im(i);
		}
	    }
	}
      if ((this->JppRe[j] != 0) || (this->JppIm[j] != 0))
	{
	  TmpCoef1 = this->JppRe[j];
	  TmpCoef2 = this->JppIm[j];
	  jinc = j + 1;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {	 	    
	      pos = this->Chain->SpiSpj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->SmiSmj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
		}
	    }
	}
      if ((this->JpzRe[j] != 0) || (this->JpzIm[j] != 0))
	{
	  TmpCoef1 = this->JpzRe[j];
	  TmpCoef2 = this->JpzIm[j];
	  jinc = j + 1;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {	 	    
	      pos = this->Chain->SpiSzj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->SpiSzj(jinc, j, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->SmiSzj(j, jinc, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->SmiSzj(jinc, j, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
		}
	    }
	}
      if ((this->SpFactorRe[j] != 0.0) || (this->SpFactorIm[j] != 0.0))
	{
	  TmpCoef1 = this->SpFactorRe[j];
	  TmpCoef2 = this->SpFactorIm[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {	 	    
	      pos = this->Chain->Spi(j, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
		}
	      pos = this->Chain->Smi(j, i, coef);
	      if (pos != dim)
		{
		  vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
		  vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
		}
	    }
	}
    }
  if (this->Jpm[j] != 0)
    {
      TmpCoef1 = this->Jpm[j];
      for (int i = firstComponent; i < LastComponent; ++i)
	{	 	    
	  pos = this->Chain->SmiSpj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      coef *= TmpCoef1;
	      vDestination.Re(pos) += coef * vSource.Re(i);
	      vDestination.Im(pos) += coef * vSource.Im(i);
	    }
	  pos = this->Chain->SmiSpj(0, j, i, coef);
	  if (pos != dim)
	    {
	      coef *= TmpCoef1;
	      vDestination.Re(pos) += coef * vSource.Re(i);
	      vDestination.Im(pos) += coef * vSource.Im(i);
	    }
	}
    }
  if ((this->JppRe[j] != 0) || (this->JppIm[j] != 0))
    {
      TmpCoef1 = this->JppRe[j];
      TmpCoef2 = this->JppIm[j];
      jinc = 0;
      for (int i = firstComponent; i < LastComponent; ++i)
	{	 	    
	  pos = this->Chain->SpiSpj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->SmiSmj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
	    }
	}
    }
  if ((this->JpzRe[j] != 0) || (this->JpzIm[j] != 0))
    {
      TmpCoef1 = this->JpzRe[j];
      TmpCoef2 = this->JpzIm[j];
      jinc = 0;
      for (int i = firstComponent; i < LastComponent; ++i)
	{	 	    
	  pos = this->Chain->SpiSzj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->SpiSzj(0, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->SmiSzj(j, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->SmiSzj(0, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
	    }
	}
    }
  if ((this->SpFactorRe[j] != 0.0) || (this->SpFactorIm[j] != 0.0))
    {
      TmpCoef1 = this->SpFactorRe[j];
      TmpCoef2 = this->SpFactorIm[j];
      for (int i = firstComponent; i < LastComponent; ++i)
	{	 	    
	  pos = this->Chain->Spi(j, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 + vSource.Re(i) * TmpCoef2);
	    }
	  pos = this->Chain->Smi(j, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
	    }
	}
    }
  return vDestination;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> PeriodicAnisotropicSpinChainHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  int Dim = this->Chain->GetHilbertSpaceDimension();
  RealSymmetricMatrix* Sx = new RealSymmetricMatrix (Dim, true);
  RealAntisymmetricMatrix* Sy = new RealAntisymmetricMatrix (Dim, true);
  RealSymmetricMatrix* Sz = new RealSymmetricMatrix (Dim, true);
  this->Chain->Sxi(this->NbrSpin - 1, *Sx);
  this->Chain->Syi(this->NbrSpin - 1, *Sy);
  this->Chain->Szi(this->NbrSpin - 1, *Sz);
  TmpList += Sx;
  TmpList += Sy;
  TmpList += Sz;
  return TmpList;
}

// evaluate diagonal matrix elements
// 

void PeriodicAnisotropicSpinChainHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coefficient;

  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      int j = 0;
      for (; j < (this->NbrSpin - 1); j++)
	{
	  this->SzSzContributions[i] += this->Jzz[j] * this->Chain->SziSzj(j, j + 1, i);
	  this->Chain->Szi(j, i, coefficient);
	  this->SzSzContributions[i] += this->SzFactor[j] * coefficient;
	}
      this->SzSzContributions[i] += this->Jzz[j] * this->Chain->SziSzj(j, 0, i);
      this->Chain->Szi(j, i, coefficient);
      this->SzSzContributions[i] += this->SzFactor[j] * coefficient;
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, PeriodicAnisotropicSpinChainHamiltonian& H) 
{
  ComplexVector TmpV2 (H.Chain->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [H.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(H.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2.Re(i - 1) = 0.0;
	}
      TmpV2.Re(i) = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Chain->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j].Re(i);
	  if (TmpV[j].Im(i) < 0.0)
	    {
	      Str << ((TmpV[j].Im(i))) << "i";
	    }
	  else
	    if (TmpV[j].Im(i) > 0.0)
	      {
		Str << "+" << TmpV[j].Im(i) << "i";
	      }
	  Str << "   ";
	}
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, PeriodicAnisotropicSpinChainHamiltonian& H) 
{
  ComplexVector TmpV2 (H.Chain->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [H.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(H.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2.Re(i - 1) = 0.0;
	}
      TmpV2.Re(i) = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Chain->GetHilbertSpaceDimension() - 1); ++i)
    {
      Str << "{";
      for (int j = 0; j < (H.Chain->GetHilbertSpaceDimension() - 1); ++j)
	{
	  Str << TmpV[j].Re(i);
	  if (TmpV[j].Im(i) < 0)
	    {
	      Str << ((TmpV[j].Im(i))) << "I";
	    }
	  else
	    if (TmpV[j].Im(i) > 0)
	      {
		Str << "+" << ((TmpV[j].Im(i))) << "I";
	      }	  
	  Str << ",";
	}
      Str << TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Re(i);
      if (TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Im(i) < 0)
	{
	  Str << ((TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Im(i))) << "I";
	}
      else
	if (TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Im(i) > 0)
	  {
	    Str << "+" << ((TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Im(i))) << "I";
	  }	  
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Chain->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j].Re(H.Chain->GetHilbertSpaceDimension() - 1);
      if (TmpV[j].Im(H.Chain->GetHilbertSpaceDimension() - 1) < 0)
	{
	  Str << ((TmpV[j].Im(H.Chain->GetHilbertSpaceDimension() - 1))) << "I";
	}
      else
	if (TmpV[j].Im(H.Chain->GetHilbertSpaceDimension() - 1) > 0)
	  {
	    Str << "+" << ((TmpV[j].Im(H.Chain->GetHilbertSpaceDimension() - 1))) << "I";
	  }	  
      Str << ",";
    }
  Str << TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Re(H.Chain->GetHilbertSpaceDimension() - 1);
  if (TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Im(H.Chain->GetHilbertSpaceDimension() - 1) < 0)
    {
      Str << ((TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Im(H.Chain->GetHilbertSpaceDimension() - 1))) << "I";
    }
  else
    if (TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Im(H.Chain->GetHilbertSpaceDimension() - 1) > 0)
      {
	Str << "+" << ((TmpV[H.Chain->GetHilbertSpaceDimension() - 1].Im(H.Chain->GetHilbertSpaceDimension() - 1))) << "I";
      }	  
  Str << "}}";
  return Str;
}

