////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                         Class author: Cecile Repellin                      //
//                                                                            //
//                                                                            //
//         class of two dimensional spin model with pseudospin                //
//                                                                            //
//                        last modification : 12/06/2016                      //
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


#include "Hamiltonian/TwoDimensionalTriangularLatticeWithPseudospinHamiltonian.h"
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

TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::TwoDimensionalTriangularLatticeWithPseudospinHamiltonian()
{
}


// constructor
//
// chain = pointer to Hilbert space of the associated system
// nbrSpinX = number of spin along the x direction
// nbrSpinY = number of spin along the y direction
// jFactor = amplitude of the Ising term
// hxFactor = amplitudes of the Zeeman term along x
// hzFactor = amplitudes of the Zeeman term along z
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::TwoDimensionalTriangularLatticeWithPseudospinHamiltonian(Spin1_2ChainWithPseudospin* chain, int nbrSpinX, int nbrSpinY, double jFactor, bool periodicBoundaryConditions, int offset)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->JFactor = jFactor;
  this->Offset = offset;
  
  // projection of kagome onto the s = 1/2 states on each triangle
  this->PseudospinCouplingElements = new double[3];    
  this->PseudospinCouplingElements[0] = 0.0;
  this->PseudospinCouplingElements[1] = 1.0/sqrt(3.0);
  this->PseudospinCouplingElements[2] = -1.0/sqrt(3.0);
  
  this-> PseudospinDiagCouplingElements = new double*[3];
  for (int i = 0; i < 3; ++i)
    this->PseudospinDiagCouplingElements[i] = new double[2];
  this->PseudospinDiagCouplingElements[0][0] = 1.0;
  this->PseudospinDiagCouplingElements[0][1] = -1.0/3.0;
  
  this->PseudospinDiagCouplingElements[1][0] = 0.0;
  this->PseudospinDiagCouplingElements[1][1] = 2.0/3.0;
  
  this->PseudospinDiagCouplingElements[2][0] = 0.0;
  this->PseudospinDiagCouplingElements[2][1] = 2.0/3.0;
  
//   //test: trivial coupling, Heisenberg model
//   this->PseudospinCouplingElements = new double[3];    
//   this->PseudospinCouplingElements[0] = 0.0;
//   this->PseudospinCouplingElements[1] = 0.0;
//   this->PseudospinCouplingElements[2] = 0.0;
//   
//   this-> PseudospinDiagCouplingElements = new double*[3];
//   for (int i = 0; i < 3; ++i)
//     this->PseudospinDiagCouplingElements[i] = new double[2];
//   this->PseudospinDiagCouplingElements[0][0] = 1.0;
//   this->PseudospinDiagCouplingElements[0][1] = 1.0;
//   
//   this->PseudospinDiagCouplingElements[1][0] = 1.0;
//   this->PseudospinDiagCouplingElements[1][1] = 1.0;
//   
//   this->PseudospinDiagCouplingElements[2][0] = 1.0;
//   this->PseudospinDiagCouplingElements[2][1] = 1.0;
  
  
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::~TwoDimensionalTriangularLatticeWithPseudospinHamiltonian() 
{
  delete[] this->SzSzContributions;
  delete[] this->PseudospinCouplingElements;
  for (int i = 0; i < 3; ++i)
    delete[] this->PseudospinDiagCouplingElements[i];
  delete[] this->PseudospinDiagCouplingElements;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (Spin1_2ChainWithPseudospin*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::ShiftHamiltonian (double shift)
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

RealVector& TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int pos3;
  int MaxPos = this->NbrSpin - 1;
  int TmpIndex1;
  int TmpIndex2;
  double TmpCoefficient;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
//       cout << (this->SzSzContributions[i]) << " " << endl;
      double TmpValue = vSource[i] * this->JFactor * 0.5;
//       TmpValue = 0.0;
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	     //AC part
	     if ((this->NbrSpinY > 1) || (this->Offset != 0))
	     {
		TmpIndex1 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
		pos = this->Chain->SpiSmj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[2]);		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		    vDestination[pos] += TmpValue * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[2];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[0]);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    vDestination[pos3] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		}
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[2]);		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		    vDestination[pos] += TmpValue * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[2];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[0]);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    vDestination[pos3] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		}	      
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2);
		TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[2];
		if (TmpCoefficient != 0.0)
		{
		  coef2 *= this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		  vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2;
		  
		  pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		  vDestination[pos2] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];		  
		}
		
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[2]);
		pos2 = this->Chain->JOffDiagonali(TmpIndex2, i, coef2);
		vDestination[pos2] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
	     }
	      

	      //AB part
	      if ((this->NbrSpinX > 1) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
		pos = this->Chain->SpiSmj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[1]);
		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		    vDestination[pos] += TmpValue * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[0]);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    vDestination[pos3] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		}
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[1]);
		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		    vDestination[pos] += TmpValue * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[0]);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    vDestination[pos3] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		}
	      
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2);
		TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		if (TmpCoefficient != 0.0)
		{
		  coef2 *= this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		  vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2;
		  
		  pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		  vDestination[pos2] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];		  
		}
		
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]);
		pos2 = this->Chain->JOffDiagonali(TmpIndex2, i, coef2);
		vDestination[pos2] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
	      }
// 	      
	      //CB part
	      if ((this->NbrSpinX > 1) && ((this->NbrSpinY > 1) || (this->Offset != 0)) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		pos = this->Chain->SpiSmj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[1]);
		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[2]);
		    vDestination[pos] += TmpValue * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[2]);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    vDestination[pos3] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
		  }
		}
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[1]);
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[2]);
		    vDestination[pos] += TmpValue * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[2]);
		    vDestination[pos2] += TmpValue * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    vDestination[pos3] += TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
		  }
		}
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2);
		TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		if (TmpCoefficient != 0.0)
		{
		  coef2 *= this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[2]);
		  vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2;
		  
		  pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		  vDestination[pos2] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];		  
		}
		
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]);
		pos2 = this->Chain->JOffDiagonali(TmpIndex2, i, coef2);
		vDestination[pos2] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
	      }
	    }
	}
    }
//   for (int i = firstComponent; i <LastComponent; ++i)
//     cout << vDestination[i] << " ";
//   cout << endl;
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

RealVector* TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
//   cout << "Error: TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::LowLevelMultipleAddMultiply is not implemented" << endl;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int LastComponent = firstComponent + nbrComponent;
  double coef;
  double coef2;
  int pos;
  int pos2;
  int pos3;
  double TmpCoefficient;
  int TmpIndex1;
  int TmpIndex2;
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
	{
	  TmpValues[k] = vSources[k][i] * this->JFactor * 0.5;
	}
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      //AC part
	     if ((this->NbrSpinY > 1) || (this->Offset != 0))
	     {
		TmpIndex1 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
		pos = this->Chain->SpiSmj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[2]);		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[2];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[0]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos3] += TmpValues[l] * coef2 * this->PseudospinCouplingElements[0];
		  }
		}
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[2]);		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[2];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[0]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos3] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		}	      
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2);
		TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[2];
		if (TmpCoefficient != 0.0)
		{
		  coef2 *= this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		  for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 2.0 * TmpValues[l] * TmpCoefficient * coef2;
		  
		  pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		  for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];		  
		}
		
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[2]);
		pos2 = this->Chain->JOffDiagonali(TmpIndex2, i, coef2);
		for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
	     }
	      

	      //AB part
	      if ((this->NbrSpinX > 1) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
		pos = this->Chain->SpiSmj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[1]);
		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[0]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos3] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		}
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[1]);
		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[0]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos3] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
		  }
		}
	      
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2);
		TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		if (TmpCoefficient != 0.0)
		{
		  coef2 *= this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[0]);
		  for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 2.0 * TmpValues[l] * TmpCoefficient * coef2;
		  
		  pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		  for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];		  
		}
		
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]);
		pos2 = this->Chain->JOffDiagonali(TmpIndex2, i, coef2);
		for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[0];
	      }
// 	      
	      //CB part
	      if ((this->NbrSpinX > 1) && ((this->NbrSpinY > 1) || (this->Offset != 0)) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		pos = this->Chain->SpiSmj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[1]);
		  
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[2]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[2]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos3] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
		  }
		}
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef);
		if (pos != dim)
		{
		  TmpCoefficient = coef * this->Chain->JDiagonali(TmpIndex1, pos, this->PseudospinDiagCouplingElements[1]);
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[2]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += TmpValues[l] * TmpCoefficient * coef2;
		  
		    pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
		  }
		
		  pos2 = this->Chain->JOffDiagonali(TmpIndex1, pos, coef2);
		  TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		  if (TmpCoefficient != 0.0)
		  {
		    coef2 = this->Chain->JDiagonali(TmpIndex2, pos2, this->PseudospinDiagCouplingElements[2]);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += TmpValues[l] * TmpCoefficient * coef2;
		    
		    pos3 = this->Chain->JOffDiagonali(TmpIndex2, pos2, coef2);
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos3] += TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
		  }
		}
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2);
		TmpCoefficient = coef * coef2 * this->PseudospinCouplingElements[1];
		if (TmpCoefficient != 0.0)
		{
		  coef2 *= this->Chain->JDiagonali(TmpIndex2, pos, this->PseudospinDiagCouplingElements[2]);
		  for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 2.0 * TmpValues[l] * TmpCoefficient * coef2;
		  
		  pos2 = this->Chain->JOffDiagonali(TmpIndex2, pos, coef2);
		  for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];		  
		}
		
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]);
		pos2 = this->Chain->JOffDiagonali(TmpIndex2, i, coef2);
		for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos2] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->PseudospinCouplingElements[2];
	      }
	    }
	  }
	}
  delete[] TmpValues;
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void TwoDimensionalTriangularLatticeWithPseudospinHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoefficient;
  int TmpIndex1;
  int TmpIndex2;
  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      //AC part
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
		TmpCoefficient = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient *= this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[2]);
		TmpCoefficient *= this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[0]);
		this->SzSzContributions[i] += TmpCoefficient;
	      }
	      
	      //AB part
	      if ((this->NbrSpinX > 1) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
		TmpCoefficient = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient *= this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]);
		TmpCoefficient *= this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[0]);
		this->SzSzContributions[i] += TmpCoefficient;
	      }

// 	      CB part
	      if ((this->NbrSpinX > 1) && ((this->NbrSpinY > 1) || (this->Offset != 0)) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		TmpCoefficient = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient *= this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]);
		TmpCoefficient *= this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[2]);
		this->SzSzContributions[i] += TmpCoefficient;
	      }
	    }
	}
    }
  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] *= this->JFactor;
    }
}

