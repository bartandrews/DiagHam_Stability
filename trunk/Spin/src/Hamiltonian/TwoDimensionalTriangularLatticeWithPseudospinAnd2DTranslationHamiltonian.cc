////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                         Class author: Cecile Repellin                      //
//                                                                            //
//                                                                            //
//       class of two dimensional spin model on the triangular lattice        //
//                   with pseudospin and translations                         //
//                                                                            //
//                     last modification : 08/12/2016                         //
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


#include "Hamiltonian/TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "Architecture/ArchitectureOperation/SUNSpinPrecalculationOperation.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;
using std::min;
using std::max;


// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// nbrSpinX = number of spin along the x direction
// nbrSpinY = number of spin along the y direction
// jFactor = amplitude of the Ising term
// hxFactor = amplitudes of the Zeeman term along x
// hzFactor = amplitudes of the Zeeman term along z
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian (Spin1_2ChainWithPseudospin* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, double jFactor, bool periodicBoundaryConditions, int offset, long memory)
{ 
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->JFactor = jFactor;
  this->JBreakD3Factor = 1.0;
  this->JBreakD3FactorDown = 1.0;
  this->Offset = offset;
  
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  
//   this->HermitianSymmetryFlag = true;
  this->HermitianSymmetryFlag = false;
    
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
  this->EvaluateExponentialFactors();
  this->EvaluateDiagonalMatrixElements();
  /*
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024l)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1l << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1l << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
      cout << endl;
      this->EnableFastMultiplication();
    }*/
}

// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// nbrSpinX = number of spin along the x direction
// nbrSpinY = number of spin along the y direction
// jFactor = amplitude of the Ising term
// hxFactor = amplitudes of the Zeeman term along x
// hzFactor = amplitudes of the Zeeman term along z
// periodicBoundaryConditions = true if periodic boundary conditions have to be used

TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian (Spin1_2ChainWithPseudospin* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, double jFactor, double jBreakD3Factor, bool periodicBoundaryConditions, int offset, long memory)
{ 
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->JFactor = jFactor;
  this->JBreakD3Factor = jBreakD3Factor;
  this->JBreakD3FactorDown = jBreakD3Factor;
  this->Offset = offset;
  
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  
//   this->HermitianSymmetryFlag = true;
  this->HermitianSymmetryFlag = false;
    
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
  
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateExponentialFactors();
  this->EvaluateDiagonalMatrixElements();

}

// destructor
//

TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::~TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian() 
{
  for (int i = 0; i < this->NbrSpinX; ++i)
    delete[] this->ExponentialFactors[i];
  delete[] this->ExponentialFactors;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (Spin1_2ChainWithPseudospinAnd2DTranslation*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::ShiftHamiltonian (double shift)
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

ComplexVector& TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
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
  int TmpIndex3;
  Complex TmpCoefficient;
  int nbrTranslationsX;
  int nbrTranslationsY;
  
  int nbrTranslationsX1;
  int nbrTranslationsY1;
  
  double* TmpOffDiagCoupling;
  if (this->JBreakD3Factor != 1.0)
  {
    TmpOffDiagCoupling = new double[3];
    TmpOffDiagCoupling[0] = sqrt(3) / 4.0;
    TmpOffDiagCoupling[1] = 0.0;
    TmpOffDiagCoupling[2] = -sqrt(3.0) / 4.0;
  }
  
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      Complex TmpValue = vSource[i] * this->JFactor * 0.5;
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      if (this->JBreakD3Factor != 1.0)
	      {
		TmpIndex1 = this->GetLinearizedIndex(j, k);
		pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		{
		  vDestination[pos] += vSource[i] * TmpOffDiagCoupling[0] * this->JBreakD3Factor * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		  vDestination[pos] += vSource[i] * TmpOffDiagCoupling[2] * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
	      }
	      
	      
	      //AC part
	     if ((this->NbrSpinY > 1) || (this->Offset != 0))
	     {
		TmpIndex1 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
		pos = this->Chain->SmiSpjJiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
// 		cout << pos << " " << dim << endl;
		if (pos != dim)
		  vDestination[pos] += TmpValue * this->JBreakD3FactorDown * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[2] * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		
		pos = this->Chain->SmiSpjJiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[2], this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * this->JBreakD3FactorDown * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] +=  TmpValue * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		

		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[0]) * this->PseudospinCouplingElements[2];
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos!= dim)
		    vDestination[pos] += 2.0 * TmpValue * this->JBreakD3FactorDown * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[2]) * this->PseudospinCouplingElements[0];
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos!= dim)
		    vDestination[pos] += 2.0 * TmpValue * this->JBreakD3FactorDown * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}  
		  
		pos2 = this->Chain->JoffiJoffj(TmpIndex1, TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		if (pos2!= dim)
		  vDestination[pos2] += 2.0 * TmpValue * this->JBreakD3FactorDown * coef * coef2 * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[2] *  this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		
	     }
	     
	      

// 	      AB part
	      if ((this->NbrSpinX > 1) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
		pos = this->Chain->SmiSpjJiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		

	      
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient = coef * this->PseudospinCouplingElements[1] * this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[0]);
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos != dim)
		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		
		TmpCoefficient = coef * this->PseudospinCouplingElements[0] *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]);
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos != dim)
		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		  
		pos2 = this->Chain->JoffiJoffj(TmpIndex1, TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		if (pos2!= dim)
		  vDestination[pos2] += 2.0 * TmpValue * coef * coef2 * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];  
		
		
	      }
	      
// 	      CB part
	      if ((this->NbrSpinX > 1) && ((this->NbrSpinY > 1) || (this->Offset != 0)) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		pos = this->Chain->SmiSpjJiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[2], this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];

		pos = this->Chain->SmiSpjJiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]) * this->PseudospinCouplingElements[2];
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos!= dim)
		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		  
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[2]) * this->PseudospinCouplingElements[1];
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos!= dim)
		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		  
		pos = this->Chain->JoffiJoffj(TmpIndex1, TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		if (pos!= dim)
		  vDestination[pos] += 2.0 * TmpValue * coef * coef2 * this->PseudospinCouplingElements[2] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];		  
		
		
		
	      }
	      
	    }
	}
    }
  if (this->JBreakD3Factor != 1.0)
    delete[] TmpOffDiagCoupling;
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

ComplexVector& TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
// //   cout << "Use Hermitian symmetry" << endl;
//   int LastComponent = firstComponent + nbrComponent;
//   int dim = this->Chain->GetHilbertSpaceDimension();
//   double coef;
//   double coef2;
//   int pos;
//   int pos2;
//   int pos3;
//   int MaxPos = this->NbrSpin - 1;
//   int TmpIndex1;
//   int TmpIndex2;
//   int TmpIndex3;
//   Complex TmpCoefficient;
//   Complex TmpSum;
//   int nbrTranslationsX;
//   int nbrTranslationsY;
//   
//   int nbrTranslationsX1;
//   int nbrTranslationsY1;
//   
//   for (int i = firstComponent; i < LastComponent; ++i)
//     {
//       TmpSum = 0.0;
//       Complex TmpValue = vSource[i] * this->JFactor * 0.5;
//       vDestination[i] += this->SzSzContributions[i] * vSource[i];
//       for (int j = 0; j < this->NbrSpinX; j++)
// 	{
// 	  for (int k = 0; k < this->NbrSpinY; k++)
// 	    {
// 	      //AC part
// 	     if ((this->NbrSpinY > 1) || (this->Offset != 0))
// 	     {
// 		TmpIndex1 = this->GetLinearizedIndex(j - this->Offset, k - 1);
// 		TmpIndex2 = this->GetLinearizedIndex(j, k);
// // 		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);		
// 		pos = this->Chain->SmiSpjJiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum += vSource[pos] * this->JFactor * 0.5 * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		if (i == 0 and pos == 0)
// 		  cout << i << " " << pos << " : " << TmpCoefficient << coef << endl;
// 		
// 		pos = this->Chain->SmiSpjJoffiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum += vSource[pos] * this->JFactor * 0.5 * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY] * this->PseudospinCouplingElements[0];
// 		}
// // 		cout << i << " " << pos << " : " << TmpCoefficient << coef << endl;
// 		
// 		pos = this->Chain->SmiSpjJiJoffj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY] * this->PseudospinCouplingElements[2];
// 		}
// 		if (i == 0 and pos == 0)
// 		  cout << i << " " << pos << " : " << coef << endl;
// 		
// 		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[2] * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY] * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[2] ;
// 		}
// // 		cout << i << " " << pos << " : " << TmpCoefficient << coef << endl;
// 		
// 		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
// 		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[0]) * this->PseudospinCouplingElements[2];
// 		if (TmpCoefficient != 0.0)
// 		{
// 		  pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
// 		  if (pos <= i)
// 		  {
// 		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		    if (pos != i)
// 		      TmpSum += vSource[pos] * this->JFactor * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  }
// 		if (i == 0 and pos == 0)
// 		  cout << i << " " << pos << " : " << TmpCoefficient << coef2 << endl;
// 		}
// 		
// 		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[2]) * this->PseudospinCouplingElements[0];
// 		if (TmpCoefficient != 0.0)
// 		{
// 		  pos = this->Chain->JOffDiagonali(TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
// 		  if (pos <= i)
// 		  {
// 		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		    if (pos != i)
// 		      TmpSum += vSource[pos] * this->JFactor * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  }
// // 		cout << i << " " << pos << endl;
// 		}  
// 		  
// 		pos = this->Chain->JoffiJoffj(TmpIndex1, TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += 2.0 * TmpValue * coef * coef2 * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[2] *  this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum += vSource[pos] * this->JFactor * coef * coef2 * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[2] *  this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// // 		cout << i << " " << pos << endl;
// 		
// 		
// 	     }
// 	     
// 	      
// 
// 	      //AB part
// 	      if ((this->NbrSpinX > 1) && (this->PeriodicBoundaryConditions || (j > 0)))
// 	      {
// 		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
// 		TmpIndex2 = this->GetLinearizedIndex(j, k);
// 		
// 		pos = this->Chain->SmiSpjJiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum += vSource[pos] * this->JFactor * 0.5 * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJoffiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum += vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJiJoffj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// // 		
// 		pos = this->Chain->SmiSpjJiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJoffiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJiJoffj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 *  coef * this->PseudospinCouplingElements[0] * this	->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 
// 		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
// 		TmpCoefficient = coef * this->PseudospinCouplingElements[1] * this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[0]);
// 		if (TmpCoefficient != 0.0)
// 		{
// 		  pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
// 		  if (pos <= i)
// 		  {
// 		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		    if (pos != i)
// 		      TmpSum +=  vSource[pos] * this->JFactor * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  }
// 		}
// 		
// 		TmpCoefficient = coef * this->PseudospinCouplingElements[0] *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]);
// 		if (TmpCoefficient != 0.0)
// 		{
// 		  pos = this->Chain->JOffDiagonali(TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
// 		  if (pos <= i)
// 		  {
// 		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		    if (pos != i)
// 		      TmpSum +=  vSource[pos] * this->JFactor * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  }
// 		}
// 		  
// 		pos = this->Chain->JoffiJoffj(TmpIndex1, TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += 2.0 * TmpValue * coef * coef2 * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];  
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * coef * coef2 * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];  
// 		}
// 		
// 		
// 	      }
// // 	      
// 	      //CB part
// 	      if ((this->NbrSpinX > 1) && ((this->NbrSpinY > 1) || (this->Offset != 0)) && (this->PeriodicBoundaryConditions || (j > 0)))
// 	      {
// 		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
// 		TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1);
// // 		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
// 		pos = this->Chain->SmiSpjJiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[2], this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJoffiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJiJoffj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 
// // 		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
// 		pos = this->Chain->SmiSpjJiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJoffiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJiJoffj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * 0.5 * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}
// 		
// 		
// 		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
// 		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]) * this->PseudospinCouplingElements[2];
// 		if (TmpCoefficient != 0.0)
// 		{
// 		  pos = this->Chain->JOffDiagonali(TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
// 		  if (pos <= i)
// 		  {
// 		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		    if (pos != i)
// 		      TmpSum +=  vSource[pos] * this->JFactor * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  }
// 		}
// 		  
// 		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[2]) * this->PseudospinCouplingElements[1];
// 		if (TmpCoefficient != 0.0)
// 		{
// 		  pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
// 		  if (pos <= i)
// 		  {
// 		    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		    if (pos != i)
// 		      TmpSum +=  vSource[pos] * this->JFactor *  TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  }
// 		}
// 		  
// 		pos = this->Chain->JoffiJoffj(TmpIndex1, TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
// 		if (pos <= i)
// 		{
// 		  vDestination[pos] += 2.0 * TmpValue * coef * coef2 * this->PseudospinCouplingElements[2] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		  if (pos != i)
// 		    TmpSum +=  vSource[pos] * this->JFactor * coef * coef2 * this->PseudospinCouplingElements[2] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 		}		
// 		
// 	      }
// 	      
// 	    }
// 	}
// 	vDestination[i] += Conj(TmpSum);
//     }
//   return vDestination;
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

ComplexVector* TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  double TmpCoefficient;
  int nbrTranslationsX;
  int nbrTranslationsY;
  Complex* TmpValues = new Complex[nbrVectors];
  
  double* TmpOffDiagCoupling;
  if (this->JBreakD3Factor != 1.0)
  {
    TmpOffDiagCoupling = new double[3];
    TmpOffDiagCoupling[0] = sqrt(3) / 4.0;
    TmpOffDiagCoupling[1] = 0.0;
    TmpOffDiagCoupling[2] = -sqrt(3.0) / 4.0;
  }
//   for (int l = 0; l < nbrVectors; ++l)
//     {
//       ComplexVector& TmpSource = vSources[l];
//       ComplexVector& TmpDestination = vDestinations[l];
//       for (int i = firstComponent; i < LastComponent; ++i)
// 	{
// 	  TmpDestination[i] += this->SzSzContributions[i] * TmpSource[i];
// 	}
//     }
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int l = 0; l < nbrVectors; ++l)
      {
	TmpValues[l] = vSources[l][i] * this->JFactor * 0.5;
	vDestinations[l][i] += this->SzSzContributions[i] * vSources[l][i];
      }
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      
	      if (this->JBreakD3Factor != 1.0)
	      {
		TmpIndex1 = this->GetLinearizedIndex(j, k);
		pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		{
		  for (int l = 0; l < nbrVectors; ++l)
		  {
		    vDestinations[l][pos] += vSources[l][i] * TmpOffDiagCoupling[0] * this->JBreakD3Factor * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		    vDestinations[l][pos] += vSources[l][i] * TmpOffDiagCoupling[2] * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		  }
		}
	      }
	      //AC part
	     if ((this->NbrSpinY > 1) || (this->Offset != 0))
	     {
		TmpIndex1 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
// 		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);		
		pos = this->Chain->SmiSpjJiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * this->JBreakD3FactorDown * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[2] * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		
// 		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
		pos = this->Chain->SmiSpjJiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[2], this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * this->JBreakD3FactorDown * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * this->JBreakD3FactorDown * coef * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		

		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[0]) * this->PseudospinCouplingElements[2];
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos!= dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 2.0 * TmpValues[l] * this->JBreakD3FactorDown * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[2]) * this->PseudospinCouplingElements[0];
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos!= dim)
		    for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += 2.0 * TmpValues[l] * this->JBreakD3FactorDown * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}  
		  
		pos2 = this->Chain->JoffiJoffj(TmpIndex1, TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		if (pos2!= dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos2] += 2.0 * TmpValues[l] * this->JBreakD3FactorDown * coef * coef2 * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[2] *  this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	     }
	     
	      

	      //AB part
	      if ((this->NbrSpinX > 1) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
// 		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
		pos = this->Chain->SmiSpjJiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
// 		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
		pos = this->Chain->SmiSpjJiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[0], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[0] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		

	      
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient = coef * this->PseudospinCouplingElements[1] * this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[0]);
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos != dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		
		TmpCoefficient = coef * this->PseudospinCouplingElements[0] *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]);
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos != dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		  
		pos2 = this->Chain->JoffiJoffj(TmpIndex1, TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		if (pos2!= dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos2] += 2.0 * TmpValues[l] * coef * coef2 * this->PseudospinCouplingElements[0] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];  
		
		
	      }
// 	      
	      //CB part
	      if ((this->NbrSpinX > 1) && ((this->NbrSpinY > 1) || (this->Offset != 0)) && (this->PeriodicBoundaryConditions || (j > 0)))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - 1, k);
		TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1);
// 		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
		pos = this->Chain->SmiSpjJiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[2], this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex2, TmpIndex1, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex2, TmpIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];

// 		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
		pos = this->Chain->SmiSpjJiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[2], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJiJoffj(TmpIndex1, TmpIndex2, i, this->PseudospinDiagCouplingElements[1], coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		pos = this->Chain->SmiSpjJoffiJoffj(TmpIndex1, TmpIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
		if (pos != dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += TmpValues[l] * coef * this->PseudospinCouplingElements[1] * this->PseudospinCouplingElements[2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
		
		coef = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[1]) * this->PseudospinCouplingElements[2];
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos!= dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		  
		TmpCoefficient = coef *  this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[2]) * this->PseudospinCouplingElements[1];
		if (TmpCoefficient != 0.0)
		{
		  pos = this->Chain->JOffDiagonali(TmpIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
		  if (pos!= dim)
		    for (int l = 0; l < nbrVectors; ++l)
		      vDestinations[l][pos] += 2.0 * TmpValues[l] * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		}
		  
		pos = this->Chain->JoffiJoffj(TmpIndex1, TmpIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
		if (pos!= dim)
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][pos] += 2.0 * TmpValues[l] * coef * coef2 * this->PseudospinCouplingElements[2] * this->PseudospinCouplingElements[1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	      }
	    }
	}
    }
  delete[] TmpValues;
  if (this->JBreakD3Factor != 1.0)
    delete[] TmpOffDiagCoupling;
  return vDestinations;
}



// evaluate diagonal matrix elements
// 

void TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoefficient;
  double** TmpDiagCoupling;
  if (this->JBreakD3Factor != 1.0)
  {
    TmpDiagCoupling = new double*[3];
    for (int l = 0; l < 3; ++l)
      TmpDiagCoupling[l] = new double[2];
    TmpDiagCoupling[0][0] = 0.0;
    TmpDiagCoupling[0][1] = -0.5;  
    TmpDiagCoupling[1][0] = -0.75;
    TmpDiagCoupling[1][1] = 0.25;
    TmpDiagCoupling[2][0] = 0.0;
    TmpDiagCoupling[2][1] = -0.5;
  }
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] = 0.0;
      
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      if (this->JBreakD3Factor != 1.0)
	      {
		TmpIndex1 = this->GetLinearizedIndex(j, k);
		
		//AB up
		TmpCoefficient = this->Chain->JDiagonali(TmpIndex1, i, TmpDiagCoupling[0]);
		this->SzSzContributions[i] += TmpCoefficient * this->JBreakD3Factor;
		
		//BC up
		TmpCoefficient = this->Chain->JDiagonali(TmpIndex1, i, TmpDiagCoupling[1]);
		this->SzSzContributions[i] += TmpCoefficient * this->JBreakD3Factor;
		
		//AC up
		TmpCoefficient = this->Chain->JDiagonali(TmpIndex1, i, TmpDiagCoupling[2]);
		this->SzSzContributions[i] += TmpCoefficient;
	      }
	      
// 	      AC part
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		TmpIndex1 = this->GetLinearizedIndex(j - this->Offset, k - 1);
		TmpIndex2 = this->GetLinearizedIndex(j, k);
		TmpCoefficient = this->Chain->SziSzj(TmpIndex1, TmpIndex2, i);
		TmpCoefficient *= this->Chain->JDiagonali(TmpIndex1, i, this->PseudospinDiagCouplingElements[2]);
		TmpCoefficient *= this->Chain->JDiagonali(TmpIndex2, i, this->PseudospinDiagCouplingElements[0]);
		this->SzSzContributions[i] += TmpCoefficient  * this->JBreakD3Factor;
	      }
	      
// 	      AB part
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
    if (this->JBreakD3Factor != 1.0)
    {
      for (int i = 0; i < 3; ++i)
	delete[] TmpDiagCoupling[i];
      delete[] TmpDiagCoupling;
    }
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

// ComplexVector& TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::HermitianLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
// 										   int firstComponent, int nbrComponent)
// {
//   int LastComponent = firstComponent + nbrComponent;
//   int Dim = this->Particles->GetHilbertSpaceDimension();
//   double Coefficient;
//   double TmpCoefficient;
//   double TmpCoef;
//   double ChargeContribution;
//   int NbrElements;
//   int TmpTotalNbrParticles;
//   AbstractSpinChain* TmpParticles = (AbstractSpinChain*) this->Particles->Clone();
//   
//   int* TmpIndexArray;
//   Complex* TmpCoefficientArray; 
//   int j;
//   int TmpNbrInteraction;
//   firstComponent -= this->PrecalculationShift;
//   LastComponent -= this->PrecalculationShift;
//   int Pos = firstComponent / this->FastMultiplicationStep; 
//   int PosMod = firstComponent % this->FastMultiplicationStep;
//   if (PosMod != 0)
//     {
//       ++Pos;
//       PosMod = this->FastMultiplicationStep - PosMod;
//     }
//   int l =  PosMod + firstComponent + this->PrecalculationShift;
//   for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
//     {
//       TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
//       TmpIndexArray = this->InteractionPerComponentIndex[Pos];
//       TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
//       Coefficient = vSource[l];
//       for (j = 0; j < TmpNbrInteraction; ++j)
// 	vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
//       vDestination[l] += this->SzSzContributions[l] * Coefficient;
//       l += this->FastMultiplicationStep;
//       ++Pos;
//     }
// 
//   int Index;  
//   Complex TmpSum;
//   int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
//   firstComponent += this->PrecalculationShift;
//   LastComponent += this->PrecalculationShift;
//   int TmpIndex1;
//   int TmpIndex2;
//   int TmpIndex3;
//   for (int k = 0; k < this->FastMultiplicationStep; ++k)
//     if (PosMod != k)
//       {
// 	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
// 	{
// 	  TmpSum = 0.0;
// 	  for (int j = 0; j < this->NbrSpinX; j++)
// 	  {
// 	    for (int k = 0; k < this->NbrSpinY; k++)
// 	      {
// 		//upward triangles
// 		TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
// 		TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
// 		TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
// 		//AB
// 		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
// 		if (pos < i)
// 		{
// 		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
// 		  TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
// 		}
// 		//AC
// 		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
// 		if (pos < i)
// 		{
// 		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
// 		  TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
// 		}
// 		//BC
// 		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
// 		if (pos < i)
// 		{
// 		  vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
// 		  TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
// 		}
// 	      
// 		//downward triangles
// 		if (this->NbrSpinX > 1)
// 		{
// 		  //BA
// 		  TmpIndex1 = this->GetLinearizedIndex(j - 1, k, 1);
// 		  TmpIndex2 = this->GetLinearizedIndex(j, k, 0);
// 		  pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
// 		  if (pos < i)
// 		  {
// 		    vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
// 		    TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
// 		  }
// 		}
// 	      
// 		if ((this->NbrSpinY > 1) || (this->Offset != 0))
// 		{
// 		  //CA
// 		  TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
// 		  TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1, 2);
// 		  pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
// 		  if (pos < i)
// 		  {
// 		    vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
// 		    TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
// 		  }
// 		  if (this->NbrSpinX > 1)
// 		  {
// 		    //BC
// 		    TmpIndex1 = this->GetLinearizedIndex(j - 1, k, 1);
// 		    TmpIndex2 = this->GetLinearizedIndex(j - this->Offset, k - 1, 2);
// 		    pos = this->Chain->SmiSpj(min(TmpIndex1, TmpIndex2), max(TmpIndex1, TmpIndex2), i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
// 		    if (pos < i)
// 		    {
// 		      vDestination[pos] += 0.5 * vSource[i] * TmpCoefficient * this->JDownEasyPlaneFactor * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY];
// 		      TmpSum += 0.5 * vSource[pos] * TmpCoefficient * this->JDownEasyPlaneFactor * Conj(this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY]);
// 		    }
// 		  }
// 		}
// 	      
// 	      }
// 	  }
// 	  vDestination[i] += (this->SzSzContributions[i] * vSource[i] + TmpSum);
// 	}
//       }
      
// }
/*
// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted
// return value = number of non-zero matrix element

long TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  double TmpCoefficient;
  long Memory = 0l;
  int NbrTranslationsX;
  int NbrTranslationsY;
  
  AbstractSpinChain* TmpParticles = (AbstractSpinChain*) this->Particles->Clone();
  int LastComponent = nbrComponent + firstComponent;
  int TmpIndex1;
  int TmpIndex2;
  int TmpIndex3;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->NbrInteractionPerComponent[i - this->PrecalculationShift] = 0;
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      //upward triangles
	      TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
	      TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
	      TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
	      //AB
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      //AC
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      //BC
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	      }
	      
	      //downward triangles
	      if (this->NbrSpinX > 1)
	      {
		//BA
		TmpIndex1 = min(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1));
		TmpIndex2 = max(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1));
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  if ((this->HermitianSymmetryFlag == false) || (pos < i))
		  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
		}
		if (this->HermitianSymmetryFlag == false)
		{
		  pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
		}
	      }
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//CA
		TmpIndex1 = min(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0));
		TmpIndex2 = max(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0));
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  if ((this->HermitianSymmetryFlag == false) || (pos < i))
		  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
		}
		if (this->HermitianSymmetryFlag == false)
		{
		  pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
		}
		if (this->NbrSpinX > 1)
		{
		  //BC
		  TmpIndex1 = min(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1));
		  TmpIndex2 = max(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1));
		  pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    if ((this->HermitianSymmetryFlag == false) || (pos < i))
		    {
		      ++Memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  }
		  if (this->HermitianSymmetryFlag == false)
		  {
		    pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		    if (pos != dim)
		    {
		      ++Memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		  }
		}
	      }
	      
	    }
	}
    }    
  delete TmpParticles;
  return Memory;
}


// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  double TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int LastComponent = nbrComponent + firstComponent;
  AbstractSpinChain* TmpParticles = (AbstractSpinChain*) this->Particles->Clone();
  
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
    
  int CurrentNbrCounting = 0;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      //upward triangles
	      TmpIndex1 = this->GetLinearizedIndex(j, k, 0);
	      TmpIndex2 = this->GetLinearizedIndex(j, k, 1);
	      TmpIndex3 = this->GetLinearizedIndex(j, k, 2);
	      //AB
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      //AC
	      pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      //BC
	      pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex3, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
	      if (pos != dim)
	      {
		if ((this->HermitianSymmetryFlag == false) || (pos < i))
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      if (this->HermitianSymmetryFlag == false)
	      {
		pos = this->Chain->SmiSpj(TmpIndex3, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  this->InteractionPerComponentIndex[Pos][position] = pos;
		  this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JEasyPlaneFactor;
		  ++position;
		}
	      }
	      
	      //downward triangles
	      if (this->NbrSpinX > 1)
	      {
		//BA
		TmpIndex1 = min(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1));
		TmpIndex2 = max(this->GetLinearizedIndex(j, k, 0), this->GetLinearizedIndex(j - 1, k, 1));
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  if ((this->HermitianSymmetryFlag == false) || (pos < i))
		  {
		    this->InteractionPerComponentIndex[Pos][position] = pos;
		    this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		    ++position;
		  }
		}
		if (this->HermitianSymmetryFlag == false)
		{
		  pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    this->InteractionPerComponentIndex[Pos][position] = pos;
		    this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		    ++position;
		  }
		}
	      }
	      if ((this->NbrSpinY > 1) || (this->Offset != 0))
	      {
		//CA
		TmpIndex1 = min(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0));
		TmpIndex2 = max(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j, k, 0));
		pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		if (pos != dim)
		{
		  if ((this->HermitianSymmetryFlag == false) || (pos < i))
		  {
		    this->InteractionPerComponentIndex[Pos][position] = pos;
		    this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		    ++position;
		  }
		}
		if (this->HermitianSymmetryFlag == false)
		{
		  pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    this->InteractionPerComponentIndex[Pos][position] = pos;
		    this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		    ++position;
		  }
		}
		if (this->NbrSpinX > 1)
		{
		  //BC
		  TmpIndex1 = min(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1));
		  TmpIndex2 = max(this->GetLinearizedIndex(j - this->Offset, k - 1, 2), this->GetLinearizedIndex(j - 1, k, 1));
		  pos = this->Chain->SmiSpj(TmpIndex1, TmpIndex2, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		  if (pos != dim)
		  {
		    if ((this->HermitianSymmetryFlag == false) || (pos < i))
		    {
		      this->InteractionPerComponentIndex[Pos][position] = pos;
		      this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		      ++position;
		    }
		  }
		  if (this->HermitianSymmetryFlag == false)
		  {
		    pos = this->Chain->SmiSpj(TmpIndex2, TmpIndex1, i, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
		    if (pos != dim)
		    {
		      this->InteractionPerComponentIndex[Pos][position] = pos;
		      this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficients * this->JDownEasyPlaneFactor;
		      ++position;
		    }
		  }
		}
	      }
	      
	    }
	}
      ++Pos;
    }
}*/

// ask if Hamiltonian implements hermitian symmetry operations
//

bool TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// evaluate all exponential factors
//   

void TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex*[this->NbrSpinX];
  for (int i = 0; i < this->NbrSpinX; ++i)
    { 
      this->ExponentialFactors[i] = new Complex[this->NbrSpinY];
      for (int j = 0; j < this->NbrSpinY; ++j)
	{ 
	  this->ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->NbrSpinX))
							       + (this->YMomentum * ((double) j) / ((double) this->NbrSpinY))));
	}
    }
}


// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

// long TwoDimensionalTriangularLatticeWithPseudospinAnd2DTranslationHamiltonian::FastMultiplicationMemory(long allowedMemory)
// {
//   this->NbrInteractionPerComponent = new int [this->Chain->GetHilbertSpaceDimension()];
//   for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); ++i)
//     this->NbrInteractionPerComponent[i] = 0;
//   timeval TotalStartingTime2;
//   timeval TotalEndingTime2;
//   double Dt2;
//   gettimeofday (&(TotalStartingTime2), 0);
//   cout << "start" << endl;
//  
//   SUNSpinPrecalculationOperation Operation(this);
//   Operation.ApplyOperation(this->Architecture);
// 
//   long Memory = 0;
//   for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); ++i)
//     Memory += this->NbrInteractionPerComponent[i];
// 
//   cout << "nbr interaction = " << Memory << endl;
//   long TmpMemory = allowedMemory - ((2 * sizeof (int*)) + sizeof (int) + sizeof(double*)) * this->Chain->GetHilbertSpaceDimension();
//   if ((TmpMemory < 0) || ((TmpMemory / ((int) ((2 * sizeof (int)) + sizeof(double)))) < Memory))
//     {
//       this->FastMultiplicationStep = 1;
//       int ReducedSpaceDimension  = this->Chain->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
//       while ((TmpMemory < 0) || ((TmpMemory / ((int) ((2 * sizeof (int)) + sizeof(double)))) < Memory))
// 	{
// 	  ++this->FastMultiplicationStep;
// 	  ReducedSpaceDimension = this->Chain->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
// 	  if (this->Chain->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
// 	    ++ReducedSpaceDimension;
// 	  TmpMemory = allowedMemory - ((2 * sizeof (int*)) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension;
// 	  Memory = 0;
// 	  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i += this->FastMultiplicationStep)
// 	    Memory += this->NbrInteractionPerComponent[i];
// 	}
//       int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
//       for (int i = 0; i < ReducedSpaceDimension; ++i)
// 	TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
//       delete[] this->NbrInteractionPerComponent;
//       this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
//       Memory = (((2 * sizeof (int*)) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension) + (Memory * ((2 * sizeof (int)) + sizeof(double)));
//     }
//   else
//     {
//       Memory = ((((2 * sizeof (int*)) + sizeof (int) + sizeof(double*)) * this->Chain->GetHilbertSpaceDimension()) + 
// 		(Memory * ((2 * sizeof (int)) + sizeof(double))));
//       this->FastMultiplicationStep = 1;
//     }
// 
//   cout << "reduction factor=" << this->FastMultiplicationStep << endl;
//   gettimeofday (&(TotalEndingTime2), 0);
//   cout << "------------------------------------------------------------------" << endl << endl;;
//   Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
//     ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
//   cout << "time = " << Dt2 << endl;
//   return Memory;
// }