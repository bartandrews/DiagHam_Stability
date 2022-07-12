////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2003 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 24/11/2003                        //
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


#include "config.h"
#include "Hamiltonian/PeriodicQuantumDots3DHamiltonian.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"

#include <iostream>
#include <math.h>

using std::ostream;
using std::cout;
using std::endl;


#define PERIODIC_HAMILTONIAN_FACTOR 150.4
#define BLOCH_FACTOR 7.644


// constructor from data
//
// space = Hilbert space
// xSize = the sample length in X direction
// ySize = the sample length in Y direction
// zSize = the sample length in Z direction
// mux = effective mass in X direction
// muy = effective mass in Y direction
// muz = effective mass in Z direction
// nbrCellX = number of steps in X direction
// nbrCellY = number of steps in Y direction
// nbrCellZ = number of steps in Z direction
// PotentialInput = pointer to a 3D potential with constant value in a cell
// waveVectorX = wave vector of Bloch function in X direction
// waveVectorY = wave vector of Bloch function in Y direction
// waveVectorZ = wave vector of Bloch function in Z direction

PeriodicQuantumDots3DHamiltonian::PeriodicQuantumDots3DHamiltonian(PeriodicThreeDOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDConstantCellPotential* PotentialInput, double waveVectorX, double waveVectorY, double waveVectorZ)
{
  this->Space = space;
  this->XSize = xSize;
  this->YSize = ySize;
  this->ZSize = zSize;
  this->Mux = mux;
  this->Muy = muy;
  this->Muz = muz;
  this->NbrCellX = nbrCellX;
  this->NbrCellY = nbrCellY;
  this->NbrCellZ = nbrCellZ;
  this->NbrStateX = this->Space->GetNbrStateX();
  this->LowerImpulsionX = this->Space->GetLowerImpulsionX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->LowerImpulsionY = this->Space->GetLowerImpulsionY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();
  this->InteractionFactors = new double** [this->NbrCellZ];  
  for (int k = 0; k < this->NbrCellZ; ++k)
    {
      this->InteractionFactors[k] = new double* [this->NbrCellY];
      for (int j = 0; j < this->NbrCellY; ++j)
	{
	  this->InteractionFactors[k][j] = new double [this->NbrCellX];
	  for (int i = 0; i < this->NbrCellX; ++i)
	    this->InteractionFactors[k][j][i] = PotentialInput->GetPotential(i, j, k);
	}
    }
  this->EvaluateInteractionFactors(waveVectorX, waveVectorY, waveVectorZ);
}


// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

PeriodicQuantumDots3DHamiltonian::PeriodicQuantumDots3DHamiltonian(const PeriodicQuantumDots3DHamiltonian& hamiltonian)
{
  this->Space = hamiltonian.Space;
  this->XSize = hamiltonian.XSize;
  this->YSize = hamiltonian.YSize;
  this->ZSize = hamiltonian.ZSize;
  this->Mux = hamiltonian.Mux;
  this->Muy = hamiltonian.Muy;
  this->Muz = hamiltonian.Muz;
  this->NbrCellX = hamiltonian.NbrCellX;
  this->NbrCellY = hamiltonian.NbrCellY;
  this->NbrCellZ = hamiltonian.NbrCellZ;
  this->NbrStateX = this->Space->GetNbrStateX();
  this->LowerImpulsionX = this->Space->GetLowerImpulsionX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->LowerImpulsionY = this->Space->GetLowerImpulsionY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();
  this->KineticElements = hamiltonian.KineticElements;
  this->NbrPrecalculatedDimension = hamiltonian.NbrPrecalculatedDimension;
  this->InteractionFactors = hamiltonian.InteractionFactors;
  this->RealWaveFunctionOverlapX = hamiltonian.RealWaveFunctionOverlapX;
  this->ImaginaryWaveFunctionOverlapX = hamiltonian.ImaginaryWaveFunctionOverlapX;
  this->RealWaveFunctionOverlapY = hamiltonian.RealWaveFunctionOverlapY;
  this->ImaginaryWaveFunctionOverlapY = hamiltonian.ImaginaryWaveFunctionOverlapY;
  this->RealWaveFunctionOverlapZ = hamiltonian.RealWaveFunctionOverlapZ;
  this->ImaginaryWaveFunctionOverlapZ = hamiltonian.ImaginaryWaveFunctionOverlapZ;
  this->RealPrecalculatedHamiltonian =  hamiltonian.RealPrecalculatedHamiltonian;
  this->ImaginaryPrecalculatedHamiltonian =  hamiltonian.ImaginaryPrecalculatedHamiltonian;
}

// destructor
//

PeriodicQuantumDots3DHamiltonian::~ PeriodicQuantumDots3DHamiltonian()
{
  delete[] this->KineticElements;
  delete[] this->InteractionFactors;
  delete[] this->RealWaveFunctionOverlapX;
  delete[] this->ImaginaryWaveFunctionOverlapX;
  delete[] this->RealWaveFunctionOverlapY;
  delete[] this->ImaginaryWaveFunctionOverlapY;
  delete[] this->RealWaveFunctionOverlapZ;
  delete[] this->ImaginaryWaveFunctionOverlapZ;
  delete[] this->RealPrecalculatedHamiltonian;
  delete[] this->ImaginaryPrecalculatedHamiltonian;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* PeriodicQuantumDots3DHamiltonian::Clone ()
{
  return new PeriodicQuantumDots3DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PeriodicQuantumDots3DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void PeriodicQuantumDots3DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->KineticElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicQuantumDots3DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
{
  double x = 0.0;
  int dim = this->Space->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicQuantumDots3DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicQuantumDots3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Space->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of idinces 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PeriodicQuantumDots3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
						       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination.Re(i) = 0.0;
      vDestination.Im(i) = 0.0;
    }
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicQuantumDots3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  int OriginX = this->NbrStateX - 1; int OriginY = this->NbrStateY - 1; int OriginZ = this->NbrStateZ - 1;

  int m1, m2, n1, n2, p1;
  int IndexX, IndexY, IndexZ;
  int** TotalIndex = new int* [this->NbrStateX]; int TmpIndex = 0;
  for (m1 = 0; m1 < this->NbrStateX; ++m1) 
    {
      TotalIndex[m1] = new int [this->NbrStateY];
      for (n1 = 0; n1 < this->NbrStateY; ++n1)	
	{
	  TotalIndex[m1][n1] = (m1 * this->NbrStateY + n1) * this->NbrStateZ;
	  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	    {	      
	      vDestination.Re(TmpIndex) += vSource.Re(TmpIndex) * this->KineticElements[TmpIndex];
	      vDestination.Im(TmpIndex) += vSource.Im(TmpIndex) * this->KineticElements[TmpIndex];
	      ++TmpIndex;
	    }
	}
    }

  int Index1, Index2;
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;
  double TmpRe = 0.0; double TmpIm = 0.0;
  int* TmpTotalIndex1; int* TmpTotalIndex2;
  int LimitZ = 0; int LengthZ = (this->NbrStateZ - 1) * 2 + 1;   
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {
      IndexX = -m1 + OriginX;
      TmpTotalIndex1 = TotalIndex[m1];
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];	  
	  for (n1 = 0; n1 < this->NbrStateY; ++n1)
	    {
	      IndexY = -n1 + OriginY;
	      for (n2 = 0; n2 < this->NbrStateY; ++n2)
		{
		  TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		  Index1 = TmpTotalIndex1[n1];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      IndexZ = -p1 + OriginZ;
		      TmpRe = 0.0; TmpIm = 0.0;
		      Index2 = TmpTotalIndex2[n2];
		      LimitZ = LengthZ - p1;
		      for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
			{
			  TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
			  TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]);  	  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		    }
   		  ++IndexY;
		}
	    }
	  ++IndexX;
	}
    }
  delete[] TotalIndex;
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

ComplexVector& PeriodicQuantumDots3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{ 
  if ((firstComponent == 0) && (nbrComponent == this->Space->GetHilbertSpaceDimension()))
    return this->LowLevelAddMultiply(vSource, vDestination);
  else
    {  
      int lastComponent = firstComponent + nbrComponent;
      int OriginX = this->NbrStateX - 1; int OriginY = this->NbrStateY - 1; int OriginZ = this->NbrStateZ - 1;
      int m1, m2, n1, n2, p1;
      int IndexX, IndexY, IndexZ;
      double* TmpRealPrecalculatedHamiltonian;
      double* TmpImaginaryPrecalculatedHamiltonian;
      double TmpRe = 0.0; double TmpIm = 0.0;
      
      int Index1 = firstComponent; int Index2 = 0;
      int ReducedIndex1 = firstComponent / this->NbrStateZ;
      int p1Begin = firstComponent - ReducedIndex1 * this->NbrStateZ;
      int m1Begin = ReducedIndex1 / this->NbrStateY;
      int n1Begin = ReducedIndex1 - m1Begin * this->NbrStateY;

      ReducedIndex1 = (lastComponent - 1) / this->NbrStateZ;
      int p1Limit = (lastComponent - 1) - ReducedIndex1 * this->NbrStateZ;
      int m1Limit = ReducedIndex1 / this->NbrStateY;
      int n1Limit = ReducedIndex1 - m1Limit * this->NbrStateY;

      int** TotalIndex = new int* [this->NbrStateX]; int TmpIndex = 0;
      for (int m = 0; m < this->NbrStateX; ++m)
        {
          TotalIndex[m] = new int [this->NbrStateY];
          for (int n = 0; n < this->NbrStateY; ++n)
            TotalIndex[m][n] = (m * this->NbrStateY + n) * this->NbrStateZ;
        }
      int* TmpTotalIndex1; int* TmpTotalIndex2;

      int LimitZ = 0; int LengthZ = (this->NbrStateZ - 1) * 2 + 1;    
      for (Index2 = firstComponent; Index2 < lastComponent; ++Index2)
	{
	  vDestination.Re(Index2) += vSource.Re(Index2) * this->KineticElements[Index2];
	  vDestination.Im(Index2) += vSource.Im(Index2) * this->KineticElements[Index2];
	}
      /*
      // p1 : p1Begin -> NbrStateZ
      for (p1 = p1Begin; p1 < this->NbrStateZ; ++p1)
	{	
	  TmpRe = 0.0; TmpIm = 0.0;
	  Index2 = 0;
	  IndexX = -m1Begin + OriginX;
	  for (m2 = 0; m2 < this->NbrStateX; ++m2)
	    {
	      IndexY = -n1Begin + OriginY;
	      for (n2 = 0; n2 < this->NbrStateY; ++n2)
		{
		  TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		  IndexZ = -p1 + OriginZ;
		  LimitZ = LengthZ - p1;
		  for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
		    {
		      TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
		      TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]); 
		    }
		  ++IndexY;
		}
	      ++IndexX;
	    }
	  vDestination.Re(Index1) += TmpRe;
	  vDestination.Im(Index1) += TmpIm;
	  ++Index1;
	}
      ++n1Begin;
      // n1 : n1Begin -> NbrStateY 
      for (n1 = n1Begin; n1 < this->NbrStateY; ++n1)	
	for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	  {	
	    TmpRe = 0.0; TmpIm = 0.0;
	    Index2 = 0;
	    IndexX = -m1Begin + OriginX;
	    for (m2 = 0; m2 < this->NbrStateX; ++m2)
	      {
		IndexY = -n1 + OriginY;
		for (n2 = 0; n2 < this->NbrStateY; ++n2)
		  {
		    TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		    TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		    IndexZ = -p1 + OriginZ;
		    LimitZ = LengthZ - p1;
		    for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
		      {
			TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
			TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]); 
		      }
		    ++IndexY;
		  }
		++IndexX;
	      }
	    vDestination.Re(Index1) += TmpRe;
	    vDestination.Im(Index1) += TmpIm;
	    ++Index1;
	  }    
      ++m1Begin;
      // m1 : m1Begin -> m1Limit - 1
      for (m1 = m1Begin; m1 < m1Limit; ++m1)		  
	for (n1 = 0; n1 < this->NbrStateY; ++n1)	  
	  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	    {	
	      TmpRe = 0.0; TmpIm = 0.0;
	      Index2 = 0;
	      IndexX = -m1 + OriginX;
	      for (m2 = 0; m2 < this->NbrStateX; ++m2)
		{
		  IndexY = -n1 + OriginY;
		  for (n2 = 0; n2 < this->NbrStateY; ++n2)
		    {
		      TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		      TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		      IndexZ = -p1 + OriginZ;
		      LimitZ = LengthZ - p1;
		      for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
			{
			  TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
			  TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]); 
			}
		      ++IndexY;
		    }
		  ++IndexX;
		}
	      vDestination.Re(Index1) += TmpRe;
	      vDestination.Im(Index1) += TmpIm;
	      ++Index1;
	    }	
      // m1 = m1Limit, n1 : 0 -> n1Limit - 1
      for (n1 = 0; n1 < n1Limit; ++n1)	
	for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	  { 
	    TmpRe = 0.0; TmpIm = 0.0;
	    Index2 = 0;
	    IndexX = -m1Limit + OriginX;
	    for (m2 = 0; m2 < this->NbrStateX; ++m2)
	      {
		IndexY = -n1 + OriginY;
		for (n2 = 0; n2 < this->NbrStateY; ++n2)
		  {
		    TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		    TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		    IndexZ = -p1 + OriginZ;
		    LimitZ = LengthZ - p1;
		    for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
		      {
			TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
			TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]); 
		      }
		    ++IndexY;
		  }
		++IndexX;
	      }
	    vDestination.Re(Index1) += TmpRe;
	      vDestination.Im(Index1) += TmpIm;
	      ++Index1;
	  }	
      // m1 = m1Limit, n1 = n1Limit, p1 : 0 -> p1Limit
      for (p1 = 0; p1 <= p1Limit; ++p1)
	{	
	  TmpRe = 0.0; TmpIm = 0.0;
	  Index2 = 0;
	  IndexX = -m1Limit + OriginX;
	  for (m2 = 0; m2 < this->NbrStateX; ++m2)
	    {
	      IndexY = -n1Limit + OriginY;
	      for (n2 = 0; n2 < this->NbrStateY; ++n2)
		{
		  TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		  IndexZ = -p1 + OriginZ;
		  LimitZ = LengthZ - p1;
		  for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
		    {
		      TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
		      TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]); 
		    }
		  ++IndexY;
		}
	      ++IndexX;
	    }
	  vDestination.Re(Index1) += TmpRe;
	  vDestination.Im(Index1) += TmpIm;
	  ++Index1;
	}
      */
      
      // p1 : p1Begin -> NbrStateZ
      IndexX = -m1Begin + OriginX;
      TmpTotalIndex1 = TotalIndex[m1Begin];
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];	  	  
	  IndexY = -n1Begin + OriginY;
	  for (n2 = 0; n2 < this->NbrStateY; ++n2)
	    {
	      TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
	      TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
	      Index1 = firstComponent;
	      for (p1 = p1Begin; p1 < this->NbrStateZ; ++p1)
		{
		  IndexZ = -p1 + OriginZ;
		  TmpRe = 0.0; TmpIm = 0.0;
		  Index2 = TmpTotalIndex2[n2];
		  LimitZ = LengthZ - p1;
		  for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
		    {
		      TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
		      TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]);  	  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}
	      ++IndexY;
	    }
	  ++IndexX;
	}      
      // n1 : n1Begin + 1 -> NbrStateY
      IndexX = -m1Begin + OriginX;
      TmpTotalIndex1 = TotalIndex[m1Begin];
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];	  
	  for (n1 = n1Begin + 1; n1 < this->NbrStateY; ++n1)
	    {
	      IndexY = -n1 + OriginY;
	      for (n2 = 0; n2 < this->NbrStateY; ++n2)
		{
		  TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		  Index1 = TmpTotalIndex1[n1];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      IndexZ = -p1 + OriginZ;
		      TmpRe = 0.0; TmpIm = 0.0;
		      Index2 = TmpTotalIndex2[n2];
		      LimitZ = LengthZ - p1;
		      for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
			{
			  TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
			  TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]);  	  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		    }
		  ++IndexY;		  
		}
	    }
	  ++IndexX;
	}      
      // m1 : m1Begin -> m1Limit - 1
      for (m1 = m1Begin + 1; m1 < m1Limit; ++m1)
	{
	  IndexX = -m1 + OriginX;
	  TmpTotalIndex1 = TotalIndex[m1];
	  for (m2 = 0; m2 < this->NbrStateX; ++m2)
	    {
	      TmpTotalIndex2 = TotalIndex[m2];	  
	      for (n1 = 0; n1 < this->NbrStateY; ++n1)
		{
		  IndexY = -n1 + OriginY;
		  for (n2 = 0; n2 < this->NbrStateY; ++n2)
		    {
		      TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		      TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		      Index1 = TmpTotalIndex1[n1];
		      for (p1 = 0; p1 < this->NbrStateZ; ++p1)
			{
			  IndexZ = -p1 + OriginZ;
			  TmpRe = 0.0; TmpIm = 0.0;
			  Index2 = TmpTotalIndex2[n2];
			  LimitZ = LengthZ - p1;
			  for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
			    {
			      TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
			      TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]);  	  
			    }
			  vDestination.Re(Index1) += TmpRe;
			  vDestination.Im(Index1) += TmpIm;
			  ++Index1;
			}
		      ++IndexY;
		    }
		}
	      ++IndexX;
	    }
	}
      // m1 = m1Limit, n1 : 0 -> n1Limit - 1
      IndexX = -m1Limit + OriginX;
      TmpTotalIndex1 = TotalIndex[m1Limit];
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];	  
	  for (n1 = 0; n1 < n1Limit; ++n1)
	    {
	      IndexY = -n1 + OriginY;
	      for (n2 = 0; n2 < this->NbrStateY; ++n2)
		{
		  TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		  Index1 = TmpTotalIndex1[n1];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      IndexZ = -p1 + OriginZ;
		      TmpRe = 0.0; TmpIm = 0.0;
		      Index2 = TmpTotalIndex2[n2];
		      LimitZ = LengthZ - p1;
		      for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
			{
			  TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
			  TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]);  	  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		    }
		  ++IndexY;		  
		}
	    }
	  ++IndexX;
	}
      // m1 = m1Limit, n1 = n1Limit, p1 : 0 -> p1Limit
      IndexX = -m1Limit + OriginX;
      TmpTotalIndex1 = TotalIndex[m1Limit];
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];	  	  
	  IndexY = -n1Limit + OriginY;
	  for (n2 = 0; n2 < this->NbrStateY; ++n2)
	    {
	      TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
	      TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
	      Index1 = TmpTotalIndex1[n1Limit];
	      for (p1 = 0; p1 <= p1Limit; ++p1)
		{
		  IndexZ = -p1 + OriginZ;
		  TmpRe = 0.0; TmpIm = 0.0;
		  Index2 = TmpTotalIndex2[n2];
		  LimitZ = LengthZ - p1;
		  for (; IndexZ < LimitZ; ++IndexZ, ++Index2)
		    {
		      TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
		      TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]);  	  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}
	      ++IndexY;
	    }
	  ++IndexX;
	}   
      delete[] TotalIndex;     
      return vDestination;
    }
}

// evaluate all interaction factors
// 
// waveVectorX = wave vector of Bloch function in X direction
// waveVectorY = wave vector of Bloch function in Y direction
// waveVectorZ = wave vector of Bloch function in Z direction  

void PeriodicQuantumDots3DHamiltonian::EvaluateInteractionFactors(double waveVectorX, double waveVectorY, double waveVectorZ)
{
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellX, this->NbrStateX, this->RealWaveFunctionOverlapX, this->ImaginaryWaveFunctionOverlapX))
    cout << "Error in evaluation of function overlap in X direction. Stop!" << endl;  
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellY, this->NbrStateY, this->RealWaveFunctionOverlapY, this->ImaginaryWaveFunctionOverlapY))
    cout << "Error in evaluation of function overlap in Y direction. Stop!" << endl;
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellZ, this->NbrStateZ, this->RealWaveFunctionOverlapZ, this->ImaginaryWaveFunctionOverlapZ))
    cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;

  double InvXFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Mux * this->XSize * this->XSize);
  double InvYFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muy * this->YSize * this->YSize);
  double InvZFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muz * this->ZSize * this->ZSize);

  double ShiftSquareK = BLOCH_FACTOR * (waveVectorX * waveVectorX / (2.0 * this->Mux) + waveVectorY * waveVectorY / (2.0 * this->Muy) + waveVectorZ * waveVectorZ / (2.0 * this->Muz));
  double ShiftKx = BLOCH_FACTOR * waveVectorX * 2.0 * M_PI/ (this->Mux * this->XSize); 
  double ShiftKy = BLOCH_FACTOR * waveVectorY * 2.0 * M_PI/ (this->Muy * this->YSize); 
  double ShiftKz = BLOCH_FACTOR * waveVectorZ * 2.0 * M_PI/ (this->Muz * this->ZSize); 
  this->KineticElements = new double[this->Space->GetHilbertSpaceDimension ()];

  double FactorX = 0.0, FactorY = 0.0;
  int TotalIndex = 0;
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      FactorX = double((i + this->LowerImpulsionX) * (i + this->LowerImpulsionX)) * InvXFactor;
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  FactorY = double((j + this->LowerImpulsionY) * (j + this->LowerImpulsionY)) * InvYFactor + FactorX;
	  for (int k = 0; k < this->NbrStateZ; ++k)
	    {	      
	      this->KineticElements[TotalIndex] = FactorY + double((k + this->LowerImpulsionZ) * (k + this->LowerImpulsionZ)) * InvZFactor;	  
	      this->KineticElements[TotalIndex] += (ShiftSquareK + ShiftKx * double(i + this->LowerImpulsionX) + ShiftKy * double(j + this->LowerImpulsionY) + ShiftKz * double(k + this->LowerImpulsionZ));
	      ++TotalIndex;
	    }
	}
    }

  int LengthX = (this->NbrStateX - 1) * 2 + 1; int LengthY = (this->NbrStateY - 1) * 2 + 1; int LengthZ = (this->NbrStateZ - 1) * 2 + 1;

  double*** TmpReal = new double** [LengthX];
  double*** TmpImaginary = new double** [LengthX];

  double TmpRe, TmpIm;
  double TmpRe2, TmpIm2;
  double* TmpRealWaveFunctionOverlapX;
  double* TmpImaginaryWaveFunctionOverlapX;
  double* TmpRealWaveFunctionOverlapY;
  double* TmpImaginaryWaveFunctionOverlapY;
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;

  for (int m = 0; m < LengthX; ++m)
    {
      TmpReal[m] = new double* [LengthY];
      TmpImaginary[m] = new double* [LengthY];
      TmpRealWaveFunctionOverlapX = this->RealWaveFunctionOverlapX[m];
      TmpImaginaryWaveFunctionOverlapX = this->ImaginaryWaveFunctionOverlapX[m];	      	  
      for (int n = 0; n < LengthY; ++n)
	{	  
	  TmpReal[m][n] = new double [this->NbrCellZ];
	  TmpImaginary[m][n] = new double [this->NbrCellZ];
	  TmpRealWaveFunctionOverlapY = this->RealWaveFunctionOverlapY[n];
	  TmpImaginaryWaveFunctionOverlapY = this->ImaginaryWaveFunctionOverlapY[n];	  
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];		  
	  for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		{
		  TmpRe2 = TmpRealWaveFunctionOverlapY[CellY];
		  TmpIm2 = TmpImaginaryWaveFunctionOverlapY[CellY];

		  for (int CellX = 0; CellX < this->NbrCellX;)
		    {		      
		      double tmpPotential = this->InteractionFactors[CellZ][CellY][CellX];
		      int tmpCellX = CellX; double TmpRe3 = 0.0, TmpIm3 = 0.0;
		      for (; this->InteractionFactors[CellZ][CellY][tmpCellX] == tmpPotential; ++tmpCellX)
			{
			  TmpRe3 += TmpRealWaveFunctionOverlapX[tmpCellX];
			  TmpIm3 += TmpImaginaryWaveFunctionOverlapX[tmpCellX];
			}
		      
		      CellX = tmpCellX;

		      TmpRe += tmpPotential * (TmpRe3 * TmpRe2 - TmpIm3 * TmpIm2);
		      TmpIm += tmpPotential * (TmpRe3 * TmpIm2 + TmpIm3 * TmpRe2);
		    }
		}
	      TmpRealPrecalculatedHamiltonian[CellZ] = TmpRe;  
	      TmpImaginaryPrecalculatedHamiltonian[CellZ] = TmpIm;  
	    }
	}
    }

  this->RealPrecalculatedHamiltonian = new double** [LengthX];
  this->ImaginaryPrecalculatedHamiltonian = new double** [LengthX];
  double* TmpRealWaveFunctionOverlapZ;
  double* TmpImaginaryWaveFunctionOverlapZ;
  for (int m = 0; m < LengthX; ++m)
    {
      this->RealPrecalculatedHamiltonian[m] = new double* [LengthY];      
      this->ImaginaryPrecalculatedHamiltonian[m] = new double* [LengthY]; 
      for (int n = 0; n < LengthY; ++n)
	{
	  this->RealPrecalculatedHamiltonian[m][n] = new double [LengthZ];      
	  this->ImaginaryPrecalculatedHamiltonian[m][n] = new double [LengthZ]; 
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];
	  for (int p = 0; p < LengthZ; ++p)
	    {
	      TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[p];
	      TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[p];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  TmpRe += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ] - TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ]);
		  TmpIm += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ] + TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ]);
		}
	      this->RealPrecalculatedHamiltonian[m][n][p] = TmpRe;
	      this->ImaginaryPrecalculatedHamiltonian[m][n][p] = TmpIm;
	    }
	}
    }
  delete[] TmpReal; delete[] TmpImaginary;
}

// evaluate the wave function overlap
//
// nbrStep = number of steps in the given direction
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool PeriodicQuantumDots3DHamiltonian::EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray)
{
  double Diff = 0.0;
  double Tmp = 0.0;
  double Tmp1 = 1.0 / double (nbrStep);
  int Length = (nbrState - 1) * 2 + 1;
  realArray = new double* [Length];
  imaginaryArray = new double* [Length];  
  int Origin = nbrState - 1;
  for (int delta = 0; delta < Length; ++delta)
    {
      realArray[delta] = new double [nbrStep];
      imaginaryArray[delta] = new double [nbrStep];
      if (delta != Origin)
	{
	  Diff = 2.0 * M_PI * double (delta - Origin);
	  Tmp = Diff / nbrStep;	
	  Diff = 1.0 / Diff;	
	  for (int i = 0; i < nbrStep; ++i)
	    {
	      realArray[delta][i] = Diff * (sin(Tmp * (i + 1)) - sin(Tmp * i));
	      imaginaryArray[delta][i] = Diff * (-cos(Tmp * (i + 1)) + cos(Tmp * i));
	    }
	}
      else
	for (int i = 0; i < nbrStep; ++i)
	  {
	    realArray[delta][i] = Tmp1;
	    imaginaryArray[delta][i] = 0.0;
	  }	
    }
  return true;
}

// determine the maximal value of partial diagonal array
//
// return = the wanted value

double PeriodicQuantumDots3DHamiltonian::MaxPartialDiagonalElement()
{
  double tmp = this->KineticElements[0];
  for (int i = 1; i < this->Space->GetHilbertSpaceDimension(); ++i)
    if (tmp < this->KineticElements[i])
      tmp = this->KineticElements[i];
  return tmp;
}
