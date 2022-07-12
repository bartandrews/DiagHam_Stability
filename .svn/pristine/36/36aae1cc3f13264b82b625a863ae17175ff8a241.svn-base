////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2004 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//     class of hamiltonian associated quantum dots in a magnetic field       //
//                                                                            //
//                      last modification : 19/04/2004                        //
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
#include "Hamiltonian/PeriodicQuantumDots3DHamiltonianInMagneticField.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"

#include <iostream>
#include <math.h>

using std::ostream;
using std::cout;
using std::endl;


#define PERIODIC_HAMILTONIAN_FACTOR 150.4
#define PARAMAGNETIC_FACTOR         3.642e-4
#define DIAMAGNETIC_FACTOR          2.198e-10
#define BLOCH_FACTOR 7.644


// constructor from data
//
// space = Hilbert space
// xSize = the sample length in X direction
// ySize = the sample length in Y direction
// zSize = the sample length in Z direction
// bx = X magnetic field component
// by = Y magnetic field component
// bz = Z magnetic field component
// nbrCellX = number of steps in X direction
// nbrCellY = number of steps in Y direction
// nbrCellZ = number of steps in Z direction
// PotentialInput = pointer to a 3D potential with constant value in a cell
// waveVectorZ = wave vector of Bloch function in Z direction

PeriodicQuantumDots3DHamiltonianInMagneticField::PeriodicQuantumDots3DHamiltonianInMagneticField(PeriodicThreeDOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, double bx, double by, double bz, ThreeDConstantCellPotential* PotentialInput, double waveVectorZ)
{
  this->Space = space;
  this->XSize = xSize;
  this->YSize = ySize;
  this->ZSize = zSize;
  this->Mux = mux;
  this->Muy = muy;
  this->Muz = muz;
  this->NbrCellX = PotentialInput->GetNbrCellX();
  this->NbrCellY = PotentialInput->GetNbrCellY();
  this->NbrCellZ = PotentialInput->GetNbrCellZ();
  //this->Bx = bx;
  //this->By = by;
  this->Bz = bz;
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
  cout << "Evaluating confinement potential ..." << endl;
  this->EvaluateConfinementPotentialFactors(waveVectorZ);
  cout << "End of confinement potential evaluation." << endl;
  cout << "Evaluating magnetic field factors ..." << endl;
  this->EvaluateMagneticFieldFactors();
  cout << "End of magnetic field evaluation." << endl;
}


// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

PeriodicQuantumDots3DHamiltonianInMagneticField::PeriodicQuantumDots3DHamiltonianInMagneticField(const PeriodicQuantumDots3DHamiltonianInMagneticField& hamiltonian)
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
  //this->Bx = hamiltonian.Bx;
  //this->By = hamiltonian.By;
  this->Bz = hamiltonian.Bz;
  this->NbrStateX = this->Space->GetNbrStateX();
  this->LowerImpulsionX = this->Space->GetLowerImpulsionX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->LowerImpulsionY = this->Space->GetLowerImpulsionY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();
  this->KineticElements = hamiltonian.KineticElements;
  this->NbrPrecalculatedDimension = hamiltonian.NbrPrecalculatedDimension;
  this->InteractionFactors = hamiltonian.InteractionFactors;
  this->RealPrecalculatedHamiltonian =  hamiltonian.RealPrecalculatedHamiltonian;
  this->ImaginaryPrecalculatedHamiltonian =  hamiltonian.ImaginaryPrecalculatedHamiltonian;
}

// destructor
//

PeriodicQuantumDots3DHamiltonianInMagneticField::~ PeriodicQuantumDots3DHamiltonianInMagneticField()
{
  delete[] this->KineticElements;
  delete[] this->InteractionFactors;
  delete[] this->RealPrecalculatedHamiltonian;
  delete[] this->ImaginaryPrecalculatedHamiltonian;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* PeriodicQuantumDots3DHamiltonianInMagneticField::Clone ()
{
  return new PeriodicQuantumDots3DHamiltonianInMagneticField(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PeriodicQuantumDots3DHamiltonianInMagneticField::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void PeriodicQuantumDots3DHamiltonianInMagneticField::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->KineticElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicQuantumDots3DHamiltonianInMagneticField::MatrixElement (RealVector& V1, RealVector& V2)
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

Complex PeriodicQuantumDots3DHamiltonianInMagneticField::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicQuantumDots3DHamiltonianInMagneticField::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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

ComplexVector& PeriodicQuantumDots3DHamiltonianInMagneticField::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
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

ComplexVector& PeriodicQuantumDots3DHamiltonianInMagneticField::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  int OriginX = this->NbrStateX - 1; int OriginY = this->NbrStateY - 1; int OriginZ = this->NbrStateZ - 1;
  int LengthX = (this->NbrStateX - 1) * 2 + 1; int LengthY = (this->NbrStateY - 1) * 2 + 1; int LengthZ = (this->NbrStateZ - 1) * 2 + 1;

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

  int* TmpTotalIndex1; int* TmpTotalIndex2;
  int Index1, Index2;
  double* TmpRealParamagnetic; double* TmpImaginaryParamagnetic;
  double TmpReal = 0.0, TmpImaginary = 0.0;
  
  for (int m = 0; m < this->NbrStateX; ++m)
    {
      TmpTotalIndex1 = TotalIndex[m];      
      TmpTotalIndex2 = TotalIndex[m]; 
      TmpRealParamagnetic = this->RealParamagneticTermPxY[m];
      TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPxY[m];
      for (n1 = 0; n1 < this->NbrStateY; ++n1)
	{
	  IndexY = -n1 + OriginY;
	  for (n2 = 0; n2 < this->NbrStateY; ++n2)
	    {
	      Index1 = TmpTotalIndex1[n1];
	      Index2 = TmpTotalIndex2[n2];	
	      TmpReal = TmpRealParamagnetic[IndexY];
	      TmpImaginary = TmpImaginaryParamagnetic[IndexY];
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  vDestination.Re(Index1) += (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
		  vDestination.Im(Index1) += (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));		  
		  ++Index1; ++Index2;
		}
	      ++IndexY;
	    }
	}
    }
  
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {
      TmpTotalIndex1 = TotalIndex[m1]; 
      IndexX = -m1 + OriginX;
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];
	  TmpRealParamagnetic = this->RealParamagneticTermPyX[IndexX];
	  TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPyX[IndexX];
	  for (int n = 0; n < this->NbrStateY; ++n)
	    {
	      Index1 = TmpTotalIndex1[n];
	      Index2 = TmpTotalIndex2[n];
	      TmpReal = TmpRealParamagnetic[n];
	      TmpImaginary = TmpImaginaryParamagnetic[n];
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  vDestination.Re(Index1) -= (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
		  vDestination.Im(Index1) -= (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));		  
		  ++Index1; ++Index2;
		}	      
	    }
	  ++IndexX;
	}
    }
  
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;
  double TmpRe = 0.0; double TmpIm = 0.0;
  int LimitZ = 0;
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

ComplexVector& PeriodicQuantumDots3DHamiltonianInMagneticField::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{ 
  if ((firstComponent == 0) && (nbrComponent == this->Space->GetHilbertSpaceDimension()))
    return this->LowLevelAddMultiply(vSource, vDestination);
  else
    {      
      int lastComponent = firstComponent + nbrComponent;
      int OriginX = this->NbrStateX - 1; int OriginY = this->NbrStateY - 1; int OriginZ = this->NbrStateZ - 1;
      int m, m1, m2, n, n1, n2, p1;
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

      double* TmpRealParamagnetic; double* TmpImaginaryParamagnetic;
      double TmpReal = 0.0, TmpImaginary = 0.0;
      
      // ParamagneticTermPxY //
      // p1 : p1Begin -> NbrStateZ
      m = m1Begin; n1 = n1Begin;
      TmpTotalIndex1 = TotalIndex[m];
      TmpTotalIndex2 = TotalIndex[m]; 
      TmpRealParamagnetic = this->RealParamagneticTermPxY[m];
      TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPxY[m];
      IndexY = -n1 + OriginY;
      for (n2 = 0; n2 < this->NbrStateY; ++n2)
	{
	  Index1 = firstComponent;
	  Index2 = TmpTotalIndex2[n2] + p1Begin;	
	  TmpReal = TmpRealParamagnetic[IndexY];
	  TmpImaginary = TmpImaginaryParamagnetic[IndexY];
	  for (int p = p1Begin; p < this->NbrStateZ; ++p)
	    {
	      vDestination.Re(Index1) += (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
	      vDestination.Im(Index1) += (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));
	      ++Index1; ++Index2;
	    }
	  ++IndexY;
	}            
      // n1 : n1Begin -> NbrStateY
      m = m1Begin;
      TmpTotalIndex1 = TotalIndex[m];
      TmpTotalIndex2 = TotalIndex[m]; 
      TmpRealParamagnetic = this->RealParamagneticTermPxY[m];
      TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPxY[m];
      for (n1 = n1Begin + 1; n1 < this->NbrStateY; ++n1)
	{      
	  IndexY = -n1 + OriginY;
	  for (n2 = 0; n2 < this->NbrStateY; ++n2)
	    {
	      Index1 = TmpTotalIndex1[n1];
	      Index2 = TmpTotalIndex2[n2];
	      TmpReal = TmpRealParamagnetic[IndexY];
	      TmpImaginary = TmpImaginaryParamagnetic[IndexY];	
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  vDestination.Re(Index1) += (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
		  vDestination.Im(Index1) += (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));		  
		  ++Index1; ++Index2;
		}
	      ++IndexY;
	    }      
	}      
      // m1 : m1Begin -> m1Limit - 1
      for (m = m1Begin + 1; m < m1Limit; ++m)
	{
	  TmpTotalIndex1 = TotalIndex[m];
	  TmpTotalIndex2 = TotalIndex[m]; 
	  TmpRealParamagnetic = this->RealParamagneticTermPxY[m];
	  TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPxY[m];
	  for (n1 = 0; n1 < this->NbrStateY; ++n1)
	    {      
	      IndexY = -n1 + OriginY;
	      for (n2 = 0; n2 < this->NbrStateY; ++n2)
		{
		  Index1 = TmpTotalIndex1[n1];
		  Index2 = TmpTotalIndex2[n2];
		  TmpReal = TmpRealParamagnetic[IndexY];
		  TmpImaginary = TmpImaginaryParamagnetic[IndexY];	
		  for (int p = 0; p < this->NbrStateZ; ++p)
		    {
		      vDestination.Re(Index1) += (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
		      vDestination.Im(Index1) += (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));	  
		      ++Index1; ++Index2;
		    }
		  ++IndexY;
		}      
	    }
	}
      // m1 = m1Limit, n1 : 0 -> n1Limit - 1
      m = m1Limit;
      TmpTotalIndex1 = TotalIndex[m];
      TmpTotalIndex2 = TotalIndex[m]; 
      TmpRealParamagnetic = this->RealParamagneticTermPxY[m];
      TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPxY[m];
      for (n1 = 0; n1 < n1Limit; ++n1)
	{      
	  IndexY = -n1 + OriginY;
	  for (n2 = 0; n2 < this->NbrStateY; ++n2)
	    {
	      Index1 = TmpTotalIndex1[n1];
	      Index2 = TmpTotalIndex2[n2];
	      TmpReal = TmpRealParamagnetic[IndexY];
	      TmpImaginary = TmpImaginaryParamagnetic[IndexY];	
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  vDestination.Re(Index1) += (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
		  vDestination.Im(Index1) += (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));		  
		  ++Index1; ++Index2;
		}
	      ++IndexY;
	    }      
	}     
      // m1 = m1Limit, n1 = n1Limit, p1 : 0 -> p1Limit
      m = m1Limit; n1 = n1Limit;
      TmpTotalIndex1 = TotalIndex[m];
      TmpTotalIndex2 = TotalIndex[m]; 
      TmpRealParamagnetic = this->RealParamagneticTermPxY[m];
      TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPxY[m];
      IndexY = -n1 + OriginY;
      for (n2 = 0; n2 < this->NbrStateY; ++n2)
	{
	  Index1 = TmpTotalIndex1[n1];
	  Index2 = TmpTotalIndex2[n2];
	  TmpReal = TmpRealParamagnetic[IndexY];
	  TmpImaginary = TmpImaginaryParamagnetic[IndexY];	
	  for (int p = 0; p <= p1Limit; ++p)
	    {
	      vDestination.Re(Index1) += (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
	      vDestination.Im(Index1) += (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));	
	      ++Index1; ++Index2;
	    }
	  ++IndexY;
	}    
      
      //  RealParamagneticTermPyX  //
      // p1 : p1Begin -> NbrStateZ
      m1 = m1Begin; n = n1Begin;
      TmpTotalIndex1 = TotalIndex[m];
      IndexX = -m1 + OriginX;
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];
	  TmpRealParamagnetic = this->RealParamagneticTermPyX[IndexX];
	  TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPyX[IndexX];	  
	  Index1 = firstComponent;
	  Index2 = TmpTotalIndex2[n] + p1Begin;
	  TmpReal = TmpRealParamagnetic[n];
	  TmpImaginary = TmpImaginaryParamagnetic[n];
	  for (int p = p1Begin; p < this->NbrStateZ; ++p)
	    {
	      vDestination.Re(Index1) -= (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
	      vDestination.Im(Index1) -= (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));		  
	      ++Index1; ++Index2;
	    }	      	    
	  ++IndexX;
	}          
      // n1 : n1Begin +1  -> NbrStateY
      m1 = m1Begin;
      TmpTotalIndex1 = TotalIndex[m1];
      IndexX = -m1 + OriginX;
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];
	  TmpRealParamagnetic = this->RealParamagneticTermPyX[IndexX];
	  TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPyX[IndexX];
	  for (int n = n1Begin + 1; n < this->NbrStateY; ++n)
	    {	  
	      Index1 = TmpTotalIndex1[n];
	      Index2 = TmpTotalIndex2[n];
	      TmpReal = TmpRealParamagnetic[n];
	      TmpImaginary = TmpImaginaryParamagnetic[n];
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  vDestination.Re(Index1) -= (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
		  vDestination.Im(Index1) -= (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));		  
		  ++Index1; ++Index2;
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
	      TmpRealParamagnetic = this->RealParamagneticTermPyX[IndexX];
	      TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPyX[IndexX];
	      for (int n = 0; n < this->NbrStateY; ++n)
		{	  
		  Index1 = TmpTotalIndex1[n];
		  Index2 = TmpTotalIndex2[n];
		  TmpReal = TmpRealParamagnetic[n];
		  TmpImaginary = TmpImaginaryParamagnetic[n];
		  for (int p = 0; p < this->NbrStateZ; ++p)
		    {
		      vDestination.Re(Index1) -= (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
		      vDestination.Im(Index1) -= (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));		  
		      ++Index1; ++Index2;
		    }	      	    
		}  
	      ++IndexX;
	    }
	}
      // m1 = m1Limit, n1 : 0 -> n1Limit - 1
      m1 = m1Limit;
      IndexX = -m1 + OriginX;
      TmpTotalIndex1 = TotalIndex[m1];
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];
	  TmpRealParamagnetic = this->RealParamagneticTermPyX[IndexX];
	  TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPyX[IndexX];
	  for (int n = 0; n < n1Limit; ++n)
	    {	  
	      Index1 = TmpTotalIndex1[n];
	      Index2 = TmpTotalIndex2[n];
	      TmpReal = TmpRealParamagnetic[n];
	      TmpImaginary = TmpImaginaryParamagnetic[n];
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  vDestination.Re(Index1) -= (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
		  vDestination.Im(Index1) -= (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));		  
		  ++Index1; ++Index2;
		}	      	    
	    }  
	  ++IndexX;
	}
      // m1 = m1Limit, n1 = n1Limit, p1 : 0 -> p1Limit
      m1 = m1Limit; n = n1Limit,
      IndexX = -m1 + OriginX;
      TmpTotalIndex1 = TotalIndex[m1];      
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];
	  TmpRealParamagnetic = this->RealParamagneticTermPyX[IndexX];
	  TmpImaginaryParamagnetic = this->ImaginaryParamagneticTermPyX[IndexX];	  
	  Index1 = TmpTotalIndex1[n];
	  Index2 = TmpTotalIndex2[n];
	  TmpReal = TmpRealParamagnetic[n];
	  TmpImaginary = TmpImaginaryParamagnetic[n];
	  for (int p = 0; p <= p1Limit; ++p)
	    {
	      vDestination.Re(Index1) -= (TmpReal * vSource.Re(Index2) - TmpImaginary * vSource.Im(Index2));
	      vDestination.Im(Index1) -= (TmpReal * vSource.Im(Index2) + TmpImaginary * vSource.Re(Index2));		  
	      ++Index1; ++Index2;
	    }	      	    
	  ++IndexX;
	}  
      

      //   PRECALCULATED HAMILTONIAN   //
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
      // n1 : n1Begin -> NbrStateY
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

// evaluate confinement potential factors
//   
// waveVectorZ = wave vector of Bloch function in Z direction 

void PeriodicQuantumDots3DHamiltonianInMagneticField::EvaluateConfinementPotentialFactors(double waveVectorZ)
{
  double** RealWaveFunctionOverlapX; double** RealWaveFunctionOverlapY; double** RealWaveFunctionOverlapZ; 
  double** ImaginaryWaveFunctionOverlapX; double** ImaginaryWaveFunctionOverlapY; double** ImaginaryWaveFunctionOverlapZ; 

  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellX, this->NbrStateX, RealWaveFunctionOverlapX, ImaginaryWaveFunctionOverlapX))
    cout << "Error in evaluation of function overlap in X direction. Stop!" << endl;  
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellY, this->NbrStateY, RealWaveFunctionOverlapY, ImaginaryWaveFunctionOverlapY))
    cout << "Error in evaluation of function overlap in Y direction. Stop!" << endl;
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellZ, this->NbrStateZ, RealWaveFunctionOverlapZ, ImaginaryWaveFunctionOverlapZ))
    cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;

  double InvXFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Mux * this->XSize * this->XSize);
  double InvYFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muy * this->YSize * this->YSize);
  double InvZFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muz * this->ZSize * this->ZSize);
  double ShiftSquareKz = BLOCH_FACTOR * waveVectorZ * waveVectorZ / (2.0 * this->Muz);
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
	      this->KineticElements[TotalIndex] += (ShiftSquareKz + ShiftKz * double(k + this->LowerImpulsionZ));
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
      TmpRealWaveFunctionOverlapX = RealWaveFunctionOverlapX[m];
      TmpImaginaryWaveFunctionOverlapX = ImaginaryWaveFunctionOverlapX[m];	      	  
      for (int n = 0; n < LengthY; ++n)
	{	  
	  TmpReal[m][n] = new double [this->NbrCellZ];
	  TmpImaginary[m][n] = new double [this->NbrCellZ];
	  TmpRealWaveFunctionOverlapY = RealWaveFunctionOverlapY[n];
	  TmpImaginaryWaveFunctionOverlapY = ImaginaryWaveFunctionOverlapY[n];	  
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];		  
	  for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		{
		  TmpRe2 = TmpRealWaveFunctionOverlapY[CellY];
		  TmpIm2 = TmpImaginaryWaveFunctionOverlapY[CellY];
		  for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		    {		      
		      TmpRe += this->InteractionFactors[CellZ][CellY][CellX] * (TmpRealWaveFunctionOverlapX[CellX] * TmpRe2 - TmpImaginaryWaveFunctionOverlapX[CellX] * TmpIm2);
		      TmpIm += this->InteractionFactors[CellZ][CellY][CellX] * (TmpRealWaveFunctionOverlapX[CellX] * TmpIm2 + TmpImaginaryWaveFunctionOverlapX[CellX] * TmpRe2);		      
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
	      TmpRealWaveFunctionOverlapZ = RealWaveFunctionOverlapZ[p];
	      TmpImaginaryWaveFunctionOverlapZ = ImaginaryWaveFunctionOverlapZ[p];
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
  delete[] RealWaveFunctionOverlapX; delete[] RealWaveFunctionOverlapY; delete[] RealWaveFunctionOverlapZ;
  delete[] ImaginaryWaveFunctionOverlapX; delete[] ImaginaryWaveFunctionOverlapY; delete[] ImaginaryWaveFunctionOverlapZ;
}

// evaluate magnetic field factors
//

void PeriodicQuantumDots3DHamiltonianInMagneticField::EvaluateMagneticFieldFactors()
{
  double* RealX; double* RealY; double* RealZ;
  double* ImaginaryX; double* ImaginaryY; double* ImaginaryZ;
  double* RealSquaredX; double* RealSquaredY; double* RealSquaredZ;
  double* ImaginarySquaredX; double* ImaginarySquaredY; double* ImaginarySquaredZ;

  if (!this->EvaluateMeanPositionOperator(this->XSize, this->NbrStateX, RealX, ImaginaryX, RealSquaredX, ImaginarySquaredX))
    cout << "Error in evaluation of mean operator in X direction" << endl;
  if (!this->EvaluateMeanPositionOperator(this->YSize, this->NbrStateY, RealY, ImaginaryY, RealSquaredY, ImaginarySquaredY))
    cout << "Error in evaluation of mean operator in Y direction" << endl;  
  if (!this->EvaluateMeanPositionOperator(this->ZSize, this->NbrStateZ, RealZ, ImaginaryZ, RealSquaredZ, ImaginarySquaredZ))
    cout << "Error in evaluation of mean operator in Z direction" << endl;

  int LengthX = (this->NbrStateX - 1) * 2 + 1; int LengthY = (this->NbrStateY - 1) * 2 + 1; int LengthZ = (this->NbrStateZ - 1) * 2 + 1;
  int OriginX = this->NbrStateX - 1; int OriginY = this->NbrStateY - 1; int OriginZ = this->NbrStateZ - 1;

  // evaluation of paramagnetic terms
  this->RealParamagneticTermPxY = new double* [this->NbrStateX];
  this->ImaginaryParamagneticTermPxY = new double* [this->NbrStateX];
  for (int m = 0; m < this->NbrStateX; ++m)
    {
      this->RealParamagneticTermPxY[m] = new double [LengthY];
      this->ImaginaryParamagneticTermPxY[m] = new double [LengthY];
      for (int delta = 0; delta < LengthY; ++delta)
	{
	  this->RealParamagneticTermPxY[m][delta] = PARAMAGNETIC_FACTOR * this->Bz * RealY[delta] * double(m + this->LowerImpulsionX) / (this->Mux * this->XSize);	
	  this->ImaginaryParamagneticTermPxY[m][delta] = PARAMAGNETIC_FACTOR * this->Bz * ImaginaryY[delta] * double(m + this->LowerImpulsionX) / (this->Mux * this->XSize);
	}
    }

  this->RealParamagneticTermPyX = new double* [LengthX];
  this->ImaginaryParamagneticTermPyX = new double* [LengthX];  
  for (int delta = 0; delta < LengthX; ++delta)
    {
      this->RealParamagneticTermPyX[delta] = new double [this->NbrStateY];
      this->ImaginaryParamagneticTermPyX[delta] = new double [this->NbrStateY];
      for (int n = 0; n < this->NbrStateY; ++n)
	{
	  this->RealParamagneticTermPyX[delta][n] = PARAMAGNETIC_FACTOR * this->Bz * RealX[delta] * double(n + this->LowerImpulsionY) / (this->Muy * this->YSize);
	  this->ImaginaryParamagneticTermPyX[delta][n] = PARAMAGNETIC_FACTOR * this->Bz * ImaginaryX[delta] * double(n + this->LowerImpulsionY) / (this->Muy * this->YSize);
	}
    }

  // evaluation of diamagnetic factor    
  for (int m = 0; m < LengthX; ++m)
    {
      this->RealPrecalculatedHamiltonian[m][OriginY][OriginZ] += (DIAMAGNETIC_FACTOR * this->Bz * this->Bz * RealSquaredX[m]) / this->Muy;
      this->ImaginaryPrecalculatedHamiltonian[m][OriginY][OriginZ] += (DIAMAGNETIC_FACTOR * this->Bz * this->Bz * ImaginarySquaredX[m]) / this->Muy;
      //this->RealPrecalculatedHamiltonian[m][OriginY][OriginZ] += (DIAMAGNETIC_FACTOR * this->By * this->By * RealSquaredX[m]) / this->Muz;
      //this->ImaginaryPrecalculatedHamiltonian[m][OriginY][OriginZ] += (DIAMAGNETIC_FACTOR * this->By * this->By * ImaginarySquaredX[m]) / this->Muz;      
    }  
   
  for (int n = 0; n < LengthY; ++n)
    {
      this->RealPrecalculatedHamiltonian[OriginX][n][OriginZ] += (DIAMAGNETIC_FACTOR * this->Bz * this->Bz * RealSquaredY[n]) / this->Mux;
      this->ImaginaryPrecalculatedHamiltonian[OriginX][n][OriginZ] += (DIAMAGNETIC_FACTOR * this->Bz * this->Bz * ImaginarySquaredY[n]) / this->Mux;    
      //this->RealPrecalculatedHamiltonian[OriginX][n][OriginZ] += (DIAMAGNETIC_FACTOR * this->Bx * this->Bx * RealSquaredY[n]) / this->Muz;
      //this->ImaginaryPrecalculatedHamiltonian[OriginX][n][OriginZ] += (DIAMAGNETIC_FACTOR * this->Bx * this->Bx * ImaginarySquaredY[n]) / this->Muz;      
    }  
  

  /*
  for (int p = 0; p < LengthZ; ++p)
    {
      this->RealPrecalculatedHamiltonian[OriginX][OriginY][p] += (DIAMAGNETIC_FACTOR * this->Bx * this->Bx * RealSquaredZ[p]) / this->Muy;
      this->ImaginaryPrecalculatedHamiltonian[OriginX][OriginY][p] += (DIAMAGNETIC_FACTOR * this->Bx * this->Bx * ImaginarySquaredZ[p]) / this->Muy;
      this->RealPrecalculatedHamiltonian[OriginX][OriginY][p] += (DIAMAGNETIC_FACTOR * this->By * this->By * RealSquaredZ[p]) / this->Mux;
      this->ImaginaryPrecalculatedHamiltonian[OriginX][OriginY][p] += (DIAMAGNETIC_FACTOR * this->By * this->By * ImaginarySquaredZ[p]) / this->Mux; 
    }
  */
}

// evaluate the wave function overlap
//
// nbrStep = number of steps in the given direction
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool PeriodicQuantumDots3DHamiltonianInMagneticField::EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray)
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

// evaluate the mean position operator in a given direction
//
// size = size of sample in the given direction
// nbrState = number of states wanted
// real = reference to the array of mean position's real component
// imaginary = reference to the array of mean position's imaginary component
// realSquared = reference to the array of X²'s real component
// imaginarySquared = reference to the array of X²'s imaginary component
// return = true if successful, otherwise false

bool PeriodicQuantumDots3DHamiltonianInMagneticField::EvaluateMeanPositionOperator(double size, int nbrState, double* &real, double* &imaginary, double* &realSquared, double* &imaginarySquared)
{
  int Length = (nbrState - 1) * 2 + 1;
  int Origin = nbrState - 1;
  
  real = new double [Length];
  imaginary = new double[Length];
  realSquared = new double [Length];
  imaginarySquared = new double [Length];

  double Tmp = 0.0; double size2 = size * size;

  for (int i = 0; i < Origin; ++i)
    {
      Tmp = 1.0 / (2.0 * M_PI * double(i - Origin));
      real[i] = 0.0;
      imaginary[i] = -(Tmp * size);
      realSquared[i] = Tmp * Tmp * 2.0 * size2;
      imaginarySquared[i] = -(Tmp * size2);
    }

  real[Origin] = 0.5 * size;
  imaginary[Origin] = 0.0;
  realSquared[Origin] = size2 / 3.0;
  imaginarySquared[Origin] = 0.0;
  
  for (int i = nbrState; i < Length; ++i)
    {
      Tmp = 1.0 / (2.0 * M_PI * double(i - Origin));
      real[i] = 0.0;
      imaginary[i] = -(Tmp * size);
      realSquared[i] = Tmp * Tmp * 2.0 * size2;
      imaginarySquared[i] = -(Tmp * size2);
    }

  return true;
}

// determine the maximal value of partial diagonal array
//
// return = the wanted value

double PeriodicQuantumDots3DHamiltonianInMagneticField::MaxPartialDiagonalElement()
{
  double tmp = this->KineticElements[0];
  for (int i = 1; i < this->Space->GetHilbertSpaceDimension(); ++i)
    if (tmp < this->KineticElements[i])
      tmp = this->KineticElements[i];
  return tmp;
}
