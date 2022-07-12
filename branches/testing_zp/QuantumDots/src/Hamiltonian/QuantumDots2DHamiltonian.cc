////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 26/02/2003                        //
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
#include "Hamiltonian/QuantumDots2DHamiltonian.h"
#include "Vector/RealVector.h"
#include "MathTools/Complex.h"
#include "Potential/ThreeDPotential.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


#define HAMILTONIAN_FACTOR 37.60



// constructor from default datas
//
// space = Hilbert space associated to the system
// xSize = system dimension in the x direction (in Angstrom unit)
// ySize = system dimension in the y direction (in Angstrom unit)
// mux = effective mass in the x direction (in electron mass unit)
// muy = effective mass in the y direction (in electron mass unit)
// nbrCellX = number of cells in the x direction
// nbrCellY = number of cells in the y direction
// overlapingFactors = tridimensionnal array where overlaping factors are stored
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

QuantumDots2DHamiltonian::QuantumDots2DHamiltonian(Confined2DOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDPotential* PotentialInput, int memory)
{
  this->Space = space;
  this->XSize = xSize;
  this->YSize = ySize;
  this->Mux = mux;
  this->Muy = muy;
  this->NbrCellX = nbrCellX;
  this->NbrCellY = nbrCellY;
  this->TotalNbrCells = this->NbrCellX * this->NbrCellY;
  this->NbrStateX = this->Space->GetNbrStateX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->InteractionFactors = PotentialInput->Potential;
  this->EvaluateInteractionFactors(memory);
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

QuantumDots2DHamiltonian::QuantumDots2DHamiltonian(const QuantumDots2DHamiltonian& hamiltonian)
{
  this->Space = hamiltonian.Space;
  this->XSize = hamiltonian.XSize;
  this->YSize = hamiltonian.YSize;
  this->Mux = hamiltonian.Mux;
  this->Muy = hamiltonian.Muy; 
  this->NbrCellX = hamiltonian.NbrCellX;
  this->NbrCellY = hamiltonian.NbrCellY;
  this->TotalNbrCells = hamiltonian.TotalNbrCells; 
  this->NbrStateX = this->Space->GetNbrStateX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->DiagonalElements = hamiltonian.DiagonalElements;
  this->NbrPrecalculatedDimension = hamiltonian.NbrPrecalculatedDimension;
  this->InteractionFactors = hamiltonian.InteractionFactors;
  this->WaveFunctionOverlapX = hamiltonian.WaveFunctionOverlapX;
  this->WaveFunctionOverlapY = hamiltonian.WaveFunctionOverlapY;
  this->FullPrecalculatedHamiltonian = hamiltonian.FullPrecalculatedHamiltonian;
}

// destructor
//

QuantumDots2DHamiltonian::~QuantumDots2DHamiltonian()
{
  delete[] this->DiagonalElements;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* QuantumDots2DHamiltonian::Clone ()
{
  return new QuantumDots2DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void QuantumDots2DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void QuantumDots2DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->DiagonalElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex QuantumDots2DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
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

Complex QuantumDots2DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& QuantumDots2DHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
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

RealVector& QuantumDots2DHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    vDestination[i] = 0.0;
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& QuantumDots2DHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Space->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& QuantumDots2DHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Space->GetHilbertSpaceDimension();

  if (this->NbrPrecalculatedDimension == 0x010)
    {
      double Tmp = 0.0;
      int Index1 = firstComponent;
      int ReducedIndex1 = Index1 / this->NbrStateZ;
      int k1 = Index1 - ReducedIndex1 * this->NbrStateZ;
      int k2 = 0;
      int CellZ;
      int ReducedDim = this->NbrStateX * this->NbrStateY;
      double* TmpWaveFunctionOverlapZ;
      double* TmpPrecalculatedHamiltonian;
      while (Index1 < LastComponent)
	{
	  vDestination[Index1] += this->DiagonalElements[Index1] * vSource[Index1];
	  int ReducedIndex2 = 0;
	  int Index2 = 0;
	  for (; ReducedIndex2 < ReducedIndex1; ++ReducedIndex2)
	    {
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	      for (k2 = 0; k2 < k1; ++k2)
		{		  
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		    {
		      Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		    }
		  vDestination[Index1] += Tmp * vSource[Index2];
		  ++Index2;
		}
	      for (; k2 < this->NbrStateZ; ++k2)
		{		  
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k2][k1];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		    {
		      Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		    }
		  vDestination[Index1] += Tmp * vSource[Index2];
		  ++Index2;
		}
	    }
	  TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	  for (k2 = 0; k2 < k1; ++k2)
	    {
	      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
	      Tmp = 0.0;
	      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		}
	      Tmp += this->PartialZPrecalculatedHamiltonian[k1][k2];		  
	      vDestination[Index1] += Tmp * vSource[Index2];
	      ++Index2;
	    }	  
	  ++k2;
	  ++Index2;	  
	  for (; k2 < this->NbrStateZ; ++k2)
	    {
	      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k2][k1];
	      Tmp = 0.0;
	      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		}
	      Tmp += this->PartialZPrecalculatedHamiltonian[k2][k1];		  
	      vDestination[Index1] += Tmp * vSource[Index2];
	      ++Index2;
	    }	  
	  ++ReducedIndex2;
	  for (; ReducedIndex2 < ReducedDim; ++ReducedIndex2)
	    {
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex2][ReducedIndex1];
	      for (k2 = 0; k2 < k1; ++k2)
		{		  
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		    {
		      Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		    }
		  vDestination[Index1] += Tmp * vSource[Index2];
		  ++Index2;
		}
	      for (; k2 < this->NbrStateZ; ++k2)
		{		  
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k2][k1];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		    {
		      Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		    }
		  vDestination[Index1] += Tmp * vSource[Index2];
		  ++Index2;
		}
	    }
	  ++k1;
	  if (k1 == this->NbrStateZ)
	    {
	      k1 = 0;
	      ++ReducedIndex1;
	    }
	  ++Index1;
	}
      return vDestination;
    }
  if (this->NbrPrecalculatedDimension == 0x100)
    {
      double Tmp = 0.0;
      for (int Index1 = firstComponent; Index1 < LastComponent; ++Index1)
	{
	  Tmp = this->DiagonalElements[Index1] * vSource[Index1];
	  int Index2 = 0;
	  for (; Index2 < Index1; ++Index2)
	    Tmp += this->FullPrecalculatedHamiltonian[Index1][Index2] * vSource[Index2];
	  ++Index2;
	  for (; Index2 < Dim; ++Index2)
	    Tmp += this->FullPrecalculatedHamiltonian[Index2][Index1] * vSource[Index2];
	  vDestination[Index1] += Tmp;
	}
      return vDestination;
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& QuantumDots2DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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

ComplexVector& QuantumDots2DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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

ComplexVector& QuantumDots2DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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

ComplexVector& QuantumDots2DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							     int firstComponent, int nbrComponent)
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> QuantumDots2DHamiltonian::LeftInteractionOperators()  
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators 
//
// return value = list of right interaction operators

List<Matrix*> QuantumDots2DHamiltonian::RightInteractionOperators()  
{
  List<Matrix*> TmpList;
  return TmpList;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, QuantumDots2DHamiltonian& H)
{
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, QuantumDots2DHamiltonian& H)
{
  return Str;
}

// evaluate all interaction factors
//   
// memory = amount of memory available to store precalculated values

void QuantumDots2DHamiltonian::EvaluateInteractionFactors(int memory)
{
  int UsedMemory = 0;
  this->WaveFunctionOverlapX = this->EvaluateWaveFunctionOverlap (this->XSize, this->NbrCellX, this->NbrStateX, UsedMemory);
  this->WaveFunctionOverlapY = this->EvaluateWaveFunctionOverlap (this->YSize, this->NbrCellY, this->NbrStateY, UsedMemory);

  UsedMemory -= this->Space->GetHilbertSpaceDimension () * sizeof(double);
  double InvXFactor = HAMILTONIAN_FACTOR / (this->Mux * this->XSize * this->XSize);
  double InvYFactor = HAMILTONIAN_FACTOR / (this->Muy * this->YSize * this->YSize);
  int TotalIndex = 0;
  double FactorX = 0.0;
  double FactorY = 0.0;
  double Factor;
  double TmpElement = 0.0;
  this->DiagonalElements = new double[ this->Space->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      FactorX = ((double) ((i + 1) * (i + 1))) * InvXFactor;
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  TmpElement = ((double) ((j + 1) * (j + 1))) * InvYFactor + FactorX;		 
	  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
	    {
	      Factor = this->WaveFunctionOverlapY[j][j][CellY];
	      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		{
		  TmpElement += this->InteractionFactors[CellY][CellX] * Factor * this->WaveFunctionOverlapX[i][i][CellX];
		}	 
	    }     
	      	    
	  this->DiagonalElements[TotalIndex] = TmpElement;
	  ++TotalIndex;	   
	}
    }
  memory -= UsedMemory;
  if (memory < 0)
    {
      cout << "Calculation from scratch is done" << endl;
      this->NbrPrecalculatedDimension = 0;
      return;
    }
  
  else
    {
      cout << "Calculation with precalculated Hamiltonian is done" << endl;
      this->NbrPrecalculatedDimension = 0x010;
      this->Partial2DPrecalculatedHamiltonian = new double** [this->NbrStateY];
      double tmp = 0.0;
      for (int n1 = 0; n1 < this->NbrStateY; ++n1)
	{
	  this->Partial2DPrecalculatedHamiltonian[n1] = new double* [this->NbrStateY];
	  for (int n2 = 0; n2 < this->NbrStateY; ++n2)
	    {
	      this->Partial2DPrecalculatedHamiltonian[n1][n2] = new double [this->NbrCellZ];
	      for (int i = 0; i < this->NbrCellX; ++i)
		{
		  tmp = 0.0;
		  for (int j = 0; j < this->NbrCellY; ++j)
		    tmp += this->WaveFunctionOverlapY[n1][n2][j] * this->InteractionFactors[j][i];
		  this->Partial2DPrecalculatedHamiltonian[n1][n2][i] = tmp;
	        }
	    }
	}      
      return;
    }

}

// evaluate wave function overlaps on a cell in a given direction
//
// size = system length in the choosen direction
// nbrStep = number of subdivision in the choosen direction
// nbrState = number of state in the choosen direction
// memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
// return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)

double*** QuantumDots2DHamiltonian::EvaluateWaveFunctionOverlap(double size, int nbrStep, int nbrState, int& memory)
{
  memory += ((((nbrState * (nbrState + 1)) >> 1) * (nbrStep * sizeof(double) + sizeof(double*)))
	     + nbrState * sizeof(double**));
  double*** TmpArray = new double** [nbrState];
  double StepInc = 1.0 / ((double) nbrStep);
  double Tmp;
  double Diff;
  for (int i = 0; i < nbrState; ++i)
    {
      TmpArray[i] = new double* [i + 1];
      for (int j = 0; j < i; ++j)
	{
	  TmpArray[i][j] = new double [nbrStep];
	  for (int k = 0; k < nbrStep; ++k)
	    {
	      Diff = (double) (i - j);
	      Tmp = M_PI * Diff * StepInc;	      
	      TmpArray[i][j][k] = M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / Diff;
	      Diff = (double) (i + j + 2);
	      Tmp = M_PI * Diff * StepInc;	      
	      TmpArray[i][j][k] -= M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / Diff;
	    }
	}
      TmpArray[i][i] = new double [nbrStep];
      for (int k = 0; k < nbrStep; ++k)
	{
	  Tmp = M_PI * (double) (2 * i + 2) * StepInc;	      
	  TmpArray[i][i][k] = StepInc - M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / ((double) (2 * i + 2));
	}     
    }
  return TmpArray;
}

// determine the maximal value of partial diagonal array
//
// return = the wanted value

double QuantumDots2DHamiltonian::MaxPartialDiagonalElement()
{
  double tmp = this->DiagonalElements[0];
  for (int i = 1; i < this->Space->GetHilbertSpaceDimension(); ++i)
    if (tmp < this->DiagonalElements[i])
      tmp = this->DiagonalElements[i];
  return tmp;
}
