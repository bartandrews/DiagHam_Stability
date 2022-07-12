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
#include "Hamiltonian/QuantumDots3DHamiltonian.h"
#include "Vector/RealVector.h"
#include "Complex.h"

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
// zSize = system dimension in the z direction (in Angstrom unit)
// preConstantRegionSize = region size in the z direction where potential is constant in every direction (region before gradiant zone)
// postConstantRegionSize = region size in the z direction where potential is constant in every direction (region after gradiant zone)
// postConstantRegionPotential = value of the potential in the region after the gradiant zone
// mux = effective mass in the x direction (in electron mass unit)
// muy = effective mass in the y direction (in electron mass unit)
// muz = effective mass in the z direction (in electron mass unit)
// nbrCellX = number of cells in the x direction
// nbrCellY = number of cells in the y direction
// nbrCellZ = number of cells in the z direction
// overlapingFactors = tridimensionnal array where overlaping factors are stored
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

QuantumDots3DHamiltonian::QuantumDots3DHamiltonian(Confined3DOneParticle* space, double xSize, double ySize, double zSize, double preConstantRegionSize,
						   double postConstantRegionSize, double postConstantRegionPotential, double mux, double muy, double muz, 
						   int nbrCellX, int nbrCellY, int nbrCellZ, double*** overlapingFactors,
						   int memory)
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
  this->TotalNbrCells = this->NbrCellX * this->NbrCellY * this->NbrCellZ;
  this->NbrStateX = this->Space->GetNbrStateX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->PreConstantRegionSize = preConstantRegionSize;
  this->PostConstantRegionSize = postConstantRegionSize;
  this->PostConstantRegionPotential = postConstantRegionPotential;
  this->InteractionFactors = overlapingFactors;  
  this->EvaluateInteractionFactors(memory);
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy 

QuantumDots3DHamiltonian::QuantumDots3DHamiltonian(const QuantumDots3DHamiltonian& hamiltonian)
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
  this->TotalNbrCells = hamiltonian.TotalNbrCells; 
  this->NbrStateX = this->Space->GetNbrStateX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->DiagonalElements = hamiltonian.DiagonalElements;
  this->NbrPrecalculatedDimension = hamiltonian.NbrPrecalculatedDimension;
  this->InteractionFactors = hamiltonian.InteractionFactors;
  this->WaveFunctionOverlapX = hamiltonian.WaveFunctionOverlapX;
  this->WaveFunctionOverlapY = hamiltonian.WaveFunctionOverlapY;
  this->WaveFunctionOverlapZ = hamiltonian.WaveFunctionOverlapZ;
  this->FullPrecalculatedHamiltonian = hamiltonian.FullPrecalculatedHamiltonian;
  this->PreConstantRegionSize = hamiltonian.PreConstantRegionSize;
  this->PostConstantRegionSize = hamiltonian.PostConstantRegionSize;
  this->PostConstantRegionPotential = hamiltonian.PostConstantRegionPotential;
}

// destructor
//

QuantumDots3DHamiltonian::~QuantumDots3DHamiltonian()
{
  delete[] this->DiagonalElements;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* QuantumDots3DHamiltonian::Clone ()
{
  return new QuantumDots3DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void QuantumDots3DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void QuantumDots3DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->DiagonalElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex QuantumDots3DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
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

Complex QuantumDots3DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& QuantumDots3DHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
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

RealVector& QuantumDots3DHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
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

RealVector& QuantumDots3DHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
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

RealVector& QuantumDots3DHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Space->GetHilbertSpaceDimension();

  // no acceleration
  if (this->NbrPrecalculatedDimension == 0)
    {
      double Tmp;      
      double Factor;
      int i1 = 0;
      int j1 = 0;
      int k1 = 0;
      for (int Index1 = firstComponent; Index1 < LastComponent; ++Index1)
	{
	  vDestination[Index1] += this->DiagonalElements[Index1] * vSource[Index1];
	  int Index2 = 0;
	  int k2 = 0;
	  for (; k2 < k1; ++k2)
	    {
	      int j2 = 0;
	      for (; j2 <= j1; ++j2)
		{
		  int i2 = 0;
		  for (; i2 <= i1; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapX[i1][i2][CellX] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
				{
				  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k2][CellZ];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapX[i2][i1][CellX] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
				{
				  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k2][CellZ];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		}
	      for (; j2 < this->NbrStateY; ++j2)
		{
		  int i2 = 0;
		  for (; i2 <= i1; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapX[i1][i2][CellX] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
				{
				  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k2][CellZ];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapX[i2][i1][CellX] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
				{
				  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k2][CellZ];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		}
	    }
	  {
	    int j2 = 0;
	    for (; j2 < j1; ++j2)
	      {
		int i2 = 0;
		for (; i2 <= i1; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapX[i1][i2][CellX] * this->WaveFunctionOverlapY[j1][j2][CellY];
			    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			      {
				Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k1][CellZ];
				}	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
		for (; i2 < this->NbrStateX; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapX[i2][i1][CellX] * this->WaveFunctionOverlapY[j1][j2][CellY];
			    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			      {
				Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k1][CellZ];
			      }	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
	      }
	    {
	      int i2 = 0;
	      for (; i2 < i1; ++i2)
		{
		  Tmp = 0.0;
		  for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		    {
		      for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			{
			  Factor = this->WaveFunctionOverlapX[i1][i2][CellX] * this->WaveFunctionOverlapY[j1][j1][CellY];
			  for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			    {
			      Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k1][CellZ];
			    }	 
			}     
		    }
		  vDestination[Index1] += Tmp * vSource[Index2];
		  ++Index2;
		}
	      ++i2;
	      ++Index2;
	      for (; i2 < this->NbrStateX; ++i2)
		{
		    Tmp = 0.0;
		    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapX[i2][i1][CellX] * this->WaveFunctionOverlapY[j1][j1][CellY];
			    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			      {
				Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k1][CellZ];
			      }	 
			  }     
		      }
		  vDestination[Index1] += Tmp * vSource[Index2];
		  ++Index2;
		}
	      
	    }
	    ++j2;
	    for (; j2 < this->NbrStateY; ++j2)
	      {
		int i2 = 0;
		for (; i2 <= i1; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapX[i1][i2][CellX] * this->WaveFunctionOverlapY[j2][j1][CellY];
			    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			      {
				Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k1][CellZ];
			      }	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
		for (; i2 < this->NbrStateX; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapX[i2][i1][CellX] * this->WaveFunctionOverlapY[j2][j1][CellY];
			    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			      {
				Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k1][k1][CellZ];
			      }	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
	      }
	  }
	  ++k2;
	  for (; k2 < this->NbrStateZ; ++k2)
	    {
	      int j2 = 0;
	      for (; j2 <= j1; ++j2)
		{
		  int i2 = 0;
		  for (; i2 <= i1; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapX[i1][i2][CellX] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
				{
				  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k2][k1][CellZ];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapX[i2][i1][CellX] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
				{
				  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k2][k1][CellZ];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		}
	      for (; j2 < this->NbrStateY; ++j2)
		{
		  int i2 = 0;
		  for (; i2 <= i1; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapX[i1][i2][CellX] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
				{
				  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k2][k1][CellZ];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapX[i2][i1][CellX] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
				{
				  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k2][k1][CellZ];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		}
	      ++i1;
	      if (i1 == this->NbrStateX)
		{
		  i1 = 0;
		  ++j1;
		  if (j1 == this->NbrStateY)
		    {
		      j1 = 0;
		      ++k1;		      
		    }
		}
	    }	      
	}
      return vDestination;
    }
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
		      Tmp += TmpWaveFunctionOverlapZ[CellZ + 1] * TmpPrecalculatedHamiltonian[CellZ];
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
		      Tmp += TmpWaveFunctionOverlapZ[CellZ + 1] * TmpPrecalculatedHamiltonian[CellZ];
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
		  Tmp += TmpWaveFunctionOverlapZ[CellZ + 1] * TmpPrecalculatedHamiltonian[CellZ];
		}
	      Tmp += TmpWaveFunctionOverlapZ[this->NbrCellZ + 1] * this->PostConstantRegionPotential;		  
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
		  Tmp += TmpWaveFunctionOverlapZ[CellZ + 1] * TmpPrecalculatedHamiltonian[CellZ];
		}
	      Tmp += TmpWaveFunctionOverlapZ[this->NbrCellZ + 1] * this->PostConstantRegionPotential;		  
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
		      Tmp += TmpWaveFunctionOverlapZ[CellZ + 1] * TmpPrecalculatedHamiltonian[CellZ];
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
		      Tmp += TmpWaveFunctionOverlapZ[CellZ + 1] * TmpPrecalculatedHamiltonian[CellZ];
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

ComplexVector& QuantumDots3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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

ComplexVector& QuantumDots3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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

ComplexVector& QuantumDots3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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

ComplexVector& QuantumDots3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							     int firstComponent, int nbrComponent)
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> QuantumDots3DHamiltonian::LeftInteractionOperators()  
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators 
//
// return value = list of right interaction operators

List<Matrix*> QuantumDots3DHamiltonian::RightInteractionOperators()  
{
  List<Matrix*> TmpList;
  return TmpList;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, QuantumDots3DHamiltonian& H)
{
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, QuantumDots3DHamiltonian& H)
{
  return Str;
}

// evaluate all interaction factors
//   
// memory = amount of memory available to store precalculated values

void QuantumDots3DHamiltonian::EvaluateInteractionFactors(int memory)
{
  int UsedMemory = 0;
  this->WaveFunctionOverlapX = this->EvaluateWaveFunctionOverlap (this->XSize, this->NbrCellX, this->NbrStateX, UsedMemory);
  this->WaveFunctionOverlapY = this->EvaluateWaveFunctionOverlap (this->YSize, this->NbrCellY, this->NbrStateY, UsedMemory);
  this->WaveFunctionOverlapZ = this->EvaluateWaveFunctionZOverlap (UsedMemory);

  UsedMemory -= this->Space->GetHilbertSpaceDimension () * sizeof(double);
  double InvXFactor = HAMILTONIAN_FACTOR / (this->Mux * this->XSize * this->XSize);
  double InvYFactor = HAMILTONIAN_FACTOR / (this->Muy * this->YSize * this->YSize);
  double InvZFactor = HAMILTONIAN_FACTOR / (this->Muz * this->ZSize * this->ZSize);
  int TotalIndex = 0;
  double FactorX = 0.0;
  double FactorY = 0.0;
  double Factor;
  double TmpElement = 0.0;
  int IncNbrCellZ = this->NbrCellZ + 1;
  this->DiagonalElements = new double[ this->Space->GetHilbertSpaceDimension ()];
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      FactorX = ((double) ((i + 1) * (i + 1))) * InvXFactor;
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  FactorY = ((double) ((j + 1) * (j + 1))) * InvYFactor + FactorX;
	  for (int k = 0; k < this->NbrStateZ; ++k)
	    {
	      TmpElement = FactorY +  ((double) ((k + 1) * (k + 1))) * InvZFactor;
	      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		{
		  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		    {
		      Factor = this->WaveFunctionOverlapX[i][i][CellX] * this->WaveFunctionOverlapY[j][j][CellY];
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  TmpElement += this->InteractionFactors[CellX][CellY][CellZ] * Factor * this->WaveFunctionOverlapZ[k][k][CellZ + 1];
			}	 
		    }     
		}
	      TmpElement += this->PostConstantRegionPotential * this->WaveFunctionOverlapZ[k][k][IncNbrCellZ];
	      this->DiagonalElements[TotalIndex] = TmpElement;
	      ++TotalIndex;
	    }
	}
    }
  memory -= UsedMemory;
  if (memory < 0)
    {
      this->NbrPrecalculatedDimension = 0;
      return;
    }

/*  if (((double) memory) > (0.5 * ((double) this->Space->GetHilbertSpaceDimension()) * 
			   ((double) this->Space->GetHilbertSpaceDimension() - 1.0) * ((double) sizeof(double))
			   + (((double) this->Space->GetHilbertSpaceDimension()) * ((double) sizeof(double*)))))
    {
      this->NbrPrecalculatedDimension = 0x100;      
      int Dim = this->Space->GetHilbertSpaceDimension();
      this->FullPrecalculatedHamiltonian = new double* [Dim];
      double Tmp;
      double* TmpWaveFunctionOverlapX;
      double* TmpWaveFunctionOverlapY;
      double* TmpWaveFunctionOverlapZ;
      int k1 = 1;
      int j1 = 0;
      int i1 = 0;
      for (int Index1 = 1; Index1 < Dim; ++Index1)
	{
	  this->FullPrecalculatedHamiltonian[Index1] = new double [Index1];
	  int Index2 = 0;
	  int i2 = 0;
	  int j2 = 0;
	  int k2 = 0;
	  while (Index2 < Index1)
	    {
	      if (i2 <= i1)
		TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i1][i2];
	      else
		TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i2][i1];
	      if (j2 <= j1)
		TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j1][j2];
	      else
		TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j2][j1];
	      if (k2 <= k1)
		TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
	      else
		TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k2][k1];
	      Tmp = 0.0;
	      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		{
		  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		    {
		      Factor = TmpWaveFunctionOverlapX[CellX] * TmpWaveFunctionOverlapY[CellY];
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * Factor * TmpWaveFunctionOverlapZ[CellZ + 1];
			}			      
		    } 
		}    
	      if ((i1 == i2) && (j1 == j2))
		{
		  Tmp += this->PostConstantRegionPotential * TmpWaveFunctionOverlapZ[IncNbrCellZ];
		}
	      this->FullPrecalculatedHamiltonian[Index1][Index2] = Tmp;  
	      ++Index2;
	      ++k2;
	      if (k2 == this->NbrStateZ)
		{
		  k2 = 0;
		  ++j2;
		  if (j2 == this->NbrStateY)
		    {
		      j2 = 0;
		      ++i2;		      
		    }
		}
	    }
	  ++k1;
	  if (k1 == this->NbrStateZ)
	    {
	      k1 = 0;
	      ++j1;
	      if (j1 == this->NbrStateY)
		{
		  j1 = 0;
		  ++i1;		      
		}
	    }

	}
      return;
    }*/

  if (((double) memory) > ((0.5 * ((double) (this->NbrStateX * this->NbrStateY)) * 
			    ((double) (this->NbrStateX * this->NbrStateY))) * 
			   (((double) sizeof(double*)) + ((double) this->NbrCellZ) * ((double) sizeof(double)))))
    {
      this->NbrPrecalculatedDimension = 0x010;      
      int Dim = this->NbrStateX * this->NbrStateY;
      this->Partial2DPrecalculatedHamiltonian = new double** [Dim];
      this->Partial2DDiagonalPrecalculatedHamiltonian = new double* [Dim];
      double Tmp;
      double* TmpWaveFunctionOverlapX;
      double* TmpWaveFunctionOverlapY;
      double* TmpPrecalculatedHamiltonian;
      int j1 = 0;
      int i1 = 0;
      for (int Index1 = 0; Index1 < Dim; ++Index1)
	{
	  this->Partial2DPrecalculatedHamiltonian[Index1] = new double* [Index1 + 1];
	  int Index2 = 0;
	  int i2 = 0;
	  int j2 = 0;
	  while (Index2 <= Index1)
	    {
	      if (i2 <= i1)
		TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i1][i2];
	      else
		TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i2][i1];
	      if (j2 <= j1)
		TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j1][j2];
	      else
		TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j2][j1];
	      this->Partial2DPrecalculatedHamiltonian[Index1][Index2] = new double [this->NbrCellZ];
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[Index1][Index2];
	      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  Tmp = 0.0;
		  for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		    {
		      for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			{
			  Tmp += this->InteractionFactors[CellX][CellY][CellZ] * TmpWaveFunctionOverlapX[CellX] * TmpWaveFunctionOverlapY[CellY];
			}			      
		    } 
		  TmpPrecalculatedHamiltonian[CellZ] = Tmp;  
		}    
	      ++Index2;
	      ++j2;
	      if (j2 == this->NbrStateY)
		{
		  j2 = 0;
		  ++i2;
		}
	    }
	  ++j1;
	  if (j1 == this->NbrStateY)
	    {
	      j1 = 0;
	      ++i1;
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

double*** QuantumDots3DHamiltonian::EvaluateWaveFunctionOverlap(double size, int nbrStep, int nbrState, int& memory)
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

// evaluate wave function overlaps on a cell in the z direction
//
// memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
// return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)

double*** QuantumDots3DHamiltonian::EvaluateWaveFunctionZOverlap(int& memory)
{
  memory += ((((this->NbrStateZ * (this->NbrStateZ + 1)) >> 1) * ((this->NbrCellZ + 2)  * sizeof(double) + sizeof(double*)))
	     + this->NbrStateZ * sizeof(double**));
  double*** TmpArray = new double** [this->NbrStateZ];
  cout << this->PreConstantRegionSize << " " << this->PostConstantRegionSize << " " << this->ZSize << " " << this->NbrCellZ << endl;
  double StepInc = (this->ZSize - this->PreConstantRegionSize - this->PostConstantRegionSize) / ((((double) this->NbrCellZ)) * this->ZSize);//(((double) this->NbrCellZ))
  double Shift = this->PreConstantRegionSize / this->ZSize;
  double PostShift = (this->ZSize -this->PostConstantRegionSize) / this->ZSize;
  double PostShift2 = this->PostConstantRegionSize / this->ZSize;
  cout << "shifts = " << StepInc << " " << Shift << " " << PostShift << endl;
  double Tmp;
  double TmpShift;
  double Diff;
  double Tmp2;
  double TmpShift2;
  double Diff2;
  for (int i = 0; i < this->NbrStateZ; ++i)
    {
      TmpArray[i] = new double* [i + 1];
      for (int j = 0; j < i; ++j)
	{
	  Diff = (double) (i - j);
	  Tmp = M_PI * Diff;
	  TmpShift = Tmp * Shift;
	  TmpArray[i][j] = new double [this->NbrCellZ + 2];
	  Diff2 = (double) (i + j + 2);
	  Tmp2 = M_PI * Diff2;	      
	  TmpShift2 = Tmp2 * Shift;
	  Diff = M_1_PI / Diff;
	  Diff2 = M_1_PI / Diff2;
	  TmpArray[i][j][0] = (sin (TmpShift) * Diff - sin (TmpShift2) * Diff2);
	  TmpArray[i][j][this->NbrCellZ + 1] = (sin (Tmp2 * PostShift) * Diff2) - (sin (Tmp * PostShift) * Diff);
	  Tmp *= StepInc;
	  Tmp2 *= StepInc;
	  for (int k = 0; k < this->NbrCellZ; ++k)
	    {
	      TmpArray[i][j][k + 1] =  ((sin ((Tmp * ((double) (k + 1))) + TmpShift) - sin ((Tmp * ((double) k)) + TmpShift)) * Diff
					- (sin ((Tmp2 * ((double) (k + 1))) + TmpShift2) - sin ((Tmp2 * ((double) k)) + TmpShift2)) * Diff2);
	    }
	}
      TmpArray[i][i] = new double [this->NbrCellZ + 2];
      Diff = (double) (2 * i + 2);
      Tmp = M_PI * Diff;	      
      TmpShift = Tmp * Shift;
      Diff = M_1_PI / Diff;
      TmpArray[i][i][0] = Shift - (sin (Tmp * Shift) * Diff);
      TmpArray[i][i][this->NbrCellZ + 1] = PostShift2 + (sin (Tmp * PostShift) * Diff);
      Tmp *= StepInc;
      for (int k = 0; k < this->NbrCellZ; ++k)
	{
	  TmpArray[i][i][k + 1] = StepInc - (sin ((Tmp * ((double) (k + 1))) + TmpShift) - sin ((Tmp * ((double) k)) + TmpShift)) * Diff;
	}     
    }
  return TmpArray;
}

