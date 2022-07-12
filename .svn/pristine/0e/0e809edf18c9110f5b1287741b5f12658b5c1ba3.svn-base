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
//                      last modification : 26/02/2004                        //
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
#include "MathTools/Complex.h"
#include "Tools/Potential/HardBoxPyramidQuantumDotThreeDConstantCellPotential.h"
#include "Tools/Potential/EllipticalDotThreeDConstantCellPotential.h"

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

QuantumDots3DHamiltonian::QuantumDots3DHamiltonian(ThreeDOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, HardBoxPyramidQuantumDotThreeDConstantCellPotential* PotentialInput, int memory)
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
  this->LeftNumber = PotentialInput->GetUnder();
  this->PreConstantRegionSize = new double [this->LeftNumber];
  this->PreConstantRegionPotential = new double [this->LeftNumber];
  for (int k = 0; k < this->LeftNumber; ++k)
    {
      this->PreConstantRegionSize[k] = PotentialInput->GetUnderSize(k);
      this->PreConstantRegionPotential[k] = PotentialInput->GetUnderPotentialValue(k);
    }
  this->RightNumber = PotentialInput->GetAbove();
  this->PostConstantRegionSize = new double [this->RightNumber];
  this->PostConstantRegionPotential = new double [this->RightNumber];
  for (int k = 0; k < this->RightNumber; ++k)
    {
      this->PostConstantRegionSize[k] = PotentialInput->GetAboveSize(k);
      this->PostConstantRegionPotential[k] = PotentialInput->GetAbovePotentialValue(k);
    }
  this->InteractionFactors = new double** [this->NbrCellZ];  
  for (int k = 0; k < this->NbrCellZ; ++k)
    {
      this->InteractionFactors[k] = new double* [this->NbrCellY];
      for (int j = 0; j < this->NbrCellY; ++j)
	{
	  this->InteractionFactors[k][j] = new double [this->NbrCellX];	
	  for (int i = 0; i < this->NbrCellX; ++i)
	    this->InteractionFactors[k][j][i] = PotentialInput->GetPotential(i, j, k + this->LeftNumber);
	}
    }
  this->EvaluateInteractionFactors(memory);
}

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

QuantumDots3DHamiltonian::QuantumDots3DHamiltonian(ThreeDOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, EllipticalDotThreeDConstantCellPotential* PotentialInput, int memory)
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
  this->LeftNumber = PotentialInput->GetUnder();
  this->PreConstantRegionSize = new double [this->LeftNumber];
  this->PreConstantRegionPotential = new double [this->LeftNumber];
  for (int k = 0; k < this->LeftNumber; ++k)
    {
      this->PreConstantRegionSize[k] = PotentialInput->GetUnderSize(k);
      this->PreConstantRegionPotential[k] = PotentialInput->GetUnderPotentialValue(k);
    }
  this->RightNumber = PotentialInput->GetAbove();
  this->PostConstantRegionSize = new double [this->RightNumber];
  this->PostConstantRegionPotential = new double [this->RightNumber];
  for (int k = 0; k < this->RightNumber; ++k)
    {
      this->PostConstantRegionSize[k] = PotentialInput->GetAboveSize(k);
      this->PostConstantRegionPotential[k] = PotentialInput->GetAbovePotentialValue(k);
    }
  this->InteractionFactors = new double** [this->NbrCellZ];  
  for (int k = 0; k < this->NbrCellZ; ++k)
    {
      this->InteractionFactors[k] = new double* [this->NbrCellY];
      for (int j = 0; j < this->NbrCellY; ++j)
	{
	  this->InteractionFactors[k][j] = new double [this->NbrCellX];	
	  for (int i = 0; i < this->NbrCellX; ++i)
	    this->InteractionFactors[k][j][i] = PotentialInput->GetPotential(i, j, k + this->LeftNumber);
	}
    }
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
  this->LeftNumber = hamiltonian.LeftNumber;
  this->PreConstantRegionSize = hamiltonian.PreConstantRegionSize;
  this->PreConstantRegionPotential = hamiltonian.PreConstantRegionPotential;
  this->RightNumber = hamiltonian.RightNumber;
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
  double Tmp = 0.0;
  double TmpTotal1 = 0.0, TmpTotal4 = 0.0, TmpTotal5 = 0.0, TmpTotal7 = 0.0;
  double TmpVSource1 = 0.0, TmpVSource4 = 0.0, TmpVSource5 = 0.0, TmpVSource7 = 0.0;
  double* TmpWaveFunctionOverlapZ;
  double* TmpPrecalculatedHamiltonian;
  int m1 = 0, n1 = 0, p1 = 0, m2 = 0, n2 = 0, p2 = 0;
  int Index1 = 0, Index2 = 0, Index3 = 0, Index4 = 0, Index5 = 0, Index6 = 0, Index7 = 0, Index8 = 0; 
  int tmpDim = this->NbrStateY * this->NbrStateZ;
  int CellZ = 0;
  int** RI = new int* [this->NbrStateX];
  
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {
      RI[m1] = new int [this->NbrStateY];
      for (n1 = 0; n1 < this->NbrStateY; ++n1)
	RI[m1][n1] = m1 * tmpDim + n1 * this->NbrStateZ;
    }
  
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {
      for (n1 = 0; n1 < this->NbrStateY; ++n1)
	{	  
	  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	    {
	      TmpTotal1 = 0.0;
	      TmpVSource1 = vSource[Index1];		  
	      for (m2 = 0; m2 < m1; ++m2)
		{ 
		  Index7 = RI[m2][n1] + p1; // Index7 = (m2 * N + n1) * H + p1
		  TmpTotal7 = 0.0;
		  TmpVSource7 = vSource[Index7];
		  for (n2 = 0; n2 < n1; ++n2)
		    {			  
		      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m1][n1][m2][n2];			  
		      Index3 = RI[m1][n1]; // Index3 = (m1 * N + n1) * H
		      Index8 = RI[m1][n2]; // Index8 = (m1 * N + n2) * H
		      Index6 = RI[m2][n1]; // Index6 = (m2 * N + n1) * H
		      Index2 = RI[m2][n2];// Index2 = (m2 * N + n2) * H
		      Index4 = Index2 + p1; // Index4 = (m2 * N + n2) * H + p1
		      Index5 = Index8 + p1; // Index5 = (m1 * N + n2) * H + p1
		      TmpTotal4 = 0.0;
		      TmpVSource4 = vSource[Index4];
		      TmpTotal5 = 0.0;
		      TmpVSource5 = vSource[Index5];
		      for (p2 = 0; p2 < p1; ++p2)
			{
			  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
			  Tmp = 0.0;
			  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
			  
			  TmpTotal1 += Tmp * vSource[Index2];
			  vDestination[Index2] += Tmp * TmpVSource1;
			  vDestination[Index3] += Tmp * TmpVSource4;
			  TmpTotal4 += Tmp * vSource[Index3];
			  TmpTotal5 += Tmp * vSource[Index6];
			  vDestination[Index6] += Tmp * TmpVSource5;			      
			  TmpTotal7 += Tmp * vSource[Index8];
			  vDestination[Index8] += Tmp * TmpVSource7;			      
			  
			  ++Index2; // Index2 = (m2 * N + n2) * H + p2			      
			  ++Index3; // Index3 = (m1 * N + n1) * H + p2
			  ++Index6; // Index6 = (m2 * N + n1) * H + p2			     
			  ++Index8; // Index8 = (m1 * N + n2) * H + p2
			}
		      // case p2 = p1
		      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p1];
		      Tmp = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		      
		      TmpTotal1 += Tmp * TmpVSource4;
		      TmpTotal4 += Tmp * TmpVSource1;			      
		      TmpTotal7 += Tmp * TmpVSource5;
		      TmpTotal5 += Tmp * TmpVSource7;
		      
		      vDestination[Index4] += TmpTotal4;
		      vDestination[Index5] += TmpTotal5;
		    } 
		  // case n2 = n1 
		  TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m1][n1][m2][n1];
		  Index3 = RI[m1][n1]; // Index3 = (m1 * N + n1) * H			  
		  Index2 = RI[m2][n1]; // Index2 = (m2 * N + n2) * H
		     		     
		  for (p2 = 0; p2 < p1; ++p2)
		    {
		      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
		      Tmp = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		      
		      TmpTotal1 += Tmp * vSource[Index2];
		      vDestination[Index2] += Tmp * TmpVSource1;
		      vDestination[Index3] += Tmp * TmpVSource7;
		      TmpTotal7 += Tmp * vSource[Index3];			      
		      
			  ++Index2; // Index2 = (m2 * N + n2) * H + p2			      
			  ++Index3; // Index3 = (m1 * N + n1) * H + p2
		    }
		  // case n2 = n1, p2 = p1
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p1];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		      
		  TmpTotal1 += Tmp * TmpVSource7;
		  TmpTotal7 += Tmp * TmpVSource1;		      
		  
		  vDestination[Index7] += TmpTotal7;
		}
	      // case m2 = m1		 		      
	      for (n2 = 0; n2 < n1; ++n2)
		{
		  TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m1][n1][m1][n2];
		  Index3 = RI[m1][n1]; // Index3 = (m1 * N + n1) * H 
		  Index2 = RI[m1][n2]; // Index2 = (m2 * N + n2) * H
		  Index4 = Index2 + p1; // Index4 = (m2 * N + n2) * H + p1
		  TmpTotal4 = 0.0;
		  TmpVSource4 = vSource[Index4];
		  for (p2 = 0; p2 < p1; ++p2)
		    {
		      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
		      Tmp = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
			  
		      TmpTotal1 += Tmp * vSource[Index2];
		      vDestination[Index2] += Tmp * TmpVSource1;
		      vDestination[Index3] += Tmp * TmpVSource4;
		      TmpTotal4 += Tmp * vSource[Index3];			      
		      
		      ++Index2; // Index2 = (m2 * N + n2) * H + p2			      
		      ++Index3; // Index3 = (m1 * N + n1) * H + p2
		    }
		  // case p2 = p1
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p1];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		      
		  TmpTotal1 += Tmp * TmpVSource4;
		  TmpTotal4 += Tmp * TmpVSource1;			      		     		      		  
		  vDestination[Index4] += TmpTotal4;		      
		} 
	      // case m2 = m1, n2 = n1 
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m1][n1][m1][n1];     
	      Index2 = RI[m1][n1]; // Index2 = (m2 * N + n2) * H
	      for (p2 = 0; p2 < p1; ++p2)
		{
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		  Tmp += this->PartialZPrecalculatedHamiltonian[p1][p2];
		  
		  TmpTotal1 += Tmp * vSource[Index2];
		  vDestination[Index2] += Tmp * TmpVSource1;			      
		  
		  ++Index2; // Index2 = (m2 * N + n2) * H + p2			      		     
		}
	      // case m2 = m1, n2 = n1, p2 = p1
	      TmpTotal1 += this->DiagonalElements[Index1] * TmpVSource1;
	      
	      vDestination[Index1] += TmpTotal1;
	      
	      ++Index1;
	    }
	}	 
    }     
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
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k1][k2][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k1][k2][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
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
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor =  this->WaveFunctionOverlapZ[k1][k2][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor =  this->WaveFunctionOverlapZ[k1][k2][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
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
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
			      }	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
		for (; i2 < this->NbrStateX; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
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
		  for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		    {
		      for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			{
			  Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j1][j1][CellY];
			  for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			    {
			      Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
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
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j1][j1][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
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
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
			      }	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
		for (; i2 < this->NbrStateX; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
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
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k2][k1][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k2][k1][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
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
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k2][k1][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k2][k1][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
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
  if (this->NbrPrecalculatedDimension == 0x001)
    {
      if ((firstComponent == 0) && (nbrComponent == this->Space->GetHilbertSpaceDimension()))
	return this->LowLevelAddMultiply(vSource, vDestination);      
      else
	{
	  int m1, m2, n1, n2, p1, p2;
	  m1 = firstComponent / (this->NbrStateY * this->NbrStateZ);
	  int tmpIndex = firstComponent - m1 * this->NbrStateY * this->NbrStateZ;
	  n1 = tmpIndex / this->NbrStateZ;
	  p1 = tmpIndex - n1 * this->NbrStateZ;
	  int Index1 = firstComponent;
	  double Tmp, TmpSum;
	  int Index2 = 0;
	  double* TmpWaveFunctionOverlapZ;
	  double* TmpPrecalculatedHamiltonian;
	  int CellZ = 0;
	  
	  for (; Index1 < LastComponent; ++Index1)
	    {
	      Index2 = 0;
	      TmpSum = 0.0;	      

	      for (m2 = 0; m2 < m1; ++m2)
		{
		  for (n2 = 0; n2 < n1; ++n2)
		    {
		      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m1][n1][m2][n2];
		      for (p2 = 0; p2 < p1; ++p2)
			{
			  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
			  Tmp = 0.0;
			  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
			  TmpSum += (Tmp * vSource[Index2]);
			  ++Index2;
			}
		      for (p2 = p1; p2 < this->NbrStateZ; ++p2)
			{
			  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p2][p1];
			  Tmp = 0.0;
			  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
			  TmpSum += (Tmp * vSource[Index2]);			  
			  ++Index2;		  
			}
		    }
		  for (n2 = n1; n2 < this->NbrStateY; ++n2)
		    {
		      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m1][n2][m2][n1];
		      for (p2 = 0; p2 < p1; ++p2)
			{
			  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
			  Tmp = 0.0;
			  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
			  TmpSum += (Tmp * vSource[Index2]);
			  ++Index2;		  
			}
		      for (p2 = p1; p2 < this->NbrStateZ; ++p2)
			{
			  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p2][p1];
			  Tmp = 0.0;
			  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
			  TmpSum += (Tmp * vSource[Index2]);
			  ++Index2;		  
			}
		    }
		}
	      // m2 = m1
	      for (n2 = 0; n2 < n1; ++n2)
		{
		  TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m2][n1][m1][n2];
		  for (p2 = 0; p2 < p1; ++p2)
		    {
		      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
		      Tmp = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
		      TmpSum += (Tmp * vSource[Index2]);
		      ++Index2;		  
		    }
		  for (p2 = p1; p2 < this->NbrStateZ; ++p2)
		    {
		      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p2][p1];
		      Tmp = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
		      TmpSum += (Tmp * vSource[Index2]);
		      ++Index2;		  
		    }
		}
	      // m2 = m1 & n2 = n1
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m1][n1][m1][n1];
	      for (p2 = 0; p2 < p1; ++p2)
		{
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
		  Tmp += this->PartialZPrecalculatedHamiltonian[p1][p2];
		  TmpSum += (Tmp * vSource[Index2]);
		  ++Index2;
		}
	      // m2 = m1 & n2 = n1 & p2 = p1
	      TmpSum += (this->DiagonalElements[Index1] * vSource[Index1]);
	      ++Index2;
	      for (p2 = p1 + 1; p2 < this->NbrStateZ; ++p2)
		{
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p2][p1];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
		  Tmp += this->PartialZPrecalculatedHamiltonian[p2][p1];
		  TmpSum += (Tmp * vSource[Index2]);
		  ++Index2;
		}	      
	      for (n2 = n1 + 1; n2 < this->NbrStateY; ++n2)
		{
		  TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m2][n2][m1][n1];
		  for (p2 = 0; p2 < p1; ++p2)
		    {
		      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
		      Tmp = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
		      TmpSum += (Tmp * vSource[Index2]);
		      ++Index2;		  
		    }
		  for (p2 = p1; p2 < this->NbrStateZ; ++p2)
		    {
		      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p2][p1];
		      Tmp = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
		      TmpSum += (Tmp * vSource[Index2]);
		      ++Index2;		  
		    }
		}	      
	      // m2 > m1
	      for (m2 = m1 + 1; m2 < this->NbrStateX; ++m2)
		{
		  for (n2 = 0; n2 < n1; ++n2)
		    {
		      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m2][n1][m1][n2];
		      for (p2 = 0; p2 < p1; ++p2)
			{
			  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
			  Tmp = 0.0;
			  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
			  TmpSum += (Tmp * vSource[Index2]);
			  ++Index2;		  
			}
		      for (p2 = p1; p2 < this->NbrStateZ; ++p2)
			{
			  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p2][p1];
			  Tmp = 0.0;
			  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
			  TmpSum += (Tmp * vSource[Index2]);
			  ++Index2;		  
			}
		    }
		  for (n2 = n1; n2 < this->NbrStateY; ++n2)
		    {
		      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m2][n2][m1][n1];
		      for (p2 = 0; p2 < p1; ++p2)
			{
			  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p1][p2];
			  Tmp = 0.0;
			  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
			  TmpSum += (Tmp * vSource[Index2]);
			  ++Index2;		  
			}
		      for (p2 = p1; p2 < this->NbrStateZ; ++p2)
			{
			  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[p2][p1];
			  Tmp = 0.0;
			  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	
			  TmpSum += (Tmp * vSource[Index2]);
			  ++Index2;		  
			}
		    }
		}
		
	      vDestination[Index1] += TmpSum; 
	      
	      ++p1;
	      if (p1 == this->NbrStateZ)
		{
		  p1 = 0;
		  ++n1;		
		  if (n1 == this->NbrStateY)
		    {
		      n1 = 0;
		      ++m1;
		    }	
		}	  
	    }	
	  return vDestination;
	}
    }
  if (this->NbrPrecalculatedDimension == 0x100)
    {
      double Tmp = 0.0; 
      for (int Index1 = firstComponent; Index1 < LastComponent; ++Index1)
	{
	  Tmp = (this->DiagonalElements[Index1] * vSource[Index1]);	  
	  for (int Index2 = 0; Index2 < Dim; ++Index2)
	    Tmp += this->FullPrecalculatedHamiltonian[Index1][Index2] * vSource[Index2];	  
	  vDestination[Index1] += Tmp;
	}
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

  this->DiagonalElements = new double[this->Space->GetHilbertSpaceDimension ()];
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      FactorX = ((double) ((i + 1) * (i + 1))) * InvXFactor;
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  FactorY = ((double) ((j + 1) * (j + 1))) * InvYFactor + FactorX;
	  for (int k = 0; k < this->NbrStateZ; ++k)
	    {
	      TmpElement = FactorY + ((double) ((k + 1) * (k + 1))) * InvZFactor  + (this->PartialZPrecalculatedHamiltonian[k][k]);
	      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		    {
		      Factor = this->WaveFunctionOverlapZ[k][k][CellZ] * this->WaveFunctionOverlapY[j][j][CellY];
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)			
			TmpElement += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i][i][CellX];			
		    }
		}
	      this->DiagonalElements[TotalIndex] = TmpElement;
	      ++TotalIndex;
	    }
	}
    }
  memory -= UsedMemory;
  if (memory < 0)
    {
      cout << "Calculation from scratch is done" << endl;
      this->NbrPrecalculatedDimension = 0;
      return;
    }

  if ( (double)memory > ((double(this->NbrStateX * this->NbrStateX) * double(this->NbrStateY * this->NbrStateY) * double(this->NbrStateZ * this->NbrStateZ)) * (double(sizeof(double)))))
    { 
      cout << "Calculation from fully calculated Hamiltonian" << endl << endl;
      this->NbrPrecalculatedDimension = 0x100;
      int Dim = this->NbrStateX * this->NbrStateY * this->NbrStateZ;
      this->FullPrecalculatedHamiltonian = new double* [Dim];
      double**** Inter1 = new double*** [this->NbrStateX];
      double tmp = 0.0;

      for (int m1 = 0; m1 < this->NbrStateX; ++m1)
	{
	  Inter1[m1] = new double** [this->NbrStateX];
	  for (int m2 = 0; m2 < m1; ++m2)
	    {
	      Inter1[m1][m2] = new double* [this->NbrCellZ];
	      Inter1[m2][m1] = new double* [this->NbrCellZ];
	      for (int k = 0; k < this->NbrCellZ; ++k)
		{
		  Inter1[m1][m2][k] = new double [this->NbrCellY];
		  Inter1[m2][m1][k] = new double [this->NbrCellY];
		  for (int j = 0; j < this->NbrCellY; ++j)
		    {
		      tmp = 0.0;
		      for (int i = 0; i < this->NbrCellX; ++i)
			tmp += this->InteractionFactors[k][j][i] * this->WaveFunctionOverlapX[m1][m2][i];
		      Inter1[m1][m2][k][j] = tmp;
		      Inter1[m2][m1][k][j] = tmp;
		    }
		}
	    }
	  Inter1[m1][m1] = new double* [this->NbrCellZ];
	  for (int k = 0; k < this->NbrCellZ; ++k)
	    {
	      Inter1[m1][m1][k] = new double [this->NbrCellY];
	      for (int j = 0; j < this->NbrCellY; ++j)
		{
		  tmp = 0.0;
		  for (int i = 0; i < this->NbrCellX; ++i)
		    tmp += this->InteractionFactors[k][j][i] * this->WaveFunctionOverlapX[m1][m1][i];
		  Inter1[m1][m1][k][j] = tmp;
		}
	    }
	}
      double***** Inter2 = new double**** [this->NbrStateY];
      for (int n1 = 0; n1 < this->NbrStateY; ++n1)
	{
	  Inter2[n1] = new double*** [this->NbrStateY];
	  for (int n2 = 0; n2 < n1; ++n2)
	    {
	      Inter2[n1][n2] = new double** [this->NbrStateX];
	      Inter2[n2][n1] = new double** [this->NbrStateX];
	      for (int m1 = 0; m1 < this->NbrStateX; ++m1)
		{
		  Inter2[n1][n2][m1] = new double* [this->NbrStateX];
		  Inter2[n2][n1][m1] = new double* [this->NbrStateX];
		  for(int m2 = 0; m2 < this->NbrStateX; ++m2)
		    {
		      Inter2[n1][n2][m1][m2] = new double [this->NbrCellZ];
		      Inter2[n2][n1][m1][m2] = new double [this->NbrCellZ];
		      for (int k = 0; k < this->NbrCellZ; ++k)
			{
			  tmp = 0.0;
			  for (int j = 0; j < this->NbrCellY; ++j)
			    tmp += this->WaveFunctionOverlapY[n1][n2][j] * Inter1[m1][m2][k][j];
			  Inter2[n1][n2][m1][m2][k] = tmp;
			  Inter2[n2][n1][m1][m2][k] = tmp;
			}
		    }
		}
	    }
	  Inter2[n1][n1] = new double** [this->NbrStateX];
	  for (int m1 = 0; m1 < this->NbrStateX; ++m1)
	    {
	      Inter2[n1][n1][m1] = new double* [this->NbrStateX];
	      for(int m2 = 0; m2 < this->NbrStateX; ++m2)
		{
		  Inter2[n1][n1][m1][m2] = new double [this->NbrCellZ];
		  for (int k = 0; k < this->NbrCellZ; ++k)
		    {
		      tmp = 0.0;
		      for (int j = 0; j < this->NbrCellY; ++j)
			tmp += this->WaveFunctionOverlapY[n1][n1][j] * Inter1[m1][m2][k][j];
		      Inter2[n1][n1][m1][m2][k] = tmp;
		    }
		}
	    }
	}

      int ReducedIndex1 = 0, ReducedIndex2 = 0, m1 = 0, m2 = 0, n1 = 0, n2 = 0, p1 = 0, p2 = 0;
     
      for (int Index1 = 0; Index1 < Dim; ++Index1)
	{
	  this->FullPrecalculatedHamiltonian[Index1] = new double [Dim];
	  ReducedIndex1 = Index1 / this->NbrStateZ;
	  p1 = Index1 - ReducedIndex1 * this->NbrStateZ;
	  m1 = ReducedIndex1 / this->NbrStateY;
	  n1 = ReducedIndex1 - m1 * this->NbrStateY;
	  for (int Index2 = 0; Index2 < Index1; ++Index2)
	    {
	      ReducedIndex2 = Index2 / this->NbrStateZ;
	      p2 = Index2 - ReducedIndex2 * this->NbrStateZ;
	      m2 = ReducedIndex2 / this->NbrStateY;
	      n2 = ReducedIndex2 - m2 * this->NbrStateY;
	      tmp = 0.0;
	      if (p1 > p2)
		for (int k = 0; k < this->NbrCellZ; ++k)
		  tmp += this->WaveFunctionOverlapZ[p1][p2][k] * Inter2[n1][n2][m1][m2][k];
	      else
		for (int k = 0; k < this->NbrCellZ; ++k)
		  tmp += this->WaveFunctionOverlapZ[p2][p1][k] * Inter2[n1][n2][m1][m2][k];
	      if (ReducedIndex1 == ReducedIndex2)
		tmp += this->PartialZPrecalculatedHamiltonian[p1][p2];
	      this->FullPrecalculatedHamiltonian[Index1][Index2] = tmp;
	      this->FullPrecalculatedHamiltonian[Index2][Index1] = tmp;	      
	    }
          this->FullPrecalculatedHamiltonian[Index1][Index1] = 0.0;
	}
      delete[] Inter1; delete[] Inter2;      
      return; 
    }
  if (((double) memory) > ((0.25 * ((double) (this->NbrStateX * this->NbrStateY)) * 
			    ((double) (this->NbrStateX * this->NbrStateY))) * 
			   (((double) sizeof(double*)) + ((double) this->NbrCellZ) * ((double) sizeof(double)))))
    {
      cout << "Calculation with symmetric precalculated Hamiltonian is done" << endl;
      this->NbrPrecalculatedDimension = 0x001;
      double* TmpPrecalculatedHamiltonian;
      double* TmpWaveFunctionOverlapX;
      double* TmpWaveFunctionOverlapY;
      int m1 = 0, m2 = 0, n1 = 0, n2 = 0;
      double Tmp = 0.0, Tmp2 = 0.0;
      int CellX = 0, CellY = 0, CellZ = 0;
      this->Partial2DPrecalculatedHamiltonian = new double**** [this->NbrStateX];
      for (m1 = 0; m1 < this->NbrStateX; ++m1)
	{	  
	  this->Partial2DPrecalculatedHamiltonian[m1] = new double*** [this->NbrStateY];	
	  for (n1 = 0; n1 < this->NbrStateY; ++n1)
	    {
	      this->Partial2DPrecalculatedHamiltonian[m1][n1] = new double** [m1 + 1];	      
	      for (m2 = 0; m2 <= m1; ++m2)
		{	      
		  this->Partial2DPrecalculatedHamiltonian[m1][n1][m2] = new double* [n1 + 1];	   
		  TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[m1][m2];  
		  for (n2 = 0; n2 <= n1; ++n2)
		    {
		      this->Partial2DPrecalculatedHamiltonian[m1][n1][m2][n2] = new double [this->NbrCellZ];
		      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[m1][n1][m2][n2];
		      TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[n1][n2];
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  Tmp = 0.0;
			  for (CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Tmp2 = TmpWaveFunctionOverlapY[CellY];
			      for (CellX = 0; CellX < this->NbrCellX; ++CellX)			     
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * TmpWaveFunctionOverlapX[CellX] * Tmp2;
			    } 
			  TmpPrecalculatedHamiltonian[CellZ] = Tmp;
			}
		    }
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
  memory += ((((this->NbrStateZ * (this->NbrStateZ + 1)) >> 1) * ((this->NbrCellZ)  * sizeof(double) + sizeof(double*)))
	     + this->NbrStateZ * sizeof(double**));
  double*** TmpArray = new double** [this->NbrStateZ];
  this->PartialZPrecalculatedHamiltonian = new double*[this->NbrStateZ];
  double LeftRegionSize = 0.0;
  for (int k = 0; k < this->LeftNumber; ++k)
    LeftRegionSize += PreConstantRegionSize[k];
  double RightRegionSize = 0.0;
  for (int k = 0; k < this->RightNumber; ++k)
    RightRegionSize += PostConstantRegionSize[k];

  double StepInc = (this->ZSize - LeftRegionSize - RightRegionSize) / ((((double) this->NbrCellZ)) * this->ZSize);
  double Shift = LeftRegionSize / this->ZSize;
  double PostShift = this->ZSize - RightRegionSize;
  double Tmp;
  double TmpShift;
  double Diff;
  double Tmp2;
  double TmpShift2;
  double Diff2;
  double TmpLength1; double TmpLength2; double Tmp3; double Tmp4; double Tmp5;
  for (int i = 0; i < this->NbrStateZ; ++i)
    {
      TmpArray[i] = new double* [i + 1];
      this->PartialZPrecalculatedHamiltonian[i] = new double [i + 1];
      for (int j = 0; j < i; ++j)
	{
	  Diff = (double) (i - j);
	  Tmp = M_PI * Diff;
	  TmpShift = Tmp * Shift;
	  TmpArray[i][j] = new double [this->NbrCellZ];
	  Diff2 = (double) (i + j + 2);
	  Tmp2 = M_PI * Diff2;
	  TmpShift2 = Tmp2 * Shift;
	  Diff = M_1_PI / Diff;
	  Diff2 = M_1_PI / Diff2;
	  TmpLength1 = 0.0; TmpLength2 = 0.0;
	  Tmp3 = Tmp/this->ZSize;
	  Tmp4 = Tmp2/this->ZSize;
	  Tmp5 = 0.0;
	  for (int k = 0; k < this->LeftNumber; ++k)
	    {
	      TmpLength1 += PreConstantRegionSize[k];
	      Tmp5 +=  ((sin(Tmp3 * TmpLength1) - sin(Tmp3 * TmpLength2)) * Diff
		                  - (sin(Tmp4 * TmpLength1) - sin(Tmp4 * TmpLength2)) * Diff2) * (this->PreConstantRegionPotential[k]);
	      TmpLength2 += PreConstantRegionSize[k];
	    }

	  Tmp *= StepInc;
	  Tmp2 *= StepInc;
	  for (int k = 0; k < this->NbrCellZ; ++k)
	    {
	      TmpArray[i][j][k] =  ((sin ((Tmp * ((double) (k + 1))) + TmpShift) - sin ((Tmp * ((double) k)) + TmpShift)) * Diff
					     - (sin ((Tmp2 * ((double) (k + 1))) + TmpShift2) - sin ((Tmp2 * ((double) k)) + TmpShift2)) * Diff2);
	    }
	  TmpLength1 = PostShift;
	  TmpLength2 = PostShift;
	  for (int k = 0; k < this->RightNumber; ++k)
	    {
	      TmpLength1 += PostConstantRegionSize[k];
	      Tmp5 +=  ((sin(Tmp3 * TmpLength1) - sin(Tmp3 * TmpLength2)) * Diff
		     - (sin(Tmp4 * TmpLength1) - sin(Tmp4 * TmpLength2)) * Diff2) * (this->PostConstantRegionPotential[k]);
	      TmpLength2 += PostConstantRegionSize[k];
	    }
	  this->PartialZPrecalculatedHamiltonian[i][j] = Tmp5;
	}

      TmpArray[i][i] = new double [this->NbrCellZ];
      Diff = (double) (2 * i + 2);
      Tmp = M_PI * Diff;
      TmpShift = Tmp * Shift;
      Diff = M_1_PI / Diff;
      Tmp3 = Tmp/this->ZSize;
      TmpLength1 = 0.0; TmpLength2 = 0.0; Tmp5 = 0.0;
      for (int k = 0; k < this->LeftNumber; ++k)
	{
	  TmpLength1 += PreConstantRegionSize[k];
	  Tmp5 += ((PreConstantRegionSize[k]/this->ZSize) - (sin(Tmp3 * TmpLength1) - sin(Tmp3 * TmpLength2)) * Diff) * (this->PreConstantRegionPotential[k]);
	  TmpLength2 += PreConstantRegionSize[k];
	}
      TmpLength1 = PostShift; TmpLength2 = PostShift;
      for (int k = 0; k < this->RightNumber; ++k)
	{
	  TmpLength1 += PostConstantRegionSize[k];
	  Tmp5 += ((PostConstantRegionSize[k]/this->ZSize) - (sin(Tmp3 * TmpLength1) - sin(Tmp3 * TmpLength2)) * Diff) * (this->PostConstantRegionPotential[k]);
	  TmpLength2 += PostConstantRegionSize[k];
	}
      this->PartialZPrecalculatedHamiltonian[i][i] = Tmp5;

      Tmp *= StepInc;
      for (int k = 0; k < this->NbrCellZ; ++k)
	{
	  TmpArray[i][i][k] = StepInc - (sin ((Tmp * ((double) (k + 1))) + TmpShift) - sin ((Tmp * ((double) k)) + TmpShift)) * Diff;
	}
    }
  return TmpArray;
}

// determine the maximal value of partial diagonal array
//
// return = the wanted value

double QuantumDots3DHamiltonian::MaxPartialDiagonalElement()
{
  double tmp = this->DiagonalElements[0];
  for (int i = 1; i < this->Space->GetHilbertSpaceDimension(); ++i)
    if (tmp < this->DiagonalElements[i])
      tmp = this->DiagonalElements[i];
  return tmp;
}
