////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of fermions on square lattice                 //
//                                  in momentum space                         //
//                                                                            //
//                        last modification : 16/02/2011                      //
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
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeNonPeriodicMomentumSpace.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnSquareLatticeMomentumSpace::FermionOnSquareLatticeMomentumSpace ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeMomentumSpace::FermionOnSquareLatticeMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateLzMax = new int [this->HilbertSpaceDimension];  
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l);
      int TmpLzMax = this->LzMax;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  this->StateLzMax[i] = TmpLzMax;
	}
      this->GenerateLookUpTable(memory);      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSquareLatticeMomentumSpace::FermionOnSquareLatticeMomentumSpace(const FermionOnSquareLatticeMomentumSpace& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnSquareLatticeMomentumSpace::~FermionOnSquareLatticeMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeMomentumSpace& FermionOnSquareLatticeMomentumSpace::operator = (const FermionOnSquareLatticeMomentumSpace& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSquareLatticeMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeMomentumSpace(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured
bool FermionOnSquareLatticeMomentumSpace::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrFermions);
  WriteLittleEndian(File, this->NbrSiteX);
  WriteLittleEndian(File, this->NbrSiteY);
  WriteLittleEndian(File, this->KxMomentum);
  WriteLittleEndian(File, this->KyMomentum);
  if (this->HilbertSpaceDimension != 0)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateLzMax[i]);
    }
  else
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateLzMax[i]);
    }
  File.close();
  return true;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSquareLatticeMomentumSpace::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> i);
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if ((Tmp & 0x1l) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ")";
    }
  Str << "]";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      this->StateDescription[pos] = 0x1ul << ((currentKx * this->NbrSiteY) + j);
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
 		  this->StateDescription[pos] = 0x1ul << ((i * this->NbrSiteY) + j);
 		  ++pos;
		}
	    }
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  unsigned long Mask = 0x1ul << ((currentKx * this->NbrSiteY) + currentKy);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos);
};


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnSquareLatticeMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    ++Count;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		++Count;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given momentum sector and fixed number of particles. 
// 
// subsytemSizeX = number of sites along the x direction that belong to the subsytem
// subsytemSizeY = number of sites along the y direction that belong to the subsytem
// subsytemStartX = x momentum marking the beginning of the rectangluar subsystem
// subsytemStartY = y momentum marking the beginning of the rectangluar subsystem
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
// kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem

HermitianMatrix FermionOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixMomentumSpace (int subsytemSizeX, int subsytemSizeY, int subsytemStartX, int subsytemStartY, int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState)
{
  subsytemSizeX = 1 + (-1 + subsytemSizeX + this->NbrSiteX) % this->NbrSiteX; // enforce subsystemSize in [1,NbrSiteX]
  subsytemSizeY = 1 + (-1 + subsytemSizeY + this->NbrSiteY) % this->NbrSiteY;
  if ((subsytemSizeX != this->NbrSiteX) && (subsytemSizeY != this->NbrSiteY))
  {
      cout << "not implemented yet!" << endl;
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
  }

  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementarySubsytemSizeX;
  int ComplementarySubsytemSizeY;
  int ComplementarySubsytemStartX;
  int ComplementarySubsytemStartY;
  if (subsytemSizeX == this->NbrSiteX)
  {
      ComplementarySubsytemSizeX = this->NbrSiteX;
      ComplementarySubsytemStartX = subsytemStartX;
      ComplementarySubsytemSizeY = this->NbrSiteY - subsytemSizeY;
      ComplementarySubsytemStartY = (subsytemStartY + subsytemSizeY) % this->NbrSiteY;
  }
  else
  {
      ComplementarySubsytemSizeX = this->NbrSiteX - subsytemSizeX;
      ComplementarySubsytemStartX = (subsytemStartX + subsytemSizeX) % this->NbrSiteX;
      ComplementarySubsytemSizeY = this->NbrSiteY;
      ComplementarySubsytemStartY = subsytemStartY;
  }

  int ComplementaryKxMomentum = (this->KxMomentum - kxSector + this->NbrSiteX) % this->NbrSiteX;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector + this->NbrSiteY) % this->NbrSiteY;

  FermionOnSquareLatticeNonPeriodicMomentumSpace SubsytemSpace (nbrParticleSector, this->NbrSiteX, this->NbrSiteY, subsytemSizeX, subsytemSizeY, subsytemStartX, subsytemStartY, kxSector, kySector);
  FermionOnSquareLatticeNonPeriodicMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, ComplementarySubsytemSizeX, ComplementarySubsytemSizeY, ComplementarySubsytemStartX, ComplementarySubsytemStartY, ComplementaryKxMomentum, ComplementaryKyMomentum);

  if (SubsytemSpace.HilbertSpaceDimension>0 && ComplementarySpace.HilbertSpaceDimension>0)
  {
      cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension
           << "; complementary Hilbert space dimension = " << ComplementarySpace.HilbertSpaceDimension << endl;
      HermitianMatrix TmpDensityMatrix (SubsytemSpace.HilbertSpaceDimension, true);
      int nbrNonZeroElements = this->EvaluatePartialDensityMatrixMomentumSpaceCore (0, ComplementarySpace.HilbertSpaceDimension, &ComplementarySpace, &SubsytemSpace, groundState, &TmpDensityMatrix);
      if (nbrNonZeroElements>0)
          return TmpDensityMatrix;
  }

  HermitianMatrix TmpDensityMatrixZero (SubsytemSpace.HilbertSpaceDimension, true);
  return TmpDensityMatrixZero;
}

// core part of the evaluation density matrix momentum space partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixMomentumSpaceCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace, ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
  FermionOnSquareLatticeMomentumSpace* TmpHilbertSpace =  (FermionOnSquareLatticeNonPeriodicMomentumSpace*) complementaryHilbertSpace;
  FermionOnSquareLatticeMomentumSpace* TmpDestinationHilbertSpace =  (FermionOnSquareLatticeNonPeriodicMomentumSpace*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
 	      int TmpLzMax = this->LzMax;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  double Coefficient = 1.0;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = TmpDestinationHilbertSpace->LzMax;
		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
		    {
		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			--Pos2;
		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState3 ^= TmpState3 >> 32;
#endif	
		      TmpState3 ^= TmpState3 >> 16;
		      TmpState3 ^= TmpState3 >> 8;
		      TmpState3 ^= TmpState3 >> 4;
		      TmpState3 ^= TmpState3 >> 2;
		      TmpState3 ^= TmpState3 >> 1;
		      Sign ^= TmpState3;
		      TmpState2 &= ~(0x1ul << Pos2);
		      --Pos2;
		    }
 		  if ((Sign & 0x1ul) == 0x0ul)		  
 		    Coefficient *= 1.0;
 		  else
 		    Coefficient *= -1.0;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = Coefficient;
		  ++Pos;
		}
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given momentum sector and fixed number of particles
// 
// subsytemSizeX = number of sites along the x direction that belong to the subsytem
// subsytemSizeY = number of sites along the y direction that belong to the subsytem
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
// kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = entanglement matrix of the subsytem
  
ComplexMatrix FermionOnSquareLatticeMomentumSpace::EvaluatePartialEntanglementMatrixMomentumSpace (int subsytemSizeX, int subsytemSizeY, int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState)
{
/* // fix later
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementarySubsytemSizeX = this->NbrSiteX - subsytemSizeX;
  int ComplementarySubsytemSizeY = this->NbrSiteY - subsytemSizeY;
  int MinComplementaryKxMomentum = (this->KxMomentum - kxSector) % ComplementarySubsytemSizeX;
  int MinComplementaryKyMomentum = (this->KyMomentum - kySector) % ComplementarySubsytemSizeY;
  int MaxComplementaryKxMomentum =  ComplementarySubsytemSizeX * ComplementaryNbrParticles - (((ComplementaryNbrParticles - 1) * ComplementaryNbrParticles) / 2);
  int MaxComplementaryKyMomentum =  ComplementarySubsytemSizeY * ComplementaryNbrParticles - (((ComplementaryNbrParticles - 1) * ComplementaryNbrParticles) / 2);
  long TotalComplementaryDimension = 0l;
  int NbrComplementarySpaces = (((MaxComplementaryKxMomentum - MinComplementaryKxMomentum + 1) / this->NbrSiteX)
				* ((MaxComplementaryKyMomentum - MinComplementaryKyMomentum + 1) / this->NbrSiteY));
  FermionOnSquareLatticeNonPeriodicMomentumSpace SubsytemSpace (nbrParticleSector, subsytemSizeX, subsytemSizeY, kxSector, kySector);
  FermionOnSquareLatticeNonPeriodicMomentumSpace** ComplementarySpaces = new FermionOnSquareLatticeNonPeriodicMomentumSpace*[NbrComplementarySpaces];
  NbrComplementarySpaces = 0;
  for (int kx = MinComplementaryKxMomentum; kx < MaxComplementaryKxMomentum; kx += this->NbrSiteX)
    for (int ky = MinComplementaryKyMomentum; ky < MaxComplementaryKyMomentum; ky += this->NbrSiteY)
      {
	ComplementarySpaces[NbrComplementarySpaces] = new FermionOnSquareLatticeNonPeriodicMomentumSpace (ComplementaryNbrParticles, ComplementarySubsytemSizeX, ComplementarySubsytemSizeY, kx, ky);
	TotalComplementaryDimension += ComplementarySpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
	++NbrComplementarySpaces;
      }
  ComplexMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), TotalComplementaryDimension, true);
  long Shift = 0;
  NbrComplementarySpaces = 0;
  for (int kx = MinComplementaryKxMomentum; kx < MaxComplementaryKxMomentum; kx += this->NbrSiteX)
    for (int ky = MinComplementaryKyMomentum; ky < MaxComplementaryKyMomentum; ky += this->NbrSiteY)
      {
	for (int i = 0; i < ComplementarySpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension(); ++i)
	  {
	    for (int j = 0; j < SubsytemSpace.GetHilbertSpaceDimension(); ++j)
	      {
	      }
	  }
	Shift += ComplementarySpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
	delete[] ComplementarySpaces[NbrComplementarySpaces];
	++NbrComplementarySpaces;	
      }
*/
  ComplexMatrix TmpEntanglementMatrix (1,2, true);
  return TmpEntanglementMatrix;
}
  
// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
// kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix FermionOnSquareLatticeMomentumSpace::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, bool removeBinomialCoefficient)
{
  ComplexMatrix TmpEntanglementMatrix;
  return TmpEntanglementMatrix;
}
 
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementaryKxMomentum = (this->KxMomentum - kxSector) % this->NbrSiteX;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector) % this->NbrSiteY;
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->NbrSiteX;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrSiteY;
  cout << "kx = " << this->KxMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
  cout << "ky = " << this->KyMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
  FermionOnSquareLatticeMomentumSpace SubsytemSpace (nbrParticleSector, this->NbrSiteX, this->NbrSiteY, kxSector, kySector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnSquareLatticeMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, ComplementaryKxMomentum, ComplementaryKyMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}


// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
											     ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
  FermionOnSquareLatticeMomentumSpace* TmpHilbertSpace =  (FermionOnSquareLatticeMomentumSpace*) complementaryHilbertSpace;
  FermionOnSquareLatticeMomentumSpace* TmpDestinationHilbertSpace =  (FermionOnSquareLatticeMomentumSpace*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  BinomialCoefficients TmpBinomial (this->NbrFermions);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
 	      int TmpLzMax = this->LzMax;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  double Coefficient = TmpInvBinomial;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = TmpDestinationHilbertSpace->LzMax;
		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
		    {
		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			--Pos2;
		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState3 ^= TmpState3 >> 32;
#endif	
		      TmpState3 ^= TmpState3 >> 16;
		      TmpState3 ^= TmpState3 >> 8;
		      TmpState3 ^= TmpState3 >> 4;
		      TmpState3 ^= TmpState3 >> 2;
		      TmpState3 ^= TmpState3 >> 1;
		      Sign ^= TmpState3;
		      TmpState2 &= ~(0x1ul << Pos2);
		      --Pos2;
		    }
 		  if ((Sign & 0x1ul) == 0x0ul)		  
 		    Coefficient *= 1.0;
 		  else
 		    Coefficient *= -1.0;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = Coefficient;
		  ++Pos;
		}
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->UnsafeAddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
}

// evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, 
												    int nbrGroundStates, ComplexVector* groundStates, double* weights, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 0.0;
	  for (int i = 0; i < nbrGroundStates; ++i)
	    TmpDensityMatrix(0, 0) += weights[i];
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 0.0;
	  for (int i = 0; i < nbrGroundStates; ++i)
	    TmpDensityMatrix(0, 0) += weights[i];
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementaryKxMomentum = (this->KxMomentum - kxSector) % this->NbrSiteX;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector) % this->NbrSiteY;
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->NbrSiteX;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrSiteY;
  cout << "kx = " << this->KxMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
  cout << "ky = " << this->KyMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
  FermionOnSquareLatticeMomentumSpace SubsytemSpace (nbrParticleSector, this->NbrSiteX, this->NbrSiteY, kxSector, kySector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnSquareLatticeMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, ComplementaryKxMomentum, ComplementaryKyMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, nbrGroundStates, groundStates, weights, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// core part of the evaluation density matrix particle partition calculation involving a sum of projetors 
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
											     int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix)
{
  FermionOnSquareLatticeMomentumSpace* TmpHilbertSpace =  (FermionOnSquareLatticeMomentumSpace*) complementaryHilbertSpace;
  FermionOnSquareLatticeMomentumSpace* TmpDestinationHilbertSpace =  (FermionOnSquareLatticeMomentumSpace*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  BinomialCoefficients TmpBinomial (this->NbrFermions);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
  Complex* TmpValues = new Complex[nbrGroundStates];
  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
 	      int TmpLzMax = this->LzMax;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  double Coefficient = TmpInvBinomial;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = TmpDestinationHilbertSpace->LzMax;
		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
		    {
		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			--Pos2;
		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState3 ^= TmpState3 >> 32;
#endif	
		      TmpState3 ^= TmpState3 >> 16;
		      TmpState3 ^= TmpState3 >> 8;
		      TmpState3 ^= TmpState3 >> 4;
		      TmpState3 ^= TmpState3 >> 2;
		      TmpState3 ^= TmpState3 >> 1;
		      Sign ^= TmpState3;
		      TmpState2 &= ~(0x1ul << Pos2);
		      --Pos2;
		    }
 		  if ((Sign & 0x1ul) == 0x0ul)		  
 		    Coefficient *= 1.0;
 		  else
 		    Coefficient *= -1.0;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = Coefficient;
		  ++Pos;
		}
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      for (int l = 0; l < nbrGroundStates; ++l)
		TmpValues[l] = weights[l] * Conj(groundStates[l][TmpStatePosition[j]]) * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    Complex Tmp = 0.0;
		    for (int l = 0; l < nbrGroundStates; ++l)
		      Tmp += TmpValues[l] * groundStates[l][TmpStatePosition[k]] * TmpStateCoefficient[k];
		    densityMatrix->UnsafeAddToMatrixElement(Pos2, TmpStatePosition2[k], Tmp);							      ;
		  }
	    }
	}
    }
  delete[] TmpValues;
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
}

// find state index from an array
//
// stateDescription = array describing the state (stored as kx1,ky1,kx2,ky2,...)
// return value = corresponding index, -1 if an error occured

int FermionOnSquareLatticeMomentumSpace::FindStateIndexFromArray(int* stateDescription)
{
  unsigned long TmpState = 0x0ul;
  for (int i = 0; i < this->NbrFermions; ++i)
    TmpState |= 0x1ul << ((stateDescription[i << 1] * this->NbrSiteY) + stateDescription[1 + (i << 1)]);
  int TmpLzMax = this->LzMax;
  while ((TmpState >> TmpLzMax) == 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}

// apply the inversion symmetry i.e (k_x,k_y)->(-k_x,-k_y) to a state 
//
// inputstate = reference on the input state
// inputSpace = pointer to the Hilbert space associated to the input state
// return value = resulting state 

ComplexVector FermionOnSquareLatticeMomentumSpace::InversionSymmetry(ComplexVector& state, FermionOnSquareLatticeMomentumSpace* inputSpace)
{
  ComplexVector OutputState (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState2 = this->StateDescription[i];
      unsigned long TmpState = 0x0ul;
      unsigned long Mask;
      unsigned long TmpState3;
      unsigned long MaskSign = 0x0ul;
      int TmpKx = 0;
      int TmpKy = 0;
      for (int j = 0; j <= this->LzMax; ++j)
	{
	  if ((TmpState2 & (0x1ul << (TmpKx * this->NbrSiteY + TmpKy))) != 0x0ul)
	    {
	      Mask = 0x1ul << (((this->NbrSiteX - TmpKx) % this->NbrSiteX) * this->NbrSiteY + ((this->NbrSiteY - TmpKy) % this->NbrSiteY));
	      TmpState3 = TmpState & (Mask - 0x1ul);
#ifdef __64_BITS__
	      TmpState3 ^= (TmpState3 >> 32);
#endif
	      TmpState3 ^= (TmpState3 >> 16);
	      TmpState3 ^= (TmpState3 >> 8);
	      TmpState3 ^= (TmpState3 >> 4);
	      TmpState3 ^= (TmpState3 >> 2);
	      MaskSign ^= (TmpState3 ^ (TmpState3 >> 1)) & 0x1ul;
	      TmpState |= Mask;
	    }
	  ++TmpKy;
	  if (TmpKy == this->NbrSiteY)
	    {
	      TmpKy = 0;
	      ++TmpKx;
	    }
	}
      int TmpLzMax = inputSpace->LzMax;
      while ((TmpState >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      int Index = inputSpace->FindStateIndex(TmpState, TmpLzMax);
      cout << i << " " << Index << " : " << hex << TmpState << " " 
	   << TmpState2 << dec << " : " ;
      this->PrintState(cout, i) << " ";
      this->PrintState(cout, Index) << " : ";
      cout << state[i] << " " << state[Index] << " " << SqrNorm(state[i]) << " " << SqrNorm(state[Index]) << endl;
      if (Index < inputSpace->HilbertSpaceDimension)
	{
	  if (MaskSign == 0ul)
	    {
	      OutputState[i] = state[Index];
	    }
	  else
	    {
	      OutputState[i] = -state[Index];
	    }
	}
      else 
	{
	  cout << "error" << endl;
	}
    }
  return OutputState;
}


Complex FermionOnSquareLatticeMomentumSpace::ComputeOverlapWaveFunctionsWithDifferentGamma (ComplexVector& firstVector, ComplexVector& secondVector, Complex * overlapMatrix)
{
   Complex Result (0,0);
   unsigned long * TmpState = new unsigned long [this->NbrFermions];
   for (int i = 0; i < this->HilbertSpaceDimension; i++)
    {
       this->GetMonomial(i, TmpState);	
       Complex Tmp =  overlapMatrix[TmpState[0]];
       for (int j = 1; j < this->NbrFermions; j++)
       {
    	 Tmp *= overlapMatrix[TmpState[j]];
       }
      Result += Conj(secondVector[i]) *  firstVector[i] * Tmp; 
    } 
 delete []  TmpState;
 return Result;
}
