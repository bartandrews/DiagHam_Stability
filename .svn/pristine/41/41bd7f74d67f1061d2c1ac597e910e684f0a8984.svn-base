////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                     class author : Cecile Repellin                         //
//                                                                            //
//                        class of fermions the CP2                           //
//                                                                            //
//                                                                            //
//                        last modification : 25/02/2013                      //
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
#include "GeneralTools/Endian.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "FermionOnCP2.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "MathTools/FactorialCoefficient.h"


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

FermionOnCP2::FermionOnCP2 ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// p = number of flux quanta (determines an irreducible representation of SU(3), along with q=0 (LLL))
// totalJz = total value of jz
// totalKz = total value of kz
// memory = amount of memory granted for precalculations

FermionOnCP2::FermionOnCP2 (int nbrFermions, int nbrFluxQuanta, int totalTz, int totalY, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->MinY = -2*this->NbrFluxQuanta;
  this->TotalTz = totalTz;
  this->TotalY = totalY;
  this->TotalR = (this->TotalY + 3*this->TotalTz + 2*this->NbrFermions*this->NbrFluxQuanta)/6;
  this->TotalS = (this->TotalY - 3*this->TotalTz + 2*this->NbrFermions*this->NbrFluxQuanta)/6;
  this->NbrLzValue = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2;
  this->LzMax = NbrLzValue - 1;  
  this->MaximumSignLookUp = 16;
  if (this->NbrFermions > 0)
    this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFluxQuanta, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0);
  else
    this->LargeHilbertSpaceDimension = 1l;
  this->quantumNumberTz = new int [this->NbrLzValue];
  this->quantumNumberY = new int [this->NbrLzValue];
  this->quantumNumberR = new int [this->NbrLzValue];
  this->quantumNumberS = new int [this->NbrLzValue];
  cout << "dim = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberTz, this->quantumNumberY, this->quantumNumberR, this->quantumNumberS); 
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateLzMax = new int [this->HilbertSpaceDimension];  
      if (this->NbrFermions > 0)
      {
	this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrFluxQuanta, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0, 0l);
	int TmpLzMax = this->LzMax;
	for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	  {
	    while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul)
	      --TmpLzMax;
	    this->StateLzMax[i] = TmpLzMax;
	  }
      }
      else
      {
	this->StateDescription[0] = 0x0ul; 
	this->StateLzMax[0] = 0;
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


// constructor for truncated CP2 geometry
// 
// nbrBosons = number of bosons
// p = number of flux quanta (determines an irreducible representation of SU(3), along with q=0 (LLL))
// totalJz = total value of jz
// totalKz = total value of kz
// memory = amount of memory granted for precalculations

FermionOnCP2::FermionOnCP2 (int nbrFermions, int nbrFluxQuanta, int minY, int totalTz, int totalY, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->MinY = minY;
  this->TotalTz = totalTz;
  this->TotalY = totalY;
  this->TotalR = (this->TotalY + 3*this->TotalTz + 2*this->NbrFermions*this->NbrFluxQuanta)/6;
  this->TotalS = (this->TotalY - 3*this->TotalTz + 2*this->NbrFermions*this->NbrFluxQuanta)/6;
  this->NbrLzValue = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2 - (2*this->NbrFluxQuanta + this->MinY) * (2*this->NbrFluxQuanta + this->MinY + 3) / 18;
  this->LzMax = this->NbrLzValue - 1;  
  this->MaximumSignLookUp = 16;
  if (this->NbrFermions > 0)
    this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFluxQuanta, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0);
  else
    this->LargeHilbertSpaceDimension = 1l;
  this->quantumNumberTz = new int [this->NbrLzValue];
  this->quantumNumberY = new int [this->NbrLzValue];
  this->quantumNumberR = new int [this->NbrLzValue];
  this->quantumNumberS = new int [this->NbrLzValue];
  cout << "dim = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberTz, this->quantumNumberY, this->quantumNumberR, this->quantumNumberS); 
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateLzMax = new int [this->HilbertSpaceDimension];  
      if (this->NbrFermions > 0)
      {
	this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrFluxQuanta, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0, 0l);
	int TmpLzMax = this->LzMax;
	for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	  {
	    while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul)
	      --TmpLzMax;
	    this->StateLzMax[i] = TmpLzMax;
	  }
      }
      else
      {
	this->StateDescription[0] = 0x0ul; 
	this->StateLzMax[0] = 0;
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

FermionOnCP2::FermionOnCP2(const FermionOnCP2& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->NbrFluxQuanta = fermions.NbrFluxQuanta;
  this->MinY = fermions.MinY;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
  this->TotalR = fermions.TotalR;
  this->TotalS = fermions.TotalS;
  this->quantumNumberTz = fermions.quantumNumberTz;
  this->quantumNumberY = fermions.quantumNumberY;  
  this->quantumNumberR = fermions.quantumNumberR;
  this->quantumNumberS = fermions.quantumNumberS; 
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

FermionOnCP2::~FermionOnCP2 ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnCP2& FermionOnCP2::operator = (const FermionOnCP2& fermions)
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
  this->NbrFluxQuanta = fermions.NbrFluxQuanta;
  this->MinY = fermions.MinY;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
  this->TotalR = fermions.TotalR;
  this->TotalS = fermions.TotalS;
  this->quantumNumberTz = fermions.quantumNumberTz;
  this->quantumNumberY = fermions.quantumNumberY;  
  this->quantumNumberR = fermions.quantumNumberR;
  this->quantumNumberS = fermions.quantumNumberS; 
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

AbstractHilbertSpace* FermionOnCP2::Clone()
{
  return new FermionOnCP2(*this);
}


// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured
bool FermionOnCP2::WriteHilbertSpace (char* fileName)
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
  WriteLittleEndian(File, this->TotalTz);
  WriteLittleEndian(File, this->TotalY);
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

ostream& FermionOnCP2::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str <<  "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> i);
      int TmpTz = this->quantumNumberTz[i];
      int TmpY = this->quantumNumberY[i];
//       cout << TmpTz << ", " << TmpY << "; " << (Tmp & 0x1l) << endl;
      if ((Tmp & 0x1l) != 0ul)
	Str << "(" << TmpTz << "," << TmpY << ")";
    }
  Str << "]";
  return Str;
}


// generate all states corresponding to the constraints
// 
// stateDescription = array that gives each state description
// nbrBosons = number of bosons
// currentJ = current value of j for a single particle
// currentJz = current value of jz for a single particle
// currentKz = current value of kz for a single particle
// currentTotalTz = current total value of Jz
// currentTotalY = current total value of Kz
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnCP2::GenerateStates(int nbrFermions, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY, long pos)
{
  if (nbrFermions < 0)
    return pos;
  if (currentTotalTz + currentTzMax*nbrFermions < this->TotalTz)
    return pos;
  if (currentTotalY + currentY*nbrFermions < this->TotalY)
    return pos;
  
  if (currentTz < -currentTzMax)
   {
     --currentTzMax;
     currentTz = currentTzMax;
     currentY = currentY - 3;
   }
    
  
  if (nbrFermions == 0)
    {
      if ((currentTotalTz == this->TotalTz) && (currentTotalY == this->TotalY))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }

  if (currentY < this->MinY)
    return pos;
  
  long TmpPos = this->GenerateStates(nbrFermions - 1, currentTz - 2, currentTzMax, currentY, currentTotalTz + currentTz, currentTotalY + currentY, pos);
  unsigned long Mask = ((0x1ul << 1) - 0x1ul) << this->GetLinearizedIndex(currentTz, currentY, 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentTz - 2, currentTzMax, currentY, currentTotalTz, currentTotalY, pos);
};



// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentJz = current value of jz for a single particle
// currentKz = current value of kz for a single particle
// currentTotalTz = current total value of Jz
// currentTotalY = current total value of Kz
// return value = Hilbert space dimension

long FermionOnCP2::EvaluateHilbertSpaceDimension(int nbrBosons, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY)
{
  if (nbrBosons < 0)
    return 0l;
  if (currentTotalTz + currentTzMax*nbrBosons < this->TotalTz)
    return 0l;
  if (currentTotalY + currentY*nbrBosons < this->TotalY)
    return 0l;
  
  if (currentTz < -currentTzMax)
   {
     --currentTzMax;
     currentTz = currentTzMax;
     currentY = currentY - 3;
   }
    
  if (nbrBosons == 0)
    {
      if ((currentTotalTz == this->TotalTz) && (currentTotalY == this->TotalY))
      {
	return 1l;
      }
      else	
	return 0l;
    }
    
  if (currentY < this->MinY)
    return 0l;
  
  long Count = 0;
  
  Count += this->EvaluateHilbertSpaceDimension(nbrBosons - 1, currentTz - 2, currentTzMax, currentY, currentTotalTz + currentTz, currentTotalY + currentY);
  Count += this->EvaluateHilbertSpaceDimension(nbrBosons, currentTz - 2, currentTzMax, currentY, currentTotalTz, currentTotalY);
  return Count;
}

// find state index from an array
//
// stateDescription = array describing the state (stored as tz1,y1,tz2,y2,...)
// return value = corresponding index, -1 if an error occured

int FermionOnCP2::FindStateIndexFromArray(int* stateDescription)
{
  unsigned long TmpState = 0x0ul;
  for (int i = 0; i < this->NbrFermions; ++i)
    TmpState |= 0x1ul << (this->GetLinearizedIndex(stateDescription[i << 1], stateDescription[1 + (i << 1)], this->NbrFermions));
  int TmpLzMax = this->LzMax;
  while ((TmpState >> TmpLzMax) == 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals
bool FermionOnCP2::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  int* rootPartition1 = new int[2*this->NbrFermions];
//   int* rootPartition2 = new int[2*this->NbrFermions];
  int TmpLzMax = this->LzMax;
  while ((this->StateDescription[index] >> TmpLzMax) == 0x0ul)
    --TmpLzMax;
  int particleIndex = 0;
  for (int ind = 0; ind <= TmpLzMax; ++ind)
    {
      if (((this->StateDescription[index] >> ind) & 0x1ul) != 0)
      {
	    rootPartition1[particleIndex] = this->quantumNumberTz[ind];
	    rootPartition1[particleIndex + 1] = this->quantumNumberY[ind];
// 	    rootPartition2[particleIndex] = this->quantumNumberR[ind];
// 	    rootPartition2[particleIndex + 1] = this->quantumNumberS[ind];
// 	      cout << particleIndex << " " << rootPartition[particleIndex] << " " << rootPartition[particleIndex + 1] << endl;
	    particleIndex += 2;
      }
    }
  int i = 0;
  while (i < this->NbrFermions - 1)
    {
      int j = i + 1;
      while  (j < this->NbrFermions)
      {
       double distance1 = (3*(rootPartition1[j << 1] - rootPartition1[i << 1]) + 2*(rootPartition1[(j << 1) + 1] - rootPartition1[(i << 1) + 1]))/6.0;
//        cout <<rootPartition[i<<1] << " " << rootPartition[j << 1]  << " " << rootPartition[(i << 1) + 1] << " " << rootPartition[(j << 1) + 1] << " " << distance << " " ;
//        this->PrintState(cout, index) << endl;
	double distance2 = (rootPartition1[(j<<1) + 1] - rootPartition1[(i<<1) + 1]) / 3;
//        cout << distance ;
//        this->PrintState(cout, index) << " " << pauliR <<  endl;
// 	  int distance = abs(rootPartition[i << 1] - rootPartition[j << 1]) + abs(rootPartition[(i << 1) + 1] - rootPartition[(j << 1) + 1]) + abs (rootPartition[i << 1] - rootPartition[j << 1] + rootPartition[(i << 1) + 1] - rootPartition[(j << 1) + 1]);
//        cout << (i<<1) << " , " << (j<<1) << " : " << rootPartition[(i << 1)] << " " << rootPartition[j << 1] << " " << rootPartition[(i << 1) + 1] << " " << rootPartition[(j << 1) + 1] << " : " << distance << endl;
       if ((distance1 < (double) pauliR) && (distance2 < (double) pauliR))
	  return false;
       else
	 j += 1;
      }
     i += 1;
    }
  return true; 
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given (Jz,Kz) sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// jzSector = Jz sector in which the density matrix has to be evaluated 
// kzSector = Kz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  FermionOnCP2::EvaluatePartialDensityMatrixParticlePartition(int nbrFermionSector, int tzSector, int ySector, RealVector& groundState, AbstractArchitecture* architecture)
{  
  if (nbrFermionSector == 0)
    {
      if ((tzSector == 0) and (ySector ==0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  if (nbrFermionSector == this->NbrFermions)
    {
      if ((tzSector == this->TotalTz) and (ySector == this->TotalY))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ComplementaryNbrFermionSector = this->NbrFermions - nbrFermionSector;
  int ComplementaryRSector = (this->TotalY - ySector + 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrFermionSector*this->NbrFluxQuanta);
  int ComplementarySSector = (this->TotalY - ySector - 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrFermionSector*this->NbrFluxQuanta);
  if ((ComplementaryRSector < 0) || (ComplementarySSector < 0) || ((ComplementaryRSector % 6) != 0) || ((ComplementarySSector % 6) != 0) || (ComplementaryRSector + ComplementarySSector > 6*this->NbrFluxQuanta*ComplementaryNbrFermionSector))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
  cout << "nbr fermions = " << nbrFermionSector << ", tz = " << tzSector << ", y = " << ySector << endl;
  FermionOnCP2 TmpDestinationHilbertSpace(nbrFermionSector, this->NbrFluxQuanta, tzSector, ySector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  if (TmpDestinationHilbertSpace.HilbertSpaceDimension == 0)
  {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
   }
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  FermionOnCP2 TmpHilbertSpace(ComplementaryNbrFermionSector, this->NbrFluxQuanta, this->TotalTz - tzSector, this->TotalY - ySector);

  
  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &TmpDestinationHilbertSpace, &TmpHilbertSpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);

  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& FermionOnCP2::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  double* SqrtCoefficients = new double [this->NbrLzValue];
  FactorialCoefficient Coef;
  if (reference >= 0l)
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0 / state[reference];
      double* InvSqrtCoefficients = new double [this->NbrLzValue];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  int r = quantumNumberR[k];
	  int s = quantumNumberS[k];
	  int t = this->NbrFluxQuanta - r - s;
	  Coef.SetToOne();
	  Coef.FactorialDivide(r);
	  Coef.FactorialDivide(s);
	  Coef.FactorialDivide(t);
	  Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
	  SqrtCoefficients[k] = sqrt(Coef.GetNumericalValue());
	  InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
	}
      unsigned long TmpState = this->StateDescription[reference];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      state[reference] = 1.0;
      delete[] InvSqrtCoefficients;
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
   }
  else
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0;
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  int r = quantumNumberR[k];
	  int s = quantumNumberS[k];
	  int t = this->NbrFluxQuanta - r - s;
	  Coef.SetToOne();
	  Coef.FactorialDivide(r);
	  Coef.FactorialDivide(s);
	  Coef.FactorialDivide(t);
	  Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
	  SqrtCoefficients[k] = sqrt(Coef.GetNumericalValue());
	  InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
	}
      unsigned long TmpState = this->StateDescription[0l];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] InvSqrtCoefficients;
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
//       for (int k = 0; k <= this->LzMax; ++k)
// 	SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k) * ((double) (this->LzMax + 1)));
//       unsigned long TmpState;
//       int Index = 0;
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	{
// 	  Index = 0;
// 	  TmpState = this->StateDescription[i];
// 	  double Coefficient = 1.0;
// 	  for (int j = this->LzMax; j >= 0; --j)
// 	    if (((TmpState >> j) & 1ul) != 0ul)
// 	      Coefficient /= SqrtCoefficients[j];
// 	  state[i] *= Coefficient;
// 	}
   }
  delete[] SqrtCoefficients;
  return state;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given (tz, y) sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// tzSector = tz sector in which the density matrix has to be evaluated
// ySector = y sector in which the density matrix has to be evaluated
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix FermionOnCP2::EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int tzSector, int ySector, RealVector& groundState, bool removeBinomialCoefficient)
{
  if (nbrFermionSector == 0)
    {
      if ((tzSector == 0) && (ySector == 0))
	{
// 	  cout <<this->TotalTz << "  " << this->TotalY << endl;
// 	  cout << this->NbrFluxQuanta << endl;
	  FermionOnCP2 TmpHilbertSpace(this->NbrFermions, this->NbrFluxQuanta, this->TotalTz - tzSector, this->TotalY - ySector);
	  RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.GetHilbertSpaceDimension(), true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      int TmpLzMax = this->LzMax;
	      unsigned long TmpState = this->StateDescription[i];
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      
	      TmpEntanglementMatrix.SetMatrixElement(0, TmpHilbertSpace.FindStateIndex(TmpState, TmpLzMax), groundState[i]);
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	
	}
    }
  
  
  if (nbrFermionSector == this->NbrFermions)
    {
      if ((tzSector == this->TotalTz) && (ySector == this->TotalY))
	{
	  FermionOnCP2 TmpDestinationHilbertSpace(nbrFermionSector, this->NbrFluxQuanta, this->TotalTz, this->TotalY);
	  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.GetHilbertSpaceDimension(), 1,true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      int TmpLzMax = this->LzMax;
	      unsigned long TmpState = this->StateDescription[i];
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      
	      TmpEntanglementMatrix.SetMatrixElement(TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax), 0, groundState[i]);
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  int ComplementaryNbrFermionSector = this->NbrFermions - nbrFermionSector;
  int ComplementaryRSector = (this->TotalY - ySector + 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrFermionSector*this->NbrFluxQuanta);
  int ComplementarySSector = (this->TotalY - ySector - 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrFermionSector*this->NbrFluxQuanta);
  if ((ComplementaryRSector < 0) || (ComplementarySSector < 0) || ((ComplementaryRSector % 6) != 0) || ((ComplementarySSector % 6) != 0) || (ComplementaryRSector + ComplementarySSector > 6*this->NbrFluxQuanta*ComplementaryNbrFermionSector))
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
  
  FermionOnCP2 TmpDestinationHilbertSpace(nbrFermionSector, this->NbrFluxQuanta, tzSector, ySector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  long TmpNbrNonZeroElements = 0;

  
  
  FermionOnCP2 TmpHilbertSpace(ComplementaryNbrFermionSector, this->NbrFluxQuanta, this->TotalLz - tzSector, this->TotalY - ySector);
  
  FactorialCoefficient Factorial;
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);

  TmpNbrNonZeroElements = this->EvaluatePartialEntanglementMatrixParticlePartitionCore(0, TmpHilbertSpace.HilbertSpaceDimension, &TmpHilbertSpace,
										       &TmpDestinationHilbertSpace, 
										       groundState, &TmpEntanglementMatrix, removeBinomialCoefficient); 
  if (TmpNbrNonZeroElements > 0)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}



// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given Lz sector and both A and B are resticted to a given number of orbitals
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// tzSector = Tz sector in which the density matrix has to be evaluated 
// ySector = Y sector in which the density matrix has to be evaluated 
// maxYA = maximum value of Y in the A partition
// minYB = minimum value of Y in the B partition
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix FermionOnCP2::EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int tzSector, int ySector,
									       int maxYA, int minYB,
									       RealVector& groundState, bool removeBinomialCoefficient)
{
  int NbrOrbitalA = (maxYA + 1)*(maxYA + 2)/2;
  int NbrOrbitalB = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2 - (2*this->NbrFluxQuanta + minYB) * (2*this->NbrFluxQuanta + minYB + 3) / 18;
  
  int ShiftYA = 2*(this->NbrFluxQuanta - maxYA);
  
  cout << NbrOrbitalA << " " << NbrOrbitalB << endl;
  int ComplementaryNbrFermionSector = this->NbrFermions - nbrFermionSector;
  
  int RSector = (ySector + 3*tzSector + 2*nbrFermionSector*this->NbrFluxQuanta);
  int SSector = (ySector - 3*tzSector + 2*nbrFermionSector*this->NbrFluxQuanta);
  
  if ((RSector < 0) || (SSector < 0) || ((RSector % 6) != 0) || ((SSector % 6) != 0) || (RSector + SSector > 6*this->NbrFluxQuanta*nbrFermionSector))
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;  
    }
    
  int ComplementaryRSector = (this->TotalY - ySector + 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrFermionSector*this->NbrFluxQuanta);
  int ComplementarySSector = (this->TotalY - ySector - 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrFermionSector*this->NbrFluxQuanta);
  
  if ((ComplementaryRSector < 0) || (ComplementarySSector < 0) || ((ComplementaryRSector % 6) != 0) || ((ComplementarySSector % 6) != 0) || (ComplementaryRSector + ComplementarySSector > 6*this->NbrFluxQuanta*ComplementaryNbrFermionSector))
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;  
    }
  
  if (nbrFermionSector == 0)
    {
      if ((tzSector == 0) && (ySector == 0))
	{
// 	  cout << this->NbrFermions << " " << this->NbrFluxQuanta << " " << minYB << " " << this->TotalTz - tzSector << " " << this->TotalY - ySector << endl;
	  FermionOnCP2 TmpHilbertSpace(this->NbrFermions, this->NbrFluxQuanta, minYB, this->TotalTz - tzSector, this->TotalY - ySector);
	  RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      int TmpLzMax = this->LzMax ;
	      unsigned long TmpState = this->StateDescription[i];
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpIndex = TmpHilbertSpace.FindStateIndex(TmpState, TmpLzMax);
	      if (TmpIndex < TmpHilbertSpace.HilbertSpaceDimension)
		TmpEntanglementMatrix.SetMatrixElement(0, TmpIndex, groundState[i]);
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	
	}
    }
  
  
  if (nbrFermionSector == this->NbrFermions)
    {
      if ((tzSector == this->TotalLz) && (ySector == this->TotalY) )
	{
	  FermionOnCP2 TmpDestinationHilbertSpace(nbrFermionSector, maxYA ,tzSector, ySector + ShiftYA);
	  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1,true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      int TmpLzMax = this->LzMax;
	      unsigned long TmpState = this->StateDescription[i];
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      
	      int TmpIndex = TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax);
	      if (TmpIndex < TmpDestinationHilbertSpace.HilbertSpaceDimension)
		TmpEntanglementMatrix.SetMatrixElement(TmpIndex, 0, groundState[i]);
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  FermionOnCP2 TmpDestinationHilbertSpace(nbrFermionSector, maxYA, tzSector, ySector + ShiftYA);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  long TmpNbrNonZeroElements = 0;

  double TmpInvBinomial = 1.0;
  if(removeBinomialCoefficient == false )
    {
      BinomialCoefficients TmpBinomial (this->NbrFermions);
      TmpInvBinomial = sqrt(1.0 / (TmpBinomial(this->NbrFermions, nbrFermionSector)));
    }
  
  
  FermionOnCP2 TmpHilbertSpace(ComplementaryNbrFermionSector, this->NbrFluxQuanta, minYB, this->TotalTz - tzSector, this->TotalY - ySector);
  
  FactorialCoefficient Factorial;
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
  
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << (this->NbrLzValue - NbrOrbitalB);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace.StateDescription[j];
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
		  int Pos2 = TmpDestinationHilbertSpace.LzMax;
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
		  
		  TmpNbrNonZeroElements++;
		  TmpEntanglementMatrix.SetMatrixElement(j, MinIndex, Coefficient*groundState[TmpPos]);
		}
	    }
	}
    }
  
  if (TmpNbrNonZeroElements > 0)
    return TmpEntanglementMatrix;
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}
// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix. The cut has to be at a given value of Y.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// tzSector = Tz sector in which the density matrix has to be evaluated 
// ySector = Y sector in which the density matrix has to be evaluated 
// maxYA = maximum value of Y in the A partition
// weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
// minYB = minimum value of Y that has to be kept in the B partition
// weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& FermionOnCP2::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int tzSector, int ySector, int maxYA, double* weightOrbitalA, 
														int minYB, double* weightOrbitalB, RealMatrix& entanglementMatrix)
{  
  int ShiftYA = 2*(this->NbrFluxQuanta - maxYA);
  
  cout << "Evaluate entanglement" << endl;
  int RSector = (ySector + 3*tzSector + 2*nbrFermionSector*this->NbrFluxQuanta);
  int SSector = (ySector - 3*tzSector + 2*nbrFermionSector*this->NbrFluxQuanta);
  
  if ((RSector < 0) || (SSector < 0) || ((RSector % 6) != 0) || ((SSector % 6) != 0) || (RSector + SSector > 6*this->NbrFluxQuanta*nbrFermionSector))
    {
      return entanglementMatrix;	  
    }
    
    
  int ComplementaryNbrFermionsSector = this->NbrFermions - nbrFermionSector;
  int ComplementaryRSector = (this->TotalY - ySector + 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrFermionsSector*this->NbrFluxQuanta);
  int ComplementarySSector = (this->TotalY - ySector - 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrFermionsSector*this->NbrFluxQuanta);
  
  if ((ComplementaryRSector < 0) || (ComplementarySSector < 0) || ((ComplementaryRSector % 6) != 0) || ((ComplementarySSector % 6) != 0) || (ComplementaryRSector + ComplementarySSector > 6*this->NbrFluxQuanta*ComplementaryNbrFermionsSector))
    {
      return entanglementMatrix;	  
    }
  
  
  FermionOnCP2 TmpDestinationHilbertSpace(nbrFermionSector, maxYA, tzSector, ySector + ShiftYA);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrFermionsSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrFermions];
  
  FermionOnCP2 TmpHilbertSpace(ComplementaryNbrFermionsSector, this->NbrFluxQuanta, minYB, this->TotalTz - tzSector, this->TotalY - ySector);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.StateDescription[i],TmpMonomial3);
      double Tmp = 1.0;
      for (int j = 0; j < nbrFermionSector; j++)
	{
	  Tmp *= weightOrbitalA[TmpMonomial3[j]];
	}
      for (int j = 0; j < TmpHilbertSpace.HilbertSpaceDimension; ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpMonomial1);
      double FormFactor = 1.0;
      for (int i = 0; i < ComplementaryNbrFermionsSector; i++)
	FormFactor *= weightOrbitalB[TmpMonomial1[i]];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
  
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;
  
  return entanglementMatrix;
}

