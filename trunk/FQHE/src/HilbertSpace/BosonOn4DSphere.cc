////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Cecile Repellin                 //
//                                                                            //
//                                                                            //
//               class of Hilbert space for bosons on 4D sphere               //
//                                                                            //
//                                                                            //
//                        last modification : 19/10/2012                      //
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
#include "HilbertSpace/BosonOn4DSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "Architecture/ArchitectureOperation/FQHESquareLatticeSymmetrizeU1U1StateOperation.h"

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

BosonOn4DSphere::BosonOn4DSphere ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// p = number of flux quanta (determines an irreducible representation of SO(5), along with q=0 (LLL))
// totalJz = total value of jz
// totalKz = total value of kz
// memory = amount of memory granted for precalculations

BosonOn4DSphere::BosonOn4DSphere (int nbrBosons, int nbrFluxQuanta, int totalJz, int totalKz, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->TotalJz = totalJz;
  this->ShiftedTotalJz = totalJz + this->NbrBosons*this->NbrFluxQuanta;
  this->TotalKz = totalKz;
  this->ShiftedTotalKz = totalKz + this->NbrBosons*this->NbrFluxQuanta;
  this->NbrLzValue = - (this->NbrFluxQuanta*(this->NbrFluxQuanta+1)*(2*this->NbrFluxQuanta+1))/6 + this->NbrFluxQuanta*(this->NbrFluxQuanta*(this->NbrFluxQuanta+1))/2 + (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 1);
  this->LzMax = NbrLzValue - 1;  
  this->Minors = 0;
  this->KeptCoordinates = 0;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0, 0);
  this->quantumNumberJ = new int [this->NbrLzValue];
  this->quantumNumberJz = new int [this->NbrLzValue];
  this->quantumNumberKz = new int [this->NbrLzValue];
  cout << "dim = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(TmpStateDescription, this->NbrBosons, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0, 0, this->NbrLzValue + this->NbrBosons, 0l);
      this->GetQuantumNumbersFromLinearizedIndex(quantumNumberJ, quantumNumberJz, quantumNumberKz);
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space : get " << TmpLargeHilbertSpaceDimension << " , should be " << this->LargeHilbertSpaceDimension << endl;
	}
      this->FermionBasis = new FermionOnSphere(this->NbrBosons, 0, this->LzMax + this->NbrBosons, TmpStateDescription, this->LargeHilbertSpaceDimension);
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
      UsedMemory += this->NbrLzValue * this->FermionBasis-> LookUpTableMemorySize * sizeof(int);
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
// bosons = reference on the hilbert space to copy to copy

BosonOn4DSphere::BosonOn4DSphere(const BosonOn4DSphere& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalJz = bosons.TotalJz;
  this->ShiftedTotalJz = bosons.TotalJz;
  this->TotalKz = bosons.TotalKz;
  this->ShiftedTotalKz = bosons.TotalKz;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->NbrLzValue = bosons.NbrLzValue;
  this->LzMax = bosons.LzMax;
  this->Minors = 0;
  this->KeptCoordinates = 0;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->quantumNumberJ = new int [this->NbrLzValue];
  this->quantumNumberJz = new int [this->NbrLzValue];
  this->quantumNumberKz = new int [this->NbrLzValue];
}

// destructor
//

BosonOn4DSphere::~BosonOn4DSphere ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOn4DSphere& BosonOn4DSphere::operator = (const BosonOn4DSphere& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->NbrLzValue = bosons.NbrLzValue;
  this->Minors = 0;
  this->KeptCoordinates = 0;
  this->TotalJz = bosons.TotalJz;
  this->ShiftedTotalJz = bosons.TotalJz;
  this->TotalKz = bosons.TotalKz;
  this->ShiftedTotalKz = bosons.TotalKz;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->quantumNumberJ = new int [this->NbrLzValue];
  this->quantumNumberJz = new int [this->NbrLzValue];
  this->quantumNumberKz = new int [this->NbrLzValue];
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOn4DSphere::Clone()
{
  return new BosonOn4DSphere(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured
bool BosonOn4DSphere::WriteHilbertSpace (char* fileName)
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
  WriteLittleEndian(File, this->NbrBosons);
  WriteLittleEndian(File, this->TotalJz);
  WriteLittleEndian(File, this->TotalKz);
  if (this->HilbertSpaceDimension != 0)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateDescription[i]);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateLzMax[i]);
    }
  else
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateDescription[i]);
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateLzMax[i]);
    }
  File.close();
  return true;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOn4DSphere::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], this->TemporaryState, this->TemporaryStateLzMax);
  for (int i = 0; i < this->NbrLzValue ; ++i)
  {
    if (i <= this->TemporaryStateLzMax)
      cout << this->TemporaryState[i];
    else
      cout << 0;
  }
  cout << " " ;
//   Str << hex << this->FermionBasis->StateDescription[state] << dec << " " <<  this->FermionBasis->StateLzMax[state] << "   ";
  Str << this->FermionBasis->StateDescription[state] <<"[";
  for (int index = 0; index <= this->TemporaryStateLzMax; ++index)
    {
      if (this->TemporaryState[index] > 0)
	{
	  for (int i = 0; i < this->TemporaryState[index]; ++i)
	    {
	      Str << "(" << this->quantumNumberJ[index] << "," << 2*this->quantumNumberJz[index] - this->quantumNumberJ[index]  << "," << 2*this->quantumNumberKz[index] + this->quantumNumberJ[index] - this->NbrFluxQuanta << ")";
	    }
	}
    }
  Str << "]";
  
//    Str << "  "  << this->FermionBasis->FindStateIndex(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state]) << " = " << state << endl;
  return Str;
}

// generate all states corresponding to the constraints
// 
// stateDescription = array that gives each state description
// nbrBosons = number of bosons
// currentJ = current value of j for a single particle
// currentJz = current value of jz for a single particle
// currentKz = current value of kz for a single particle
// currentTotalJz = current total value of Jz
// currentTotalKz = current total value of Kz
// currentFermionicPosition = current fermionic position within the state description
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
  
long BosonOn4DSphere::GenerateStates(unsigned long* stateDescription, int nbrBosons, int currentJ, int currentJz, int currentKz, int currentTotalJz, int currentTotalKz, int currentFermionicPosition, long pos)
{

  if (nbrBosons < 0)
    return pos;
  if (currentTotalKz > this->ShiftedTotalKz)
    return pos;
  if (currentTotalJz > this->ShiftedTotalJz)
    return pos;
  if (currentTotalJz + nbrBosons*(currentJ + this->NbrFluxQuanta) < this->ShiftedTotalJz)
    return pos;
  if (currentTotalKz + nbrBosons*(2*this->NbrFluxQuanta + currentJ) < this->ShiftedTotalKz)
    return pos;
  
  if (currentKz < 0)
   {
     currentKz = this->NbrFluxQuanta - currentJ;
     currentJz--;
   }
    
  if (currentJz < 0)
  {
   currentJ--;
   currentJz = currentJ;
   currentKz = this->NbrFluxQuanta - currentJ;
  }
  if (nbrBosons == 0)
    {
      if ((currentTotalJz == this->ShiftedTotalJz) && (currentTotalKz == this->ShiftedTotalKz))
	{
	  stateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }

  if (currentJ < 0)
    return pos;

  for (int k = nbrBosons; k > 0; --k)
    {
      long TmpPos = this->GenerateStates(stateDescription, nbrBosons - k, currentJ, currentJz, currentKz - 1, currentTotalJz + (k * (2*currentJz - currentJ + this->NbrFluxQuanta)), currentTotalKz + (k * (2*currentKz + currentJ)), currentFermionicPosition - k - 1, pos);
      unsigned long Mask = ((0x1ul << k) - 0x1ul) << (currentFermionicPosition - k - 1);
      for (; pos < TmpPos; ++pos)
	stateDescription[pos] |= Mask;
    }
  return this->GenerateStates(stateDescription, nbrBosons, currentJ, currentJz, currentKz - 1, currentTotalJz, currentTotalKz, currentFermionicPosition - 1, pos);
};


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentJ = current value of j for a single particle
// currentJz = current value of jz for a single particle
// currentKz = current value of kz for a single particle
// currentTotalJz = current total value of Jz
// currentTotalKz = current total value of Kz
// return value = Hilbert space dimension

long BosonOn4DSphere::EvaluateHilbertSpaceDimension(int nbrBosons, int currentJ, int currentJz, int currentKz, int currentTotalJz, int currentTotalKz)
{
  if (nbrBosons < 0)
    return 0l;
  if (currentTotalKz > this->ShiftedTotalKz)
    return 0l;
  if (currentTotalJz > this->ShiftedTotalJz)
    return 0l;
  if (currentTotalJz + nbrBosons*(currentJ + this->NbrFluxQuanta) < this->ShiftedTotalJz)
    return 0l;
  if (currentTotalKz + nbrBosons*(2*this->NbrFluxQuanta + currentJ) < this->ShiftedTotalKz)
    return 0l;
  
  if (currentKz < 0)
   {
     currentKz = this->NbrFluxQuanta - currentJ;
     currentJz--;
   }
    
  if (currentJz < 0)
  {
   currentJ--;
   currentJz = currentJ;
   currentKz = this->NbrFluxQuanta - currentJ;
  }
  
  if (nbrBosons == 0)
    {
      if ((currentTotalJz == this->ShiftedTotalJz) && (currentTotalKz == this->ShiftedTotalKz))
      {
	return 1l;
      }
      else	
	return 0l;
    }
    
  if (currentJ < 0)
    return 0l;
  
  long Count = 0;
  for (int k = nbrBosons; k >= 0; --k)
    Count += this->EvaluateHilbertSpaceDimension(nbrBosons - k, currentJ, currentJz, currentKz - 1, currentTotalJz + (k * (2*currentJz - currentJ + this->NbrFluxQuanta)), currentTotalKz + (k * (2*currentKz + currentJ)));
  return Count;
}



// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given (Jz,Kz) sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// jzSector = Jz sector in which the density matrix has to be evaluated 
// kzSector = Kz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOn4DSphere::EvaluatePartialDensityMatrixParticlePartition(int nbrBosonSector, int jzSector, int kzSector, RealVector& groundState, AbstractArchitecture* architecture)
{  
  if (nbrBosonSector == 0)
    {
      if ((jzSector == 0) and (kzSector ==0))
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

  if (nbrBosonSector == this->NbrBosons)
    {
      if ((jzSector == this->TotalJz) and (kzSector == this->TotalKz))
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

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;

  if ((abs(this->TotalJz - jzSector) + abs(this->TotalKz - kzSector)) > (ComplementaryNbrBosonSector * this->NbrFluxQuanta))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
  cout << "nbr boson = " << nbrBosonSector << ", jz = " << jzSector << ", kz = " << kzSector << endl;
  BosonOn4DSphere TmpDestinationHilbertSpace(nbrBosonSector, this->NbrFluxQuanta, jzSector, kzSector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  BosonOn4DSphere TmpHilbertSpace(ComplementaryNbrBosonSector, this->NbrFluxQuanta, this->TotalJz - jzSector, this->TotalKz - kzSector);

  
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

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to be stored
// return value = number of components that have been added to the density matrix

long BosonOn4DSphere::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace, RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
   BosonOn4DSphere* TmpHilbertSpace =  (BosonOn4DSphere*) complementaryHilbertSpace;
   BosonOn4DSphere* TmpDestinationHilbertSpace =  (BosonOn4DSphere*) destinationHilbertSpace;
   int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
   int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
   unsigned long* TmpMonomial2 = new unsigned long [NbrBosonSector];
   unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
   unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
   int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int MaxIndex = minIndex + nbrIndex;
   long TmpNbrNonZeroElements = 0l;

   double* LogFactorials = new double[this->NbrBosons + 1];
   LogFactorials[0] = 0.0;
   LogFactorials[1] = 0.0;
   for (int i = 2 ; i <= this->NbrBosons; ++i)
     LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
   double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];

   for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
	TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
   for (; minIndex < MaxIndex; ++minIndex)    
     {
      int Pos = 0;
      TmpHilbertSpace->ConvertToMonomial(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpMonomial1);
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->FermionBasis->StateDescription[j], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[j], TmpMonomial2);
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < NbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < NbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < NbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 
	      
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      TmpStateCoefficient[Pos] = exp(TmpFactorial);
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		    ++TmpNbrNonZeroElements;
		  }
	    }
	}
     }
   delete[] TmpMonomial2;
   delete[] TmpMonomial1;
   delete[] TmpMonomial3;
   delete[] TmpStatePosition;
   delete[] TmpStatePosition2;
   delete[] TmpStateCoefficient;
   delete[] TmpDestinationLogFactorials;
   return TmpNbrNonZeroElements;
}
