////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of fermions on lattice in real space              //
//                                                                            //
//                        last modification : 09/09/2014                      //
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
#include "HilbertSpace/FermionOnLatticeRealSpace.h"
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
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>
#include <bitset>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnLatticeRealSpace::FermionOnLatticeRealSpace ()
{
  this->NbrFermions = 0;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrSite = 0;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateLzMax = 0;  
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// memory = amount of memory granted for precalculations

FermionOnLatticeRealSpace::FermionOnLatticeRealSpace (int nbrFermions, int nbrSite, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
      this->StateLzMax = new int [this->LargeHilbertSpaceDimension];  
      int CurrentLzMax = this->NbrLzValue;
      while ((((this->StateDescription[0] >> CurrentLzMax) & 0x1ul) == 0x0ul)&&(CurrentLzMax>=0))
	--CurrentLzMax;
      this->StateLzMax[0] = CurrentLzMax;
      for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while (((this->StateDescription[i] >> CurrentLzMax) & 0x1ul) == 0x0ul)
	    --CurrentLzMax;
	  this->StateLzMax[i] = CurrentLzMax;
 	}
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
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

FermionOnLatticeRealSpace::FermionOnLatticeRealSpace(const FermionOnLatticeRealSpace& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSite = fermions.NbrSite;
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

FermionOnLatticeRealSpace::~FermionOnLatticeRealSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeRealSpace& FermionOnLatticeRealSpace::operator = (const FermionOnLatticeRealSpace& fermions)
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
  this->NbrSite = fermions.NbrSite;
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

AbstractHilbertSpace* FermionOnLatticeRealSpace::Clone()
{
  return new FermionOnLatticeRealSpace(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeRealSpace::GenerateStates(int nbrFermions, int currentSite, long pos)
{
  if (nbrFermions == 0)
    {
      this->StateDescription[pos] = 0x0ul;	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentSite; j >= 0; --j)
	{
	  this->StateDescription[pos] = 0x1ul << j;
	  ++pos;
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, pos);
  unsigned long Mask = 0x1ul << currentSite;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentSite - 1, pos);
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// return value = Hilbert space dimension
long FermionOnLatticeRealSpace::EvaluateHilbertSpaceDimension(int nbrFermions)
{
  BinomialCoefficients binomials(this->NbrSite);
  long dimension = binomials(this->NbrSite, nbrFermions);
  return dimension;
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = kx sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeRealSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  FermionOnLatticeRealSpace SubsytemSpace (nbrParticleSector, this->NbrSite);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite);
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

// evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
// nbrKeptOrbitals = array of orbitals that have to be kept
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem

ComplexMatrix FermionOnLatticeRealSpace::EvaluatePartialEntanglementMatrix (int nbrParticleSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  if ((nbrParticleSector >  nbrKeptOrbitals) ||  (ComplementaryNbrParticles > (this->NbrSite - nbrKeptOrbitals) ))
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
  if (nbrKeptOrbitals == 0)
    {
      if (nbrParticleSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, this->HilbertSpaceDimension, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      TmpEntanglementMatrix[i][0] = groundState[i];
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  
  this->KeptOrbitals = new int [nbrKeptOrbitals];
  for (int i = 0 ; i < nbrKeptOrbitals; ++i) 
    {
      this->KeptOrbitals[i] = keptOrbitals[i];
    }
  
  if (nbrParticleSector == 0)
    {
      FermionOnLatticeRealSpace ComplementarySpace (this->NbrFermions, this->NbrSite - nbrKeptOrbitals);
      ComplexMatrix TmpEntanglementMatrix(1, ComplementarySpace.HilbertSpaceDimension, true);
      long TmpEntanglementMatrixZero = this->EvaluatePartialEntanglementMatrixZeroParticuleCase (0, ComplementarySpace.HilbertSpaceDimension, &ComplementarySpace, groundState, &TmpEntanglementMatrix);
      if (TmpEntanglementMatrixZero > 0)
	{
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrixZero;
	  return TmpEntanglementMatrixZero;
	}    
    }

  if (ComplementaryNbrParticles == 0)
    {
      FermionOnLatticeRealSpace SubsytemSpace (this->NbrFermions, nbrKeptOrbitals);
      ComplexMatrix TmpEntanglementMatrix(SubsytemSpace.HilbertSpaceDimension, 1, true);
      long TmpEntanglementMatrixZero = this->EvaluatePartialEntanglementMatrixMaxParticuleCase (&SubsytemSpace, groundState, &TmpEntanglementMatrix); 
      if (TmpEntanglementMatrixZero > 0)
	{
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrixZero;
	  return TmpEntanglementMatrixZero;
	}    
    }

  if (nbrKeptOrbitals == this->NbrSite)
    {
      if (nbrParticleSector == this->NbrFermions)
	{
	  ComplexMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      TmpEntanglementMatrix[0][i] = groundState[i];
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }

  FermionOnLatticeRealSpace SubsytemSpace (nbrParticleSector, nbrKeptOrbitals);
  FermionOnLatticeRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite - nbrKeptOrbitals);
  ComplexMatrix TmpEntanglementMatrix (SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySpace.HilbertSpaceDimension, true);
  
  long TmpEntanglementMatrixZero = this->EvaluatePartialEntanglementMatrixCore(0, ComplementarySpace.HilbertSpaceDimension, &ComplementarySpace, &SubsytemSpace, groundState, &TmpEntanglementMatrix);
  if (TmpEntanglementMatrixZero > 0)
    {
      return TmpEntanglementMatrix;
    }
   else
     {
       ComplexMatrix TmpEntanglementMatrixZero;
       return TmpEntanglementMatrixZero;
     }    
}
  
// core part of the evaluation orbital cut entanglement matrix calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnLatticeRealSpace::EvaluatePartialEntanglementMatrixCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace, ComplexVector& groundState, ComplexMatrix* entanglementMatrix)
{
//cout <<"Entering core"<<endl;
  FermionOnLatticeRealSpace* TmpHilbertSpace = (FermionOnLatticeRealSpace*) complementaryHilbertSpace;
  FermionOnLatticeRealSpace* TmpDestinationHilbertSpace = (FermionOnLatticeRealSpace*) destinationHilbertSpace;
  long TmpNbrNonZeroElements = 0;
  int* TraceOutOrbitals = new int [this->LzMax - TmpDestinationHilbertSpace->LzMax];
  int MaxIndex = minIndex + nbrIndex;
  int TmpIndex = 0;
  for (int i = 0; i < this->LzMax; ++i)
    {
      if (SearchInArray<int>(i, this->KeptOrbitals, TmpDestinationHilbertSpace->LzMax) < 0)
	TraceOutOrbitals[TmpIndex++] = i;
    }
/*  for(int i = 0; i <this->LzMax - TmpDestinationHilbertSpace->LzMax;i ++)
    cout <<TraceOutOrbitals[i]<<" ";
  cout <<endl;*/
  for (; minIndex < MaxIndex; ++minIndex) 
    {
      int Pos = 0;
      unsigned long TmpStateCompact = TmpHilbertSpace->StateDescription[minIndex];
      unsigned long TmpState = 0x0ul;
//      cout <<"TmpStateCompact = "<<std::bitset<16>(TmpStateCompact)<<endl;
      for (int i = 0 ; i < TmpHilbertSpace->LzMax; ++i)
	TmpState |= (TmpStateCompact >> i) << TraceOutOrbitals[i];
//      cout <<"TmpState = "<<std::bitset<16>(TmpState)<<endl;
      
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpStateCompact2 = TmpDestinationHilbertSpace->StateDescription[j];
	  unsigned long TmpState2 = 0x0ul;
	  for (int i = 0 ; i < TmpDestinationHilbertSpace->LzMax; ++i)
	    TmpState2 |= ((TmpStateCompact2 >> i) & 0x1ul) << this->KeptOrbitals[i];
	  
	//  cout <<"TmpState2 = "<<std::bitset<16>(TmpState2)<<endl;
      
	  unsigned long TmpState3 = TmpState | TmpState2;

	  int TmpLzMax = this->LzMax; 
	  while ((TmpState3 >> TmpLzMax) == 0x0ul)
	    --TmpLzMax; 
//	  cout <<std::bitset<16>(TmpState3)<<" "<<dec<<TmpLzMax<<endl;
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
	      
	      entanglementMatrix->AddToMatrixElement(j, minIndex, Coefficient * groundState[TmpPos]);
	      ++TmpNbrNonZeroElements;
	    }
	}
    }
  delete[] TraceOutOrbitals;
  return TmpNbrNonZeroElements;
}


long FermionOnLatticeRealSpace::EvaluatePartialEntanglementMatrixZeroParticuleCase (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace, ComplexVector& groundState, ComplexMatrix* entanglementMatrix)
{
  FermionOnLatticeRealSpace* TmpHilbertSpace = (FermionOnLatticeRealSpace*) complementaryHilbertSpace;
  long TmpNbrNonZeroElements = 0;
  int* TraceOutOrbitals = new int [TmpHilbertSpace->LzMax];
  int MaxIndex = minIndex + nbrIndex;
  int TmpIndex = 0;
  for (int i = 0; i < this->LzMax; ++i)
    {
      if (SearchInArray<int>(i, this->KeptOrbitals, this->LzMax - TmpHilbertSpace->LzMax) < 0)
	TraceOutOrbitals[TmpIndex++] = i;
    }
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpStateCompact = TmpHilbertSpace->StateDescription[minIndex];
      unsigned long TmpState = 0x0ul;
      for (int i = 0 ; i < TmpHilbertSpace->LzMax; ++i)
	TmpState |= (TmpStateCompact >> i) << TraceOutOrbitals[i];
      
      int TmpLzMax = this->LzMax; 
      while ((TmpState >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
      if (TmpPos != this->HilbertSpaceDimension)
	{
	  entanglementMatrix->AddToMatrixElement(0, minIndex,  groundState[TmpPos]);
	  ++TmpNbrNonZeroElements;
	}
    }
  delete[] TraceOutOrbitals;
  return TmpNbrNonZeroElements;
}


// core part of the evaluation orbital cut entanglement matrix calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnLatticeRealSpace::EvaluatePartialEntanglementMatrixMaxParticuleCase (ParticleOnSphere* destinationHilbertSpace, ComplexVector& groundState, ComplexMatrix* entanglementMatrix)
{
  FermionOnLatticeRealSpace* TmpDestinationHilbertSpace = (FermionOnLatticeRealSpace*) destinationHilbertSpace;
  long TmpNbrNonZeroElements = 0;
  int* TraceOutOrbitals = new int [this->LzMax - TmpDestinationHilbertSpace->LzMax];
  int TmpIndex = 0;
  for (int i = 0; i < this->LzMax; ++i)
    {
      if (SearchInArray<int>(i, this->KeptOrbitals, TmpDestinationHilbertSpace->LzMax) < 0)
	TraceOutOrbitals[TmpIndex++] = i;
    }
  
  for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
    {
      unsigned long TmpStateCompact2 = TmpDestinationHilbertSpace->StateDescription[j];
      unsigned long TmpState2 = 0x0ul;
      for (int i = 0 ; i < TmpDestinationHilbertSpace->LzMax; ++i)
	TmpState2 |= ((TmpStateCompact2 >> i) & 0x1ul) << this->KeptOrbitals[i];
	  
      int TmpLzMax = this->LzMax; 
      while ((TmpState2 >> TmpLzMax) == 0x0ul)
	--TmpLzMax; 
      int TmpPos = this->FindStateIndex(TmpState2, TmpLzMax);
      if (TmpPos != this->HilbertSpaceDimension)
	{
	  entanglementMatrix->AddToMatrixElement(j, 0, groundState[TmpPos]);
	  ++TmpNbrNonZeroElements;
	}
    }
  delete[] TraceOutOrbitals;
  return TmpNbrNonZeroElements;
}
