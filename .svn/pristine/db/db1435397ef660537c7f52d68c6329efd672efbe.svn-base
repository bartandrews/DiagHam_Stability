////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2014 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on the 4D manifold T2 x S2               //
// forbidding mutliple orbital occupations and nearest orbital occupations    //
//                       and using magnetic translations                      //
//                                                                            //
//                        last modification : 06/03/2017                      //
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
#include "HilbertSpace/BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "MathTools/IntegerAlgebraTools.h"

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

BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// nbrFluxQuantumTorus = number of flux quanta piercing the torus
// kxMomentum = momentum in the x direction for the torus
// kyMomentum = momentum in the y direction for the torus
// nbrFluxQuantumSphere = number of flux quanta piercing the sphere
// totalLz = projection of the total angular momentum along the z axis for the sphere
// memory = amount of memory granted for precalculations

BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations (int nbrBosons, int nbrFluxQuantumTorus, 
																int kxMomentum, int kyMomentum,
																int nbrFluxQuantumSphere, int totalLz, unsigned long memory)
{  
  this->NbrFermions = nbrBosons;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->NbrFluxQuantumTorus = nbrFluxQuantumTorus;
  this->NbrFluxQuantumSphere = nbrFluxQuantumSphere;
  this->NbrOrbitalSphere = this->NbrFluxQuantumSphere + 1;
  this->ShiftedTotalLz = (this->TotalLz + (this->NbrFermions * this->NbrFluxQuantumSphere)) / 2;
  this->XMomentum = kxMomentum;
  this->YMomentum = kyMomentum;
  this->MaxMomentum = this->NbrFluxQuantumTorus * this->NbrOrbitalSphere;
  this->MomentumModulo = FindGCD(this->NbrFermions, this->NbrFluxQuantumTorus);

  this->NbrMomentum = this->MaxMomentum + 1;
  this->StateShift = this->MaxMomentum / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;
#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->MaxMomentum + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->MaxMomentum + 1) >> 1);
#endif
  if ((this->MaxMomentum & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;
  this->MaximumSignLookUp = 16;
  this->GenerateSignLookUpTable();

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFluxQuantumTorus - 1, this->NbrOrbitalSphere -1 , 0, 0);
  cout << "dim = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrFluxQuantumTorus - 1, this->NbrOrbitalSphere - 1, 0, 0, this->MaxMomentum - 1, 0l);
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space : get " << TmpLargeHilbertSpaceDimension << " , should be " << this->LargeHilbertSpaceDimension << endl;
	}
      this->InvertShift = 0;
      this->InvertUnshift = 0;
      if (this->NbrFermions == 0)
	{
	  this->StateMaxMomentum = new int [this->LargeHilbertSpaceDimension];
	  this->StateDescription[0] = 0x0ul; 
	  this->StateMaxMomentum[0] = 0;
	}
      else
	{
	  if (this->NbrFluxQuantumTorus > 1)
	    {
	      unsigned long TmpMask = (0x1ul << this->NbrOrbitalSphere) - 0x1ul;
	      unsigned long TmpMask2 = (0x1ul << this->NbrOrbitalSphere) | 0x1ul;
	      int TmpShift = (this->NbrFluxQuantumTorus - 1) * this->NbrOrbitalSphere;
	      TmpLargeHilbertSpaceDimension = 0l;
	      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
		{
		  if (((this->StateDescription[i] & (this->StateDescription[i] >> TmpShift)) & TmpMask) != 0x0ul)
		    {
		      this->StateDescription[i] = 0x0ul;
		    }
		  else
		    {		      
		      int TmpIndex = 0;
		      while ((TmpIndex < this->MaxMomentum) && (((this->StateDescription[i] >> TmpIndex) & TmpMask2) != TmpMask2))
			{
			  ++TmpIndex;		  
			}
		      if (TmpIndex < this->MaxMomentum)
			{
			  this->StateDescription[i] = 0x0ul;
			}
		      else
			{
			  ++TmpLargeHilbertSpaceDimension;
			}
		    }
		}
	    }
	  else
	    {
	      TmpLargeHilbertSpaceDimension = this->LargeHilbertSpaceDimension;
	    }
	  if (TmpLargeHilbertSpaceDimension > 0l)
	    {
	      unsigned long* TmpStateDescriptions = new unsigned long[TmpLargeHilbertSpaceDimension];
	      TmpLargeHilbertSpaceDimension = 0l;
	      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
		{
		  if (this->StateDescription[i] != 0x0ul)
		    {
		      TmpStateDescriptions[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
		      ++TmpLargeHilbertSpaceDimension;
		    }
		}
	      delete[] this->StateDescription;
	      this->StateDescription = TmpStateDescriptions;
	      this->LargeHilbertSpaceDimension = TmpLargeHilbertSpaceDimension;	      
	      if (this->LargeHilbertSpaceDimension >= (1l << 30))
		this->HilbertSpaceDimension = 0;
	      else
		this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
	      this->LargeHilbertSpaceDimension  = this->GenerateStates();
	      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
	      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
	      if (this->LargeHilbertSpaceDimension > 0l)
		{
		  this->StateMaxMomentum = new int [this->LargeHilbertSpaceDimension];
		  int TmpLzMax = this->MaxMomentum - 1;
		  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
		    {
		      while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul)
			--TmpLzMax;
		      this->StateMaxMomentum[i] = TmpLzMax;
		    }
		}
	    }
	  else
	    {
	      this->LargeHilbertSpaceDimension = 0l;
	      this->HilbertSpaceDimension = 0;
	    }
	}
//       for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
// 	{
// 	  this->PrintState(cout, (int) i) << endl;
// 	}
      this->MaximumSignLookUp = 16;
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
      UsedMemory = this->MaxMomentum * sizeof(int);
//      UsedMemory += this->MaxMomentum * this->FermionBasis->LookUpTableMemorySize * sizeof(int);
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

BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations(const BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrFermions = bosons.NbrFermions;
  this->IncNbrFermions = bosons.IncNbrFermions;
  this->TotalLz = bosons.TotalLz; 
  this->NbrFluxQuantumTorus = bosons.NbrFluxQuantumTorus;
  this->NbrFluxQuantumSphere = bosons.NbrFluxQuantumSphere;
  this->ShiftedTotalLz = bosons.ShiftedTotalLz;
  this->NbrFluxQuantumTorus = bosons.NbrFluxQuantumTorus;
  this->NbrOrbitalSphere = bosons.NbrOrbitalSphere;
  this->XMomentum = bosons.XMomentum;
  this->YMomentum = bosons.YMomentum;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->MomentumModulo = bosons.MomentumModulo;
  this->InvertShift = bosons.InvertShift;
  this->InvertUnshift = bosons.InvertUnshift;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
  this->InvertShift = bosons.InvertShift;
  this->InvertUnshift = bosons.InvertUnshift;
}

// destructor
//

BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::~BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations & BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::operator = (const BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations & bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrFermions = bosons.NbrFermions;
  this->IncNbrFermions = bosons.IncNbrFermions;
  this->TotalLz = bosons.TotalLz;
  this->NbrFluxQuantumTorus = bosons.NbrFluxQuantumTorus;
  this->NbrFluxQuantumSphere = bosons.NbrFluxQuantumSphere;
  this->ShiftedTotalLz = bosons.ShiftedTotalLz;
  this->NbrFluxQuantumTorus = bosons.NbrFluxQuantumTorus;
  this->NbrOrbitalSphere = bosons.NbrOrbitalSphere;
  this->XMomentum = bosons.XMomentum;
  this->YMomentum = bosons.YMomentum;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->MomentumModulo = bosons.MomentumModulo;
  this->InvertShift = bosons.InvertShift;
  this->InvertUnshift = bosons.InvertUnshift;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
  this->InvertShift = bosons.InvertShift;
  this->InvertUnshift = bosons.InvertUnshift;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::Clone()
{
  return new BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations(*this);
}

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::GenerateStates()
{
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslations;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslations) == this->StateDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = 0x0ul;
	    }
	}
      else
	{
	  this->StateDescription[i] = 0x0ul;
	}
    }
  unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	  this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);
	  ++TmpLargeHilbertSpaceDimension;
	}	
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  return TmpLargeHilbertSpaceDimension;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentLz = current angular momentum along z for a single particle
// currentTotalKy = current total momentum along y
// currentTotalLz = current total angular momentum along z
// currentFermionicPosition = current fermionic position within the state description
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
  
long BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::GenerateStates(int nbrBosons, int currentKy, int currentLz, int currentTotalKy, int currentTotalLz, int currentFermionicPosition, long pos)
{

  if (nbrBosons < 0)
    return pos;
  if (currentLz < 0)
    {
      currentFermionicPosition -= currentLz + 1;
      currentLz = this->NbrOrbitalSphere - 1;
      currentKy--;
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKy % this->NbrFluxQuantumTorus) == this->YMomentum) && (currentTotalLz == this->ShiftedTotalLz))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKy < 0)
    return pos;

  long TmpPos = this->GenerateStates(nbrBosons - 1, currentKy, currentLz - 2, currentTotalKy + currentKy, currentTotalLz + currentLz, currentFermionicPosition - 2, pos);
  unsigned long Mask = 0x1ul << currentFermionicPosition ;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  return this->GenerateStates(nbrBosons, currentKy, currentLz - 1, currentTotalKy, currentTotalLz, currentFermionicPosition - 1, pos);
};


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentLz = current angular momentum along z for a single particle
// currentTotalKy = current total momentum along y
// currentTotalLz = current total angular momentum along z
// return value = Hilbert space dimension

long BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentLz, int currentTotalKy, int currentTotalLz)
{
  if (nbrBosons < 0)
    return 0l;
  if (currentLz < 0)
    {
      currentLz = this->NbrOrbitalSphere - 1;
      currentKy--;
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKy % this->NbrFluxQuantumTorus) == this->YMomentum) && (currentTotalLz == this->ShiftedTotalLz))
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long Count = 0;
  Count += this->EvaluateHilbertSpaceDimension(nbrBosons - 1, currentKy, currentLz - 2, currentTotalKy + currentKy, currentTotalLz + currentLz);
  Count += this->EvaluateHilbertSpaceDimension(nbrBosons, currentKy, currentLz - 1, currentTotalKy, currentTotalLz);
  return Count;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  Str << "[";
  int TmpKx;
  int TmpKy;
  for (int i = 0; i <= this->MaxMomentum; ++i)
    {
      if (((TmpState >> i) & 0x1ul) != 0x0ul)
	{
	  this->GetLinearizedIndex(i, TmpKx, TmpKy);
	  Str << "(" << TmpKx << "," << TmpKy << ")";
	}
    }
  Str << "]";
  return Str;
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals

bool BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  unsigned long TmpState = this->StateDescription[index];
  unsigned long TmpMask = (0x3ul << this->NbrOrbitalSphere) | 0x3ul;
  unsigned long TmpPattern0110 = (0x1ul << this->NbrOrbitalSphere) | 0x2ul;
  unsigned long TmpPattern1001 = (0x2ul << this->NbrOrbitalSphere) | 0x1ul;
  int TmpNbrPatterns = 0;
  if (this->NbrFluxQuantumTorus > 2)
    {
      for (int i = 1; i < this->NbrFluxQuantumTorus; ++i)
	{
	  for (int j = 1;  j < this->NbrOrbitalSphere; ++j)
	    {
	      if (((TmpState & TmpMask) == TmpPattern1001) || ((TmpState & TmpMask) == TmpPattern0110))
		{
		  if ((this->GetSafeOccupationWithPBC(index, i - 2, j - 1) != 0x0ul) || (this->GetSafeOccupationWithPBC(index, i - 2, j) != 0x0ul) ||
		      (this->GetSafeOccupationWithPBC(index, i + 1, j - 1) != 0x0ul) || (this->GetSafeOccupationWithPBC(index, i + 1, j) != 0x0ul) ||
		      (this->GetSafeOccupationWithPBC(index, i - 1, j - 2) != 0x0ul) || (this->GetSafeOccupationWithPBC(index, i, j - 2) != 0x0ul) ||
		      (this->GetSafeOccupationWithPBC(index, i - 1, j + 1) != 0x0ul) || (this->GetSafeOccupationWithPBC(index, i, j + 1) != 0x0ul))
		    {
		      return false;
		    }
		  ++TmpNbrPatterns;
		}
	      TmpMask <<= 1;
	      TmpPattern0110 <<= 1;
	      TmpPattern1001 <<= 1;
	    }
	  TmpMask <<= 1;
	  TmpPattern0110 <<= 1;      
	  TmpPattern1001 <<= 1;
	}
      TmpMask = 0x3ul | (0x3ul << ((this->NbrFluxQuantumTorus - 1) *  this->NbrOrbitalSphere));
      TmpPattern0110 = (0x1ul << ((this->NbrFluxQuantumTorus - 1) *  this->NbrOrbitalSphere)) | 0x2ul;
      TmpPattern1001 = (0x2ul << ((this->NbrFluxQuantumTorus - 1) *  this->NbrOrbitalSphere)) | 0x1ul;
      for (int j = 1;  j < this->NbrOrbitalSphere; ++j)
	{
	  if (((TmpState & TmpMask) == TmpPattern1001) || ((TmpState & TmpMask) == TmpPattern0110))
	    {
	      if ((this->GetSafeOccupationWithPBC(index, -2, j - 1) != 0x0ul) || (this->GetSafeOccupationWithPBC(index, -2, j) != 0x0ul) ||
		  (this->GetSafeOccupationWithPBC(index, 1, j - 1) != 0x0ul) || (this->GetSafeOccupationWithPBC(index, 1, j) != 0x0ul) ||
		  (this->GetSafeOccupationWithPBC(index, -1, j - 2) != 0x0ul) || (this->GetSafeOccupationWithPBC(index, 0, j - 2) != 0x0ul) ||
		  (this->GetSafeOccupationWithPBC(index, -1, j + 1) != 0x0ul) || (this->GetSafeOccupationWithPBC(index, 0, j + 1) != 0x0ul))
		{
		  return false;
		}
	      ++TmpNbrPatterns;
	    }
	  TmpMask <<= 1;
	  TmpPattern0110 <<= 1;
	  TmpPattern1001 <<= 1;
	}      
    }
  else
    {
      for (int i = 1; i < this->NbrFluxQuantumTorus; ++i)
	{
	  for (int j = 1;  j < this->NbrOrbitalSphere; ++j)
	    {
	      if (((TmpState & TmpMask) == TmpPattern1001) || ((TmpState & TmpMask) == TmpPattern0110))
		{
		  if ((this->GetSafeOccupation(index, i - 2, j - 1) != 0x0ul) || (this->GetSafeOccupation(index, i - 2, j) != 0x0ul) ||
		      (this->GetSafeOccupation(index, i + 1, j - 1) != 0x0ul) || (this->GetSafeOccupation(index, i + 1, j) != 0x0ul) ||
		      (this->GetSafeOccupation(index, i - 1, j - 2) != 0x0ul) || (this->GetSafeOccupation(index, i, j - 2) != 0x0ul) ||
		      (this->GetSafeOccupation(index, i - 1, j + 1) != 0x0ul) || (this->GetSafeOccupation(index, i, j + 1) != 0x0ul))
		    {
		      return false;
		    }
		  ++TmpNbrPatterns;
		}
	      TmpMask <<= 1;
	      TmpPattern0110 <<= 1;
	      TmpPattern1001 <<= 1;
	    }
	  TmpMask <<= 1;
	  TmpPattern0110 <<= 1;      
	  TmpPattern1001 <<= 1;
	}
    }
  if (TmpNbrPatterns > 0)
    {
      int TmpNbrTranslations;
      if (TmpState != this->FindCanonicalForm(this->FindCanonical(TmpState, 0, 0), TmpNbrTranslations))
	{
	  return false;
	}
    }
  return true;
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals

unsigned long BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::FindCanonical(unsigned long state, int xPosition, int yPosition)
{
  if (yPosition >= (this->NbrOrbitalSphere - 1))
    {
      yPosition = 0;
      xPosition++;
    }
  if (xPosition >= this->NbrFluxQuantumTorus)
    {
      return state;
    }
  int TmpShift;
  unsigned long TmpMask;
  unsigned long TmpPattern0110;
  unsigned long TmpPattern1001;
  if (xPosition == (this->NbrFluxQuantumTorus - 1))
    {
      TmpShift = yPosition;
      TmpMask = 0x3ul | (0x3ul << ((this->NbrFluxQuantumTorus - 1) *  this->NbrOrbitalSphere));
      TmpPattern0110 = (0x1ul << ((this->NbrFluxQuantumTorus - 1) *  this->NbrOrbitalSphere)) | 0x2ul;
      TmpPattern1001 = (0x2ul << ((this->NbrFluxQuantumTorus - 1) *  this->NbrOrbitalSphere)) | 0x1ul;
    }
  else
    {
      TmpShift = (xPosition *  this->NbrOrbitalSphere) + yPosition;
      TmpMask = (0x3ul << this->NbrOrbitalSphere) | 0x3ul;
      TmpPattern0110 = (0x1ul << this->NbrOrbitalSphere) | 0x2ul;
      TmpPattern1001 = (0x2ul << this->NbrOrbitalSphere) | 0x1ul;
    }
  TmpMask <<= TmpShift;
  TmpPattern0110 <<= TmpShift;
  TmpPattern1001 <<= TmpShift;
  if ((state & TmpMask) == TmpPattern0110)
    {
      unsigned long TmpState1 = this->FindCanonical(state, xPosition, yPosition + 1);
      unsigned long TmpState2 = this->FindCanonical((state & ~TmpMask) | TmpPattern1001, xPosition, yPosition + 1);
      if (TmpState1 < TmpState2)
	return TmpState2;
      else
	return TmpState1;
    }
  if ((state & TmpMask) == TmpPattern1001)
    {
      unsigned long TmpState1 = this->FindCanonical(state, xPosition, yPosition + 1);
      unsigned long TmpState2 = this->FindCanonical((state & ~TmpMask) | TmpPattern0110, xPosition, yPosition + 1);
      if (TmpState1 < TmpState2)
	return TmpState2;
      else
	return TmpState1;      
    }
  return this->FindCanonical(state, xPosition, yPosition + 1);
}
