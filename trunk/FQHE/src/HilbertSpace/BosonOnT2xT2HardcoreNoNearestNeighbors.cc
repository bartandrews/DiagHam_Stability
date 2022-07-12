////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2014 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on the 4D manifold T2 x T2               //
// forbidding mutliple orbital occupations and nearest orbital occupations    //
//                                                                            //
//                        last modification : 16/02/2016                      //
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
#include "HilbertSpace/BosonOnT2xT2HardcoreNoNearestNeighbors.h"
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

BosonOnT2xT2HardcoreNoNearestNeighbors::BosonOnT2xT2HardcoreNoNearestNeighbors ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// nbrFluxQuanta1 = number of flux quanta for the first torus
// nbrFluxQuanta2 = number of flux quanta for the second torus
// totalKy1 = total momentum along y for the first torus
// totalKy2 = total momentum along y for the second torus

BosonOnT2xT2HardcoreNoNearestNeighbors::BosonOnT2xT2HardcoreNoNearestNeighbors (int nbrBosons, int nbrFluxQuanta1, int nbrFluxQuanta2, int totalKy1, int totalKy2, unsigned long memory)
{  
  this->NbrFermions = nbrBosons;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrSiteX =  nbrFluxQuanta1;
  this->NbrSiteY =  nbrFluxQuanta2;
  this->KxMomentum = totalKy1;
  this->KyMomentum = totalKy2;
  this->LzMax = this->NbrSiteX * this->NbrSiteY - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0);
  cout << "dim = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->LzMax, 0, 0l);
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space : get " << TmpLargeHilbertSpaceDimension << " , should be " << this->LargeHilbertSpaceDimension << endl;
	}
      this->InvertShift = 0;
      this->InvertUnshift = 0;
      if (this->NbrFermions == 0)
	{
	  this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
	  this->StateDescription[0] = 0x0ul; 
	  this->StateLzMax[0] = 0;
	}
      else
	{
	  if (this->NbrSiteX > 1)
	    {
	      unsigned long TmpMask = (0x1ul << this->NbrSiteY) - 0x1ul;
	      unsigned long TmpMask2 = (0x1ul << this->NbrSiteY) | 0x1ul;
	      int TmpShift = (this->NbrSiteX - 1) * this->NbrSiteY;
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
		      while ((TmpIndex < this->NbrLzValue) && (((this->StateDescription[i] >> TmpIndex) & TmpMask2) != TmpMask2))
			{
			  ++TmpIndex;		  
			}
		      if (TmpIndex < this->NbrLzValue)
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
	      this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
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
      UsedMemory = this->NbrLzValue * sizeof(int);
//      UsedMemory += this->NbrLzValue * this->FermionBasis->LookUpTableMemorySize * sizeof(int);
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

BosonOnT2xT2HardcoreNoNearestNeighbors::BosonOnT2xT2HardcoreNoNearestNeighbors(const BosonOnT2xT2HardcoreNoNearestNeighbors& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrFermions = bosons.NbrFermions;
  this->IncNbrFermions = bosons.IncNbrFermions;
  this->TotalLz = bosons.TotalLz;
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
  this->InvertShift = bosons.InvertShift;
  this->InvertUnshift = bosons.InvertUnshift;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnT2xT2HardcoreNoNearestNeighbors::~BosonOnT2xT2HardcoreNoNearestNeighbors ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnT2xT2HardcoreNoNearestNeighbors & BosonOnT2xT2HardcoreNoNearestNeighbors::operator = (const BosonOnT2xT2HardcoreNoNearestNeighbors & bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrFermions = bosons.NbrFermions;
  this->IncNbrFermions = bosons.IncNbrFermions;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
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

AbstractHilbertSpace* BosonOnT2xT2HardcoreNoNearestNeighbors::Clone()
{
  return new BosonOnT2xT2HardcoreNoNearestNeighbors(*this);
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentFermionicPosition = current fermionic position within the state description
// pos = position in StateDescription array where to store states
// occupationLastOrbital = occupation of the last orbital at the current kx value
// return value = position from which new states have to be stored
  
long BosonOnT2xT2HardcoreNoNearestNeighbors::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentFermionicPosition, int occupationLastOrbital, long pos)
{

  if (nbrBosons < 0)
    return pos;
  if (currentKy < 0)
    {
      currentFermionicPosition -= currentKy + 1;
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrBosons == 0)
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

  long TmpPos = 0;
  if ((currentKy != 0) || (occupationLastOrbital == 0))
    {
      if (currentKy == (this->NbrSiteY - 1))
	{
	  TmpPos += this->GenerateStates(nbrBosons - 1, currentKx, currentKy - 2, currentTotalKx + currentKx, currentTotalKy + currentKy, currentFermionicPosition - 2, 
					 1, pos);
	}
      else
	{
	  TmpPos += this->GenerateStates(nbrBosons - 1, currentKx, currentKy - 2, currentTotalKx + currentKx, currentTotalKy + currentKy, currentFermionicPosition - 2, 
					 occupationLastOrbital, pos);
	}
    }
  unsigned long Mask = 0x1ul << currentFermionicPosition ;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  if (currentKy == (this->NbrSiteY - 1))
    {
      return this->GenerateStates(nbrBosons, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, currentFermionicPosition - 1, 0, pos);
    }
  else
    {
      return this->GenerateStates(nbrBosons, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, currentFermionicPosition - 1, occupationLastOrbital, pos);
    }
};


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// occupationLastOrbital = occupation of the last orbital at the current kx value
// return value = Hilbert space dimension

long BosonOnT2xT2HardcoreNoNearestNeighbors::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int occupationLastOrbital)
{
  if (nbrBosons < 0)
    return 0l;
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if ((currentKy != 0) || (occupationLastOrbital == 0))
    {
      if (currentKy == (this->NbrSiteY - 1))
	{
	  Count += this->EvaluateHilbertSpaceDimension(nbrBosons - 1, currentKx, currentKy - 2, currentTotalKx + currentKx, currentTotalKy + currentKy, 1);
	}
      else
	{
	  Count += this->EvaluateHilbertSpaceDimension(nbrBosons - 1, currentKx, currentKy - 2, currentTotalKx + currentKx, currentTotalKy + currentKy, occupationLastOrbital);
	}
    }
  if (currentKy == (this->NbrSiteY - 1))
    {
      Count += this->EvaluateHilbertSpaceDimension(nbrBosons, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, 0);
    }
  else
    {
      Count += this->EvaluateHilbertSpaceDimension(nbrBosons, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, occupationLastOrbital);
    }
  return Count;
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals

bool BosonOnT2xT2HardcoreNoNearestNeighbors::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  unsigned long TmpState = this->StateDescription[index];
  unsigned long TmpMask = (0x3ul << this->NbrSiteY) | 0x3ul;
  unsigned long TmpPattern0110 = (0x1ul << this->NbrSiteY) | 0x2ul;
  unsigned long TmpPattern1001 = (0x2ul << this->NbrSiteY) | 0x1ul;
  int TmpNbrPatterns = 0;
  if (this->NbrSiteX > 2)
    {
      if (this->NbrSiteY > 2)
	{
	  for (int i = 1; i <= this->NbrSiteX; ++i)
	    {
	      for (int j = 1;  j <= this->NbrSiteY; ++j)
		{
		  unsigned long TmpPattern = ((this->GetSafeOccupationWithPBC(index, i - 1, j - 1)) |
					      (this->GetSafeOccupationWithPBC(index, i, j - 1) << 1) |
					      (this->GetSafeOccupationWithPBC(index, i - 1, j) << 2) |
					      (this->GetSafeOccupationWithPBC(index, i, j) << 3));
		  if ((TmpPattern == 0x9ul) || (TmpPattern == 0x6ul))
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
		}
	    }
	}
      else
	{
	}
    }
  else
    {
      if (this->NbrSiteY > 2)
	{
	  for (int i = 1; i < this->NbrSiteX; ++i)
	    {
	      for (int j = 1;  j < this->NbrSiteY; ++j)
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
      else
	{
	}
    }
  if (TmpNbrPatterns > 0)
    {
      if (TmpState != this->FindCanonical(TmpState, 0, 0))
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

unsigned long BosonOnT2xT2HardcoreNoNearestNeighbors::FindCanonical(unsigned long state, int xPosition, int yPosition)
{
  if (yPosition >= this->NbrSiteY)
    {
      yPosition = 0;
      xPosition++;
    }
  if (xPosition >= this->NbrSiteX)
    {
      return state;
    }
  int TmpShift;
  unsigned long TmpMask;
  unsigned long TmpPattern0110;
  unsigned long TmpPattern1001;
  if (xPosition == (this->NbrSiteX - 1))
    {
      if (yPosition == (this->NbrSiteY - 1))
	{
	  TmpShift = 0;
	  TmpMask = ((0x1ul | (0x1ul << (this->NbrSiteY - 1))) 
		     | ((0x1ul | (0x1ul << (this->NbrSiteY - 1))) << ((this->NbrSiteX - 1) *  this->NbrSiteY)));
	  TmpPattern0110 = (0x1ul << (this->NbrSiteY - 1)) | (0x1ul << ((this->NbrSiteX - 1) *  this->NbrSiteY));
	  TmpPattern1001 = 0x1ul | (0x1ul << ((this->NbrSiteX * this->NbrSiteY) - 1));
	}
      else
	{
	  TmpShift = yPosition;
	  TmpMask = 0x3ul | (0x3ul << ((this->NbrSiteX - 1) *  this->NbrSiteY));
	  TmpPattern0110 = (0x1ul << ((this->NbrSiteX - 1) *  this->NbrSiteY)) | 0x2ul;
	  TmpPattern1001 = (0x2ul << ((this->NbrSiteX - 1) *  this->NbrSiteY)) | 0x1ul;
	}
    }
  else
    {
      if (yPosition == (this->NbrSiteY - 1))
	{
	  TmpShift = (xPosition *  this->NbrSiteY);
	  TmpMask = 0x1ul | (0x3ul << (this->NbrSiteY - 1)) | (0x1ul << ((2 * this->NbrSiteY) - 1));
	  TmpPattern0110 = 0x1ul | (0x1ul << ((2 * this->NbrSiteY) - 1));
	  TmpPattern1001 = 0x3ul << (this->NbrSiteY - 1);
	}
      else
	{
	  TmpShift = (xPosition *  this->NbrSiteY) + yPosition;
	  TmpMask = (0x3ul << this->NbrSiteY) | 0x3ul;
	  TmpPattern0110 = (0x1ul << this->NbrSiteY) | 0x2ul;
	  TmpPattern1001 = (0x2ul << this->NbrSiteY) | 0x1ul;
	}
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
