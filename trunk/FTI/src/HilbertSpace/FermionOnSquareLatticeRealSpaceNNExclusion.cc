////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice in real space           //
//                        with nearest neighbor exclusion                     //
//                                                                            //
//                        last modification : 11/02/2018                      //
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
#include "HilbertSpace/FermionOnSquareLatticeRealSpaceNNExclusion.h"
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

FermionOnSquareLatticeRealSpaceNNExclusion::FermionOnSquareLatticeRealSpaceNNExclusion ()
{
  this->NbrSitesX = 0;
  this->NbrSitesY = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeRealSpaceNNExclusion::FermionOnSquareLatticeRealSpaceNNExclusion (int nbrFermions, int nbrSitesX, int nbrSitesY, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;

  this->NbrSitesX = nbrSitesX;
  this->NbrSitesY = nbrSitesY;  
  this->NbrSite = this->NbrSitesX * this->NbrSitesY;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSite - 1, this->NbrSitesY - 1, 0x0ul, 0x0ul);
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
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrSitesY - 1, 0x0ul, 0x0ul, 0l);
      
//       for (long i = 1l; i < TmpLargeHilbertSpaceDimension; ++i)
// 	{
// 	  if (this->StateDescription[i - 1] < this->StateDescription[i])
// 	    {
// 	      cout << i << " " << this->StateDescription[i - 1] << " " << this->StateDescription[i] << endl;
// 	    }
// 	}
//       unsigned long TmpMask = (0x1ul << this->NbrSitesY) - 0x1ul;
//       long TmpDim = 0;
//       for (long i = 0l; i < TmpLargeHilbertSpaceDimension; ++i)
// 	{
// 	  bool Flag = false;
// 	  unsigned long Tmp = this->StateDescription[i];
// 	  for (int j = 1; (j < this->NbrSitesX) && (Flag == false); ++j)
// 	    {
// 	      if (((Tmp & (Tmp >> this->NbrSitesY)) & TmpMask) != 0x0ul)
// 		{
// 		  Flag = true;
// 		}
// 	      Tmp >>= this->NbrSitesY;
// 	    }
// 	  if (Flag == false)
// 	    {
// 	      ++TmpDim;
// 	    }
// 	}
//       cout << "reduced dim = " << TmpDim << endl;

      this->StateLzMax = new int [this->LargeHilbertSpaceDimension];  
      int CurrentLzMax = this->NbrLzValue;
      while ((((this->StateDescription[0] >> CurrentLzMax) & 0x1ul) == 0x0ul) && (CurrentLzMax >= 0))
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
      
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  if (this->FindStateIndex(this->StateDescription[i],  this->StateLzMax[i]) != i)
	    {
	      cout << i << " " << this->FindStateIndex(this->StateDescription[i],  this->StateLzMax[i]) << endl;
	    }
	}


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

FermionOnSquareLatticeRealSpaceNNExclusion::FermionOnSquareLatticeRealSpaceNNExclusion(const FermionOnSquareLatticeRealSpaceNNExclusion& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSitesX = fermions.NbrSitesX;
  this->NbrSitesY = fermions.NbrSitesY;  
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

FermionOnSquareLatticeRealSpaceNNExclusion::~FermionOnSquareLatticeRealSpaceNNExclusion ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeRealSpaceNNExclusion& FermionOnSquareLatticeRealSpaceNNExclusion::operator = (const FermionOnSquareLatticeRealSpaceNNExclusion& fermions)
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
  this->NbrSitesX = fermions.NbrSitesX;
  this->NbrSitesY = fermions.NbrSitesY;  
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

AbstractHilbertSpace* FermionOnSquareLatticeRealSpaceNNExclusion::Clone()
{
  return new FermionOnSquareLatticeRealSpaceNNExclusion(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site linearized index
// currentSiteY = y coordinate of the current site
// previousXConfiguration = configuraton at the previous x coordinate
// currentXConfiguration = configuraton at the current x coordinate
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeRealSpaceNNExclusion::GenerateStates(int nbrFermions, int currentSite, int currentSiteY, unsigned long previousXConfiguration, unsigned long currentXConfiguration, long pos)
{
  if (nbrFermions == 0)
    {
      this->StateDescription[pos] = 0x0ul;	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    {
      return pos;
    }
  if (currentSiteY < 0)
    {
      currentSiteY = this->NbrSitesY - 1;
      previousXConfiguration = currentXConfiguration;
      currentXConfiguration = 0x0ul;
    }
  long TmpPos;
  if ((previousXConfiguration & (0x1ul << currentSiteY)) == 0x0ul)
    {
      if (currentSiteY == 0)
	{
	  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, currentSiteY - 1, previousXConfiguration, currentXConfiguration | (0x1ul << currentSiteY), pos);
	}
      else
	{
	  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 2, currentSiteY - 2, previousXConfiguration, currentXConfiguration | (0x1ul << currentSiteY), pos);      
	}
      unsigned long Mask = 0x1ul << currentSite;
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return this->GenerateStates(nbrFermions, currentSite - 1, currentSiteY - 1, previousXConfiguration, currentXConfiguration, pos);
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentSite = current site linearized index
// currentSiteY = y coordinate of the current site
// previousXConfiguration = configuraton at the previous x coordinate
// currentXConfiguration = configuraton at the current x coordinate
// return value = Hilbert space dimension

long FermionOnSquareLatticeRealSpaceNNExclusion::EvaluateHilbertSpaceDimension(int nbrFermions, int currentSite, int currentSiteY, unsigned long previousXConfiguration, unsigned long currentXConfiguration)
{
  if (nbrFermions == 0)
    {
      return 1l;
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    {
      return 0l;
    }
  if (currentSiteY < 0)
    {
      currentSiteY = this->NbrSitesY - 1;
      previousXConfiguration = currentXConfiguration;
      currentXConfiguration = 0x0ul;
    }
  long TmpDimension = 0l;
  if ((previousXConfiguration & (0x1ul << currentSiteY)) == 0x0ul)
    {
      if (currentSiteY == 0)
	{
	  TmpDimension = this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentSite - 1, currentSiteY - 1, previousXConfiguration, currentXConfiguration | (0x1ul << currentSiteY));
	}
      else
	{
	  TmpDimension = this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentSite - 2, currentSiteY - 2, previousXConfiguration, currentXConfiguration | (0x1ul << currentSiteY));
	}
    }
  TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions, currentSite - 1, currentSiteY - 1, previousXConfiguration, currentXConfiguration);
  return TmpDimension;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSquareLatticeRealSpaceNNExclusion::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  Str << "| ";
  for (int i = 0; i < this->NbrSite; i += this->NbrSitesY)
    {
      for (int j = 0; j < this->NbrSitesY; ++j)
	{
	  Str << ((TmpState >> (i + j)) & 0x1ul) << " ";
	}
      Str << "| ";
    }
  return Str;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSquareLatticeRealSpaceNNExclusion::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    {
      return this->HilbertSpaceDimension;
    }
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}


  
