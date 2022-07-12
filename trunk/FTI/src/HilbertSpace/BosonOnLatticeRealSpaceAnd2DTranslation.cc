////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on lattice			              //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 12/09/2014                      //
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
#include "HilbertSpace/BosonOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
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
#include "Architecture/ArchitectureOperation/FQHETorusParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

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

BosonOnLatticeRealSpaceAnd2DTranslation::BosonOnLatticeRealSpaceAnd2DTranslation ()
{
  this->NbrBosons = 0;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = 0;
  this->NbrSite = 0;
  this->MomentumModulo = 1;
  this->KxMomentum = 0; 
  this->KyMomentum = 0;
  this->StateShift = 2 * (this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrBosons * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->LastMomentumMask = ((unsigned long) 1);
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrBosons = number of bosons
// nbrSite = number of sites
// xMomentum = momentum sector in the x direction
// maxYMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction 
// memory = amount of memory granted for precalculations

BosonOnLatticeRealSpaceAnd2DTranslation::BosonOnLatticeRealSpaceAnd2DTranslation (int nbrBosons, int nbrSite, int xMomentum, int  maxXMomentum,
										      int yMomentum, int maxYMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->FermionicMaxMomentum = this->MaxMomentum + this->NbrBosons - 1;
  this->NbrMomentum = this->MaxMomentum;
  this->MaxXMomentum = maxXMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->TemporaryState = new unsigned long [this->NbrMomentum];
  this->ProdATemporaryState = new unsigned long [this->NbrMomentum];
  this->TemporaryStateOperators = new unsigned long [this->NbrMomentum];

  this->StateXShift = 1;
  this->MomentumIncrement = (this->NbrBosons * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->LastMomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->KxMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = this->NbrSite / this->MaxXMomentum;
  this->ComplementaryStateXShift = this->MaxMomentum - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->MaxYMomentum = maxYMomentum;
  this->NbrSitePerUnitCell = this->NbrSite / (this->MaxXMomentum * this->MaxYMomentum);
  this->KyMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = this->MaxXMomentum;
  this->StateYShift = (this->NbrSite / (this->MaxXMomentum * this->MaxYMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
	
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->LargeHilbertSpaceDimension  = this->GenerateStates();
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  /*this->StateMaxMomentum = new int [this->LargeHilbertSpaceDimension];  
	  int CurrentMaxMomentum = this->MaxMomentum;
	  while (((this->StateDescription[0] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
	    --CurrentMaxMomentum;
	  this->StateMaxMomentum[0] = CurrentMaxMomentum;
	  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      while (((this->StateDescription[i] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
		--CurrentMaxMomentum;
	      this->StateMaxMomentum[i] = CurrentMaxMomentum;
	    }
*/
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
	  UsedMemory = this->NbrMomentum * sizeof(int);
	  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
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
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnLatticeRealSpaceAnd2DTranslation::BosonOnLatticeRealSpaceAnd2DTranslation(const BosonOnLatticeRealSpaceAnd2DTranslation& bosons)
{
  this->NbrBosons = bosons.NbrBosons;  
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->NbrSite = bosons.NbrSite;

  this->MaxMomentum = bosons.MaxMomentum;
  this->MaxXMomentum = bosons.MaxXMomentum;
  this->MaxYMomentum = bosons.MaxYMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->MomentumModulo = bosons.MomentumModulo;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->FermionicMaxMomentum = bosons.FermionicMaxMomentum;

  this->StateXShift = bosons.StateXShift;
  this->StateYShift = bosons.StateYShift;
  this->NbrYMomentumBlocks = bosons.NbrYMomentumBlocks;
 
  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->LastMomentumMask = bosons.LastMomentumMask;
  this->NbrSitePerUnitCell = bosons.NbrSitePerUnitCell;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;

  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;

  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->Flag = bosons.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  this->TemporaryState = new unsigned long [this->NbrMomentum];
  this->ProdATemporaryState = new unsigned long [this->NbrMomentum];	
  this->TemporaryStateOperators = new unsigned long [this->NbrMomentum];
}

// destructor
//

BosonOnLatticeRealSpaceAnd2DTranslation::~BosonOnLatticeRealSpaceAnd2DTranslation ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeRealSpaceAnd2DTranslation& BosonOnLatticeRealSpaceAnd2DTranslation::operator = (const BosonOnLatticeRealSpaceAnd2DTranslation& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
      delete[] this->TemporaryState;
      delete[] this->ProdATemporaryState;
      delete[] this->TemporaryStateOperators;
    }
  this->NbrBosons = bosons.NbrBosons;  
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->NbrSite = bosons.NbrSite;

  this->MaxMomentum = bosons.MaxMomentum;
  this->MaxXMomentum = bosons.MaxXMomentum;
  this->MaxYMomentum = bosons.MaxYMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->MomentumModulo = bosons.MomentumModulo;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->FermionicMaxMomentum = bosons.FermionicMaxMomentum;

  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->LastMomentumMask = bosons.LastMomentumMask;
  this->NbrSitePerUnitCell = bosons.NbrSitePerUnitCell;
  
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateXShift = bosons.StateXShift;
  this->StateYShift = bosons.StateYShift;
  this->NbrYMomentumBlocks = bosons.NbrYMomentumBlocks;
 

  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;

  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->Flag = bosons.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->TemporaryState = new unsigned long [this->NbrMomentum];
  this->ProdATemporaryState = new unsigned long [this->NbrMomentum];	
  this->TemporaryStateOperators = new unsigned long [this->NbrMomentum];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeRealSpaceAnd2DTranslation::Clone()
{
  return new BosonOnLatticeRealSpaceAnd2DTranslation(*this);
}


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// return value = Hilbert space dimension
 long BosonOnLatticeRealSpaceAnd2DTranslation::EvaluateHilbertSpaceDimension(int nbrBosons)
{
  BinomialCoefficients binomials(this->NbrSite);
  long dimension = binomials(this->NbrSite+this->NbrBosons - 1, this->NbrBosons);
  return dimension;
}

// generate all states 
// 
// nbrBosons = number of bosons
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnLatticeRealSpaceAnd2DTranslation::RawGenerateStates(int nbrBosons, int currentSite, long pos)
{
  if (nbrBosons == 0)
    {
      this->StateDescription[pos] = 0x0ul;	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrBosons < 0))
    return pos;
  if (nbrBosons == 1)
    {
      for (int j = currentSite; j >= 0; --j)
	{
	  this->StateDescription[pos] = 0x1ul << j;
	  ++pos;
	}
      return pos;
    }

  long TmpPos = pos;
  for (int i = nbrBosons; i > 0; --i)
    {
       TmpPos = this->RawGenerateStates(nbrBosons - i,  currentSite - 1, pos);
       unsigned long Mask = ((0x1ul << i) - 0x1ul) << (currentSite + nbrBosons - i);
	for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
  return this->RawGenerateStates(nbrBosons, currentSite - 1, pos);
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long BosonOnLatticeRealSpaceAnd2DTranslation::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(this->NbrBosons, this->NbrSite - 1, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY) == this->StateDescription[i]))
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
  //cout << "new dim = " << TmpLargeHilbertSpaceDimension << endl;
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

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnLatticeRealSpaceAnd2DTranslation::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;
  if (this->LookUpTableShift[lzmax] < 0)
    return this->HilbertSpaceDimension;

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
 
// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnLatticeRealSpaceAnd2DTranslation::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescription[state], this->MaxMomentum + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  int i = 0;
  for (; i <= this->TemporaryStateKyMax; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i < this->NbrMomentum; ++i)
    Str << "0 ";
 return Str;
}



// apply a^+_m_sigma a_n_sigma operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator including the orbital and the spin index
// n = index of the annihilation operator including the orbital and the spin index
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int BosonOnLatticeRealSpaceAnd2DTranslation::AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{

  this->FermionToBoson(this->StateDescription[index], this->FermionicMaxMomentum,this->TemporaryState,this->TemporaryStateKyMax); 
  if ((n > this->TemporaryStateKyMax)||(this->TemporaryState[n]==0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }	

  coefficient = (double) this->TemporaryState[n];
  --this->TemporaryState[n];

  if (m >  this->TemporaryStateKyMax)
	{
		for(int i = this->TemporaryStateKyMax +1; i <= m; i++)
	   	 this->TemporaryState[i] =0;
                this->TemporaryStateKyMax = m;
	}
  this->TemporaryState[m]++;
  coefficient *= (double) this->TemporaryState[m];
  coefficient = sqrt(coefficient); 
  this->ProdATemporaryStateNbrStateInOrbit  =  this->NbrStateInOrbit[index];
  unsigned long TmpState = this->BosonToFermion(this->TemporaryState,this->TemporaryStateKyMax);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}
  
// apply a^+_m1_sigma a^+_m2_sigma operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int BosonOnLatticeRealSpaceAnd2DTranslation::AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  for(int i =0; i < this->NbrMomentum ;i++)
      this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  coefficient = this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);

  this->TemporaryStateKyMax = this->MaxMomentum - 1;
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0)
    --this->TemporaryStateKyMax; 
  unsigned long TmpState = this->BosonToFermion(this->TemporaryState,this->TemporaryStateKyMax);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}


// apply a^+_m operator to the state produced using AuAu method (without destroying it)
//
// m = first index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int BosonOnLatticeRealSpaceAnd2DTranslation::Ad (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  for(int i =0; i < this->NbrMomentum ;i++)
      this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m];
  coefficient *= this->TemporaryState[m];
  coefficient = sqrt(coefficient);
  
  this->TemporaryStateKyMax = this->MaxMomentum - 1;
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0)
    --this->TemporaryStateKyMax; 
  unsigned long TmpState = this->BosonToFermion(this->TemporaryState,this->TemporaryStateKyMax);  
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);  
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = kx sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnLatticeRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementaryKxMomentum = (this->KxMomentum - kxSector);
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->MaxXMomentum;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector);
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->MaxYMomentum;
  BosonOnLatticeRealSpaceAnd2DTranslation SubsytemSpace (nbrParticleSector, this->NbrSite, kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnLatticeRealSpaceAnd2DTranslation ComplementarySpace (ComplementaryNbrParticles, this->NbrSite, ComplementaryKxMomentum, this->MaxXMomentum, 
							      ComplementaryKyMomentum, this->MaxYMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  FQHETorusParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, (ParticleOnTorusWithMagneticTranslations*) &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  cout << "nbr matrix elements non zero = " << Operation.GetNbrNonZeroMatrixElements() << endl;
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

long BosonOnLatticeRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, 
												 ParticleOnTorusWithMagneticTranslations* complementaryHilbertSpace,  
												 ParticleOnTorusWithMagneticTranslations* destinationHilbertSpace,
												 ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
  BosonOnLatticeRealSpaceAnd2DTranslation* TmpDestinationHilbertSpace =  (BosonOnLatticeRealSpaceAnd2DTranslation*) destinationHilbertSpace;
  BosonOnLatticeRealSpaceAnd2DTranslation* TmpHilbertSpace = (BosonOnLatticeRealSpaceAnd2DTranslation*) complementaryHilbertSpace;
  BosonOnLatticeRealSpace* TmpDestinationFullHilbertSpace = new BosonOnLatticeRealSpace(TmpDestinationHilbertSpace->NbrBosons,
											TmpDestinationHilbertSpace->NbrSite);
  unsigned long* TmpMonomial2 = new unsigned long [TmpDestinationHilbertSpace->NbrBosons];
  unsigned long* TmpMonomial1 = new unsigned long [TmpHilbertSpace->NbrBosons];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
  Complex** FourrierCoefficientsDestination = new Complex* [this->MomentumModulo];
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  long TmpNbrNonZeroElements = 0l;
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficientsDestination[i] = new Complex [this->MaxYMomentum];
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficientsDestination[i][j] = Phase (2.0 * M_PI * ((double) (i * TmpDestinationHilbertSpace->KxMomentum) / ((double) this->MaxXMomentum) + (double) (j * TmpDestinationHilbertSpace->KyMomentum) / ((double) this->MaxYMomentum)));
	  FourrierCoefficients[i][j] = Phase (2.0 * M_PI * ((double) (i * this->KxMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->KyMomentum) / ((double) this->MaxYMomentum)));
	}
    }
  
  ComplexMatrix TmpEntanglementMatrix (nbrIndex, TmpDestinationHilbertSpace->GetHilbertSpaceDimension(), true);
  int TmpMinIndex = minIndex;
  int MaxIndex = minIndex + nbrIndex;
  this->ProdATemporaryStateNbrStateInOrbit = 1;
  TmpDestinationHilbertSpace->ProdATemporaryStateNbrStateInOrbit = 1;  

  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
  double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[TmpHilbertSpace->NbrBosons] - LogFactorials[TmpDestinationHilbertSpace->NbrBosons];

  double* TmpDestinationLogFactorials = new double [TmpDestinationFullHilbertSpace->HilbertSpaceDimension];
  Complex* TmpDestinationFactors = new Complex[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
  int* TmpDestinationRealIndices = new int[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
  for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
    {
      unsigned long TmpCanonicalState2 = TmpDestinationFullHilbertSpace->FermionBasis->StateDescription[j];
      int TmpDestinationNbrTranslationX;
      int TmpDestinationNbrTranslationY;
      double TmpDestinationCoefficient = 1.0;
      int RealDestinationIndex = TmpDestinationHilbertSpace->SymmetrizeAdAdResult(TmpCanonicalState2, TmpDestinationCoefficient, 
										  TmpDestinationNbrTranslationX, TmpDestinationNbrTranslationY);
      TmpDestinationRealIndices[j] = RealDestinationIndex;
      if (RealDestinationIndex < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
	{
	  TmpDestinationFactors[j] = TmpDestinationCoefficient * FourrierCoefficientsDestination[TmpDestinationNbrTranslationX][TmpDestinationNbrTranslationY];
	  TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationFullHilbertSpace->FermionBasis->StateDescription[j],
						     TmpDestinationHilbertSpace->FermionicMaxMomentum, 
						     TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateKyMax);
	  double TmpFactor = 0.0;
	  for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateKyMax; ++k)
	    TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
	  TmpDestinationLogFactorials[j] = TmpFactor;
	}
    }
  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->ConvertToMonomial(TmpHilbertSpace->StateDescription[minIndex], TmpHilbertSpace->FermionicMaxMomentum, TmpMonomial1);
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->StateDescription[minIndex], TmpHilbertSpace->FermionicMaxMomentum, 
				      TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateKyMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateKyMax; ++k)
 	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      double TmpRescalingFactor = sqrt((double) TmpHilbertSpace->NbrStateInOrbit[minIndex]);
      for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  if (TmpDestinationRealIndices[j] < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
	    {
	      TmpDestinationFullHilbertSpace->ConvertToMonomial(TmpDestinationFullHilbertSpace->FermionBasis->StateDescription[j], TmpDestinationFullHilbertSpace->FermionBasis->StateLzMax[j], TmpMonomial2);
	      int TmpIndex2 = 0;
	      int TmpIndex3 = 0;
	      int TmpIndex4 = 0;
	      while ((TmpIndex2 < TmpHilbertSpace->NbrBosons) && (TmpIndex3 < TmpDestinationFullHilbertSpace->NbrBosons)) 
		{
		  while ((TmpIndex2 < TmpHilbertSpace->NbrBosons) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		      ++TmpIndex2;
		      ++TmpIndex4;		  
		    }
		  if (TmpIndex2 < TmpHilbertSpace->NbrBosons)
		    {
		      while ((TmpIndex3 < TmpDestinationFullHilbertSpace->NbrBosons) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
			{
			  TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
			  ++TmpIndex3;
			  ++TmpIndex4;		  
			}
		    }
		}
	      while (TmpIndex2 < TmpHilbertSpace->NbrBosons)
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      while (TmpIndex3 < TmpDestinationFullHilbertSpace->NbrBosons)
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		  ++TmpIndex3;
		  ++TmpIndex4;		  
		}
	      int TmpNbrTranslationX;
	      int TmpNbrTranslationY;
	      unsigned long TmpState3 = this->ConvertFromMonomial(TmpMonomial3);
	      double TmpCoefficient = TmpRescalingFactor;
	      int TmpPos = this->SymmetrizeAdAdResult(TmpState3, TmpCoefficient, TmpNbrTranslationX, TmpNbrTranslationY);
	      if (TmpPos < this->HilbertSpaceDimension)
		{      
		  this->FermionToBoson(TmpState3, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
		  double TmpFactorial = 0.0;	      
		  for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
		    TmpFactorial += LogFactorials[this->TemporaryState[k]];
		  TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
		  TmpFactorial *= 0.5; 
		  TmpCoefficient *= exp(TmpFactorial);
		  TmpEntanglementMatrix[TmpDestinationRealIndices[j]][minIndex - TmpMinIndex] += (TmpCoefficient * groundState[TmpPos] * TmpDestinationFactors[j]
												  * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY]); 
		}

	    }
	}
    }
  
  for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      for (int j = i; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  densityMatrix->UnsafeAddToMatrixElement(i, j, TmpEntanglementMatrix[j] * TmpEntanglementMatrix[i]);
	  ++TmpNbrNonZeroElements;
	}
    }
  delete TmpDestinationFullHilbertSpace;
  delete[] TmpDestinationFactors;
  delete[] TmpDestinationRealIndices;
  return TmpNbrNonZeroElements;
}

