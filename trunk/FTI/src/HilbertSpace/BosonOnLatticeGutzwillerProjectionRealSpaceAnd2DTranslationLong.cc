////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons hardcore on lattice                      //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 30/06/2016                      //
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
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong.h"
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

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong ()
{
  this->NbrBosons = 0;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = 0;
  this->NbrSite = 0;
  this->MomentumModulo = 0;
  this->XMomentum = 0;
  this->YMomentum = 0;
  this->StateShift = 2 * (this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrBosons * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((ULONGLONG) 0x1ul);
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateMaxMomentum = 0;
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrBosons = number of fermions
// nbrSite = number of sites
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction 
// memory = amount of memory granted for precalculations

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong (int nbrBosons, int nbrSite, int xMomentum, int  maxXMomentum,  int yMomentum, int maxYMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;

  this->MomentumModulo = this->MaxXMomentum;
  this->StateShift = this->MaxMomentum / this->MomentumModulo;

  this->StateXShift = 1;
  this->MomentumIncrement = (this->NbrBosons * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (((ULONGLONG) 0x1ul) << this->StateShift) - ((ULONGLONG) 0x1ul);

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = this->NbrSite / this->MaxXMomentum;
  this->ComplementaryStateXShift = this->MaxMomentum - this->StateXShift;
  this->XMomentumMask = (((ULONGLONG) 0x1ul) << this->StateXShift) - ((ULONGLONG) 0x1ul);

  this->MaxYMomentum = maxYMomentum;
  this->NbrSitePerUnitCell = this->NbrSite / (this->MaxXMomentum * this->MaxYMomentum);
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = this->MaxXMomentum;
  this->StateYShift = (this->NbrSite / (this->MaxXMomentum * this->MaxYMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (((ULONGLONG) 0x1ul) << this->StateYShift) - ((ULONGLONG) 0x1ul);
  this->YMomentumBlockMask = (((ULONGLONG) 0x1ul) << this->YMomentumBlockSize) - ((ULONGLONG) 0x1ul);  
  this->YMomentumFullMask = ((ULONGLONG) 0x0ul);
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 
	
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
      //      cout <<" After generate states"<<endl;
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  this->StateMaxMomentum = new int [this->LargeHilbertSpaceDimension];  
	  int CurrentMaxMomentum = this->MaxMomentum;
	  while (((this->StateDescription[0] >> CurrentMaxMomentum) & ((ULONGLONG) 0x1ul) ) ==  ((ULONGLONG) 0x0ul))
	    --CurrentMaxMomentum;
	  this->StateMaxMomentum[0] = CurrentMaxMomentum;
	  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      while (((this->StateDescription[i] >> CurrentMaxMomentum) & ((ULONGLONG) 0x1ul) ) ==  ((ULONGLONG) 0x0ul))
		--CurrentMaxMomentum;
	      this->StateMaxMomentum[i] = CurrentMaxMomentum;
	    }
	    
	  this->GenerateLookUpTable(memory);
	  //	  cout <<" After Look up table generation"<<endl;	  
#ifdef __DEBUG__
	  long UsedMemory = 0;
	  UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
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
  this->SignLookUpTable = 0;
  this->NbrParticleLookUpTable = 0;
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong(const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong& bosons)
{
  this->NbrBosons = bosons.NbrBosons;  
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->NbrSite = bosons.NbrSite;

  this->MaxXMomentum = bosons.MaxXMomentum;
  this->XMomentum = bosons.XMomentum;
  this->StateXShift = bosons.StateXShift;
  this->ComplementaryStateXShift = bosons.ComplementaryStateXShift;
  this->XMomentumMask = bosons.XMomentumMask;
  this->MaxYMomentum = bosons.MaxYMomentum;
  this->YMomentum = bosons.YMomentum;
  this->NbrYMomentumBlocks = bosons.NbrYMomentumBlocks;
  this->StateYShift = bosons.StateYShift;
  this->YMomentumBlockSize = bosons.YMomentumBlockSize;
  this->ComplementaryStateYShift = bosons.ComplementaryStateYShift;
  this->YMomentumMask = bosons.YMomentumMask;
  this->YMomentumBlockMask = bosons.YMomentumBlockMask;  
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->ComplementaryYMomentumFullMask = bosons.ComplementaryYMomentumFullMask;
  this->YMomentumFullMask = bosons.YMomentumFullMask;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->MomentumModulo = bosons.MomentumModulo;
  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->NbrSitePerUnitCell = bosons.NbrSitePerUnitCell;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;

  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;

  this->NbrParticleLookUpTable = bosons.NbrParticleLookUpTable;

  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->Flag = bosons.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::~BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong& BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::operator = (const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrBosons = bosons.NbrBosons;  
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->NbrSite = bosons.NbrSite;

  this->MaxXMomentum = bosons.MaxXMomentum;
  this->XMomentum = bosons.XMomentum;
  this->StateXShift = bosons.StateXShift;
  this->ComplementaryStateXShift = bosons.ComplementaryStateXShift;
  this->XMomentumMask = bosons.XMomentumMask;
  this->MaxYMomentum = bosons.MaxYMomentum;
  this->YMomentum = bosons.YMomentum;
  this->NbrYMomentumBlocks = bosons.NbrYMomentumBlocks;
  this->StateYShift = bosons.StateYShift;
  this->YMomentumBlockSize = bosons.YMomentumBlockSize;
  this->ComplementaryStateYShift = bosons.ComplementaryStateYShift;
  this->YMomentumMask = bosons.YMomentumMask;
  this->YMomentumBlockMask = bosons.YMomentumBlockMask;  
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;

  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->MomentumModulo = bosons.MomentumModulo;
  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->NbrSitePerUnitCell = bosons.NbrSitePerUnitCell;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;

  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;

  this->NbrParticleLookUpTable = bosons.NbrParticleLookUpTable;

  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->Flag = bosons.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::Clone()
{
  return new BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong(*this);
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::GenerateStates()
{
  this->StateDescription = new  ULONGLONG [this->LargeHilbertSpaceDimension];
  long TmpHilbertSpaceDimenion = this->RawGenerateStates(this->NbrBosons, this->NbrSite - 1, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      //      this->PrintState(cout,i)<<endl;
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY) == this->StateDescription[i]))
	{
	  //	  cout <<"After FindCanonicalForm"<<endl;
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      //  cout <<"After TestMomentum"<<endl;
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = ((ULONGLONG) 0x0ul);
	    }
	}
      else
	{
	  this->StateDescription[i] = ((ULONGLONG) 0x0ul);
	}
    }
  if (TmpLargeHilbertSpaceDimension == 0l)
    return 0l;
  
  ULONGLONG * TmpStateDescription = new ULONGLONG [TmpLargeHilbertSpaceDimension];  
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != ((ULONGLONG) 0x0ul))
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
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::RawGenerateStates(int nbrBosons, int currentSite, long pos)
{
  if (nbrBosons == 0)
    {
      this->StateDescription[pos] = ((ULONGLONG) 0x0ul);	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrBosons < 0))
    return pos;

  if (nbrBosons == 1)
    {
      for (int j = currentSite; j >= 0; --j)
	{
	  this->StateDescription[pos] = (((ULONGLONG) 0x1ul) << j);
	  ++pos;
	}
      return pos;
    }
  long TmpPos = this->RawGenerateStates(nbrBosons - 1, currentSite - 1, pos);
  ULONGLONG Mask = (((ULONGLONG) 0x1ul) << currentSite);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->RawGenerateStates(nbrBosons, currentSite - 1, pos);
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

int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  ULONGLONG State = this->StateDescription[index];
  
  if ((State & (((ULONGLONG) 0x1ul) << n)) == ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }

  coefficient = 1.0;

  State &= ~(((ULONGLONG) 0x1ul) << n);
  if ((State & (((ULONGLONG) 0x1ul) << m))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  State |= (((ULONGLONG) 0x1ul) << m);
  this->ProdATemporaryNbrStateInOrbit =  this->NbrStateInOrbit[index];
  return this->SymmetrizeAdAdResult(State, coefficient, nbrTranslationX, nbrTranslationY);
}
  
// apply a^+_m1_sigma a^+_m2_sigma operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != ((ULONGLONG) 0x0ul)) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != ((ULONGLONG) 0x0ul)) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  
  coefficient = 1.0;
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  for (int i = 0; i < this->MaxMomentum; ++i)
    Str << (unsigned long)  ((TmpState >> i) & ((ULONGLONG) 0x1ul)) << " ";

  /*  if (this->FindStateIndex(this->StateDescription[state], this->StateMaxMomentum[state]) != state)
    {
      Str << "  error";
      }*/
  return Str;
}

// apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::A (int index, int n)
{
  this->ProdATemporaryStateMaxMomentum = this->StateMaxMomentum[index];
  this->ProdATemporaryState = this->StateDescription[index];
  if ((n >  this->ProdATemporaryStateMaxMomentum) || ((this->ProdATemporaryState  & (((ULONGLONG) 0x1ul) << n)) == ((ULONGLONG) 0x0ul)))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = 1.0;
  this->ProdATemporaryState  &= ~(((ULONGLONG) 0x1ul) << n);
  if (this->ProdATemporaryState == ((ULONGLONG) 0x0ul))
    {
      this->ProdATemporaryStateMaxMomentum = 0;
    }
  else
    {
      if (this->ProdATemporaryStateMaxMomentum == n)
	while ((this->ProdATemporaryState >> this->ProdATemporaryStateMaxMomentum) ==  ((ULONGLONG) 0x0ul))
	  --this->ProdATemporaryStateMaxMomentum;
    }
  return Coefficient;
}


// apply a^+_m operator to the state produced using AuAu method (without destroying it)
//
// m = first index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::Ad (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  if ((TmpState  & (((ULONGLONG) 0x1ul) << m))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  
  TmpState &= ~(((ULONGLONG) 0x1ul) << m);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}


// apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::AA (int index, int n1, int n2)
{
//  cout <<"using BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::A "<<endl;
  this->ProdATemporaryStateMaxMomentum = this->StateMaxMomentum[index];
  this->ProdATemporaryState = this->StateDescription[index];
  if ((n1 >  this->ProdATemporaryStateMaxMomentum) || (n2 >  this->ProdATemporaryStateMaxMomentum) ||  ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n1)) == ((ULONGLONG) 0x0ul)) || 
      ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n2)) == ((ULONGLONG) 0x0ul)) || (n1 == n2))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = 1.0;
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);

  
  if (this->ProdATemporaryState == ((ULONGLONG) 0x0ul))
    {
      this->ProdATemporaryStateMaxMomentum = 0;
    }
  else
    {
      if (this->ProdATemporaryStateMaxMomentum == n1)
	while ((this->ProdATemporaryState >> this->ProdATemporaryStateMaxMomentum) ==  ((ULONGLONG) 0x0ul))
	  --this->ProdATemporaryStateMaxMomentum;
    }
  return Coefficient;
}



// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// return value = Hilbert space dimension

long BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::EvaluateHilbertSpaceDimension(int nbrBosons)
{
  BinomialCoefficients binomials(this->NbrSite);
  long dimension = binomials(this->NbrSite, nbrBosons);
  return dimension;
}



// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrMomentum);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrMomentum)
    this->MaximumLookUpShift = this->NbrMomentum;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrMomentum];
  this->LookUpTableShift = new int [this->NbrMomentum];
  for (int i = 0; i < this->NbrMomentum; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentMaxMomentum = this->StateMaxMomentum[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
  if (CurrentMaxMomentum < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentMaxMomentum] = 0;
  else
    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  --CurrentMaxMomentum;
	  while (CurrentMaxMomentum > this->StateMaxMomentum[i])
	    {
	      this->LookUpTableShift[CurrentMaxMomentum] = -1;
	      --CurrentMaxMomentum;
	    }
 	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
	  if (CurrentMaxMomentum < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentMaxMomentum] = 0;
	  else
	    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;
  this->RescalingFactors = new double* [this->NbrMomentum];
  for (int i = 1; i <= this->MaxMomentum; ++i)
    {
      this->RescalingFactors[i] = new double [this->NbrMomentum];
      for (int j = 1; j <= this->MaxMomentum; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}
