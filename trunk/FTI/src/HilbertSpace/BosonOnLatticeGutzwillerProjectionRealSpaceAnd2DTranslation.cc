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
//                        last modification : 10/09/2014                      //
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
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
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

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation ()
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
  this->MomentumMask = ((unsigned long) 1);
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

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation (int nbrBosons, int nbrSite, int xMomentum, int  maxXMomentum,  int yMomentum, int maxYMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->StateXShift = 1;
  this->MomentumIncrement = (this->NbrBosons * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = this->NbrSite / this->MaxXMomentum;
  this->ComplementaryStateXShift = this->MaxMomentum - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->MaxYMomentum = maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrSitePerUnitCell = this->NbrSite /  (this->MaxYMomentum * this->MaxXMomentum);
  this->NbrYMomentumBlocks = this->MaxXMomentum;
  this->StateYShift = (this->NbrSite / (this->MaxXMomentum * this->MaxYMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
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
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  this->StateMaxMomentum = new int [this->LargeHilbertSpaceDimension];  
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
 this->SignLookUpTable = 0;
 this->NbrParticleLookUpTable = 0;
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation(const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation& bosons)
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

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::~BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation& BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::operator = (const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation& bosons)
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

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;
  this->NbrSitePerUnitCell = bosons.NbrSitePerUnitCell;

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

AbstractHilbertSpace* BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::Clone()
{
  return new BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation(*this);
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::GenerateStates()
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
	//cout <<  this->StateDescription[i] <<endl;
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

// apply a^+_m_sigma a_n_sigma operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator including the orbital and the spin index
// n = index of the annihilation operator including the orbital and the spin index
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{

//  cout <<"inside int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)"<<endl;
//  cout <<"creation operator = "<< m  <<endl;
//  cout <<"annihilation operator = "<< n  <<endl;
  unsigned long State = this->StateDescription[index];
//  cout <<State<<endl;
  if ((State & (0x1ul << n)) == 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  State &= ~(0x1ul << n);
//  cout <<State<<endl;
  if ((State & (0x1ul << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  State |= (0x1ul << m);
//  cout <<" m = " << m << " n = " <<n << " State = "<< this->StateDescription[index] <<" -> "<<"State = " << State<<endl;
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

int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  TmpState |= (0x1ul << m2);
  TmpState |= (0x1ul << m1);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->MaxMomentum; ++i)
    Str << ((TmpState >> i) & 0x1ul) << " ";
  cout <<TmpState<<endl;
   if (this->FindStateIndex(this->StateDescription[state], this->StateMaxMomentum[state]) != state)
     {
       Str << "  error";
     }
  return Str;
}

// apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::A (int index, int n)
{
//  cout <<"using BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::A "<<endl;
  this->ProdATemporaryStateMaxMomentum = this->StateMaxMomentum[index];
  this->ProdATemporaryState = this->StateDescription[index];
  if ((n >  this->ProdATemporaryStateMaxMomentum) || ((this->ProdATemporaryState & (0x1ul << n)) == 0x0ul))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = 1.0;
  this->ProdATemporaryState &= ~(0x1ul << n);
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdATemporaryStateMaxMomentum = 0;
    }
  else
    {
      if (this->ProdATemporaryStateMaxMomentum == n)
	while ((this->ProdATemporaryState >> this->ProdATemporaryStateMaxMomentum) == 0x0ul)
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

int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::Ad (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
//  cout <<"using BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::Ad "<<endl;
  unsigned long TmpState = this->ProdATemporaryState;
  if ((TmpState & (0x1ul << m)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  
  TmpState |= (0x1ul << m);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}


// apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::AA (int index, int n1, int n2)
{
//  cout <<"using BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::A "<<endl;
  this->ProdATemporaryStateMaxMomentum = this->StateMaxMomentum[index];
  this->ProdATemporaryState = this->StateDescription[index];
  if ((n1 >  this->ProdATemporaryStateMaxMomentum) || (n2 >  this->ProdATemporaryStateMaxMomentum) || ((this->ProdATemporaryState & (0x1ul << n1)) == 0x0ul) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0x0ul))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = 1.0;
  this->ProdATemporaryState &= ~(0x1ul << n2);
  this->ProdATemporaryState &= ~(0x1ul << n1);
  
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdATemporaryStateMaxMomentum = 0;
    }
  else
    {
      if (this->ProdATemporaryStateMaxMomentum == n1)
	while ((this->ProdATemporaryState >> this->ProdATemporaryStateMaxMomentum) == 0x0ul)
	  --this->ProdATemporaryStateMaxMomentum;
    }
  return Coefficient;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = kx sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
      if ((kxSector == this->XMomentum) && (kySector == this->YMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementaryKxMomentum = (this->XMomentum - kxSector);
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->MaxXMomentum;
  int ComplementaryKyMomentum = (this->YMomentum - kySector);
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->MaxYMomentum;
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation SubsytemSpace (nbrParticleSector, this->NbrSite, kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation ComplementarySpace (ComplementaryNbrParticles, this->NbrSite, ComplementaryKxMomentum, this->MaxXMomentum, 
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
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, 
														     ParticleOnTorusWithMagneticTranslations* complementaryHilbertSpace,  
														     ParticleOnTorusWithMagneticTranslations* destinationHilbertSpace,
														     ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  Complex** FourrierCoefficientsDestination = new Complex* [this->MomentumModulo];
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  long TmpNbrNonZeroElements = 0l;
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation* TmpDestinationHilbertSpace =  (BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation*) destinationHilbertSpace;
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation* TmpHilbertSpace = (BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation*) complementaryHilbertSpace;
  BosonOnLatticeGutzwillerProjectionRealSpace* TmpDestinationFullHilbertSpace = new BosonOnLatticeGutzwillerProjectionRealSpace(TmpDestinationHilbertSpace->NbrBosons,
																TmpDestinationHilbertSpace->NbrSite);
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficientsDestination[i] = new Complex [this->MaxYMomentum];
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficientsDestination[i][j] = Phase (2.0 * M_PI * ((double) (i * TmpDestinationHilbertSpace->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * TmpDestinationHilbertSpace->YMomentum) / ((double) this->MaxYMomentum)));
	  FourrierCoefficients[i][j] = Phase (2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
	}
    }
  
  ComplexMatrix TmpEntanglementMatrix (nbrIndex, TmpDestinationHilbertSpace->GetHilbertSpaceDimension(), true);
  int TmpMinIndex = minIndex;
  int MaxIndex = minIndex + nbrIndex;
  BinomialCoefficients TmpBinomial (this->NbrBosons);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrBosons, TmpDestinationHilbertSpace->NbrBosons));
  this->ProdATemporaryNbrStateInOrbit = 1;
  TmpDestinationHilbertSpace->ProdATemporaryNbrStateInOrbit = 1;
  
  Complex* TmpDestinationFactors = new Complex[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
  int* TmpDestinationRealIndices = new int[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
  for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
    {
      unsigned long TmpCanonicalState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
      int TmpDestinationNbrTranslationX;
      int TmpDestinationNbrTranslationY;
      double TmpDestinationCoefficient = 1.0;
      int RealDestinationIndex = TmpDestinationHilbertSpace->SymmetrizeAdAdResult(TmpCanonicalState2, TmpDestinationCoefficient, TmpDestinationNbrTranslationX, TmpDestinationNbrTranslationY);
      TmpDestinationRealIndices[j] = RealDestinationIndex;
      if (RealDestinationIndex < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
	{
	  TmpDestinationFactors[j] = (TmpDestinationCoefficient * TmpInvBinomial) * FourrierCoefficientsDestination[TmpDestinationNbrTranslationX][TmpDestinationNbrTranslationY];
	}
    }
  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
      double TmpRescalingFactor = sqrt((double) TmpHilbertSpace->NbrStateInOrbit[minIndex]);
      for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
	      if (TmpDestinationRealIndices[j] < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
		{
		  int TmpNbrTranslationX;
		  int TmpNbrTranslationY;
		  double TmpCoefficient = TmpRescalingFactor;
		  unsigned long TmpState3 = TmpState | TmpState2;
		  int TmpPos = this->SymmetrizeAdAdResult(TmpState3, TmpCoefficient, TmpNbrTranslationX, TmpNbrTranslationY);		  
		  if (TmpPos < this->HilbertSpaceDimension)
		    {      
		      Complex Coefficient = TmpCoefficient * TmpDestinationFactors[j] * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY];
		      TmpEntanglementMatrix[TmpDestinationRealIndices[j]][minIndex - TmpMinIndex] += groundState[TmpPos] * Coefficient;
		    }
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
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      delete[] FourrierCoefficientsDestination[i];
      delete[] FourrierCoefficients[i];
    }
  delete[] FourrierCoefficientsDestination;
  delete[] FourrierCoefficients;
  return TmpNbrNonZeroElements;
}

