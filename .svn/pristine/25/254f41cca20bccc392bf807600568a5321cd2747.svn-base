////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                    class of bosons on sphere with SU(4) spin               //
//                                                                            //
//                        last modification : 19/12/2011                      //
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
#include "HilbertSpace/BosonOnSphereWithSU4SpinAllEntanglement.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/StringTools.h"

#include <cmath>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;

// flag for switching testing output
//#define TEST_SU4_ALL_E

// default constructor
// 

BosonOnSphereWithSU4SpinAllEntanglement::BosonOnSphereWithSU4SpinAllEntanglement ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// totalSz = twice the total sz projection
// totalIsospin = twice the total isospin value (number imbalance between Plus and Minus)
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU4SpinAllEntanglement::BosonOnSphereWithSU4SpinAllEntanglement (int nbrBosons, int totalLz, int lzMax, int totalSz, int totalIsospin,
						    unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSz;
  this->TotalIsospin = totalIsospin;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();

  this->TemporaryStatePlus = new unsigned long[2*this->NbrLzValue];
  this->TemporaryStateMinus = new unsigned long[2*this->NbrLzValue];
  
  this->TemporaryStateDownMinus = this->TemporaryStateMinus;
  this->TemporaryStateUpMinus = this->TemporaryStateMinus + this->NbrLzValue;
  this->TemporaryStateDownPlus = this->TemporaryStatePlus;
  this->TemporaryStateUpPlus = this->TemporaryStatePlus + this->NbrLzValue;

  this->ProdATemporaryStatePlus = new unsigned long[2*this->NbrLzValue];
  this->ProdATemporaryStateMinus = new unsigned long[2*this->NbrLzValue];
  
  this->ProdATemporaryStateUpPlus = this->ProdATemporaryStatePlus + this->NbrLzValue;
  this->ProdATemporaryStateUpMinus = this->ProdATemporaryStateMinus + this->NbrLzValue;
  this->ProdATemporaryStateDownPlus = this->ProdATemporaryStatePlus;
  this->ProdATemporaryStateDownMinus = this->ProdATemporaryStateMinus;

  int NUp = this->NbrBosons + this->TotalSpin;
  int NDown = this->NbrBosons - this->TotalSpin;
  int NPlus = this->NbrBosons + this->TotalIsospin;
  int NMinus = this->NbrBosons - this->TotalIsospin;
  this->NPlusLzMax = 2*this->LzMax + NPlus;
  this->NMinusLzMax = 2*this->LzMax + NMinus;

  //cout << "NUp="<<NUp<<" NDown="<<NDown<<" NPlus="<<NPlus<<" NMinus="<<NMinus<<endl;

  if ((NUp < 0) || ((NUp & 0x1) != 0) || (NMinus < 0) || ((NMinus & 0x1) != 0) ||
      (NDown < 0) || ((NDown & 0x1) != 0) || (NPlus < 0) || ((NPlus & 0x1) != 0))
    {
      cout << "Error: TotalSz and TotalIsoSpin need to have the same parity as nbrBosons !"<<endl;
      this->LargeHilbertSpaceDimension = 0l;
      exit(-1);
    }
  else
    {
      NUp >>= 1;
      NMinus >>= 1;
      NDown >>= 1;
      NPlus >>= 1;
      this->LargeHilbertSpaceDimension = EvaluateHilbertSpaceDimension(NbrBosons, LzMax, NPlus, TotalLz, NUp);
    }
  this->StateDescriptionPlus =  new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionMinus =  new unsigned long [this->LargeHilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(NbrBosons, 4*LzMax+3, NPlus, (TotalLz + LzMax * NbrBosons) >> 1, NUp, 0l);
#ifdef TEST_SU4_ALL_E
  for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
    cout << i << " : " << hex << this->StateDescriptionPlus[i] << " " << this->StateDescriptionMinus[i] << dec << endl;  
#endif
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU4SpinAllEntanglement!" << endl;
      exit(1);
    }
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;
  this->GenerateLookUpTable(memory);
#ifdef TEST_SU4_ALL_E
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
    {
      cout << i << " : ";
      this->PrintState(cout, i);
      int Searched = this->FindStateIndex(this->StateDescriptionPlus[i], this->StateDescriptionMinus[i]);
      cout << Searched << endl;
      if (Searched!=i)
	cout << "Error finding state "<<i<<endl;
    }
#endif

#ifdef __DEBUG__
  long UsedMemory = LookUpMemorySize;
  UsedMemory += this->HilbertSpaceDimension * (2 * sizeof(unsigned long));
  cout << "memory requested for Hilbert space = ";
  PrintMemorySize (cout,UsedMemory);
  cout << ", including search tables of ";
  PrintMemorySize (cout,LookUpMemorySize)<<endl;
#endif
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereWithSU4SpinAllEntanglement::BosonOnSphereWithSU4SpinAllEntanglement(const BosonOnSphereWithSU4SpinAllEntanglement& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpPlusLzMax = bosons.LzMax;
  this->NUpMinusLzMax = bosons.LzMax;
  this->NDownPlusLzMax = bosons.LzMax;
  this->NDownMinusLzMax = bosons.LzMax;
  this->NPlusLzMax = bosons.NPlusLzMax;
  this->NMinusLzMax = bosons.NMinusLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->TotalSpin = bosons.TotalSpin;
  this->TotalIsospin = bosons.TotalIsospin;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;

  this->TemporaryStatePlus = new unsigned long[2*this->NbrLzValue];
  this->TemporaryStateMinus = new unsigned long[2*this->NbrLzValue];
  
  this->TemporaryStateUpPlus = this->TemporaryStatePlus + this->NbrLzValue;
  this->TemporaryStateUpMinus = this->TemporaryStateMinus + this->NbrLzValue;
  this->TemporaryStateDownPlus = this->TemporaryStatePlus;
  this->TemporaryStateDownMinus = this->TemporaryStateMinus;

  this->ProdATemporaryStatePlus = new unsigned long[2*this->NbrLzValue];
  this->ProdATemporaryStateMinus = new unsigned long[2*this->NbrLzValue];
  
  this->ProdATemporaryStateUpPlus = this->ProdATemporaryStatePlus + this->NbrLzValue;
  this->ProdATemporaryStateUpMinus = this->ProdATemporaryStateMinus + this->NbrLzValue;
  this->ProdATemporaryStateDownPlus = this->ProdATemporaryStatePlus;
  this->ProdATemporaryStateDownMinus = this->ProdATemporaryStateMinus;
    
  this->StateDescriptionPlus = bosons.StateDescriptionPlus;
  this->StateDescriptionMinus = bosons.StateDescriptionMinus;

  this->FullLookUp = bosons.FullLookUp;
  this->LookUpTablePlus = bosons.LookUpTablePlus;
  this->LookUpTableMinus = bosons.LookUpTableMinus;
  this->NbrUniqueStateDescriptionPlus = bosons.NbrUniqueStateDescriptionPlus;
  this->UniqueStateDescriptionPlus = bosons.UniqueStateDescriptionPlus;
  this->UniqueStateDescriptionSubArrayIndicesPlus = bosons.UniqueStateDescriptionSubArrayIndicesPlus;
  this->LookUpMemorySize = bosons.LookUpMemorySize;
}

// destructor
//

BosonOnSphereWithSU4SpinAllEntanglement::~BosonOnSphereWithSU4SpinAllEntanglement ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionPlus;
      delete[] this->StateDescriptionMinus;
      if (FullLookUp)
	{
	  delete [] LookUpTablePlus;
	  delete [] LookUpTableMinus;
	}
      else
	{
	  delete[] this->UniqueStateDescriptionPlus;
	  delete[] this->UniqueStateDescriptionSubArrayIndicesPlus;
	}
    }
  delete[] this->TemporaryStatePlus;
  delete[] this->TemporaryStateMinus;
  delete[] this->ProdATemporaryStatePlus;
  delete[] this->ProdATemporaryStateMinus;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU4SpinAllEntanglement& BosonOnSphereWithSU4SpinAllEntanglement::operator = (const BosonOnSphereWithSU4SpinAllEntanglement& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionPlus;
      delete[] this->StateDescriptionMinus;
      if (FullLookUp)
	{
	  delete [] LookUpTablePlus;
	  delete [] LookUpTableMinus;
	}
      else
	{
	  delete[] this->UniqueStateDescriptionPlus;
	  delete[] this->UniqueStateDescriptionSubArrayIndicesPlus;
	}
    }
  delete[] this->TemporaryStatePlus;
  delete[] this->TemporaryStateMinus;
  delete[] this->ProdATemporaryStatePlus;
  delete[] this->ProdATemporaryStateMinus;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->TotalIsospin = bosons.TotalIsospin;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpPlusLzMax = bosons.LzMax;
  this->NUpMinusLzMax = bosons.LzMax;
  this->NDownPlusLzMax = bosons.LzMax;
  this->NDownMinusLzMax = bosons.LzMax;
  this->NPlusLzMax = bosons.NPlusLzMax;
  this->NMinusLzMax = bosons.NMinusLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  
  this->TemporaryStatePlus = new unsigned long[2*this->NbrLzValue];
  this->TemporaryStateMinus = new unsigned long[2*this->NbrLzValue];
  
  this->TemporaryStateUpPlus = this->TemporaryStatePlus + this->NbrLzValue;
  this->TemporaryStateUpMinus = this->TemporaryStateMinus + this->NbrLzValue;
  this->TemporaryStateDownPlus = this->TemporaryStatePlus;
  this->TemporaryStateDownMinus = this->TemporaryStateMinus;

  this->ProdATemporaryStatePlus = new unsigned long[2*this->NbrLzValue];
  this->ProdATemporaryStateMinus = new unsigned long[2*this->NbrLzValue];
  
  this->ProdATemporaryStateUpPlus = this->ProdATemporaryStatePlus + this->NbrLzValue;
  this->ProdATemporaryStateUpMinus = this->ProdATemporaryStateMinus + this->NbrLzValue;
  this->ProdATemporaryStateDownPlus = this->ProdATemporaryStatePlus;
  this->ProdATemporaryStateDownMinus = this->ProdATemporaryStateMinus;

  this->StateDescriptionPlus = bosons.StateDescriptionPlus;
  this->StateDescriptionMinus = bosons.StateDescriptionMinus;
  this->FullLookUp = bosons.FullLookUp;
  this->LookUpTablePlus = bosons.LookUpTablePlus;
  this->LookUpTableMinus = bosons.LookUpTableMinus;
  this->NbrUniqueStateDescriptionPlus = bosons.NbrUniqueStateDescriptionPlus;
  this->UniqueStateDescriptionPlus = bosons.UniqueStateDescriptionPlus;
  this->UniqueStateDescriptionSubArrayIndicesPlus = bosons.UniqueStateDescriptionSubArrayIndicesPlus;
  this->LookUpMemorySize = bosons.LookUpMemorySize;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereWithSU4SpinAllEntanglement::Clone()
{
  return new BosonOnSphereWithSU4SpinAllEntanglement(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereWithSU4SpinAllEntanglement::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereWithSU4SpinAllEntanglement::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereWithSU4SpinAllEntanglement::ExtractSubspace (AbstractQuantumNumber& q, 
								 SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m_up a_m_up operator to a given state  (only spin up isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_um a_m_um

double  BosonOnSphereWithSU4SpinAllEntanglement::AdupAup (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  return (double) (this->TemporaryStateUpPlus[m]);
}

// apply a^+_m_um a_m_um operator to a given state  (only spin up isospin minus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_um a_m_um

double BosonOnSphereWithSU4SpinAllEntanglement::AdumAum (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  return (double) (this->TemporaryStateUpMinus[m]);  
}
 
// apply a^+_m_dp a_m_dp operator to a given state (only spin down isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dm a_m_dm

double BosonOnSphereWithSU4SpinAllEntanglement::AddpAdp (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  return (double) (this->TemporaryStateDownPlus[m]);
}

// apply a^+_m_dm a_m_dm operator to a given state (only spin down isospin minus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dm a_m_dm

double BosonOnSphereWithSU4SpinAllEntanglement::AddmAdm (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  return (double) (this->TemporaryStateDownMinus[m]);
}

/*
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  cout << "in AddmAdm [ -(" <<this->TemporaryStateMinus[0];
  for (int i = 1; i <2*NbrLzValue ; ++i)
    cout << ", "<<this->TemporaryStateMinus[i];
  cout << ") +("<<this->TemporaryStatePlus[0];
  for (int i = 1; i <2*NbrLzValue ; ++i)
    cout << ", "<<this->TemporaryStatePlus[i];
  cout <<")]"<<endl;

  for (int i = 0; i <NbrLzValue ; ++i)
    {
      if (this->TemporaryStateMinus[i]!=this->TemporaryStateDownMinus[i])
	cout << "Discrepancy in Down Minus "<<i<<endl;
      if (this->TemporaryStateMinus[i+NbrLzValue]!=this->TemporaryStateUpMinus[i])
	cout << "Discrepancy in Up Minus "<<i<<endl;
    }

  double ReturnVal = (double) (this->TemporaryStateDownMinus[m]);
  this->PrintState(cout,index)<< ": dm["<<m<<"]=" <<ReturnVal<<endl;
  return ReturnVal;
  //return (double) (this->TemporaryStateDownMinus[m]);  
}
*/

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSU4SpinAllEntanglement::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionPlus[state], this->StateDescriptionMinus[state], 
		       this->TemporaryStatePlus, this->TemporaryStateMinus);
  //  Str <<hex << this->StateDescriptionPlus[state] <<" "<< this->StateDescriptionMinus[state]<<dec;
  Str << " | ";
  for (int i = this->LzMax; i >=0 ; --i)
    {
      Str << "(" << this->TemporaryStateUpPlus[i] << "," << this->TemporaryStateUpMinus[i] << "," 
	  << this->TemporaryStateDownPlus[i] << "," << this->TemporaryStateDownMinus[i] << ") | ";
    }
  return Str;
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSU4SpinAllEntanglement::GenerateLookUpTable(unsigned long memory)
{

#ifdef TEST_SU4_ALL_E
  unsigned long LastPlusPart=StateDescriptionPlus[0];
  int BlockStartIndex=0;
  int BlockLength=1;
  int NbrBlocks=1;
  unsigned long LastMinusPart=StateDescriptionMinus[0];
  for (long i=1; i<LargeHilbertSpaceDimension; ++i)
    {
      if (StateDescriptionPlus[i]!=LastPlusPart)
	{
	  if (LastPlusPart < StateDescriptionPlus[i])
	    {
	      cout << "States not monotonously ordered in index 1 at block "<<NbrBlocks<<", line "<<i<<endl;
	    }
	  BlockStartIndex=i;
	  BlockLength=1;
	  LastPlusPart=StateDescriptionPlus[i];
	  ++NbrBlocks;
	}
      else
	{
	  if (LastMinusPart <= StateDescriptionMinus[i])
	    {
	      cout << "States not monotonously ordered in index 2 for block "<<NbrBlocks<<" starting at "<<BlockStartIndex<<endl;
	    }
	  ++BlockLength;
	}
      LastMinusPart=StateDescriptionMinus[i];
    }
#endif

  // set shift for higher bits
  int MaxLookUpSize = (NPlusLzMax>NMinusLzMax ? NPlusLzMax+1 : NMinusLzMax+1);
  if (MaxLookUpSize<std::log((double)memory/sizeof(unsigned long))/std::log(2.0))
    {
      cout << "Using full look-up table"<<endl;
      this->FullLookUp=true;
      int LookUpSizePlus = NPlusLzMax+1;
      int LookUpSizeMinus = NMinusLzMax+1;
      // look-up table for higher bits
      this->LookUpTablePlus = new unsigned long[0x1l<<LookUpSizePlus];
      // look-up table with two entries : the first one used lzmax value of the state an the second 
      this->LookUpTableMinus = new unsigned long[0x1l<<LookUpSizeMinus];
      this->LookUpMemorySize = ((0x1l<<LookUpSizePlus)+(0x1l<<LookUpSizeMinus))*sizeof(unsigned long);
      unsigned long LastPlusValue = this->StateDescriptionPlus[0];
      this->LookUpTablePlus[LastPlusValue]=0;
      this->LookUpTableMinus[this->StateDescriptionMinus[0]]=0;
      unsigned SectorCount=1;
      for (long i=1; i<this->LargeHilbertSpaceDimension; ++i)
	{
	  if (LastPlusValue!=this->StateDescriptionPlus[i])
	    {
	      LastPlusValue=this->StateDescriptionPlus[i];
	      this->LookUpTablePlus[LastPlusValue]=i;
	      this->LookUpTableMinus[this->StateDescriptionMinus[i]]=0;
	      SectorCount=1;
	    }
	  else
	    {
	      this->LookUpTableMinus[this->StateDescriptionMinus[i]]=SectorCount;
	      ++SectorCount;
	    }
	}
    }
  else
    {
      cout << "partial look-up not fully implemented, yet: using no look-up table"<<endl;
      this->FullLookUp=false;
      long TmpUniquePartition = 1l;
      for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescriptionPlus[i - 1] == this->StateDescriptionPlus[i]))
	    {
	      ++i;
	    }
	  if (i < this->LargeHilbertSpaceDimension)
	    ++TmpUniquePartition;
	}
      
      this->NbrUniqueStateDescriptionPlus = TmpUniquePartition;
      this->UniqueStateDescriptionPlus = new unsigned long [this->NbrUniqueStateDescriptionPlus];
      this->UniqueStateDescriptionSubArrayIndicesPlus = new int [this->NbrUniqueStateDescriptionPlus+1];
      this->LookUpMemorySize = this->NbrUniqueStateDescriptionPlus*sizeof(unsigned long ) + (this->NbrUniqueStateDescriptionPlus+1) * sizeof(int );
      TmpUniquePartition = 0l;
      this->UniqueStateDescriptionPlus[0l] = this->StateDescriptionPlus[0l];
      this->UniqueStateDescriptionSubArrayIndicesPlus[0l] = 0;
      for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescriptionPlus[i - 1] == this->StateDescriptionPlus[i]))
	    {
	      ++i;
	    }
	  if (i < this->LargeHilbertSpaceDimension)
	    {
	      ++TmpUniquePartition;
	      this->UniqueStateDescriptionPlus[TmpUniquePartition] = this->StateDescriptionPlus[i];
	      this->UniqueStateDescriptionSubArrayIndicesPlus[TmpUniquePartition] = i; 
	    }
	}
      this->UniqueStateDescriptionSubArrayIndicesPlus[this->NbrUniqueStateDescriptionPlus] = this->HilbertSpaceDimension;
    }
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons remaining to be placed
// lzMax = momentum maximum value for a boson in the state 
// currentLzMax = momentum maximum value for bosons that are still to be placed
//                (counting from 0 to lzMax for down and lzMax+1 to 2LzMax+1 for up)
// totalLz = momentum total value
// nPlus = remaining number of particles in plus / block I
// nUp = remaining number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU4SpinAllEntanglement::GenerateStates(int nbrBosons, int currentLzMax, int nPlus, int totalLz, int nUp, long pos)
{
  int CurrentBlock = currentLzMax/(this->NbrLzValue<<1);
  int CurrentMomentum = currentLzMax%(this->NbrLzValue<<1);
  int CurrentSpinUp = CurrentMomentum/(this->NbrLzValue);
  CurrentMomentum = CurrentMomentum%this->NbrLzValue;
  //cout << "BosonOnSphereWithSU4SpinAllEntanglement::GenerateStates(nBos="<< nbrBosons<<", currentLzMax="<<currentLzMax<<", nPlus="<<nPlus<<", totalLz="<< totalLz<<", nUp="<< nUp<<") - Block = "<<CurrentBlock<<"Momentum="<<CurrentMomentum<<", Spin = "<<CurrentSpinUp<<endl;
  if ((nbrBosons < 0) || (nPlus < 0 ) || (nUp<0) || (currentLzMax < 0) || (totalLz < 0) || ((currentLzMax<NbrLzValue)&&((nbrBosons * CurrentMomentum) < totalLz)))
    return pos;
  if (nbrBosons == 0)
    {
      if ((totalLz == 0)&&(nPlus==0)&&(nUp==0))
	{
	  this->StateDescriptionPlus[pos] = 0x0ul;
	  this->StateDescriptionMinus[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else return pos;
    }
  if (currentLzMax<NbrLzValue)
    {
      if (((nbrBosons * currentLzMax) == totalLz)&&(nPlus==0)&&(nUp==0))
	{
	  this->StateDescriptionPlus[pos] = 0x0ul;
	  this->StateDescriptionMinus[pos] = ((0x1ul<<nbrBosons)-1ul) << CurrentMomentum;
	  return (pos + 1l);
	}     
      if ((nbrBosons == 1)&&(nUp == 0)&&(nPlus == 0))
	{
	  if (currentLzMax >= totalLz)
	    {
	      this->StateDescriptionPlus[pos] = 0x0ul;
	      this->StateDescriptionMinus[pos] = 0x1ul<<totalLz;
	      return (pos + 1l);
	    }
	  else return pos;
	}
    }
  int ReducedCurrentLzMax = currentLzMax - 1;
  long TmpPos;
  int MaxToPlace = nbrBosons;
  int NumberShift = nbrBosons;
  int SpinShift = 0;
  unsigned long* CurrentStateDescription=this->StateDescriptionMinus;
  if (CurrentBlock==1)
    {
      MaxToPlace = nPlus;
      NumberShift = nPlus;
      CurrentStateDescription=this->StateDescriptionPlus;
    }
  if (CurrentSpinUp==1)
    SpinShift = NbrLzValue;
  unsigned long Mask;
  for (int ToPlace=MaxToPlace; ToPlace>=0; --ToPlace)
    {
      TmpPos = this->GenerateStates(nbrBosons-ToPlace, ReducedCurrentLzMax, nPlus-ToPlace*CurrentBlock,
				    totalLz-ToPlace*CurrentMomentum, nUp-ToPlace*CurrentSpinUp, pos);
      
      Mask = ((0x1ul << ToPlace) - 1ul) << (CurrentMomentum + NumberShift - ToPlace + SpinShift);
      for (; pos < TmpPos; ++pos)
	CurrentStateDescription[pos] |= Mask;
    }
  return pos;
}


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// totalNPlus = total number of particles in block I
// totalNUp = total number of particles with up spin
// return value = Hilbert space dimension

long BosonOnSphereWithSU4SpinAllEntanglement::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalNPlus, int totalLz, int totalNUp)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, 4*lzMax+3, totalNPlus, (totalLz + lzMax * nbrBosons) >> 1,totalNUp, 0);
}


// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons remaining to be placed
// lzMax = momentum maximum value for a boson in the state 
// currentLzMax = momentum maximum value for bosons that are still to be placed
//                (counting from 0 to lzMax for down and lzMax+1 to 2LzMax+1 for up)
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU4SpinAllEntanglement::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int currentLzMax, int nPlus, int totalLz, int nUp, int level)
{
  int CurrentBlock = currentLzMax/(this->NbrLzValue<<1);
  int CurrentMomentum = currentLzMax%(this->NbrLzValue<<1);
  int CurrentSpinUp = CurrentMomentum/(this->NbrLzValue);
  CurrentMomentum = CurrentMomentum%this->NbrLzValue;
  
//   for (int i=0; i<level; ++i) cout << "  ";
//   cout << "SEV(n="<<nbrBosons<<", lz="<<currentLzMax<<", nP="<<nPlus<<", totalLz="<<totalLz<<", nUp="<<nUp<<") - Block="<<CurrentBlock<<", Lz="<<CurrentMomentum<<" Sz="<<CurrentSpinUp<<endl;
  if ((nbrBosons < 0) || (nPlus < 0 ) || (nUp<0) || (currentLzMax < 0) || (totalLz < 0) || ((currentLzMax<NbrLzValue)&&((nbrBosons * CurrentMomentum) < totalLz)))
    {
//       for (int i=0; i<level; ++i) cout << "  ";
//       cout << "c0 - add 0"<<endl;
      return 0;
    }
  if (nbrBosons == 0)
    {
      if ((totalLz == 0)&&(nPlus==0)&&(nUp==0))
	{
//   	  for (int i=0; i<level; ++i) cout << "  ";
//   	  cout << "c1 - add 1"<<endl;
	  return 1;
	}
      else
	{
//   	  for (int i=0; i<level; ++i) cout << "  ";	  
//   	  cout << "c1 - add 0"<<endl;
	  return 0;
	}
    }
  if (currentLzMax<NbrLzValue)
    {
      if (((nbrBosons * currentLzMax) == totalLz)&&(nPlus==0)&&(nUp==0))
	{
// 	  for (int i=0; i<level; ++i) cout << "  ";
// 	  cout << "c2a - add 1"<<endl;
	  return 1;
	}
      if ((nbrBosons == 1)&&(nPlus==0)&&(nUp==0))
	{
	  if (currentLzMax >= totalLz)
	    {
//  	      for (int i=0; i<level; ++i) cout << "  ";
//  	      cout << "c2b - add 1"<<endl;
	      return 1;
	    }
	  else
	    {
//  	      for (int i=0; i<level; ++i) cout << "  ";
//  	      cout << "c2b - add 0"<<endl;
	      return 0;
	    }
	}
      else
	{
	  
	}
    }
  int ReducedCurrentLzMax = currentLzMax - 1;
  long TmpDim = 0;
  int MaxToPlace = nbrBosons;
  if (CurrentBlock==1) MaxToPlace = nPlus;
  for (int ToPlace=MaxToPlace; ToPlace>=0; --ToPlace)
    TmpDim += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons-ToPlace, ReducedCurrentLzMax, nPlus-ToPlace*CurrentBlock,
							 totalLz-ToPlace*CurrentMomentum, nUp-ToPlace*CurrentSpinUp, level+1);
//    for (int i=0; i<level; ++i) cout << "  ";
//    cout << "add: "<<TmpDim<<endl;
  return TmpDim;
}




// apply a^+_m_up a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AdupAup (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  if (this->TemporaryStateUpPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpPlus[n];
  --this->TemporaryStateUpPlus[n];
  ++this->TemporaryStateUpPlus[m];
  coefficient *= (double) this->TemporaryStateUpPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->StateDescriptionMinus[index]);
}

// apply a^+_m_up a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AdupAum (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  if (this->TemporaryStateUpMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpMinus[n];
  --this->TemporaryStateUpMinus[n];
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  ++this->TemporaryStateUpPlus[m];
  coefficient *= (double) this->TemporaryStateUpPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->BosonToFermion(this->TemporaryStateMinus)); 
}

// apply a^+_m_up a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AdupAdp (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  if (this->TemporaryStateDownPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownPlus[n];
  --this->TemporaryStateDownPlus[n];
  ++this->TemporaryStateUpPlus[m];
  coefficient *= (double) this->TemporaryStateUpPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->StateDescriptionMinus[index]);
}

// apply a^+_m_up a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AdupAdm (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  if (this->TemporaryStateDownMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownMinus[n];
  --this->TemporaryStateDownMinus[n];
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  ++this->TemporaryStateUpPlus[m];
  coefficient *= (double) this->TemporaryStateUpPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a^+_m_um a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AdumAup (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  if (this->TemporaryStateUpPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpPlus[n];
  --this->TemporaryStateUpPlus[n];
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  ++this->TemporaryStateUpMinus[m];
  coefficient *= (double) this->TemporaryStateUpMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a^+_m_um a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AdumAum (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  if (this->TemporaryStateUpMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpMinus[n];
  --this->TemporaryStateUpMinus[n];
  ++this->TemporaryStateUpMinus[m];
  coefficient *= (double) this->TemporaryStateUpMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionPlus[index], this->BosonToFermion(this->TemporaryStateMinus)); 
}

// apply a^+_m_um a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AdumAdp (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  if (this->TemporaryStateDownPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownPlus[n];
  --this->TemporaryStateDownPlus[n];
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  ++this->TemporaryStateUpMinus[m];
  coefficient *= (double) this->TemporaryStateUpMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a^+_m_um a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AdumAdm (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  if (this->TemporaryStateDownMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownMinus[n];
  --this->TemporaryStateDownMinus[n];
  ++this->TemporaryStateUpMinus[m];
  coefficient *= (double) this->TemporaryStateUpMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionPlus[index], this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a^+_m_dp a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AddpAup (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  if (this->TemporaryStateUpPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpPlus[n];
  --this->TemporaryStateUpPlus[n];
  ++this->TemporaryStateDownPlus[m];
  coefficient *= (double) this->TemporaryStateDownPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->StateDescriptionMinus[index]);
}

// apply a^+_m_dp a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AddpAum (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  if (this->TemporaryStateUpMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpMinus[n];
  --this->TemporaryStateUpMinus[n];
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  ++this->TemporaryStateDownPlus[m];
  coefficient *= (double) this->TemporaryStateDownPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a^+_m_dp a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AddpAdp (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  if (this->TemporaryStateDownPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownPlus[n];
  --this->TemporaryStateDownPlus[n];
  ++this->TemporaryStateDownPlus[m];
  coefficient *= (double) this->TemporaryStateDownPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->StateDescriptionMinus[index]);  
}

// apply a^+_m_dp a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AddpAdm (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  if (this->TemporaryStateDownMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownMinus[n];
  --this->TemporaryStateDownMinus[n];
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  ++this->TemporaryStateDownPlus[m];
  coefficient *= (double) this->TemporaryStateDownPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a^+_m_dm a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AddmAup (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  if (this->TemporaryStateUpPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpPlus[n];
  --this->TemporaryStateUpPlus[n];
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  ++this->TemporaryStateDownMinus[m];
  coefficient *= (double) this->TemporaryStateDownMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a^+_m_dm a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AddmAum (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  if (this->TemporaryStateUpMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpMinus[n];
  --this->TemporaryStateUpMinus[n];
  ++this->TemporaryStateDownMinus[m];
  coefficient *= (double) this->TemporaryStateDownMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionPlus[index], this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a^+_m_dm a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AddmAdp (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->TemporaryStatePlus);
  if (this->TemporaryStateDownPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  coefficient = (double) this->TemporaryStateDownPlus[n];
  --this->TemporaryStateDownPlus[n];
  ++this->TemporaryStateDownMinus[m];
  coefficient *= (double) this->TemporaryStateDownMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStatePlus), this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a^+_m_dm a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinAllEntanglement::AddmAdm (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->TemporaryStateMinus);
  if (this->TemporaryStateDownMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownMinus[n];
  --this->TemporaryStateDownPlus[n];
  ++this->TemporaryStateDownPlus[m];
  coefficient *= (double) this->TemporaryStateDownMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionPlus[index], this->BosonToFermion(this->TemporaryStateMinus));
}

// apply a_n1_up a_n2_up operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AupAup (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateUpPlus[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUpPlus[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  double Coefficient = this->ProdATemporaryStateUpPlus[n2];
  --this->ProdATemporaryStateUpPlus[n2];
  Coefficient *= this->ProdATemporaryStateUpPlus[n1];
  --this->ProdATemporaryStateUpPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_up a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AupAum (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  if (this->ProdATemporaryStateUpPlus[n1] == 0)
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  if (this->ProdATemporaryStateUpMinus[n2] == 0)
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryStateUpMinus[n2];
  --this->ProdATemporaryStateUpMinus[n2];
  Coefficient *= this->ProdATemporaryStateUpPlus[n1];
  --this->ProdATemporaryStateUpPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_up a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AupAdp (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateDownPlus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  double Coefficient = this->ProdATemporaryStateDownPlus[n2];
  --this->ProdATemporaryStateDownPlus[n2];
  Coefficient *= this->ProdATemporaryStateUpPlus[n1];
  --this->ProdATemporaryStateUpPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_up a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AupAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  if (this->ProdATemporaryStateUpPlus[n1] == 0)
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  if (this->ProdATemporaryStateDownMinus[n2] == 0)
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownMinus[n2];
  Coefficient *= this->ProdATemporaryStateUpPlus[n1];
  --this->ProdATemporaryStateUpPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_um a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AumAum (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  if ((this->ProdATemporaryStateUpMinus[n1] == 0) || (this->ProdATemporaryStateUpMinus[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateUpMinus[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  double Coefficient = this->ProdATemporaryStateUpMinus[n2];
  --this->ProdATemporaryStateUpMinus[n2];
  Coefficient *= this->ProdATemporaryStateUpMinus[n1];
  --this->ProdATemporaryStateUpMinus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_um a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AumAdp (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  if (this->ProdATemporaryStateUpMinus[n1] == 0)
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  if (this->ProdATemporaryStateDownPlus[n2] == 0)
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryStateDownPlus[n2];
  --this->ProdATemporaryStateDownPlus[n2];
  Coefficient *= this->ProdATemporaryStateUpMinus[n1];
  --this->ProdATemporaryStateUpMinus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_um a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AumAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  if ((this->ProdATemporaryStateUpMinus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownMinus[n2];
  Coefficient *= this->ProdATemporaryStateUpMinus[n1];
  --this->ProdATemporaryStateUpMinus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_dp a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AdpAdp (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  if ((this->ProdATemporaryStateDownPlus[n1] == 0) || (this->ProdATemporaryStateDownPlus[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDownPlus[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  double Coefficient = this->ProdATemporaryStateDownPlus[n2];
  --this->ProdATemporaryStateDownPlus[n2];
  Coefficient *= this->ProdATemporaryStateDownPlus[n1];
  --this->ProdATemporaryStateDownPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_dp a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AdpAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  if (this->ProdATemporaryStateDownPlus[n1] == 0) 
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  if (this->ProdATemporaryStateDownMinus[n2] == 0)
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownMinus[n2];
  Coefficient *= this->ProdATemporaryStateDownPlus[n1];
  --this->ProdATemporaryStateDownPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_dm a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinAllEntanglement::AdmAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionMinus[index], this->NMinusLzMax, this->ProdATemporaryStateMinus);
  if ((this->ProdATemporaryStateDownMinus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDownMinus[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionPlus[index], this->NPlusLzMax, this->ProdATemporaryStatePlus);
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownMinus[n2];
  Coefficient *= this->ProdATemporaryStateDownMinus[n1];
  --this->ProdATemporaryStateDownMinus[n1];
  return sqrt(Coefficient);
}

// apply a^+_m1_up a^+_m2_up operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AdupAdup (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateUpPlus, coefficient);
}

// apply a^+_m1_up a^+_m2_um operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AdupAdum (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, coefficient);
}

// apply a^+_m1_up a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AdupAddp (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateDownPlus, coefficient);
}

// apply a^+_m1_up a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AdupAddm (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateDownMinus, coefficient);
}

// apply a^+_m1_um a^+_m2_um operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AdumAdum (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpMinus, this->TemporaryStateUpMinus, coefficient);
}

// apply a^+_m1_um a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AdumAddp (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpMinus, this->TemporaryStateDownPlus, coefficient);
}

// apply a^+_m1_um a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AdumAddm (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpMinus, this->TemporaryStateDownMinus, coefficient);
}

// apply a^+_m1_dp a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AddpAddp (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDownPlus, this->TemporaryStateDownPlus, coefficient);
}

// apply a^+_m1_dp a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AddpAddm (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDownPlus, this->TemporaryStateDownMinus, coefficient);
}

// apply a^+_m1_dm a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinAllEntanglement::AddmAddm (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDownMinus, this->TemporaryStateDownMinus, coefficient);
}


/*
// convert a state from one SU(4) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector

void BosonOnSphereWithSU4SpinAllEntanglement::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU4Indices = new int [this->NbrBosons];
  int* TmpSU4Indices2 = new int [this->NbrBosons];
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  targetState.ClearVector();
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUpPlus[i], this->StateDescriptionUpMinus[i], 
			   this->StateDescriptionDownPlus[i], this->StateDescriptionDownMinus[i],
			   this->TemporaryStateUpPlus, this->TemporaryStateUpMinus,
			   this->TemporaryStateDownPlus, this->TemporaryStateDownMinus); 
      double OccupationCoefficient = 0.0;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (unsigned l = 0; l < this->TemporaryStateDownMinus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 3;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateDownPlus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 2;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateUpMinus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateUpPlus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUpPlus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUpMinus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDownPlus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDownMinus[j]];
	}
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSU4Indices, TmpSU4Indices2, oneBodyBasis, OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] OccupationCoefficientArray;
  delete[] TmpMomentumIndices;
  delete[] TmpSU4Indices;
  delete[] TmpSU4Indices2;
}

// compute the transformation matrix from one SU(4) basis to another, transforming the one body basis in each momentum sector
//
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// return value = transformation matrix

ComplexMatrix BosonOnSphereWithSU4SpinAllEntanglement::TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU4Indices = new int [this->NbrBosons];
  int* TmpSU4Indices2 = new int [this->NbrBosons];
  ComplexMatrix TmpMatrix(this->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUpPlus[i], this->StateDescriptionUpMinus[i], 
			   this->StateDescriptionDownPlus[i], this->StateDescriptionDownMinus[i],
			   this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, 
			   this->TemporaryStateDownPlus, this->TemporaryStateDownMinus); 
      int TmpIndex = 0;
      double OccupationCoefficient = 0.0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (unsigned l = 0; l < this->TemporaryStateDownMinus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 3;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateDownPlus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 2;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateUpMinus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateUpPlus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUpPlus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUpMinus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDownPlus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDownMinus[j]];
	}
      this->TransformOneBodyBasisRecursive(TmpMatrix[i], 1.0, 0, TmpMomentumIndices, TmpSU4Indices, TmpSU4Indices2, oneBodyBasis,
					   OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSU4Indices;
  delete[] TmpSU4Indices2;
  delete[] OccupationCoefficientArray;
  return TmpMatrix;
}

// recursive part of the convertion from a state from one SU(4) basis to another, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSU4Indices = array that gives the spin dressing the initial n-body state
// currentSU4Indices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
// occupationCoefficientArray = array that provides 1/2 ln (N!)

void BosonOnSphereWithSU4SpinAllEntanglement::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
							      int position, int* momentumIndices, int* initialSU4Indices, int* currentSU4Indices, ComplexMatrix* oneBodyBasis,
							      double occupationCoefficient, double* occupationCoefficientArray) 
{
  if (position == this->NbrBosons)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->TemporaryStateUpPlus[i] = 0ul;
	  this->TemporaryStateUpMinus[i] = 0ul; 
	  this->TemporaryStateDownPlus[i] = 0ul;
	  this->TemporaryStateDownMinus[i] = 0ul;
	}
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  switch (currentSU4Indices[i])
	    {
	    case 0:
	      this->TemporaryStateUpPlus[momentumIndices[i]]++;
	      break;
	    case 1:
	      this->TemporaryStateUpMinus[momentumIndices[i]]++;
	      break;
	    case 2:
	      this->TemporaryStateDownPlus[momentumIndices[i]]++;
	      break;
	    case 3:
	      this->TemporaryStateDownMinus[momentumIndices[i]]++;
	      break;
	    }
	}
      int Index = this->FindStateIndex(this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, 
				       this->TemporaryStateDownPlus, this->TemporaryStateDownMinus);
      if (Index < this->HilbertSpaceDimension)
	{
	  
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateUpPlus[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateUpMinus[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateDownPlus[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateDownMinus[i]];
	    }
	  targetState[Index] += coefficient * exp (occupationCoefficient);
	}
      return;      
    }
  else
    {
      currentSU4Indices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][3 - initialSU4Indices[position]][2]), position + 1, momentumIndices, initialSU4Indices, currentSU4Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU4Indices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][3 - initialSU4Indices[position]][1]), position + 1, momentumIndices, initialSU4Indices, currentSU4Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU4Indices[position] = 2;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][3 - initialSU4Indices[position]][0]), position + 1, momentumIndices, initialSU4Indices, currentSU4Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU4Indices[position] = 3;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][3 - initialSU4Indices[position]][0]), position + 1, momentumIndices, initialSU4Indices, currentSU4Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
    }
}

// compute the projection matrix from the SU(4) Hilbert space to an U(1) Hilbert space
// 
// targetSpace = pointer to the U(1) Hilbert space
// type = type of particles that has to be kept (0 for type up-plus, 1 for type up-minus, 2 for type down-plus, 3 for type down-minus)
// return value = projection matrix

ComplexMatrix BosonOnSphereWithSU4SpinAllEntanglement::TransformationMatrixSU4ToU1(BosonOnSphereShort* targetSpace, int type)
{
  ComplexMatrix TmpMatrix (targetSpace->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  unsigned long* TmpStateDescription=NULL;
  unsigned long* TmpStateDescriptionOther1=NULL;
  unsigned long* TmpStateDescriptionOther2=NULL;
  unsigned long* TmpStateDescriptionOther3=NULL;
  switch (type)
    {
    case 0:
      {
	TmpStateDescription = this->StateDescriptionUpPlus;
	TmpStateDescriptionOther1 = this->StateDescriptionUpMinus;
	TmpStateDescriptionOther2 = this->StateDescriptionDownPlus;
	TmpStateDescriptionOther3 = this->StateDescriptionDownMinus;
      }
      break;
    case 1:
      {
	TmpStateDescription = this->StateDescriptionUpMinus;
	TmpStateDescriptionOther1 = this->StateDescriptionUpPlus;
	TmpStateDescriptionOther2 = this->StateDescriptionDownPlus;
	TmpStateDescriptionOther3 = this->StateDescriptionDownMinus;
      }
      break;
    case 2:
      {
	TmpStateDescription = this->StateDescriptionDownPlus;
	TmpStateDescriptionOther1 = this->StateDescriptionUpMinus;
	TmpStateDescriptionOther2 = this->StateDescriptionUpPlus;
	TmpStateDescriptionOther3 = this->StateDescriptionDownMinus;
      }
      break;
    case 3:
      {
	TmpStateDescription = this->StateDescriptionDownMinus;
	TmpStateDescriptionOther1 = this->StateDescriptionUpMinus;
	TmpStateDescriptionOther2 = this->StateDescriptionUpPlus;
	TmpStateDescriptionOther3 = this->StateDescriptionDownPlus;
      }
      break;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if ((TmpStateDescriptionOther1[i] == 0x0ul) && (TmpStateDescriptionOther2[i] == 0x0ul) && 
	  (TmpStateDescriptionOther3[i] == 0x0ul))
	{
	  unsigned long TmpState = TmpStateDescription[i];
	  int TmpLzMax = this->FermionicLzMax;
	  while ((TmpState >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  int Index = targetSpace->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (Index < targetSpace->HilbertSpaceDimension)
	    {
	      TmpMatrix[i][Index] = 1.0;
	    }
	}
    }
  return TmpMatrix;
}
*/
