////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere including two                //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 19/05/2009                      //
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
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/StringTools.h"


#include <cmath>
#include <bitset>
#include <cstdlib>
#include <algorithm>
#include <map>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;
using std::map;
using std::pair;

#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif


// default constructor
//

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels()
{
}

// basic constructor with contraint on the number of particles per Landau level component
// 
// nbrFermionsUp = number of fermions in level N=1
// nbrFermionsDown = number of fermions in level N=0
// totalLz = twice the momentum total value
// lzMaxUp = twice the maximum Lz value reached by a fermion with a spin up
// lzMaxDown = twice the maximum Lz value reached by a fermion with a spin down
// memory = amount of memory granted for precalculations

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels (int nbrFermionsUp, int nbrFermionsDown, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory)
{
  cout << "Constructor with fixed nbr of particles per LL"<<endl;
  this->NbrFermions = nbrFermionsUp + nbrFermionsDown;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->NbrFermionsUp = nbrFermionsUp;
  this->NbrFermionsDown = nbrFermionsDown;
  this->LzMaxUp = lzMaxUp;
  this->LzMaxDown = lzMaxDown;
  if (this->LzMaxUp >= this->LzMaxDown)
    {
      this->LzMax = this->LzMaxUp;
      this->LzShiftUp = 0;
      this->LzShiftDown = (this->LzMaxUp - this->LzMaxDown) >> 1;
    }
  else
    {
      this->LzMax = this->LzMaxDown;
      this->LzShiftDown = 0;
      this->LzShiftUp = (this->LzMaxDown - this->LzMaxUp) >> 1;
    }
  this->LzTotalShift = this->LzMaxDown + this->LzMaxUp;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  /*
  //For testing purposes: this constructs the full Hilbert space and then picks configurations with a given number of particles in N=0 and N=1

  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateFullHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  int TmpDimension = this->GenerateFullStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0);
  if (TmpDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count and State Generation in FermionOnSphereTwoLandauLevels! " << this->HilbertSpaceDimension << " " << TmpDimension  << endl;
      for (int i = 0; i < TmpDimension; ++i)
           this->PrintState(cout, i) << endl;
       exit(1);
    }

   cout << "Full HS dim= "<<this->HilbertSpaceDimension<<endl;
   int NewHilbertSpaceDimension = 0;
   int TmpNbrUp, TmpNbrDown;
   unsigned long Tmp;
   unsigned long* FlagArray = new unsigned long [this->HilbertSpaceDimension];
   unsigned long* NewStateDescription = new unsigned long [this->HilbertSpaceDimension];
   int* NewStateHighestBit = new int [this->HilbertSpaceDimension];
   
   for (int i = 0; i < this->HilbertSpaceDimension; i++)
      {
         TmpNbrUp = 0;
         TmpNbrDown = 0;
         for (int j = 0; j <= this->LzMax; j++)
            {
              Tmp = ((this->StateDescription[i] >> (j << 1)) & ((unsigned long) 0x3));
              if ((Tmp == 0x1l) || (Tmp == 0x3l))
                TmpNbrDown++;
               if ((Tmp == 0x2l) || (Tmp == 0x3l))
                TmpNbrUp++;
            }

         if ((TmpNbrUp == this->NbrFermionsUp) && (TmpNbrDown == this->NbrFermionsDown))
            {
               FlagArray[i] = 1;
               NewStateDescription[i] = StateDescription[i];
               NewStateHighestBit[i] = StateHighestBit[i];
               NewHilbertSpaceDimension++;
            }  
          else
            {
               FlagArray[i] = 0; 
               NewStateDescription[i] = 0;
               NewStateHighestBit[i] = 0;
            } 
      }

    delete[] this->StateHighestBit;
    delete[] this->StateDescription;


    this->StateDescription = new unsigned long [NewHilbertSpaceDimension];
    this->StateHighestBit = new int [NewHilbertSpaceDimension];
    int counter = 0;  
    for (int i = 0; i < this->HilbertSpaceDimension; i++)
      if (FlagArray[i] > 0)
          {
            this->StateDescription[counter] = NewStateDescription[i];
            this->StateHighestBit[counter] = NewStateHighestBit[i];
            counter++;  
          } 

    this->HilbertSpaceDimension = NewHilbertSpaceDimension;
    cout << "Reduced HilbertSpaceDimension= "<<this->HilbertSpaceDimension<<endl;
*/
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermionsUp, this->NbrFermionsDown, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  cout << "HS dim= "<<this->HilbertSpaceDimension<<endl;
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermionsUp, this->NbrFermionsDown, this->LzMaxUp, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0);
  cout << "HS dim check "<<this->HilbertSpaceDimension<<endl;
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

/*
  this->NbrFermions = nbrFermionsUp + nbrFermionsDown;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->NbrFermionsUp = nbrFermionsUp;
  this->NbrFermionsDown = nbrFermionsDown;
  this->LzMaxUp = lzMaxUp;
  this->LzMaxDown = lzMaxDown;
  if (this->LzMaxUp >= this->LzMaxDown)
    {
      this->LzMax = this->LzMaxUp;
      this->LzShiftUp = 0;
      this->LzShiftDown = (this->LzMaxUp - this->LzMaxDown) >> 1;
    }
  else
    {
      this->LzMax = this->LzMaxDown;
      this->LzShiftDown = 0;
      this->LzShiftUp = (this->LzMaxDown - this->LzMaxUp) >> 1;
    }
  this->LzTotalShift = this->LzMaxDown + this->LzMaxUp;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermionsUp, this->NbrFermionsDown, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  cout << this->LargeHilbertSpaceDimension<<endl;

  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  int TmpDimension = this->GenerateStates(this->NbrFermionsUp, this->NbrFermionsDown, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0);
  cout << "TmpDim= " << TmpDimension << endl;
  for (int i = 0; i < TmpDimension; ++i)
    this->PrintState(cout, i) << endl;
  exit(2);

  if (TmpDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count and State Generation in FermionOnSphereTwoLandauLevels! " << this->HilbertSpaceDimension << " " << TmpDimension  << endl;
  for (int i = 0; i < TmpDimension; ++i)
    this->PrintState(cout, i) << endl;
       exit(1);
    }

  //  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermionsUp, this->NbrFermionsDown, this->LzMaxUp, this->LzMaxDown, );
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
*/
}

// basic constructor with no contraint on the number of particles per spin component
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMaxUp = twice the maximum Lz value reached by a fermion with a spin up
// lzMaxDown = twice the maximum Lz value reached by a fermion with a spin down
// memory = amount of memory granted for precalculations

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels (int nbrFermions, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->LzMaxUp = lzMaxUp;
  this->LzMaxDown = lzMaxDown;
  if (this->LzMaxUp >= this->LzMaxDown)
    {
      this->LzMax = this->LzMaxUp;
      this->LzShiftUp = 0;
      this->LzShiftDown = (this->LzMaxUp - this->LzMaxDown) >> 1;
    }
  else
    {
      this->LzMax = this->LzMaxDown;
      this->LzShiftDown = 0;
      this->LzShiftUp = (this->LzMaxDown - this->LzMaxUp) >> 1;
    }
  this->LzTotalShift = this->LzMaxDown + this->LzMaxUp;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateFullHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  int TmpDimension = this->GenerateFullStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0);
  if (TmpDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count and State Generation in FermionOnSphereTwoLandauLevels! " << this->HilbertSpaceDimension << " " << TmpDimension  << endl;
  for (int i = 0; i < TmpDimension; ++i)
    this->PrintState(cout, i) << endl;
       exit(1);
    }

  //  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermionsUp, this->NbrFermionsDown, this->LzMaxUp, this->LzMaxDown, );
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels(const FermionOnSphereTwoLandauLevels& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->LzMaxUp = fermions.LzMaxUp;
  this->LzMaxDown = fermions.LzMaxDown;
  this->LzShiftUp = fermions.LzShiftUp;
  this->LzShiftDown = fermions.LzShiftDown;
  this->LzTotalShift = fermions.LzTotalShift;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
}

// destructor
//

FermionOnSphereTwoLandauLevels::~FermionOnSphereTwoLandauLevels ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereTwoLandauLevels& FermionOnSphereTwoLandauLevels::operator = (const FermionOnSphereTwoLandauLevels& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->LzMaxUp = fermions.LzMaxUp;
  this->LzMaxDown = fermions.LzMaxDown;
  this->LzShiftUp = fermions.LzShiftUp;
  this->LzShiftDown = fermions.LzShiftDown;
  this->LzTotalShift = fermions.LzTotalShift;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereTwoLandauLevels::Clone()
{
  return new FermionOnSphereTwoLandauLevels(*this);
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereTwoLandauLevels::AduAu (int index, int m)
{
  if ((this->StateDescription[index] & (0x2l << ((m + this->LzShiftUp) << 1))) != 0)
    return 1.0;
  else
    return 0.0;
}

// apply a^+_d_m a_d_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_d_m a_d_m

double FermionOnSphereTwoLandauLevels::AddAd (int index, int m)
{
  if ((this->StateDescription[index] & (0x1l << ((m + this->LzShiftDown) << 1))) != 0)
    return 1.0;
  else
    return 0.0;
}



// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereTwoLandauLevels::AduAu (int index, int m, int n, double& coefficient)
{  
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m += this->LzShiftUp;
  n += this->LzShiftUp;
  m = (m<<1) + 1;
  n = (n<<1) + 1;  
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(((unsigned long) (0x1ul)) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (((unsigned long) (0x1ul)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (((unsigned long) (0x1ul)) << m);
  return this->FindStateIndex(State, NewLargestBit);
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereTwoLandauLevels::AddAd (int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m += this->LzShiftDown;
  n += this->LzShiftDown;
  m <<= 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(((unsigned long) (0x1ul)) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (((unsigned long) (0x1ul)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (((unsigned long) (0x1ul)) << m);
  return this->FindStateIndex(State, NewLargestBit);
}


// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereTwoLandauLevels::AduAd (int index, int m, int n, double& coefficient)
  {
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m += this->LzShiftUp;
  n += this->LzShiftDown;
  m = (m<<1) + 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(((unsigned long) (0x1ul)) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (((unsigned long) (0x1ul)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (((unsigned long) (0x1ul)) << m);
  return this->FindStateIndex(State, NewLargestBit);
}



// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereTwoLandauLevels::AddAu (int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m += this->LzShiftDown;
  n += this->LzShiftUp;
  m <<= 1;
  n = (n<<1) + 1;  
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(((unsigned long) (0x1ul)) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (((unsigned long) (0x1ul)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (((unsigned long) (0x1ul)) << m);
  return this->FindStateIndex(State, NewLargestBit);
}


// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::AuAu (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 += this->LzShiftUp;
  n2 += this->LzShiftUp;
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  ++n2;
  if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::AdAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 += this->LzShiftDown;
  n2 += this->LzShiftDown;
  n1 <<= 1;
  n2 <<= 1;
  if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::AuAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 += this->LzShiftUp;
  n2 += this->LzShiftDown;
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ( ((this->ProdATemporaryState >> this->ProdALzMax) == 0) && (this->ProdALzMax>0))
    --this->ProdALzMax;
  return Coefficient;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::AduAdu (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 += this->LzShiftUp;
  m2 += this->LzShiftUp;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  ++m2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::AddAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 += this->LzShiftDown;
  m2 += this->LzShiftDown;
  m1 <<= 1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::AduAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 += this->LzShiftUp;
  m2 += this->LzShiftDown;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	Index = ((n[i] + this->LzShiftDown) << 1) + spinIndices[i];
      else
 	Index = ((n[i] + this->LzShiftUp) << 1) + spinIndices[i];
     if ((this->ProdATemporaryState & (0x1l << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      this->ProdATemporaryState &= ~(0x1l << Index);
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;

  return Coefficient;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::ProdA (int index, int* n, int spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (((spinIndices >> i) & 0x1) == 0)
	Index = ((n[i] + this->LzShiftDown) << 1) + ((spinIndices >> i) & 0x1);
      else
	Index = ((n[i] + this->LzShiftUp) << 1) + ((spinIndices >> i) & 0x1);
      if ((this->ProdATemporaryState & (0x1l << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      this->ProdATemporaryState &= ~(0x1l << Index);
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;

  return Coefficient;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	Index = ((m[i] + this->LzShiftDown) << 1) + spinIndices[i];
      else
 	Index = ((m[i] + this->LzShiftUp) << 1) + spinIndices[i];
      if ((TmpState & (0x1l << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (((spinIndices >> i) & 0x1) == 0)
	Index = ((m[i] + this->LzShiftDown) << 1) + ((spinIndices >> i) & 0x1);
      else
	Index = ((m[i] + this->LzShiftUp) << 1) + ((spinIndices >> i) & 0x1);
      if ((TmpState & (0x1l << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// generate all states corresponding to the constraints
// 
// nbrFermionsUp = number of fermions with spin up
// nbrFermionsDown = number of fermions with spin down
// lzMaxUp = momentum maximum value for a fermion
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereTwoLandauLevels::GenerateStates(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz, long pos)
{
  if ((nbrFermionsUp < 0) || (nbrFermionsDown < 0) || (totalLz < 0))
    return pos;
  if ((nbrFermionsUp == 0) && (totalLz == 0) && (nbrFermionsDown == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  
  if (lzMax < 0)
    return pos;
  
  if ((nbrFermionsUp == 1) && (nbrFermionsDown == 0)) 
    {
      if (lzMax >= totalLz)
        if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
         {
          this->StateDescription[pos] = 0x2ul << (totalLz << 1);
          return (pos + 1l);
         }
      return pos;
    }

  if ((nbrFermionsUp == 0) && (nbrFermionsDown == 1)) 
    {
      if (lzMax >= totalLz)
        if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
         {
          this->StateDescription[pos] = 0x1ul << (totalLz << 1);
          return (pos + 1l);
         }
      return pos;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return pos;
  
  if (((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp)) &&
      ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown)))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp - 1, nbrFermionsDown - 1, lzMax - 1, totalLz - (2 * lzMax), pos);
      unsigned long Mask = 0x3ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
         this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp - 1, nbrFermionsDown, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x2ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
        this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp, nbrFermionsDown - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x1ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
        this->StateDescription[pos] |= Mask;
    }
  
  return this->GenerateStates(nbrFermionsUp, nbrFermionsDown, lzMax - 1, totalLz, pos); 

  /*
  cout << "warning : untested code" << endl;
  if ((nbrFermionsUp < 0) || (nbrFermionsDown < 0) || (totalLz < 0))
    return pos;
  if ((nbrFermionsUp == 0) && (totalLz == 0) && (nbrFermionsDown == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  
  if (lzMax < 0)
    return pos;
  
  if ((nbrFermionsUp == 1) && (nbrFermionsDown == 0)) 
    {
      if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	{
	  this->StateDescription[pos] = 0x2ul << (totalLz << 1);
	  return (pos + 1l);
	}
      else
	return pos;
    }
  if ((nbrFermionsDown == 1) && (nbrFermionsUp == 0)) 
    {
      if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	{
	  this->StateDescription[pos] = 0x1ul << (totalLz << 1);
	  return (pos + 1l);
	}
      else
	return pos;
    }
  
  if ((lzMax == 0) && (totalLz != 0))
    return pos;
  
  
  if (((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp)) &&
      ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown)))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp - 1, nbrFermionsDown - 1, lzMax - 1, totalLz - (2 * lzMax), pos);
      unsigned long Mask = 0x3ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp - 1, nbrFermionsDown, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x2ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp, nbrFermionsDown  - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x1ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return this->GenerateStates(nbrFermionsUp, nbrFermionsDown, lzMax - 1, totalLz, pos);
  */
};

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereTwoLandauLevels::GenerateFullStates(int nbrFermions, int lzMax, int totalLz, long pos)
{
  if (nbrFermions < 0)
    return pos;
  if ((nbrFermions == 0) && (totalLz == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  
  if (lzMax < 0)
    return pos;
  
  if (nbrFermions == 1) 
    {
      if (lzMax >= totalLz)
	{
	  if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	    {
	      this->StateDescription[pos] = 0x2ul << (totalLz << 1);
	      ++pos;
	    }
	  if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	    {
	      this->StateDescription[pos] = 0x1ul << (totalLz << 1);
	      ++pos;
	    }
	}
      return pos;
    }
  
  if ((lzMax == 0) && (totalLz != 0))
    return pos;
  
  if (((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp)) &&
      ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown)))
    {
      long TmpPos = this->GenerateFullStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), pos);
      unsigned long Mask = 0x3ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      long TmpPos = this->GenerateFullStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x2ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    {
      long TmpPos = this->GenerateFullStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x1ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  
  return this->GenerateFullStates(nbrFermions, lzMax - 1, totalLz, pos);
}

// evaluate Hilbert space dimension
//
// nbrFermionsUp = number of fermions with spin up
// nbrFermionsDown = number of fermions with spin down
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereTwoLandauLevels::ShiftedEvaluateHilbertSpaceDimension(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz)
{
 if ((nbrFermionsUp < 0) || (nbrFermionsDown < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermionsUp == 0) && (nbrFermionsDown == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if ((nbrFermionsUp == 1) && (nbrFermionsDown == 0))
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
        {
           if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
             ++Tmp;
         }
      return Tmp;
    }

  if ((nbrFermionsUp == 0) && (nbrFermionsDown == 1))
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
        {
           if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
             ++Tmp;
         }
      return Tmp;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;
  
  long Tmp = 0l;
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
          Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp - 1, nbrFermionsDown - 1, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp -1, nbrFermionsDown, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))    
    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp, nbrFermionsDown - 1, lzMax - 1, totalLz - lzMax);
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp, nbrFermionsDown, lzMax - 1, totalLz);
  return Tmp;

/*
  cout << "warning : untested code" << endl;
  if ((nbrFermionsUp < 0) || (nbrFermionsDown < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermionsUp == 0) && (nbrFermionsDown == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if ((nbrFermionsUp == 1) && (nbrFermionsDown == 0)) 
    {
      if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	return 1l;
      else
	return 0l;
    }

  if ((nbrFermionsUp == 0) && (nbrFermionsDown == 1)) 
    {
      if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	return 1l;
      else
	return 0l;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
	 Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp - 1, nbrFermionsDown - 1, lzMax - 1, totalLz - (2 * lzMax));

      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp - 1, nbrFermionsDown, lzMax - 1, totalLz - lzMax);
    }

  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp, nbrFermionsDown - 1, lzMax - 1, totalLz - lzMax);

  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp, nbrFermionsDown, lzMax- 1, totalLz);

  return Tmp;
  */
}

// evaluate Hilbert space dimension without constraint on the number of particles per level
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereTwoLandauLevels::ShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if (nbrFermions == 1) 
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
	{
	  if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	    ++Tmp;
	  if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	    ++Tmp;
	}
      return Tmp;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;
  
  long Tmp = 0l;
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
	Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions -1, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))    
    Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax);
  Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz);
  return Tmp;
}


// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector FermionOnSphereTwoLandauLevels::ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace)
{
  RealVector FinalState(this->HilbertSpaceDimension, true);
  for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
    {
      unsigned long TmpUpState = upStateSpace.StateDescription[j] << this->LzShiftUp;
      int TmpPos = upStateSpace.LzMax + this->LzShiftUp;
      while (TmpPos > 0)
	{
	  unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
	  TmpUpState |= Tmp << TmpPos;
	  TmpUpState ^= Tmp;
	  --TmpPos;
	}
      TmpUpState <<= 1;
      double TmpComponent = upState[j];
      int Max = 63;
      while ((TmpUpState & (0x1ul << Max)) == 0x0ul)
	--Max;
      int Min = 0;
      while ((TmpUpState & (0x1ul << Min)) == 0x0ul)
	++Min;
      unsigned long TmpUpStateMask = (0x1ul << Max) - 1;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if ((this->StateDescription[i] & TmpUpState) == TmpUpState)
	  {	    
	    unsigned long TmpUpState3 = this->StateDescription[i] & TmpUpStateMask;
	    unsigned long TmpUpState2 = TmpUpState3;
#ifdef  __64_BITS__
	    TmpUpState3 &= 0x5555555555555555ul;
	    TmpUpState2 &= 0xaaaaaaaaaaaaaaaaul;
#else
	    TmpUpState3 &= 0x55555555ul;
	    TmpUpState2 &= 0xaaaaaaaaul;
#endif	    
	    unsigned long Sign = 0x0;
	    int Pos = this->LzMax << 1;
	    while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
	      Pos -= 2;
	    while (Pos > 0)
	      {
		unsigned long TmpUpState4 = TmpUpState2 & ((0x1ul << Pos) - 1ul);
#ifdef  __64_BITS__
		TmpUpState4 ^= TmpUpState4 >> 32;
#endif	
		TmpUpState4 ^= TmpUpState4 >> 16;
		TmpUpState4 ^= TmpUpState4 >> 8;
		TmpUpState4 ^= TmpUpState4 >> 4;
		TmpUpState4 ^= TmpUpState4 >> 2;
		TmpUpState4 ^= TmpUpState4 >> 1;
		Sign ^= TmpUpState4;
		Pos -= 2;
		while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
		  Pos -= 2;
	      }
	    if ((Sign & 0x1ul) == 0x0ul)
	      FinalState[i] = TmpComponent;
	    else
	      FinalState[i] = -TmpComponent;
	  }
    }

  for (int j = 0; j < downStateSpace.HilbertSpaceDimension; ++j)
    {
      unsigned long TmpDownState = downStateSpace.StateDescription[j] << this->LzShiftDown;
      int TmpPos = downStateSpace.LzMax + this->LzShiftDown;
      while (TmpPos > 0)
	{
	  unsigned long Tmp = TmpDownState & (0x1ul << TmpPos);
	  TmpDownState |= Tmp << TmpPos;
	  TmpDownState ^= Tmp;
	  --TmpPos;
	}
      double TmpComponent = downState[j];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if ((this->StateDescription[i] & TmpDownState) == TmpDownState)
	  {
	    FinalState[i] *= TmpComponent;
	  }
    }

  return FinalState;
}

// compute the product of a bosonic state and a fermionic state belonging in two Landau levels
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, 
								    BosonOnSphereShort* bosonSpace, FermionOnSphereTwoLandauLevels* finalSpace, int firstComponent, int nbrComponent)
{
  map <unsigned long, double> SortingMap;
  map <unsigned long, double>::iterator It;
  
  unsigned long* Monomial = new unsigned long [this->NbrFermions];
  unsigned long* Slater = new unsigned long [this->NbrFermions];
  int MaxComponent = firstComponent + nbrComponent;
  
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j] != 0)
	{
	  this->ConvertToMonomial(j, Slater);		
	  for (int i = firstComponent; i < MaxComponent; i++)
	    {
	      if(bosonState[i] != 0)
		{
		  bosonSpace->GetMonomial(i, Monomial);
		  
		  for (int Index=0; Index < this->NbrFermions;Index++)
		    Monomial[Index]*=2;
		  
		  this->MonomialsTimesSlater(Slater, Monomial, SortingMap, finalSpace);
		 
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++ )
 		    {
		      int TmpLzMax = 2 * finalSpace->LzMaxUp + 1;
		      while ( ( (*It).first >> TmpLzMax ) == 0x0ul)
			--TmpLzMax;
		      
		      outputVector[finalSpace->FindStateIndex( (*It).first, TmpLzMax )] += bosonState[i] * fermionState[j] * (*It).second;
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}

// compute the product of a monomial and a Slater determinant belonging in two Landau levels
// 
// slater = array where the Slater determinant is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereTwoLandauLevels::MonomialsTimesSlater(unsigned long* slater,unsigned long* monomial, map <unsigned long, double>& sortingMap, 
							  FermionOnSphereTwoLandauLevels* finalSpace)
{
  unsigned long * State =new unsigned long[this->NbrFermions];
  unsigned long * TmpSlater = new unsigned long [this->NbrFermions];
  int NbrPermutation = 0;
  bool Bool = true;
  double Coef = 1.0;
  int k = 1;
  
  for (int Index = 0; Index < this->NbrFermions; Index++)
    State[Index] = slater[Index] + monomial[Index];
  
  finalSpace->GeneratesDifferentState(sortingMap, slater, State, this, 0, Coef);

  while (std::prev_permutation(monomial, monomial + this->NbrFermions))
    {		
      NbrPermutation = 0;
      Bool = true;
      for (int Index = 0; Index < this->NbrFermions;Index++)
	{
	  State[Index] = slater[Index] + monomial[Index];
	  TmpSlater[Index] = slater[Index];
	}
      
      SortArrayDownOrdering(State, TmpSlater, this->NbrFermions, NbrPermutation);
      k = 1;
      while( ( k < this->NbrFermions ) && ( Bool == true))
	{
	  if(State[k-1] == State[k])
	    {
	      Bool = false;
	    }
	  k++;
	}
      if(Bool == true)
	{			
	  if((NbrPermutation & 1) == 0)
	    Coef = 1.0;
	  else
	    Coef = -1.0;
	  
	  finalSpace->GeneratesDifferentState(sortingMap, TmpSlater, State, this, 0, Coef);
	}
    }
  delete [] State;
  delete [] TmpSlater;
}



// compute the product of a bosonic state and a fermionic state belonging in two Landau levels
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::BosonicStateTimeFermionicState(LongRationalVector& bosonState, LongRationalVector& fermionState, LongRationalVector& outputVector, 
								    BosonOnSphereShort* bosonSpace, FermionOnSphereTwoLandauLevels* finalSpace, int firstComponent, int nbrComponent)
{
  map <unsigned long, LongRational> SortingMap;
  map <unsigned long, LongRational>::iterator It;
  
  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  
  int MaxComponent = firstComponent + nbrComponent;
  
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j].IsZero() == false)
	{
	  this->ConvertToMonomial(j, Slater);		
	  for (int i = firstComponent; i < MaxComponent; i++)
	    {
	      if( bosonState[i].IsZero() == false )
		{
		  bosonSpace->GetMonomial(i, Monomial);
		  for (int Index=0; Index < this->NbrFermions;Index++)
		    Monomial[Index] *= 2;
		  
		  this->MonomialsTimesSlater(Slater, Monomial, SortingMap, finalSpace);
		  
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++ )
		    {
		      int TmpLzMax = 2 * finalSpace->LzMaxUp + 1;
		      
		      while (( (*It).first >> TmpLzMax) == 0x0ul)
			--TmpLzMax;
		      
		      outputVector[finalSpace->FindStateIndex( (*It).first, TmpLzMax)] += bosonState[i] * fermionState[j] * (*It).second;
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}

// compute the product of a monomial and a Slater determinant belonging in two Landau levels
// 
// slater = array where the Slater determinant is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereTwoLandauLevels::MonomialsTimesSlater(unsigned long* slater,unsigned long* monomial, map <unsigned long, LongRational> & sortingMap, 
							  FermionOnSphereTwoLandauLevels* finalSpace)
{
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long * TmpSlater = new unsigned long [this->NbrFermions];
  
  int NbrPermutation = 0;
  bool Bool = true;
  LongRational Coef = 1l;
  int k = 1;

  for (int Index = 0; Index < this->NbrFermions; Index++)
    State[Index] = slater[Index] + monomial[Index];
  
  finalSpace->GeneratesDifferentState( sortingMap, slater, State, this, 0, Coef);
  while (std::prev_permutation(monomial, monomial + this->NbrFermions))
    {		
      NbrPermutation = 0;
      Bool = true;
      for (int Index = 0; Index < this->NbrFermions;Index++)
	{
	  State[Index] = slater[Index] + monomial[Index];
	  TmpSlater[Index] = slater[Index];
	}
      
      SortArrayDownOrdering(State, TmpSlater, this->NbrFermions, NbrPermutation);
      k = 1;
      while( ( k < this->NbrFermions ) && ( Bool == true ))
	{
	  if(State[k-1] == State[k])
	    {
	      Bool = false;
	    }
	  k++;
	}
      if(Bool == true) 
	{	
	  if ((NbrPermutation & 1) == 0)
	    Coef = 1l;
	  else
	    Coef=-1l;
	  
	  finalSpace->GeneratesDifferentState(sortingMap, TmpSlater, State, this, 0, Coef);
	}
    }
  delete [] State;
  delete [] TmpSlater;
}



// generate the different states that appear in the product of a monomial and a Slater determinant in the two Landau levels
//
// sortingMap = map in which the generated states and their coefficient will be stored
// slater = array where the Slater determinant is stored in its monomial representation
// state = array where the obtained state is stored in its monomial representation
// slaterSpace = pointer to the Hilbert Space which the Slater determinant belongs to
// index = index of the particle being examinate
// coef = coefficient of the state being generate

void FermionOnSphereTwoLandauLevels::GeneratesDifferentState( map <unsigned long, double>& sortingMap, unsigned long* slater, unsigned long* state, 
							      FermionOnSphereTwoLandauLevels* slaterSpace, int index, double coef)
{
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;

  while( ( ( ( state[index] & 0x01ul ) == 0x0ul ) || ( state[index] == ( (this->LzMaxUp << 1) + 0x1ul ) ) ) && ( index < (this->NbrFermions - 1) ) )
    index++;
  
  if( ( index == this->NbrFermions - 1 ) || ( (state[index] >> 1) == 0 ))
    {

      InsertionResult = sortingMap.insert (pair <unsigned long, double> ( this->ConvertFromMonomial(state), coef));
      
      if (InsertionResult.second == false)
	{
	  InsertionResult.first->second += coef;
	}

      if(( ( state[index] & 0x01ul) == 1) && ( (state[index] >> 1) != 0 ) )
	{
	  coef *= -((double)((double)(slater[index]>>1)/(double)slaterSpace->LzMaxUp)-((double)(state[index]>>1)/(double)this->LzMaxUp));
	  state[index]--;

	  InsertionResult = sortingMap.insert (pair <unsigned long, double> ( this->ConvertFromMonomial(state), coef));
	  
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += coef;
	    }
	  state[index]++;
	}
      return;
    }
  
  this->GeneratesDifferentState( sortingMap, slater, state, slaterSpace, index + 1, coef);
  
  if( state[index] != ( state[index + 1] + 1 ) )
	{
	  coef *= -((double)((double)(slater[index]>>1)/(double)slaterSpace->LzMaxUp)-((double)(state[index]>>1)/(double)this->LzMaxUp));
	  state[index]--;
	  this->GeneratesDifferentState( sortingMap, slater, state, slaterSpace, index + 1, coef);
	  state[index]++;
	}
}

// generate the different states that appear in the product of a monomial and a Slater determinant in the two Landau levels
//
// sortingMap = map in which the generated states and their coefficient will be stored
// slater = array where the Slater determinant is stored in its monomial representation
// state = array where the obtained state is stored in its monomial representation
// slaterSpace = pointer to the Hilbert Space which the Slater determinant belongs to
// index = index of the particle being examinate
// nbrStates = number of different obtained states
// coef = coefficient of the state being generate

void FermionOnSphereTwoLandauLevels::GeneratesDifferentState( map <unsigned long, LongRational>& sortingMap, unsigned long* slater, unsigned long * state, 
							      FermionOnSphereTwoLandauLevels* slaterSpace, int index, LongRational coef)
{
  pair <map <unsigned long, LongRational>::iterator, bool> InsertionResult;
  
  while(( ( (state[index]&0x01ul) == 0x0ul ) || ( state[index] == ( (this->LzMaxUp << 1) + 0x1ul) ) ) && ( index < (this->NbrFermions - 1 ) ) )
    index++;
  
  if( ( index == this->NbrFermions - 1 ) || (( state[index] >> 1 ) == 0 ) )
    {
      InsertionResult = sortingMap.insert (pair <unsigned long, LongRational> ( this->ConvertFromMonomial(state), coef));
      
      if (InsertionResult.second == false)
	{
	  InsertionResult.first->second += coef;
	}
      
      if( ( ( state[index] & 0x01ul ) == 1 ) && ( ( state[index] >> 1 ) != 0 ) )
	{
	  coef *= -( ((long) slater[index] >> 1 ) * ((long) this->LzMaxUp) - ((long)state[index]>>1) * ((long)slaterSpace->LzMaxUp));
	  coef /= ((long) (slaterSpace->LzMaxUp * this->LzMaxUp));
	  state[index]--;
	  
	  InsertionResult = sortingMap.insert (pair <unsigned long, LongRational> ( this->ConvertFromMonomial(state), coef));
	  
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += coef;
	    }

	  state[index]++;
	}
      return;
    }
  
  this->GeneratesDifferentState( sortingMap, slater, state, slaterSpace, index + 1, coef);
  
  if( state[index] != ( state[index + 1] + 1 ) )
    {
      coef *= -( ((long) slater[index] >> 1 ) * ((long) this->LzMaxUp) - ((long)state[index]>>1) * ((long)slaterSpace->LzMaxUp));
      coef /= ((long) (slaterSpace->LzMaxUp * this->LzMaxUp));
      state[index]--;
      this->GeneratesDifferentState(sortingMap, slater, state, slaterSpace, index + 1, coef);
      state[index]++;
    }
}

// compute the projection of the product of a bosonic state and a fermionic state
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed
// reverseFluxFlag = true if it a reverse flux attachment

void FermionOnSphereTwoLandauLevels::BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, 
								    BosonOnSphereShort* bosonSpace,FermionOnSphere* finalSpace, 
								    int firstComponent, int nbrComponent, bool reverseFluxFlag)
{
  map<unsigned long , double> SortingMap;
  map<unsigned long , double>::iterator It;
  
  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent + nbrComponent;
  int NbrVariable = 0;
  unsigned long* Variable = 0;
	
  BinomialCoefficients * TmpBinomials = 0;
  
  if (reverseFluxFlag == true)
    TmpBinomials = new BinomialCoefficients(finalSpace->LzMax + this->LzMaxDown);
  
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j] != 0)
	{
	  
	  if (reverseFluxFlag == false)
	    {
	      Variable = new unsigned long[this->NbrFermions];
	      this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	    }
	  else
	    this->ConvertToMonomial(j,Slater);
	  
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(bosonState[i] != 0)
		{
		  bosonSpace->GetMonomial(i, Monomial);
		  
		  if (reverseFluxFlag == false)
		    this->MonomialsTimesSlaterProjection(Slater, Monomial, Variable, NbrVariable, SortingMap, finalSpace);
		  else
		    {
		      this->MonomialsTimesSlaterProjectionReverse(Slater, Monomial, SortingMap, (*TmpBinomials), finalSpace );
		    }
		  for ( It = SortingMap.begin() ; It != SortingMap.end(); It++)
		    {
		      int TmpLzMax = finalSpace->LzMax;
		      while ((( (*It).first >> TmpLzMax) & 0x1ul) == 0x0ul)
			--TmpLzMax;
		      outputVector[finalSpace->FindStateIndex((*It).first, TmpLzMax)] += bosonState[i] * fermionState[j] *  (*It).second;
		    }
		}
	      SortingMap.clear();
	    }
	}
    }
}

// compute the projection of the product of a bosonic state and a fermionic state
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::BosonicStateTimeFermionicState( LongRationalVector& bosonState, LongRationalVector& fermionState, LongRationalVector& outputVector, 
								     BosonOnSphereShort* bosonSpace, FermionOnSphere* finalSpace, int firstComponent, int nbrComponent)
{
  map<unsigned long , LongRational> SortingMap;
  map<unsigned long , LongRational>::iterator It;
  
  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent + nbrComponent;
  int NbrVariable = 0;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j].IsZero() == false)
	{
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(bosonState[i].IsZero() == false)
		{
		  bosonSpace->GetMonomial(i, Monomial);
		  this->MonomialsTimesSlaterProjection(Slater, Monomial, Variable, NbrVariable, SortingMap, finalSpace);
		  
		  for ( It = SortingMap.begin() ; It != SortingMap.end(); It++ )
		    {
		      int TmpLzMax = finalSpace->LzMax;
		      while ((( (*It).first >> TmpLzMax) & 0x1ul) == 0x0ul)
			--TmpLzMax;
		      outputVector[finalSpace->FindStateIndex( (*It).first, TmpLzMax)] += bosonState[i] * fermionState[j] *  (*It).second;
		    }
		    SortingMap.clear();	
		}
	    }
	}
    }
}

// compute the projection of the product of a bosonic state and a fermionic state using the lz->-lz symmetry
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed
// reverseFluxFlag = true if it a reverse flux attachment

void FermionOnSphereTwoLandauLevels::BosonicStateTimeFermionicStateSymmetric(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, 
									     BosonOnSphereShort* bosonSpace,FermionOnSphere* finalSpace, int firstComponent,
									     int nbrComponent, bool reverseFluxFlag)
{ 
  map<unsigned long , double> SortingMap;
  map<unsigned long , double>::iterator It;

  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent + nbrComponent;
  int NbrVariable = 0;
  unsigned long* Variable = 0;
  BinomialCoefficients * TmpBinomials = 0;
  
  if (reverseFluxFlag == true)
    TmpBinomials = new BinomialCoefficients(finalSpace->LzMax + this->LzMaxDown);
  
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j] != 0)
	{
	  if (reverseFluxFlag == false)
	    {
	      Variable = new unsigned long[this->NbrFermions];
	      this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	    }
	  else
	    this->ConvertToMonomial(j,Slater);
	  
	  for (int i = firstComponent ; i < NbrMax ; i++)
	    {
	      if(bosonState[i] != 0)
		{
		  unsigned long TmpState = bosonSpace->FermionBasis->GetSymmetricState(bosonSpace->FermionBasis->StateDescription[i]);
		  int BTmpLzMax = bosonSpace->LzMax + this->NbrFermions - 1;
		  while (((TmpState >> BTmpLzMax) & 0x1ul) == 0x0ul)
		    --BTmpLzMax;
		  int SymmetricIndex = bosonSpace->FermionBasis->FindStateIndex(TmpState,BTmpLzMax);
		  if( SymmetricIndex > i)
		    {
		      bosonSpace->GetMonomial(i,Monomial);
		      
		      if (reverseFluxFlag == false)
			this->MonomialsTimesSlaterProjection(Slater, Monomial, Variable, NbrVariable, SortingMap, finalSpace);
		      else
			this->MonomialsTimesSlaterProjectionReverse(Slater, Monomial, SortingMap,*TmpBinomials, finalSpace);
		      
		      for ( It = SortingMap.begin() ; It != SortingMap.end(); It++ )
			{
			  int TmpLzMax = finalSpace->LzMax;
			  while ((( (*It).first >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  outputVector[finalSpace->FindStateIndex( (*It).first, TmpLzMax)] += bosonState[i] * fermionState[j] * (*It).second;
			  unsigned long TmpState = finalSpace->GetSymmetricState((*It).first);
			  TmpLzMax = finalSpace->LzMax;
			  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  outputVector[finalSpace->FindStateIndex(TmpState,TmpLzMax)] += bosonState[i] * fermionState[j] * (*It).second;
			}
		    }
		  else if( SymmetricIndex == i )
		    {
		      bosonSpace->GetMonomial(i, Monomial);
		      
		      if (reverseFluxFlag == false)
			this->MonomialsTimesSlaterProjection(Slater, Monomial, Variable, NbrVariable, SortingMap, finalSpace);
		      else
			this->MonomialsTimesSlaterProjectionReverse(Slater, Monomial, SortingMap,*TmpBinomials, finalSpace);
		      
		      
		      for ( It = SortingMap.begin(); It != SortingMap.end(); It++ )    		      
			{
			  int TmpLzMax = finalSpace->LzMax;
			  while ((( (*It).first >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  outputVector[finalSpace->FindStateIndex( (*It).first, TmpLzMax)] += bosonState[i] * fermionState[j] * (*It).second;
			}
		    }
		  SortingMap.clear();
		}
	    }
	}
    }
}

// compute the projection of the product of a bosonic state and a fermionic state using the lz->-lz symmetry
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::BosonicStateTimeFermionicStateSymmetric(LongRationalVector& bosonState, LongRationalVector& fermionState, LongRationalVector& outputVector, 
									     BosonOnSphereShort* bosonSpace,FermionOnSphere* finalSpace, int firstComponent,int nbrComponent)
{
  map<unsigned long , LongRational> SortingMap;
  map<unsigned long , LongRational>::iterator It;

  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent + nbrComponent;
  int NbrVariable = 0;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if( fermionState[j].IsZero() == false )
	{
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if( bosonState[i].IsZero() == false)
		{
		  unsigned long TmpState=bosonSpace->FermionBasis->GetSymmetricState(bosonSpace->FermionBasis->StateDescription[i]);
		  int BTmpLzMax = bosonSpace->LzMax + this->NbrFermions - 1;
		  
		  while (((TmpState >> BTmpLzMax) & 0x1ul) == 0x0ul)
		    --BTmpLzMax;
		  
		  if( bosonSpace->FermionBasis->FindStateIndex(TmpState, BTmpLzMax) > i )
		    {
		      bosonSpace->GetMonomial(i,Monomial);
		      this->MonomialsTimesSlaterProjection(Slater, Monomial, Variable, NbrVariable, SortingMap, finalSpace);
		      
		      for ( It = SortingMap.begin() ; It != SortingMap.end(); It++ )
			{
			  int TmpLzMax = finalSpace->LzMax;
			  
			  while ((( (*It).first >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  
			  outputVector[finalSpace->FindStateIndex( (*It).first, TmpLzMax)] += bosonState[i]*fermionState[j]* (*It).second;
			  unsigned long TmpState = finalSpace->GetSymmetricState((*It).first);
			  TmpLzMax = finalSpace->LzMax;
			  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  outputVector[finalSpace->FindStateIndex(TmpState,TmpLzMax)] += bosonState[i] * fermionState[j] * (*It).second;
			}
		    }
		  if( bosonSpace->FermionBasis->FindStateIndex(TmpState,BTmpLzMax) == i)
		    {
		      bosonSpace->GetMonomial(i, Monomial);
		      this->MonomialsTimesSlaterProjection(Slater, Monomial, Variable, NbrVariable, SortingMap, finalSpace);
		      for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
			{
			  int TmpLzMax = finalSpace->LzMax;
			  while (( ( (*It).first  >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  outputVector[finalSpace->FindStateIndex((*It).first, TmpLzMax)] += bosonState[i] * fermionState[j] * (*It).second;
			}
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}

// compute the product and the projection of a Slater determinant and a monomial 
// 
// slater = array where the slater is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalStates = array where the obtained states are stored in their fermionic representation

void FermionOnSphereTwoLandauLevels::MonomialsTimesSlaterProjection(unsigned long* slater, unsigned long* monomial, unsigned long * variable, 
								    int nbrVariable, map <unsigned long, double> & sortingMap, FermionOnSphere* finalSpace)
{
  unsigned long* State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0;

  pair <map <unsigned long, double>::iterator, bool> InsertionResult;
  
  bool Bool = true;
  double Coef = 1.0;

  long TmpLzMaxUp = this->LzMaxUp;
  long TmpFinalLzMaxUp = 2l + finalSpace->LzMax;
  double InverseFactor = 1.0 / (((double) TmpLzMaxUp) * ((double) TmpFinalLzMaxUp));

  for (int i = 0; i < this->NbrFermions ; i++)
    State[i] = slater[i] + monomial[i];

  for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
    {
      long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
      if (Numerator == 0l)
	Coef = 0.0;
      else
	Coef *= ((double) Numerator) * InverseFactor;
    }
  
  unsigned long Mask;
  unsigned long Sign = 0ul;
  if(Coef != 0.0)
    {
      for(int i = 0; ( i < this->NbrFermions ) && ( Bool == true ); i++)
	{
	  Mask= (1ul << ( State[i] - 1 ));
	  if ( (TmpState & Mask) != 0ul)
	    Bool = false;
	  unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef _64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= TmpState2 >> 16;
	  TmpState2 ^= TmpState2 >> 8;
	  TmpState2 ^= TmpState2 >> 4;
	  TmpState2 ^= TmpState2 >> 2;
	  TmpState2 ^= TmpState2 >> 1;
	  Sign ^= TmpState2;
	  TmpState |= Mask;
	}
      if (Bool == true)
	{
	  if ((Sign & 0x1ul) != 0ul)
	    Coef *= -1.0;
	  
	  InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState, Coef));
	  
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  while (std::prev_permutation(monomial,monomial+this->NbrFermions))
    {
      Coef = 1.0;
      for(int i = 0; i < this->NbrFermions; i++)
	{
	  State[i] = slater[i] + monomial[i];
	}
      for(int k = 0; (k < nbrVariable) && (Coef != 0.0); k++)
	{
	  long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
	  if (Numerator == 0l)
	    Coef = 0.0;
	  else
	    Coef *= ((double) Numerator) * InverseFactor;
	}
      if (Coef != 0.0)
	{
	  Bool = true;
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i = 0; (i < this->NbrFermions) && (Bool == true); i++)
	    {
	      Mask = (1ul << (State[i] - 1));
	      if((TmpState&Mask) != 0)
					Bool=false;
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif	
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  if (Bool == true)
	    {
	      if ((Sign & 0x1ul) != 0ul)
		Coef *= -1.0;	
	  
	      InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState, Coef));
	      
	      if (InsertionResult.second == false)
		{
		  InsertionResult.first->second += Coef;
		}
	    }
	}
    }
  delete [] State;
}

// compute the projection of the product of a Slater determinant and a monomial 
// 
// slater = array where the slater is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereTwoLandauLevels::MonomialsTimesSlaterProjection(unsigned long* slater,unsigned long* monomial,unsigned long* variable,
								    int nbrVariable, map <unsigned long, LongRational> & sortingMap, FermionOnSphere* finalSpace)
{
  unsigned long* State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0x0ul;
  
  bool Bool = true;
  pair <map <unsigned long, LongRational>::iterator, bool> InsertionResult;
  
  LongRational Coef (1l);
  long TmpLzMaxUp = this->LzMaxUp;
  long TmpFinalLzMaxUp = 2l + finalSpace->LzMax;
  long InverseFactor = (TmpLzMaxUp * TmpFinalLzMaxUp);
  
  for (int i = 0; i < this->NbrFermions ; i++)
    State[i] = slater[i] + monomial[i];
  
  for(int k = 0 ; (k < nbrVariable) && (Coef.IsZero() == false); k++)
    {
      long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
      if (Numerator == 0l)
	Coef = 0l;
      else
	{
	  Coef *=  Numerator;
	  Coef /= InverseFactor;
	}
    }
  
  unsigned long Mask;
  unsigned long Sign = 0ul;
  if( Coef.IsZero() == false )
    {
      for(int i = 0; (i < this->NbrFermions) && (Bool == true); i++)
	{
	  Mask = (1ul << (State[i]-1));
	  if ( (TmpState & Mask) != 0ul)
	    Bool = false;
	  unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef _64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= TmpState2 >> 16;
	  TmpState2 ^= TmpState2 >> 8;
	  TmpState2 ^= TmpState2 >> 4;
	  TmpState2 ^= TmpState2 >> 2;
	  TmpState2 ^= TmpState2 >> 1;
	  Sign ^= TmpState2;
	  TmpState |= Mask;
	}
      if(Bool == true)
	{
	  if ((Sign & 0x1ul) != 0ul)
	    Coef *= -1l;
	  	 
	  InsertionResult = sortingMap.insert (pair <unsigned long, LongRational> (TmpState, Coef));
	  
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  while (std::prev_permutation(monomial,monomial+this->NbrFermions))
    {
      Coef = 1l;
      for(int i = 0; i < this->NbrFermions; i++)
	{
	  State[i] = slater[i] + monomial[i];
	}
      for(int k = 0; (k < nbrVariable) && (Coef.IsZero() == false); k++)
	{
	  long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
	  if (Numerator == 0l)
	    Coef = 0l;
	  else
	    {
	      Coef *=  Numerator;
	      Coef /= InverseFactor;
	    }
	}
      if( Coef.IsZero() == false )
	{
	  Bool = true;
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i=0; (i < this->NbrFermions)&&(Bool);i++)
	    {
	      Mask = (1ul << (State[i] - 1));
	      if( ( TmpState & Mask ) != 0)
		Bool = false;
	      
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif	
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  if(Bool == true)
	    {
	      if ((Sign & 0x1ul) != 0ul)
		Coef *= -1l;
	      
	      InsertionResult = sortingMap.insert (pair <unsigned long, LongRational> (TmpState, Coef));
	  
	      if (InsertionResult.second == false)
		{
		  InsertionResult.first->second += Coef;
		}
	      
	    }
	}
    }
  delete [] State;
}


// compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed
// reverseFluxFlag = true if it a reverse flux attachment

void FermionOnSphereTwoLandauLevels::LLLFermionicStateTimeFermionicState(RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, 
									 FermionOnSphere* lllFermionSpace,BosonOnSphereShort* finalSpace, 
									 int firstComponent,int nbrComponent, bool reverseFluxFlag)
{
  map<unsigned long , double> SortingMap;
  map<unsigned long , double>::iterator It;
  
  unsigned long* LLLSlater = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];

  int NbrMax = firstComponent + nbrComponent;
  int NbrVariable = 0;

  FactorialCoefficient Coefficient;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j] != 0)
	{
	  NbrVariable = 0;
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(lllFermionState[i] != 0)
		{
		  lllFermionSpace->GetMonomial(i, LLLSlater);
		  if (reverseFluxFlag == false)
		    this->SlaterTimesSlaterProjection(Slater,LLLSlater,Variable,NbrVariable, SortingMap, finalSpace);
		  else
		    this->SlaterTimesSlaterProjectionReverse(Slater,LLLSlater,Variable,NbrVariable, SortingMap, finalSpace);
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
		    {
		      int FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
		      while ((( ((*It).first) >> FTmpLzMax) & 0x1ul) == 0x0ul)
			--FTmpLzMax;
		      finalSpace->FermionToBoson((*It).first, FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
		      Coefficient.SetToOne();
		      for(int p = 0; p < finalSpace->TemporaryStateLzMax + 1; p++)
			{
			  Coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			}
		      outputVector[finalSpace->FermionBasis->FindStateIndex((*It).first, FTmpLzMax)] += lllFermionState[i] * fermionState[j] * (*It).second * Coefficient.GetIntegerValue();
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}

// compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels
//
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed
// reverseFluxFlag = true if it a reverse flux attachment

void FermionOnSphereTwoLandauLevels::LLLFermionicStateTimeFermionicState(LongRationalVector& lllFermionState, LongRationalVector& fermionState, LongRationalVector& outputVector, 
									 FermionOnSphere* lllFermionSpace,BosonOnSphereShort* finalSpace, 
									 int firstComponent,int nbrComponent, bool reverseFluxFlag)
{
  map<unsigned long , LongRational> SortingMap;
  map<unsigned long , LongRational>::iterator It;

  unsigned long* LLLSlater = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent + nbrComponent;
  int NbrVariable = 0;
  FactorialCoefficient Coefficient;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j].IsZero()== false)
	{
	  NbrVariable = 0;
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(lllFermionState[i].IsZero() == false)
		{
		  lllFermionSpace->GetMonomial(i, LLLSlater);
		  if (reverseFluxFlag == false)
		    this->SlaterTimesSlaterProjection(Slater,LLLSlater,Variable,NbrVariable, SortingMap, finalSpace);
		  else
		    this->SlaterTimesSlaterProjectionReverse(Slater,LLLSlater,Variable,NbrVariable, SortingMap, finalSpace);
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
		    {
		      int FTmpLzMax = finalSpace->LzMax+this->NbrFermions-1;
		      
		      while ( ( ( (*It).first >> FTmpLzMax ) & 0x1ul) == 0x0ul)
						--FTmpLzMax;
		      finalSpace->FermionToBoson((*It).first, FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
		      Coefficient.SetToOne();
		      for(int p = 0; p < finalSpace->TemporaryStateLzMax + 1; p++)
			{
			  Coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			}
		      outputVector[finalSpace->FermionBasis->FindStateIndex((*It).first, FTmpLzMax)] += lllFermionState[i] * fermionState[j]* (*It).second * Coefficient.GetIntegerValue();
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}


// compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels using lz->-lz symmetry
//
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed
// reverseFluxFlag = true if it a reverse flux attachment

void FermionOnSphereTwoLandauLevels::LLLFermionicStateTimeFermionicStateSymmetric(RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, 
										  FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, 
										  int firstComponent, int nbrComponent, bool reverseFluxFlag)
{
  unsigned long* LLLSlater = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];

  map<unsigned long , double> SortingMap;
  map<unsigned long , double>::iterator It;

  int NbrMax = firstComponent + nbrComponent;
  int NbrVariable = 0;

  unsigned long* Variable = new unsigned long[this->NbrFermions];
  FactorialCoefficient coefficient;

  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j] != 0)
	{
	  NbrVariable = 0;
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent ; i < NbrMax ; i++)
	    {
	      if(lllFermionState[i] != 0)
		{
		  unsigned long TmpState = lllFermionSpace->GetSymmetricState(lllFermionSpace->StateDescription[i]);
		  int TmpLzMax = lllFermionSpace->LzMax;
		  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		    --TmpLzMax;
		  if(lllFermionSpace->FindStateIndex(TmpState,TmpLzMax) > i)
		    {
		      lllFermionSpace->GetMonomial(i,LLLSlater);
		      if (reverseFluxFlag == false)
			this->SlaterTimesSlaterProjection(Slater, LLLSlater, Variable, NbrVariable, SortingMap, finalSpace);
		      else
			this->SlaterTimesSlaterProjectionReverse(Slater, LLLSlater, Variable, NbrVariable, SortingMap, finalSpace);
		      for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
			{
			  int FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
			  while ((( (*It).first >> FTmpLzMax ) & 0x1ul) == 0x0ul)
			    --FTmpLzMax;
			  finalSpace->FermionToBoson((*It).first, FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
			  coefficient.SetToOne();
			  for(int p = 0; p <= finalSpace->TemporaryStateLzMax; p++)
			    {
			      coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			    }
			  
			  double TmpCoefficient = lllFermionState[i] * fermionState[j] * (*It).second * coefficient.GetIntegerValue();
			  outputVector[finalSpace->FermionBasis->FindStateIndex((*It).first, FTmpLzMax)] += TmpCoefficient;
			  
			  unsigned long TmpState = finalSpace->FermionBasis->GetSymmetricState((*It).first);
			  FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
			  while (((TmpState >> FTmpLzMax) & 0x1ul) == 0x0ul)
			    --FTmpLzMax;
			  
			  outputVector[finalSpace->FermionBasis->FindStateIndex(TmpState,FTmpLzMax)] += TmpCoefficient;
			}
		    }
		  if(lllFermionSpace->FindStateIndex(TmpState,TmpLzMax) == i)
		    {
		      lllFermionSpace->GetMonomial(i, LLLSlater);
		      this->SlaterTimesSlaterProjection(Slater, LLLSlater, Variable, NbrVariable, SortingMap, finalSpace);
		      for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
			{
			  int FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
			  while( ( ( (*It).first >> FTmpLzMax ) & 0x1ul ) == 0x0ul )
			    --FTmpLzMax;
			  finalSpace->FermionToBoson((*It).first, FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
			  coefficient.SetToOne();
			  for(int p = 0; p <= finalSpace->TemporaryStateLzMax; p++)
			    {
			      coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			    }			  
			  outputVector[finalSpace->FermionBasis->FindStateIndex((*It).first, FTmpLzMax)] += lllFermionState[i] * fermionState[j] * (*It).second * coefficient.GetIntegerValue();
			}
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}

// compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels using lz->-lz symmetry
//
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed
// reverseFluxFlag = true if it a reverse flux attachment

void FermionOnSphereTwoLandauLevels::LLLFermionicStateTimeFermionicStateSymmetric(LongRationalVector& lllFermionState, LongRationalVector& fermionState, LongRationalVector& outputVector, 
										  FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, 
										  int firstComponent, int nbrComponent, bool reverseFluxFlag)
{
  unsigned long* LLLSlater = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];

  map<unsigned long , LongRational> SortingMap;
  map<unsigned long , LongRational>::iterator It;
  
  int NbrMax = firstComponent+nbrComponent;
  int NbrVariable = 0;

  unsigned long* Variable = new unsigned long[this->NbrFermions];
  FactorialCoefficient coefficient;
  
  for (int j = 0 ; j < this->HilbertSpaceDimension ; j++)
    {
      if(fermionState[j].IsZero() == false)
	{
	  NbrVariable = 0;
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent ; i < NbrMax ; i++)
	    {
	      if( lllFermionState[i].IsZero() == false )
		{
		  unsigned long TmpState = lllFermionSpace->GetSymmetricState(lllFermionSpace->StateDescription[i]);
		  int TmpLzMax = lllFermionSpace->LzMax;
		  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		    --TmpLzMax;
		  if(lllFermionSpace->FindStateIndex(TmpState,TmpLzMax) > i)
		    {
		      lllFermionSpace->GetMonomial(i,LLLSlater);
		      if (reverseFluxFlag == false)
			this->SlaterTimesSlaterProjection(Slater, LLLSlater, Variable, NbrVariable,SortingMap, finalSpace);
		      else
			this->SlaterTimesSlaterProjectionReverse(Slater, LLLSlater, Variable, NbrVariable,SortingMap, finalSpace);
		      for ( It = SortingMap.begin(); It != SortingMap.end(); It++)		      
			{
			  int FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
			  while ((( (*It).first >> FTmpLzMax) & 0x1ul) == 0x0ul)
			    --FTmpLzMax;
			  finalSpace->FermionToBoson((*It).first, FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
			  coefficient.SetToOne();
			  for(int p = 0; p <= finalSpace->TemporaryStateLzMax + 1; p++)
			    {
			      coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			    }
			  LongRational TmpCoefficient = lllFermionState[i] * fermionState[j] * (*It).second * coefficient.GetIntegerValue();
			  outputVector[finalSpace->FermionBasis->FindStateIndex((*It).first, FTmpLzMax)] += TmpCoefficient;
			  
			  unsigned long TmpState = finalSpace->FermionBasis->GetSymmetricState((*It).first);
			  FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
			  while (((TmpState >> FTmpLzMax) & 0x1ul) == 0x0ul)
			    --FTmpLzMax;
			  
			  outputVector[finalSpace->FermionBasis->FindStateIndex(TmpState,FTmpLzMax)] += TmpCoefficient;
			}
		    }
		  if( lllFermionSpace->FindStateIndex(TmpState, TmpLzMax) == i )
		    {
		      lllFermionSpace->GetMonomial(i, LLLSlater);
		      this->SlaterTimesSlaterProjection(Slater, LLLSlater, Variable, NbrVariable, SortingMap, finalSpace);
	
		      for ( It = SortingMap.begin(); It != SortingMap.end(); It++)		      
			{
			  int FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
			  while (( (  (*It).first >> FTmpLzMax) & 0x1ul ) == 0x0ul )
			    --FTmpLzMax;
			  finalSpace->FermionToBoson( (*It).first, FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
			  coefficient.SetToOne();
			  for(int p = 0; p <= finalSpace->TemporaryStateLzMax; p++)
			    {
			      coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			    }			  
			  outputVector[finalSpace->FermionBasis->FindStateIndex((*It).first,FTmpLzMax)] += lllFermionState[i] * fermionState[j] * (*It).second * coefficient.GetIntegerValue();
			}
		    }
		}
		SortingMap.clear();
	    }
	}
    }
}

// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace


void FermionOnSphereTwoLandauLevels::SlaterTimesSlaterProjection(unsigned long* slater, unsigned long* lllslater, unsigned long * variable, int nbrVariable, 
								 map <unsigned long, double> & sortingMap, BosonOnSphereShort* finalSpace)
{
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0ul;
  double Coef = 1.0;

  long TmpLzMaxUp = this->LzMaxUp;
  long TmpFinalLzMaxUp = 2l + finalSpace->LzMax;
  double InverseFactor = 1.0 / (((double) TmpLzMaxUp) * ((double) TmpFinalLzMaxUp));

  pair <map <unsigned long, double>::iterator, bool> InsertionResult;

  for (int i = 0; i < this->NbrFermions ; i++)
    State[i] = slater[i] + lllslater[i];
  for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
    {
      long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
      if (Numerator == 0l)
	Coef = 0.0;
      else
	Coef *= ((double) Numerator) * InverseFactor;
    }
  
  unsigned long Mask = 0ul;
  unsigned long Sign = 0ul;
  if( Coef != 0.0 )
    {
      for (int i = 0; i < this->NbrFermions; i++)
	{
	  State[i]--;
	}
      InsertionResult = sortingMap.insert (pair<unsigned long,double> (finalSpace->ConvertFromMonomial(State), Coef));
      if (InsertionResult.second == false)
	{
	  InsertionResult.first->second += Coef;
	}
    }
  while (std::prev_permutation(lllslater, lllslater + this->NbrFermions))
    {
      Coef = 1.0;
      for(int i = 0 ; i < this->NbrFermions; i++)
	State[i] = slater[i] + lllslater[i];
      for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
	{
	  long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
	  if (Numerator == 0l)
	    Coef = 0.0;
	  else
	    Coef *= ((double) Numerator) * InverseFactor;
	}
      if (Coef != 0.0)
	{
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i=0; (i < this->NbrFermions);i++)
	    {
	      State[i]--;
	      Mask = (1ul << lllslater[i]);
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  SortArrayDownOrdering(State,this->NbrFermions);
	 
	  if ((Sign & 0x1ul) != 0ul)
	    Coef *= -1.0;
	 
	  InsertionResult = sortingMap.insert ( pair <unsigned long, double> (finalSpace->ConvertFromMonomial(State), Coef));
	  
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  delete [] State;
}

// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels with reverse flux attachment
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereTwoLandauLevels::SlaterTimesSlaterProjectionReverse(unsigned long* slater, unsigned long* lllslater, unsigned long * variable, int nbrVariable, 
									map<unsigned long, double> & sortingMap, BosonOnSphereShort* finalSpace)
{
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0ul;
  double Coef = 1.0;

  long TmpLzMaxUp = this->LzMaxUp;
  long TmpFinalLzMaxUp = 2l + finalSpace->LzMax;
  double InverseFactor = 1.0 / (((double) TmpLzMaxUp) * ((double) TmpFinalLzMaxUp));

  pair <map <unsigned long, double>::iterator, bool> InsertionResult;

  for (int i = 0; i < this->NbrFermions ; i++)
    State[i] = slater[i] + lllslater[i];
  for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
    {
      long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
      if (Numerator == 0l)
	Coef = 0.0;
      else
	Coef *= ((double) Numerator) * InverseFactor;
    }
  
  unsigned long Mask = 0ul;
  unsigned long Sign = 0ul;
  if( Coef != 0.0 )
    {
      for (int i = 0; i < this->NbrFermions; i++)
	{
	  State[i]--;
	}
      InsertionResult = sortingMap.insert (pair<unsigned long,double> (finalSpace->ConvertFromMonomial(State), Coef));
      if (InsertionResult.second == false)
	{
	  InsertionResult.first->second += Coef;
	}
    }
  while (std::prev_permutation(lllslater, lllslater + this->NbrFermions))
    {
      Coef = 1.0;
      for(int i = 0 ; i < this->NbrFermions; i++)
	State[i] = slater[i] + lllslater[i];
      for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
	{
	  long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
	  if (Numerator == 0l)
	    Coef = 0.0;
	  else
	    Coef *= ((double) Numerator) * InverseFactor;
	}
      if( Coef != 0.0 )
	{
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i=0; (i < this->NbrFermions);i++)
	    {
	      State[i]--;
	      Mask = (1ul << lllslater[i]);
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  SortArrayDownOrdering(State,this->NbrFermions);
	 
	  if ((Sign & 0x1ul) != 0ul)
	    Coef *= -1.0;
	 
	  InsertionResult = sortingMap.insert ( pair <unsigned long, double> (finalSpace->ConvertFromMonomial(State), Coef));
	  
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  delete [] State;
}

// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace
// return value = number of different obtained states

void FermionOnSphereTwoLandauLevels::SlaterTimesSlaterProjection(unsigned long* slater,unsigned long* lllslater,unsigned long * variable,int nbrVariable, 
								 map <unsigned long, LongRational> & sortingMap, BosonOnSphereShort* finalSpace)
{
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0ul;
  LongRational Coef(1l);

  long TmpLzMaxUp = this->LzMaxUp;
  long TmpFinalLzMaxUp = 2l + finalSpace->LzMax;
  long InverseFactor = TmpFinalLzMaxUp * TmpLzMaxUp;

  pair <map <unsigned long, LongRational>::iterator, bool> InsertionResult;

  for (int i = 0; i < this->NbrFermions; i++)
    State[i] = slater[i] + lllslater[i];

  for(int k = 0; (k < nbrVariable) && (Coef.IsZero() == false); k++)
    {
      long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
      if (Numerator == 0l)
	Coef = 0l;
      else
	{
	  Coef *= Numerator;
	  Coef /= InverseFactor;
	}
    }
	
  unsigned long Mask = 0ul;
  unsigned long Sign = 0ul;
  if(Coef.IsZero() == false)
    {
      for (int i=0; (i < this->NbrFermions);i++)
	{
	  State[i]--;
	}
      
      InsertionResult = sortingMap.insert ( pair <unsigned long, LongRational> (finalSpace->ConvertFromMonomial(State), Coef));
      
      if (InsertionResult.second == false)
	{
	  InsertionResult.first->second += Coef;
	}
    }
  while (std::prev_permutation(lllslater, lllslater + this->NbrFermions))
    {
      Coef = 1l;

      for(int i = 0 ; i < this->NbrFermions; i++)
	State[i] = slater[i] + lllslater[i];

      for(int k = 0 ; (k < nbrVariable) && (Coef.IsZero() == false); k++)
	{
	  long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
	  if (Numerator == 0l)
	    Coef = 0l;
	  else
	    {
	      Coef *= Numerator;
	      Coef /= InverseFactor;
	    }
	}
      if( Coef.IsZero() == false )
	{
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i=0; (i < this->NbrFermions);i++)
	    {
	      State[i]--;
	      Mask = (1ul << lllslater[i]);
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  SortArrayDownOrdering(State,this->NbrFermions);
	  
	  if ((Sign & 0x1ul) != 0ul)
	    Coef *= -1l;

	  InsertionResult = sortingMap.insert ( pair <unsigned long, LongRational> (finalSpace->ConvertFromMonomial(State), Coef));
	  
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  delete [] State;
}

// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels with reverse flux attachment
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace
// return value = number of different obtained states

void FermionOnSphereTwoLandauLevels::SlaterTimesSlaterProjectionReverse(unsigned long* slater,unsigned long* lllslater,unsigned long * variable,int nbrVariable, 
									map <unsigned long, LongRational> & sortingMap, BosonOnSphereShort* finalSpace)
{
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0ul;
  LongRational Coef(1l);

  long TmpLzMaxUp = this->LzMaxUp;
  long TmpFinalLzMaxUp = 2l + finalSpace->LzMax;
  long InverseFactor = TmpFinalLzMaxUp * TmpLzMaxUp;

  pair <map <unsigned long, LongRational>::iterator, bool> InsertionResult;

  for (int i = 0; i < this->NbrFermions; i++)
    State[i] = slater[i] + lllslater[i];

  for(int k = 0; (k < nbrVariable) && (Coef.IsZero() == false); k++)
    {
      long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
      if (Numerator == 0l)
	Coef = 0l;
      else
	{
	  Coef *= Numerator;
	  Coef /= InverseFactor;
	}
    }
	
  unsigned long Mask = 0ul;
  unsigned long Sign = 0ul;
  if(Coef.IsZero() == false)
    {
      for (int i=0; (i < this->NbrFermions);i++)
	{
	  State[i]--;
	}
      
      InsertionResult = sortingMap.insert ( pair <unsigned long, LongRational> (finalSpace->ConvertFromMonomial(State), Coef));
      
      if (InsertionResult.second == false)
	{
	  InsertionResult.first->second += Coef;
	}
    }
  while (std::prev_permutation(lllslater, lllslater + this->NbrFermions))
    {
      Coef = 1l;

      for(int i = 0 ; i < this->NbrFermions; i++)
	State[i] = slater[i] + lllslater[i];

      for(int k = 0 ; (k < nbrVariable) && (Coef.IsZero() == false); k++)
	{
	  long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
	  if (Numerator == 0l)
	    Coef = 0l;
	  else
	    {
	      Coef *= Numerator;
	      Coef /= InverseFactor;
	    }
	}
      if( Coef.IsZero() == false )
	{
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i=0; (i < this->NbrFermions);i++)
	    {
	      State[i]--;
	      Mask = (1ul << lllslater[i]);
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  SortArrayDownOrdering(State,this->NbrFermions);
	  
	  if ((Sign & 0x1ul) != 0ul)
	    Coef *= -1l;

	  InsertionResult = sortingMap.insert ( pair <unsigned long, LongRational> (finalSpace->ConvertFromMonomial(State), Coef));
	  
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  delete [] State;
}

// compute the number of particles in each Landau level
//
// state = ID of the state to handle
// LLOccupationConfiguration = array where the decomposition will be store

void  FermionOnSphereTwoLandauLevels::LandauLevelOccupationNumber(int state, int* lLOccupationConfiguration)
{
  unsigned long State = this->StateDescription[state];
  for(int i = this->LzMax; i>=0; i--)
    {
      switch ((State >> (i << 1)) & 0x3ul)
	{
	case 0x1ul:
	  {
	    lLOccupationConfiguration[0]++;
	  }
	  break;
	case 0x2ul:
	  {
	    lLOccupationConfiguration[1]++;
	  }
	  break;
	case 0x3ul:
	  {
	    lLOccupationConfiguration[0]++;
	    lLOccupationConfiguration[1]++;
	  }
	  break;
	default :
	  {
	    break;
	  }
	}
    }
}


// project out any configurations that have particles on levels other than lll
//
// inputVector = vector to apply the projection to
// outputVector = projected vector
// finalSpace = reference to space of output vector space

void FermionOnSphereTwoLandauLevels::ProjectionInTheLowestLevel(RealVector &inputVector, RealVector & outputVector, FermionOnSphere * finalSpace)
{
  unsigned long * Landau = new unsigned long [this->NbrFermions];
  unsigned long * State = new unsigned long [this->NbrFermions];
  for(int i = 0; i < this->NbrFermions ; i++)
      Landau[i] = 0;
  for(int i = 0 ; i < finalSpace->GetHilbertSpaceDimension() ; i++)
    {
      finalSpace->GetMonomial(i,State);
      unsigned long Etat = this->ConvertFromPowerLandauRepresentation(State,Landau);
      int TmpLzMax = 2*this->LzMaxUp+1;
      while ((Etat >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
	outputVector[i] = inputVector[this->CarefulFindStateIndex(Etat,TmpLzMax)];
    }
}

// compute the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::LLLFermionicStateTimeFermionicState(RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, 
									 FermionOnSphere* lllFermionSpace, BosonOnSphereTwoLandauLevels * finalSpace, 
									 int firstComponent,int nbrComponent)
{
  map<unsigned long , double> SortingMap;
  map<unsigned long , double>::iterator It;
  
  unsigned long* LLLSlater = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  
  int NbrMax = firstComponent + nbrComponent;
  
  FactorialCoefficient Coefficient;
  
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j] != 0)
	{
	  this->ConvertToMonomial(j, Slater);
			
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(lllFermionState[i] != 0)
		{
		  lllFermionSpace->GetMonomial(i, LLLSlater);
		  
		  for (int Index = 0; Index < this->NbrFermions;Index++)
		    {
		      LLLSlater[Index] *= 2;
		    }
		  
		  this->SlaterTimesSlater(Slater, LLLSlater , SortingMap, finalSpace); 
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
		    {
		      long Index = finalSpace->FindStateIndex(((*It).first));
		      finalSpace->FermionToBoson(finalSpace->StateDescription[Index],finalSpace->StateLzMax[Index], finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
		      Coefficient.SetToOne();
		      for(int p = 0; p < finalSpace->TemporaryStateLzMax + 1; p++)
			{
			  Coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			}
		      outputVector[Index] += lllFermionState[i] * fermionState[j] * (*It).second * Coefficient.GetIntegerValue();
		    }
		}
		SortingMap.clear();
	    }
	}
    }
}

// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereTwoLandauLevels::SlaterTimesSlater(unsigned long* slater, unsigned long* lllslater, map <unsigned long, double> & sortingMap, 
						       BosonOnSphereTwoLandauLevels* finalSpace)
{
  unsigned long * State = new unsigned long[this->NbrFermions];
  double Coef = 1.0;
  unsigned long TmpState = 0ul;
  unsigned long Mask = 0ul;
  unsigned long Sign = 0ul;
  
  for (int Index = 0; Index < this->NbrFermions; Index++)
    State[Index] = slater[Index] + lllslater[Index];
  
  finalSpace->GeneratesDifferentState(sortingMap , slater , State, this, 0, Coef);
  
  while (std::prev_permutation(lllslater, lllslater + this->NbrFermions))
    {	
      Coef = 1.0;
      for (int Index = 0; Index < this->NbrFermions; Index++)
	{
	  State[Index] = slater[Index] + lllslater[Index];
	}
      
      TmpState = 0ul;
      Sign = 0ul;
      for (int i = 0; i < this->NbrFermions; i++)
	{
	  Mask = (1ul << (lllslater[i] >> 1) );
	  unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= TmpState2 >> 16;
	  TmpState2 ^= TmpState2 >> 8;
	  TmpState2 ^= TmpState2 >> 4;
	  TmpState2 ^= TmpState2 >> 2;
	  TmpState2 ^= TmpState2 >> 1;
	  Sign ^= TmpState2;
	  TmpState |= Mask;
	}
      
      if ((Sign & 0x1ul) != 0ul)
	Coef *= -1l;
      
      finalSpace->GeneratesDifferentState(sortingMap , slater , State, this, 0, Coef);
    }
  delete [] State;
}

// find state index. not using lookup table at the moment
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereTwoLandauLevels::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  int start, end, mid;
  
  start = 0;					//index of start of range
  end = this->HilbertSpaceDimension;		//index of end of range + 1 
  
  while ( (end - start) > 0 ) 
    {
      mid = (start + end) >> 1 ; 		//work out the mid-point
      if ( stateDescription > this->StateDescription[mid] )	
	{
	  end = mid;
	}
      else if ( stateDescription < this->StateDescription[mid] )
	{
	  start = mid + 1;	  	 
	}
      else
	{
	  return mid;
	}
    }	
  return this->HilbertSpaceDimension;
}

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int FermionOnSphereTwoLandauLevels::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != (this->LzMax + 1))
    return -1;
  unsigned long TmpState = 0x0ul;
  int TmpNbrParticles = 0;
  int TmpTotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      if (TmpDescription[i][0] == 'u')
	{
	  TmpState |= 0x2ul << (2 * i);
	  TmpTotalLz += i;
	  ++TmpNbrParticles;	  
	}
      else
	{
	  if (TmpDescription[i][0] == 'd')
	    {
	      TmpState |= 0x1ul << (2 * i);
	      TmpTotalLz += i;
	      ++TmpNbrParticles;	  
	    }
	  else
	    {
	      if (TmpDescription[i][0] == 'X')
		{
		  TmpState |= 0x3ul << (2 * i);
		  TmpTotalLz += 2 * i;
		  TmpNbrParticles += 2;	  
		}
	      else
		{
		  if (TmpDescription[i][0] != '0')
		    {
		      return -1;
		    }
		}
	    }
	}
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if ((TmpNbrParticles != this->NbrFermions) 
      || (TmpTotalLz != ((this->TotalLz + this->NbrFermions * this->LzMax) >> 1)))
    return -1;
  int TmpLzMax = 2 * this->LzMax + 1;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}

// compute the product and the projection of a Slater determinant and a monomial with reverse flux attachment
// 
// slater = array where the slater is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// landau =  array where the landau level of fermions is stored
// sortingMap = map in which the generated states and their coefficient will be stored
// binomialsCoefficient = binomials coefficient needed in the computation
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereTwoLandauLevels::MonomialsTimesSlaterProjectionReverse(unsigned long* slater, unsigned long* monomial, map <unsigned long, double> & sortingMap, 
									   BinomialCoefficients& binomialsCoefficient, FermionOnSphere* finalSpace)
{
  unsigned long* State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0ul;
  
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;
  
  bool Bool = true;
  double Coef;
  double TmpCoef;
  const long OtherSpaceNphi =  finalSpace->LzMax + this->LzMaxDown;
  
  do
    {
      Coef = 1.0;
      for(int k = 0 ; (k < this->NbrFermions) && (Coef != 0.0); k++)
	{
	  const int NbrLL = (slater[k] & 0x1ul);
	  
	  if((monomial[k] + ((long) (slater[k]>>1)) >= this->LzMaxDown + 1l)&&( monomial[k] + ((long) (slater[k]>>1)) <= (this->LzMaxDown + 1l +finalSpace->LzMax)))
	    {
	      State[k] = monomial[k] + ((long) (slater[k]>>1)) - this->LzMaxDown - 1l;
	      
	      if(NbrLL == 1)
		{
		  TmpCoef = ((long)(((long) monomial[k] + 1l) * finalSpace->LzMax) - ((long)State[k])*( NbrLL + 1l + OtherSpaceNphi));
		  TmpCoef *= binomialsCoefficient.GetNumericalCoefficient(this->LzMaxDown + NbrLL + 1l , (slater[k]>>1));
		  TmpCoef /= (this->LzMaxDown + NbrLL +1);
		  TmpCoef /= binomialsCoefficient.GetNumericalCoefficient(OtherSpaceNphi, monomial[k]);
		  TmpCoef /= ((OtherSpaceNphi + 1 + NbrLL)*(OtherSpaceNphi + NbrLL));
		}
	      else
		{
		  TmpCoef = binomialsCoefficient.GetNumericalCoefficient(this->LzMaxDown , (slater[k]>>1) - 1l);
		  TmpCoef /= binomialsCoefficient.GetNumericalCoefficient(OtherSpaceNphi , monomial[k]);
		  TmpCoef /= (OtherSpaceNphi + 1);
		}
	      TmpCoef *= binomialsCoefficient.GetNumericalCoefficient(finalSpace->LzMax , State[k]);
	      TmpCoef *= finalSpace->LzMax + 1;
	      Coef *= TmpCoef;
	      
	    }
	  else
	    {
	      Coef = 0.0;
	    }
	}
      if (Coef != 0.0)
	{				
	  TmpState = 0ul;
	  unsigned long Mask;
	  unsigned long Sign = 0ul;
	  Bool = true;
	  for(int i = 0; ( i < this->NbrFermions ) && ( Bool == true ); i++)
	    {
	      Mask= 1ul << State[i];
	      if ( (TmpState & Mask) != 0ul)
		Bool = false;
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef _64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  if (Bool == true)
	    {
	      if ((Sign & 0x1ul) != 0ul)
		Coef *= -1.0;
	      
	      InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState, Coef));
	      
	      if (InsertionResult.second == false)
		{
		  InsertionResult.first->second += Coef;
		}
	    }
	}
    }
  while (std::prev_permutation(monomial,monomial+this->NbrFermions));
  delete [] State;
}



// compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::FermionicStateTimeFermionicState(RealVector& fermionState1, RealVector& fermionState2, RealVector& outputVector, 
								      FermionOnSphereTwoLandauLevels*  fermionSpace2 , BosonOnSphereShort* finalSpace, int firstComponent,int nbrComponent)
{
  map<unsigned long , double> SortingMap;
  map<unsigned long , double>::iterator It;
  
  unsigned long* Slater1 = new unsigned long[this->NbrFermions];
  unsigned long* Slater2 = new unsigned long[this->NbrFermions];

  int NbrMax = firstComponent + nbrComponent;
		
  FactorialCoefficient Coefficient;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState1[j] != 0)
	{
	  this->ConvertToMonomial(this->StateDescription[j], Slater1);
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(fermionState2[i] != 0)
		{
			fermionSpace2->ConvertToMonomial(fermionSpace2->StateDescription[j], Slater2);
		  this->SecondLandauLevelSlaterTimesSlaterProjection(Slater1,Slater2, SortingMap, finalSpace);
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
		    {
		      int FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
		      while ((( ((*It).first) >> FTmpLzMax) & 0x1ul) == 0x0ul)
			--FTmpLzMax;
		      finalSpace->FermionToBoson((*It).first, FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
		      Coefficient.SetToOne();
		      for(int p = 0; p < finalSpace->TemporaryStateLzMax + 1; p++)
			{
			  Coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			}
		      outputVector[finalSpace->FermionBasis->FindStateIndex((*It).first, FTmpLzMax)] += fermionState2[i] * fermionState1[j] * (*It).second * Coefficient.GetIntegerValue();
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}

// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in three Landau levels
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereTwoLandauLevels::SecondLandauLevelSlaterTimesSlaterProjection(unsigned long* slater1, unsigned long* slater2, map <unsigned long, double> & sortingMap, 
										  BosonOnSphereShort* finalSpace)
{
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0ul;
  
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;
  int LzMax2 = finalSpace->LzMax + 4 - this->LzMax;
  double Coef = 1.0;
  long PowerIn1;
  long PowerIn2;
  long PowerOut;
  long Numerator;
  long AlphaIn = this->LzMax * (this->LzMax-1);
  
  long FinalDenominator = (finalSpace->LzMax+4) * (finalSpace->LzMax+3);
  long InitialDenominator = LzMax2 * this->LzMax;
  long BetaDenominator = InitialDenominator * (finalSpace->LzMax+4);
  
  
  long TmpLzMaxUp = this->LzMaxUp;
  long TmpFinalLzMaxUp = 2l + finalSpace->LzMax;
  long InverseFactor = TmpFinalLzMaxUp * TmpLzMaxUp;
  
  
  for (int i = 0; (i < this->NbrFermions) && (Coef != 0.0); i++)
    {
      State[i] = (slater1[i]>>1) + (slater2[i]>>1);
      
      PowerIn1 = (slater1[i]>>1);
      PowerIn2 = (slater2[i]>>1);
      PowerOut = State[i];
      if((slater1[i] & 0x1ul) != 0ul)
	{
	  if((slater2[i] & 0x1ul) != 0ul)
	    {
	      long Beta = PowerIn1 * (LzMax2) * (finalSpace->LzMax + 0x4l) + PowerIn2 * (this->LzMax)  * (finalSpace->LzMax + 0x4l) - 0x2ul * PowerOut * (LzMax2) * (this->LzMax);
	      Numerator = PowerIn1 * PowerIn2 * FinalDenominator * BetaDenominator * (finalSpace->LzMax + 2) -  (PowerOut-0x1ul)* PowerOut * InitialDenominator*BetaDenominator * (finalSpace->LzMax + 2) + Beta * (PowerOut-0x1ul) * FinalDenominator*InitialDenominator;
	      if (Numerator == 0x0l)
		Coef = 0.0;
	      else
		Coef *= ((double)Numerator/((double) FinalDenominator*InitialDenominator*(finalSpace->LzMax + 2)*BetaDenominator));
	      
	    }
	  else
	    {
	      
	      Numerator = -(PowerIn1 * (2l + finalSpace->LzMax)) + ((PowerOut - 1) * this->LzMaxUp);
	      if (Numerator == 0l)
		Coef = 0.0;
	      else
		{
		  Coef *= ((double) Numerator);
		  Coef /= ((double) (2l + finalSpace->LzMax)* this->LzMaxUp);
		}
	    }
	}
      else
	{
	  if ((slater2[i] & 0x1ul) != 0ul)
	    {
	      Numerator = -(PowerIn2 * (2l + finalSpace->LzMax)) + ((PowerOut - 1) * LzMax2);
	      if (Numerator == 0l)
		Coef = 0.0;
	      else
		{
		  Coef *= ((double) Numerator);
		  Coef /= ((double) (2l + finalSpace->LzMax)* LzMax2);
		}
	    }
	}
      
    }
  
  unsigned long Mask=0ul;
  unsigned long Sign = 0ul;
  if(Coef != 0.0)
    {
      
      for (int i = 0; (i < this->NbrFermions) ; i++)
	{
	  State[i] -= 2;
	}
      for (int i = 0; (i < this->NbrFermions) ; i++)
	
	InsertionResult = sortingMap.insert (pair<unsigned long,double> (finalSpace->ConvertFromMonomial(State), Coef));
      if (InsertionResult.second == false)
	{
	  InsertionResult.first->second += Coef;
	}
    }
  while (std::prev_permutation(slater2, slater2 + this->NbrFermions))
    {
      Coef = 1.0;
      for (int i = 0; (i < this->NbrFermions) && (Coef != 0.0); i++)
	{
	  State[i] = (slater1[i]>>1) + (slater2[i]>>1);
	  PowerIn1 = (slater1[i]>>1);
	  PowerIn2 = (slater2[i]>>1);
	  PowerOut = State[i];
	  
	  if ((slater1[i] & 0x1ul) != 0ul)
	    {
	      if ((slater2[i] & 0x1ul) != 0ul)
		{
		  long Beta = PowerIn1 * (LzMax2) * (finalSpace->LzMax + 4) + PowerIn2 * (this->LzMax)  * (finalSpace->LzMax + 4) - 2 * PowerOut * (LzMax2) * (this->LzMax);
					Numerator = PowerIn1 * PowerIn2 * FinalDenominator * BetaDenominator * (finalSpace->LzMax + 2) -  (PowerOut-0x1ul)* PowerOut * InitialDenominator*BetaDenominator * (finalSpace->LzMax + 2) + Beta * (PowerOut-0x1ul) * FinalDenominator*InitialDenominator;
					if(Numerator == 0x0l)
					  Coef = 0.0;
					else
					  Coef *= ((double)Numerator/((double) FinalDenominator*InitialDenominator*(finalSpace->LzMax + 2)*BetaDenominator));
					
		}
	      else
		{
		  
		  Numerator = -(PowerIn1 * (2l + finalSpace->LzMax)) + ((PowerOut - 1) * this->LzMaxUp);
		  if (Numerator == 0l)
		    Coef = 0.0;
		  else
		    {
		      Coef *= ((double) Numerator);
		      Coef /= ((double) (2l + finalSpace->LzMax)* this->LzMaxUp);
		    }
		}
	    }
	  else
	    {
	      if ((slater2[i] & 0x1ul) != 0ul)
		{
		  Numerator = -(PowerIn2 * (2l + finalSpace->LzMax)) + ((PowerOut - 1) * LzMax2);
		  if (Numerator == 0l)
		    Coef = 0.0;
		  else
		    {
		      Coef *= ((double) Numerator);
		      Coef /= ((double) (2l + finalSpace->LzMax)* LzMax2);
		    }
		}
	    }
	  
	}
      
      if( Coef != 0.0 )
	{
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i = 0; i < this->NbrFermions;i++)
	    {
	      State[i] -= 2;
				
	      Mask = (1ul << slater2[i]);
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  SortArrayDownOrdering(State,this->NbrFermions);
	  if((Sign & 0x1ul) != 0ul)
	    Coef *= -1.0;
	  InsertionResult = sortingMap.insert (pair<unsigned long,double> (finalSpace->ConvertFromMonomial(State), Coef));
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  delete [] State;
}
