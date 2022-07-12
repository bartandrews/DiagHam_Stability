////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                    class of bosons on sphere with SU(3) spin               //
//                                                                            //
//                        last modification : 03/12/2011                      //
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
#include "HilbertSpace/BosonOnSphereWithSU3Spin.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"

#include <math.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

BosonOnSphereWithSU3Spin::BosonOnSphereWithSU3Spin ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// totalTz = twice the total Tz value
// totalY = three time the total Y value
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU3Spin::BosonOnSphereWithSU3Spin (int nbrBosons, int totalLz, int lzMax, int totalTz, int totalY,
						    unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalY = totalY;
  this->TotalTz = totalTz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryState1 = new unsigned long[this->NbrLzValue];
  this->TemporaryState2 = new unsigned long[this->NbrLzValue];
  this->TemporaryState3 = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryState1;
  this->TemporaryStateSigma[1] = this->TemporaryState2;
  this->TemporaryStateSigma[2] = this->TemporaryState3;
  this->ProdATemporaryState1 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState2 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState3 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryState1;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryState2;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryState3;

  int N1 = (2 * nbrBosons) + totalY + (3 * totalTz);
  int N2 = (2 * nbrBosons) + totalY - (3 * totalTz);
  int N3 = nbrBosons - totalY;
  this->N1LzMax = this->LzMax + N1 - 1;
  this->N2LzMax = this->LzMax + N2 - 1;
  this->N3LzMax = this->LzMax + N3 - 1;
  this->FermionicLzMax = this->N1LzMax;
  if (this->N2LzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->N2LzMax;
  if (this->N3LzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->N3LzMax;
  if ((N1 < 0) || (N2 < 0) || (N3 < 0) || ((N1 % 6) != 0) || ((N2 % 6) != 0) || ((N3 % 3) != 0))
    this->LargeHilbertSpaceDimension = 0l;
  else
    {
      N1 /= 6;
      N2 /= 6;
      N3 /= 3;
      this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, N1, N2, N3);
    }

  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescription1 = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescription2 = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescription3 = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionSigma[0] = this->StateDescription1;
      this->StateDescriptionSigma[1] = this->StateDescription2;
      this->StateDescriptionSigma[2] = this->StateDescription3;
      long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 
							   N1, N2, N3, 0l);
      if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
	  cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU3Spin!" << endl;
	  exit(1);
	}
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;  
     SortTripleElementArrayDownOrdering<unsigned long>(this->StateDescription1, this->StateDescription2, this->StateDescription3, this->LargeHilbertSpaceDimension);
     this->GenerateLookUpTable(memory);
//      for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
//        {
// 	 cout << i << " : ";
// 	 this->PrintState(cout, i);
// 	 cout << this->FindStateIndex(this->StateDescription1[i], this->StateDescription2[i], this->StateDescription3[i]);
// 	 cout << endl;
//        }
#ifdef __DEBUG__
      int UsedMemory = 0;
      UsedMemory += this->HilbertSpaceDimension * (3 * sizeof(unsigned long));
      cout << "memory requested for Hilbert space = ";
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

BosonOnSphereWithSU3Spin::BosonOnSphereWithSU3Spin(const BosonOnSphereWithSU3Spin& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->N1LzMax = bosons.N1LzMax;
  this->N2LzMax = bosons.N2LzMax;
  this->N3LzMax = bosons.N3LzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryState1 = new unsigned long[this->NbrLzValue];
  this->TemporaryState2 = new unsigned long[this->NbrLzValue];
  this->TemporaryState3 = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryState1;
  this->TemporaryStateSigma[1] = this->TemporaryState2;
  this->TemporaryStateSigma[2] = this->TemporaryState3;
  this->ProdATemporaryState1 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState2 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState3 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryState1;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryState2;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryState3;
  this->StateDescription1 = bosons.StateDescription1;
  this->StateDescription2 = bosons.StateDescription2;
  this->StateDescription3 = bosons.StateDescription3; 
  this->StateDescriptionSigma[0] = this->StateDescription1;
  this->StateDescriptionSigma[1] = this->StateDescription2;
  this->StateDescriptionSigma[2] = this->StateDescription3;
  this->NbrUniqueStateDescription1 = bosons.NbrUniqueStateDescription1;
  this->UniqueStateDescription1 = bosons.UniqueStateDescription1;
  this->UniqueStateDescriptionSubArraySize1 = bosons.UniqueStateDescriptionSubArraySize1;
  this->NbrUniqueStateDescription2 = bosons.NbrUniqueStateDescription2;
  this->UniqueStateDescription2 = bosons.UniqueStateDescription2;
  this->UniqueStateDescriptionSubArraySize2 = bosons.UniqueStateDescriptionSubArraySize2;
  this->FirstIndexUniqueStateDescription2 = bosons.FirstIndexUniqueStateDescription2;
}

// destructor
//

BosonOnSphereWithSU3Spin::~BosonOnSphereWithSU3Spin ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription1;
      delete[] this->StateDescription2;
      delete[] this->StateDescription3;
      delete[] this->UniqueStateDescription1;
      delete[] this->UniqueStateDescriptionSubArraySize1;
      delete[] this->NbrUniqueStateDescription2;
      for (long i = 0l; i < this->NbrUniqueStateDescription1; ++i)
	{
	  delete[] this->UniqueStateDescription2[i];
	  delete[] this->UniqueStateDescriptionSubArraySize2[i];
	  delete[] this->FirstIndexUniqueStateDescription2[i];
	}
      delete[] this->UniqueStateDescription2;
      delete[] this->UniqueStateDescriptionSubArraySize2;
      delete[] this->FirstIndexUniqueStateDescription2;
    }
  delete[] this->TemporaryState1;
  delete[] this->TemporaryState2;
  delete[] this->TemporaryState3;
  delete[] this->ProdATemporaryState1;
  delete[] this->ProdATemporaryState2;
  delete[] this->ProdATemporaryState3;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU3Spin& BosonOnSphereWithSU3Spin::operator = (const BosonOnSphereWithSU3Spin& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription1;
      delete[] this->StateDescription2;
      delete[] this->StateDescription3;
      delete[] this->UniqueStateDescription1;
      delete[] this->UniqueStateDescriptionSubArraySize1;
      delete[] this->NbrUniqueStateDescription2;
      for (long i = 0l; i < this->NbrUniqueStateDescription1; ++i)
	{
	  delete[] this->UniqueStateDescription2[i];
	  delete[] this->UniqueStateDescriptionSubArraySize2[i];
	  delete[] this->FirstIndexUniqueStateDescription2[i];
	}
      delete[] this->UniqueStateDescription2;
      delete[] this->UniqueStateDescriptionSubArraySize2;
      delete[] this->FirstIndexUniqueStateDescription2;
    }
  delete[] this->TemporaryState1;
  delete[] this->TemporaryState2;
  delete[] this->TemporaryState3;
  delete[] this->ProdATemporaryState1;
  delete[] this->ProdATemporaryState2;
  delete[] this->ProdATemporaryState3;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->N1LzMax = bosons.N1LzMax;
  this->N2LzMax = bosons.N2LzMax;
  this->N3LzMax = bosons.N3LzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryState1 = new unsigned long[this->NbrLzValue];
  this->TemporaryState2 = new unsigned long[this->NbrLzValue];
  this->TemporaryState3 = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryState1;
  this->TemporaryStateSigma[1] = this->TemporaryState2;
  this->TemporaryStateSigma[2] = this->TemporaryState3;
  this->ProdATemporaryState1 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState2 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState3 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryState1;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryState2;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryState3;
  this->StateDescription1 = bosons.StateDescription1;
  this->StateDescription2 = bosons.StateDescription2;
  this->StateDescription3 = bosons.StateDescription3;
  this->StateDescriptionSigma[0] = this->StateDescription1;
  this->StateDescriptionSigma[1] = this->StateDescription2;
  this->StateDescriptionSigma[2] = this->StateDescription3;
  this->NbrUniqueStateDescription1 = bosons.NbrUniqueStateDescription1;
  this->UniqueStateDescription1 = bosons.UniqueStateDescription1;
  this->UniqueStateDescriptionSubArraySize1 = bosons.UniqueStateDescriptionSubArraySize1;
  this->NbrUniqueStateDescription2 = bosons.NbrUniqueStateDescription2;
  this->UniqueStateDescription2 = bosons.UniqueStateDescription2;
  this->UniqueStateDescriptionSubArraySize2 = bosons.UniqueStateDescriptionSubArraySize2;
  this->FirstIndexUniqueStateDescription2 = bosons.FirstIndexUniqueStateDescription2;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereWithSU3Spin::Clone()
{
  return new BosonOnSphereWithSU3Spin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereWithSU3Spin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereWithSU3Spin::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereWithSU3Spin::ExtractSubspace (AbstractQuantumNumber& q, 
								 SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_1 a_m_1

double  BosonOnSphereWithSU3Spin::Ad1A1 (int index, int m)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->TemporaryState1);
  return (double) (this->TemporaryState1[m]);  
}

// apply a^+_m_2 a_m_2 operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_2 a_m_2

double BosonOnSphereWithSU3Spin::Ad2A2 (int index, int m)
{
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->TemporaryState2);
  return (double) (this->TemporaryState2[m]);  
}
 
// apply a^+_m_3 a_m_3 operator to a given state (only state 3 Tz=0, Y=-2/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_3 a_m_3

double BosonOnSphereWithSU3Spin::Ad3A3 (int index, int m)
{
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->TemporaryState3);
  return (double) (this->TemporaryState3[m]);  
}

// find state index
//
// stateDescription1 = unsigned integer describing the fermionic state for type 1 particles
// stateDescription2 = unsigned integer describing the fermionic state for type 2 particles
// stateDescription3 = unsigned integer describing the fermionic state for type 3 particles
// return value = corresponding index

int BosonOnSphereWithSU3Spin::FindStateIndex(unsigned long stateDescription1, unsigned long stateDescription2, unsigned long stateDescription3)
{
//  cout << hex << stateDescription1 << " " << stateDescription2 << " " << stateDescription3 << dec << endl;
  int PosMin = 0;
  int PosMax = this->NbrUniqueStateDescription1 - 1;
  int PosMid = (PosMin + PosMax) >> 1;
//  cout << "entering " << PosMin << " " << PosMax << endl;
  unsigned long CurrentState = this->UniqueStateDescription1[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescription1))
    {
       if (CurrentState > stateDescription1)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->UniqueStateDescription1[PosMid];
    }
  if (CurrentState != stateDescription1)
    PosMid = PosMax;
//  cout << "pass 1 : " << PosMid << endl;
  unsigned long* TmpStateDescriptionArray = this->UniqueStateDescription2[PosMid];
  int* TmpFirstIndexUniqueStateDescription2 = this->FirstIndexUniqueStateDescription2[PosMid];
  int* TmpUniqueStateDescriptionSubArraySize2 = this->UniqueStateDescriptionSubArraySize2[PosMid];
  PosMin = 0;
  PosMax = this->NbrUniqueStateDescription2[PosMid] - 1;
  PosMid = (PosMin + PosMax) >> 1;
//  cout << "entring pass 2 : " << PosMin << " " << PosMax << endl;
  CurrentState = TmpStateDescriptionArray[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescription2))
    {
       if (CurrentState > stateDescription2)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = TmpStateDescriptionArray[PosMid];
    }
  if (CurrentState != stateDescription2)
    PosMid = PosMax;
  PosMin = TmpFirstIndexUniqueStateDescription2[PosMid];
  PosMax = PosMin + TmpUniqueStateDescriptionSubArraySize2[PosMid] - 1;
//  cout << "pass2 : " << PosMin << " " << PosMax << endl;
  PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescription3[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescription3))
    {
       if (CurrentState > stateDescription3)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->StateDescription3[PosMid];
    }
  if (CurrentState != stateDescription3)
    return PosMax;
  else
    return PosMid;
}



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSU3Spin::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescription1[state], this->StateDescription2[state], this->StateDescription3[state],
		       this->TemporaryState1, this->TemporaryState2, this->TemporaryState3); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = this->LzMax; i >=0 ; --i)
    {
      Str << "(" << this->TemporaryState1[i] << "," << this->TemporaryState2[i] << "," << this->TemporaryState3[i] << ") | ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU3Spin::GenerateStates(int nbrBosons, int lzMax1, int lzMax2, int lzMax3, int totalLz, int nbrN1, int nbrN2, int nbrN3, long pos)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0))
    return pos;
  if ((nbrBosons == 0) && (totalLz == 0))
    {
      this->StateDescription1[pos] = 0x0ul;
      this->StateDescription2[pos] = 0x0ul;
      this->StateDescription3[pos] = 0x0ul;
      return (pos + 1l);
    }
  if ((lzMax1 < 0) || (lzMax2 < 0) || (lzMax3 < 0))
    return pos;

  if (nbrBosons == 1) 
    {
      if ((nbrN3 == 1) && (lzMax3 >= totalLz))
	{
	  this->StateDescription1[pos] = 0x0ul;
	  this->StateDescription2[pos] = 0x0ul;
	  this->StateDescription3[pos] = 0x1ul << totalLz;
	  return (pos + 1l);
	}
      return pos;
    }

  long TmpPos;
  unsigned long Mask1;
  unsigned long Mask2;
  unsigned long Mask3;

  if (nbrN1 == 0)
    {
      if (nbrN2 == 0)
	{
	  for (int k = nbrN3; k > 0; --k)
	    {
	      TmpPos = this->GenerateStates(nbrBosons - k, 0, 0, lzMax3 - 1, totalLz - (lzMax3 * k), 
					    0, 0, nbrN3 - k, pos); 
	      Mask3 = ((0x1ul << k) - 1ul) << (lzMax3 + nbrN3 - k);
	      for (; pos < TmpPos; ++pos)
		{
		  this->StateDescription3[pos] |= Mask3;
		}
	    }
	  pos = this->GenerateStates(nbrBosons, 0, 0, lzMax3 - 1, totalLz, 0, 0, nbrN3, pos);
	  return pos;
	}
      TmpPos = this->GenerateStates(nbrBosons - nbrN2, 0, 0, lzMax3, totalLz - (lzMax2 * nbrN2), 
				    0, 0, nbrN3, pos); 
      Mask2 = ((0x1ul << nbrN2) - 1ul) << lzMax2;
      for (; pos < TmpPos; ++pos)
	this->StateDescription2[pos] |= Mask2;
      for (int j = nbrN2 - 1; j > 0; --j)
	{
	  TmpPos = this->GenerateStates(nbrBosons - j, 0, lzMax2 - 1, lzMax3, totalLz - (lzMax2 * j), 
					0, nbrN2 - j, nbrN3, pos); 
	  Mask2 = ((0x1ul << j) - 1ul) << (lzMax2 + nbrN2 - j);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription2[pos] |= Mask2;
	}
      pos = this->GenerateStates(nbrBosons, 0, lzMax2 - 1, lzMax3, totalLz, 0, nbrN2, nbrN3, pos);
      return pos;
    }
  
  TmpPos = this->GenerateStates(nbrBosons - nbrN1, 0, lzMax2, lzMax3, totalLz - (lzMax1 * nbrN1), 
				0, nbrN2, nbrN3, pos); 
  Mask1 = ((0x1ul << nbrN1) - 1ul) << lzMax1;
  for (; pos < TmpPos; ++pos)
    this->StateDescription1[pos] |= Mask1;
  for (int i = nbrN1 - 1; i > 0; --i)
    {
      TmpPos = this->GenerateStates(nbrBosons - i, lzMax1 - 1, lzMax2, lzMax3, totalLz - (lzMax1 * i), 
				    nbrN1 - i, nbrN2, nbrN3, pos); 
      Mask1 = ((0x1ul << i) - 1ul) << (lzMax1 + nbrN1 - i);
      for (; pos < TmpPos; ++pos)
	{
	  this->StateDescription1[pos] |= Mask1;
	}
    }
  pos = this->GenerateStates(nbrBosons, lzMax1 - 1, lzMax2, lzMax3, totalLz, nbrN1, nbrN2, nbrN3, pos);
  return pos;
};


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSU3Spin::GenerateLookUpTable(unsigned long memory)
{  
  long TmpUniquePartition = 1l;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescription1[i - 1] == this->StateDescription1[i]))
	{
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	++TmpUniquePartition;
    }

  this->NbrUniqueStateDescription1 = TmpUniquePartition;
  this->UniqueStateDescription1 = new unsigned long [this->NbrUniqueStateDescription1];
  this->UniqueStateDescriptionSubArraySize1 = new int [this->NbrUniqueStateDescription1];
  TmpUniquePartition = 0l;
  this->UniqueStateDescription1[0l] = this->StateDescription1[0l];
  this->UniqueStateDescriptionSubArraySize1[0l] = 1;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescription1[i - 1] == this->StateDescription1[i]))
	{
	  ++this->UniqueStateDescriptionSubArraySize1[TmpUniquePartition];
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	{
	  ++TmpUniquePartition;
	  this->UniqueStateDescription1[TmpUniquePartition] = this->StateDescription1[i];
	  this->UniqueStateDescriptionSubArraySize1[TmpUniquePartition] = 1; 
	}
    }

  this->NbrUniqueStateDescription2 = new int [this->NbrUniqueStateDescription1];
  TmpUniquePartition = 0;
  long TmpIndex = 0l;
  while (TmpIndex < this->LargeHilbertSpaceDimension)
    {
      long Lim = TmpIndex + this->UniqueStateDescriptionSubArraySize1[TmpUniquePartition];
      this->NbrUniqueStateDescription2[TmpUniquePartition] = 1;
      ++TmpIndex;
      while (TmpIndex < Lim)
	{
	  while ((TmpIndex < Lim) && (this->StateDescription2[TmpIndex - 1] == this->StateDescription2[TmpIndex]))
	    ++TmpIndex;
	  if (TmpIndex < Lim)
	    {
	      ++this->NbrUniqueStateDescription2[TmpUniquePartition];
	      ++TmpIndex;
	    }
	}
      ++TmpUniquePartition;
    }
  this->UniqueStateDescription2 = new unsigned long* [this->NbrUniqueStateDescription1];
  this->UniqueStateDescriptionSubArraySize2 = new int* [this->NbrUniqueStateDescription1];
  this->FirstIndexUniqueStateDescription2 = new int* [this->NbrUniqueStateDescription1];
  for (long i = 0l; i < this->NbrUniqueStateDescription1; ++i)
    {
      this->UniqueStateDescription2[i] = new unsigned long [this->NbrUniqueStateDescription2[i]];
      this->UniqueStateDescriptionSubArraySize2[i] = new int [this->NbrUniqueStateDescription2[i]];
      this->FirstIndexUniqueStateDescription2[i] = new int [this->NbrUniqueStateDescription2[i]];
    }

  TmpUniquePartition = 0;
  TmpIndex = 0l;
  while (TmpIndex < this->LargeHilbertSpaceDimension)
    {
      long Lim = TmpIndex + this->UniqueStateDescriptionSubArraySize1[TmpUniquePartition];
      int TmpUniquePartition2 = 0;
      this->UniqueStateDescription2[TmpUniquePartition][TmpUniquePartition2] = this->StateDescription2[TmpIndex];
      this->UniqueStateDescriptionSubArraySize2[TmpUniquePartition][TmpUniquePartition2] = 1;
      this->FirstIndexUniqueStateDescription2[TmpUniquePartition][TmpUniquePartition2] = TmpIndex;
      ++TmpIndex;
      while (TmpIndex < Lim)
	{
	  while ((TmpIndex < Lim) && (this->StateDescription2[TmpIndex - 1] == this->StateDescription2[TmpIndex]))
	    {
	      ++this->UniqueStateDescriptionSubArraySize2[TmpUniquePartition][TmpUniquePartition2];	      
	      ++TmpIndex;
	    }
	  if (TmpIndex < Lim)
	    {
	      ++TmpUniquePartition2;
	      this->UniqueStateDescription2[TmpUniquePartition][TmpUniquePartition2] = this->StateDescription2[TmpIndex];
	      this->UniqueStateDescriptionSubArraySize2[TmpUniquePartition][TmpUniquePartition2] = 1;
	      this->FirstIndexUniqueStateDescription2[TmpUniquePartition][TmpUniquePartition2] = TmpIndex;
	      ++TmpIndex;
	    }
	}
      ++TmpUniquePartition;
    }
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// return value = Hilbert space dimension

long BosonOnSphereWithSU3Spin::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0))
    return 0l;
  if ((nbrBosons == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0)
    return 0l;
  if (nbrBosons == 1)
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }
  long Tmp = 0l;
  for (int i = nbrN1; i >= 0; --i)
    for (int j = nbrN2; j >= 0; --j)
      for (int k = nbrN3; k >= 0; --k)
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j + k), lzMax - 1, totalLz - (lzMax * (i + j + k)), 
							  nbrN1 - i, nbrN2 - j, nbrN3 - k);
  return  Tmp;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// totalTz = twice the total Tz value
// totalY = three time the total Y value
// return value = Hilbert space dimension

long BosonOnSphereWithSU3Spin::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalTz, int totalY)
{
  int N1 = (2 * nbrBosons) + totalY + (3 * totalTz);
  int N2 = (2 * nbrBosons) + totalY - (3 * totalTz);
  int N3 = nbrBosons - totalY;
  if ((N1 >= 0) && (N2 >= 0) && (N3 >= 0) && ((N1 % 6) == 0) && ((N2 % 6) == 0) && ((N3 % 3) == 0))
    {
      N1 /= 6;
      N2 /= 6;
      N3 /= 3;
      return this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, (totalLz + (nbrBosons * lzMax)) >> 1, N1, N2, N3);
    }
  else
    return 0l;
}

// apply a^+_m_1 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU3Spin::Ad1A1 (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->TemporaryState1);
  if (this->TemporaryState1[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState1[n];
  --this->TemporaryState1[n];
  ++this->TemporaryState1[m];
  coefficient *= (double) this->TemporaryState1[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState1), this->StateDescription2[index], 
			      this->StateDescription3[index]);  
}

// apply a^+_m_1 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU3Spin::Ad1A2 (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->TemporaryState1);
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->TemporaryState2);
  if (this->TemporaryState2[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState2[n];
  --this->TemporaryState2[n];
  ++this->TemporaryState1[m];
  coefficient *= (double) this->TemporaryState1[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState1), this->BosonToFermion(this->TemporaryState2), 
			      this->StateDescription3[index]);  
}

// apply a^+_m_1 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU3Spin::Ad1A3 (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->TemporaryState1);
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->TemporaryState3);
  if (this->TemporaryState3[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState3[n];
  --this->TemporaryState3[n];
  ++this->TemporaryState1[m];
  coefficient *= (double) this->TemporaryState1[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState1), this->StateDescription2[index],
			      this->BosonToFermion(this->TemporaryState3));  
}

// apply a^+_m_2 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU3Spin::Ad2A1 (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->TemporaryState1);
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->TemporaryState2);
  if (this->TemporaryState1[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState1[n];
  --this->TemporaryState1[n];
  ++this->TemporaryState2[m];
  coefficient *= (double) this->TemporaryState2[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState1), this->BosonToFermion(this->TemporaryState2), 
			      this->StateDescription3[index]);  
}

// apply a^+_m_2 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU3Spin::Ad2A2 (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->TemporaryState2);
  if (this->TemporaryState2[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState2[n];
  --this->TemporaryState2[n];
  ++this->TemporaryState2[m];
  coefficient *= (double) this->TemporaryState2[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescription1[index], this->BosonToFermion(this->TemporaryState2), 
			      this->StateDescription3[index]);  
}

// apply a^+_m_2 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU3Spin::Ad2A3 (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->TemporaryState3);
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->TemporaryState2);
  if (this->TemporaryState3[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState3[n];
  --this->TemporaryState3[n];
  ++this->TemporaryState2[m];
  coefficient *= (double) this->TemporaryState2[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescription1[index], this->BosonToFermion(this->TemporaryState2), 
			      this->BosonToFermion(this->TemporaryState3));  
}

// apply a^+_m_3 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU3Spin::Ad3A1 (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->TemporaryState1);
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->TemporaryState3);
  if (this->TemporaryState1[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState1[n];
  --this->TemporaryState1[n];
  ++this->TemporaryState3[m];
  coefficient *= (double) this->TemporaryState3[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState1), this->StateDescription2[index],
			      this->BosonToFermion(this->TemporaryState3));  
}

// apply a^+_m_3 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU3Spin::Ad3A2 (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->TemporaryState3);
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->TemporaryState2);
  if (this->TemporaryState2[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState2[n];
  --this->TemporaryState2[n];
  ++this->TemporaryState3[m];
  coefficient *= (double) this->TemporaryState3[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescription1[index], this->BosonToFermion(this->TemporaryState2), 
			      this->BosonToFermion(this->TemporaryState3));  
}

// apply a^+_m_3 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU3Spin::Ad3A3 (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->TemporaryState3);
  if (this->TemporaryState3[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState3[n];
  --this->TemporaryState3[n];
  ++this->TemporaryState3[m];
  coefficient *= (double) this->TemporaryState3[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescription1[index], this->StateDescription2[index], 
			      this->BosonToFermion(this->TemporaryState3));  
}

// apply a_n1_1 a_n2_1 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A1A1 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->ProdATemporaryState1);
  if ((this->ProdATemporaryState1[n1] == 0) || (this->ProdATemporaryState1[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState1[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->ProdATemporaryState2);
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->ProdATemporaryState3);
  double Coefficient = this->ProdATemporaryState1[n2];
  --this->ProdATemporaryState1[n2];
  Coefficient *= this->ProdATemporaryState1[n1];
  --this->ProdATemporaryState1[n1];
  return sqrt(Coefficient);
}

// apply a_n1_1 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A1A2 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->ProdATemporaryState1);
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->ProdATemporaryState2);
  if ((this->ProdATemporaryState1[n1] == 0) || (this->ProdATemporaryState2[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->ProdATemporaryState3);
  double Coefficient = this->ProdATemporaryState2[n2];
  --this->ProdATemporaryState2[n2];
  Coefficient *= this->ProdATemporaryState1[n1];
  --this->ProdATemporaryState1[n1];
  return sqrt(Coefficient);
}

// apply a_n1_1 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A1A3 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->ProdATemporaryState1);
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->ProdATemporaryState3);
  if ((this->ProdATemporaryState1[n1] == 0) || (this->ProdATemporaryState3[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->ProdATemporaryState2);
  double Coefficient = this->ProdATemporaryState3[n2];
  --this->ProdATemporaryState3[n2];
  Coefficient *= this->ProdATemporaryState1[n1];
  --this->ProdATemporaryState1[n1];
  return sqrt(Coefficient);
}

// apply a_n1_2 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A2A2 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->ProdATemporaryState2);
  if ((this->ProdATemporaryState2[n1] == 0) || (this->ProdATemporaryState2[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryState2[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->ProdATemporaryState1);
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->ProdATemporaryState3);
  double Coefficient = this->ProdATemporaryState2[n2];
  --this->ProdATemporaryState2[n2];
  Coefficient *= this->ProdATemporaryState2[n1];
  --this->ProdATemporaryState2[n1];
  return sqrt(Coefficient);
}

// apply a_n1_2 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A2A3 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->ProdATemporaryState2);
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->ProdATemporaryState3);
  if ((this->ProdATemporaryState2[n1] == 0) || (this->ProdATemporaryState3[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->ProdATemporaryState1);
  double Coefficient = this->ProdATemporaryState3[n2];
  --this->ProdATemporaryState3[n2];
  Coefficient *= this->ProdATemporaryState2[n1];
  --this->ProdATemporaryState2[n1];
  return sqrt(Coefficient);
}

// apply a_n1_3 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A3A3 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->ProdATemporaryState3);
  if ((this->ProdATemporaryState3[n1] == 0) || (this->ProdATemporaryState3[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryState3[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->ProdATemporaryState1);
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->ProdATemporaryState2);
  double Coefficient = this->ProdATemporaryState3[n2];
  --this->ProdATemporaryState3[n2];
  Coefficient *= this->ProdATemporaryState3[n1];
  --this->ProdATemporaryState3[n1];
  return sqrt(Coefficient);
}

// apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad1Ad1 (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryState1, this->TemporaryState1, coefficient);
}

// apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad1Ad2 (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryState1, this->TemporaryState2, coefficient);
}

// apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad1Ad3 (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryState1, this->TemporaryState3, coefficient);
}

// apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad2Ad2 (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryState2, this->TemporaryState2, coefficient);
}

// apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad2Ad3 (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryState2, this->TemporaryState3, coefficient);
}

// apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad3Ad3 (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryState3, this->TemporaryState3, coefficient);
}

// convert a state from one SU(3) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void BosonOnSphereWithSU3Spin::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, 
						     long firstComponent, long nbrComponents)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU3Indices = new int [this->NbrBosons];
  int* TmpSU3Indices2 = new int [this->NbrBosons];
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  targetState.ClearVector();
  long LastComponent = firstComponent + nbrComponents;
  if (nbrComponents == 0)
    LastComponent = this->LargeHilbertSpaceDimension;
  for (long i = firstComponent; i < LastComponent; ++i)
    {
      this->FermionToBoson(this->StateDescription1[i], this->StateDescription2[i], this->StateDescription3[i],
			   this->TemporaryState1, this->TemporaryState2, this->TemporaryState3); 
      double OccupationCoefficient = 0.0;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (int l = 0; l < this->TemporaryState3[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU3Indices[TmpIndex] = 2;
	      ++TmpIndex;	      
	    }
	  for (int l = 0; l < this->TemporaryState2[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU3Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (int l = 0; l < this->TemporaryState1[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU3Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryState1[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryState2[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryState3[j]];
	}
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSU3Indices, TmpSU3Indices2, oneBodyBasis, OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] OccupationCoefficientArray;
  delete[] TmpMomentumIndices;
  delete[] TmpSU3Indices;
  delete[] TmpSU3Indices2;
}

// compute the transformation matrix from one SU(3) basis to another, transforming the one body basis in each momentum sector
//
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// return value = transformation matrix

ComplexMatrix BosonOnSphereWithSU3Spin::TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU3Indices = new int [this->NbrBosons];
  int* TmpSU3Indices2 = new int [this->NbrBosons];
  ComplexMatrix TmpMatrix(this->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->StateDescription1[i], this->StateDescription2[i], this->StateDescription3[i],
			   this->TemporaryState1, this->TemporaryState2, this->TemporaryState3); 
      int TmpIndex = 0;
      double OccupationCoefficient = 0.0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (int l = 0; l < this->TemporaryState3[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU3Indices[TmpIndex] = 2;
	      ++TmpIndex;	      
	    }
	  for (int l = 0; l < this->TemporaryState2[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU3Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (int l = 0; l < this->TemporaryState1[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU3Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryState1[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryState2[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryState3[j]];
	}
      this->TransformOneBodyBasisRecursive(TmpMatrix[i], 1.0, 0, TmpMomentumIndices, TmpSU3Indices, TmpSU3Indices2, oneBodyBasis,
					   OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSU3Indices;
  delete[] TmpSU3Indices2;
  delete[] OccupationCoefficientArray;
  return TmpMatrix;
}

// recursive part of the convertion from a state from one SU(3) basis to another, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSU3Indices = array that gives the spin dressing the initial n-body state
// currentSU3Indices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
// occupationCoefficientArray = array that provides 1/2 ln (N!)

void BosonOnSphereWithSU3Spin::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
							      int position, int* momentumIndices, int* initialSU3Indices, int* currentSU3Indices, ComplexMatrix* oneBodyBasis,
							      double occupationCoefficient, double* occupationCoefficientArray) 
{
  if (position == this->NbrBosons)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->TemporaryState1[i] = 0ul;
	  this->TemporaryState2[i] = 0ul; 
	  this->TemporaryState3[i] = 0ul;
	}
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  switch (currentSU3Indices[i])
	    {
	    case 0:
	      this->TemporaryState1[momentumIndices[i]]++;
	      break;
	    case 1:
	      this->TemporaryState2[momentumIndices[i]]++;
	      break;
	    case 2:
	      this->TemporaryState3[momentumIndices[i]]++;
	      break;
	    }
	}
      int Index = this->FindStateIndex(this->TemporaryState1, this->TemporaryState2, this->TemporaryState3);
      if (Index < this->HilbertSpaceDimension)
	{
	  
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryState1[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryState2[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryState3[i]];
	    }
	  targetState[Index] += coefficient * exp (occupationCoefficient);
	}
      return;      
    }
  else
    {
      currentSU3Indices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][2 - initialSU3Indices[position]][2]), position + 1, momentumIndices, initialSU3Indices, currentSU3Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU3Indices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][2 - initialSU3Indices[position]][1]), position + 1, momentumIndices, initialSU3Indices, currentSU3Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU3Indices[position] = 2;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][2 - initialSU3Indices[position]][0]), position + 1, momentumIndices, initialSU3Indices, currentSU3Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
    }
}

// compute the projection matrix from the SU(3) Hilbert space to an U(1) Hilbert space
// 
// targetSpace = pointer to the U(1) Hilbert space
// type = type of particles that has to be kept (0 for type 1, 1 for type 2, 2 for type 3
// return value = projection matrix

ComplexMatrix BosonOnSphereWithSU3Spin::TransformationMatrixSU3ToU1(BosonOnSphereShort* targetSpace, int type)
{
  ComplexMatrix TmpMatrix (targetSpace->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  unsigned long* TmpStateDescription;
  unsigned long* TmpStateDescriptionOther1;
  unsigned long* TmpStateDescriptionOther2;
  switch (type)
    {
    case 0:
      {
	TmpStateDescription = this->StateDescription1;
	TmpStateDescriptionOther1 = this->StateDescription2;
	TmpStateDescriptionOther2 = this->StateDescription3;
      }
      break;
    case 1:
      {
	TmpStateDescription = this->StateDescription2;
	TmpStateDescriptionOther1 = this->StateDescription1;
	TmpStateDescriptionOther2 = this->StateDescription3;
      }
      break;
    case 2:
      {
	TmpStateDescription = this->StateDescription3;
	TmpStateDescriptionOther1 = this->StateDescription2;
	TmpStateDescriptionOther2 = this->StateDescription1;
      }
      break;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if ((TmpStateDescriptionOther1[i] == 0x0ul) && (TmpStateDescriptionOther2[i] == 0x0ul))
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
