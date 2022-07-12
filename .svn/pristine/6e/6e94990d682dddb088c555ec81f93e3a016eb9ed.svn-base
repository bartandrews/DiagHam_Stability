////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                   class for bosons on sphere with SU(2) spin               //
//                     without contraint on the total momentum                //
//                                                                            //
//                        last modification : 30/03/2018                      //
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
#include "HilbertSpace/BosonOnSphereWithSU2SpinAllLz.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinTwoLandauLevels.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/Permutations.h"
#include "GeneralTools/Endian.h"
#include "Architecture/ArchitectureOperation/FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <map>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

BosonOnSphereWithSU2SpinAllLz::BosonOnSphereWithSU2SpinAllLz ()
{
}

// basic constructor without any constraint on Sz
// 
// nbrBosons = number of bosons
// lzMax = twice the maximum Lz value reached by a boson

BosonOnSphereWithSU2SpinAllLz::BosonOnSphereWithSU2SpinAllLz (int nbrBosons, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NbrBosonsUp = 0;
  this->NbrBosonsDown = 0;
  this->NUpLzMax = this->LzMax + this->NbrBosons - 1;
  this->NDownLzMax = this->LzMax + this->NbrBosons - 1;
  this->FermionicLzMax = this->NUpLzMax;
  if (this->NDownLzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->NDownLzMax;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax);
  
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, this->FermionicLzMax, this->FermionicLzMax, 0l);

  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2SpinAllLz!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  



  this->GenerateLookUpTable(10000000);

#ifdef __DEBUG__
   long UsedMemory = 0l;
   UsedMemory += this->LargeHilbertSpaceDimension * (2l * sizeof(unsigned long));
   cout << "memory requested for Hilbert space = ";
   if (UsedMemory >= 1024l)
    if (UsedMemory >= 1048576l)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// basic constructor
// 
// nbrBosons = number of bosons
// lzMax = twice the maximum Lz value reached by a boson
// totalSpin = twice the total spin value
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU2SpinAllLz::BosonOnSphereWithSU2SpinAllLz (int nbrBosons, int lzMax, int totalSpin, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NbrBosonsUp = this->NbrBosons + this->TotalSpin;
  this->NbrBosonsDown = this->NbrBosons - this->TotalSpin;
  if ((this->NbrBosonsUp < 0) || ((this->NbrBosonsUp & 0x1) != 0) ||
      (this->NbrBosonsDown < 0) || ((this->NbrBosonsDown & 0x1) != 0))
    this->LargeHilbertSpaceDimension = 0l;
  else
    {
      this->NbrBosonsUp >>= 1;
      this->NbrBosonsDown >>= 1;
      this->NUpLzMax = this->LzMax + this->NbrBosonsUp - 1;
      this->NDownLzMax = this->LzMax + this->NbrBosonsDown - 1;
      this->FermionicLzMax = this->NUpLzMax;
      if (this->NDownLzMax > this->FermionicLzMax)
	this->FermionicLzMax = this->NDownLzMax;
      this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, this->NbrBosonsUp, this->NbrBosonsDown);
    }
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, this->NbrBosonsUp, this->NbrBosonsDown, 0l);
//   for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
//     cout << i << " : " << hex << this->StateDescriptionUp[i] << " " << this->StateDescriptionDown[i] << dec << endl;
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2SpinAllLz!" << endl;
      exit(1);
    }
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  this->TargetSpace = this;

  this->GenerateLookUpTable(memory);
  
#ifdef __DEBUG__
  long UsedMemory = 0;
  UsedMemory += this->LargeHilbertSpaceDimension * (2l * sizeof(unsigned long));
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024l)
    if (UsedMemory >= 1048576l)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// constructor from a binary file that describes the Hilbert space
// 
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU2SpinAllLz::BosonOnSphereWithSU2SpinAllLz (char* fileName, unsigned long memory)
{
  this->ReadHilbertSpace(fileName);

  this->NbrLzValue = this->LzMax + 1;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->Flag.Initialize();
 
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NbrBosonsUp = this->NbrBosons + this->TotalSpin;
  this->NbrBosonsDown = this->NbrBosons - this->TotalSpin;
  this->NbrBosonsUp >>= 1;
  this->NbrBosonsDown >>= 1;
  this->NUpLzMax = this->LzMax + this->NbrBosonsUp - 1;
  this->NDownLzMax = this->LzMax + this->NbrBosonsDown - 1;
  this->FermionicLzMax = this->NUpLzMax;
  if (this->NDownLzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->NDownLzMax;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->TargetSpace = this;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += this->LargeHilbertSpaceDimension * (4l * sizeof(unsigned long));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024l)
	if (UsedMemory >= 1048576l)
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

BosonOnSphereWithSU2SpinAllLz::BosonOnSphereWithSU2SpinAllLz(const BosonOnSphereWithSU2SpinAllLz& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->TotalSpin = bosons.TotalSpin;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnSphereWithSU2SpinAllLz::~BosonOnSphereWithSU2SpinAllLz ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU2SpinAllLz& BosonOnSphereWithSU2SpinAllLz::operator = (const BosonOnSphereWithSU2SpinAllLz& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      delete[] this->UniqueStateDescriptionUp;
      delete[] this->UniqueStateDescriptionSubArraySizeUp;
      delete[] this->FirstIndexUniqueStateDescriptionUp;
    }
  delete[] this->TemporaryStateUp;
  delete[] this->TemporaryStateDown;
  delete[] this->ProdATemporaryStateUp;
  delete[] this->ProdATemporaryStateDown;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereWithSU2SpinAllLz::Clone()
{
  return new BosonOnSphereWithSU2SpinAllLz(*this);
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMaxUp = momentum maximum value for a boson in the state with up
// lzMaxDown = momentum maximum value for a boson in the state with down
// nbrNUp = number of particles with quantum up
// nbrNDown = number of particles with quantum down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU2SpinAllLz::GenerateStates(int nbrBosons, int lzMaxUp, int lzMaxDown, 
						   int nbrNUp, int nbrNDown, long pos)
{
  if ((nbrBosons < 0) || (nbrNUp < 0) || (nbrNDown < 0))
    return pos;
  if (nbrBosons == 0)
    {
      this->StateDescriptionUp[pos] = 0x0ul;
      this->StateDescriptionDown[pos] = 0x0ul;
      return (pos + 1l);
    }
  if ((lzMaxUp < 0) || (lzMaxDown < 0))
    return pos;

  long TmpPos;
  unsigned long MaskUp;
  unsigned long MaskDown;

  if (nbrNUp == 0)
    {
      for (int l = nbrNDown; l > 0; --l)
	{
	  TmpPos = this->GenerateStates(nbrBosons - l, 0, lzMaxDown - 1, 
					0, nbrNDown - l, pos); 
	  MaskDown = ((0x1ul << l) - 1ul) << (lzMaxDown + nbrNDown - l);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
      pos = this->GenerateStates(nbrBosons, 0, lzMaxDown - 1, 0, nbrNDown, pos);
      return pos;
    }
  
  TmpPos = this->GenerateStates(nbrBosons - nbrNUp, 0, lzMaxDown, 0, nbrNDown, pos); 
  MaskUp = ((0x1ul << nbrNUp) - 1ul) << lzMaxUp;
  for (; pos < TmpPos; ++pos)
    this->StateDescriptionUp[pos] |= MaskUp;
  for (int i = nbrNUp - 1; i > 0; --i)
    {
      TmpPos = this->GenerateStates(nbrBosons - i, lzMaxUp - 1, lzMaxDown, nbrNUp - i, nbrNDown, pos); 
      MaskUp = ((0x1ul << i) - 1ul) << (lzMaxUp + nbrNUp - i);
      for (; pos < TmpPos; ++pos)
	{
	  this->StateDescriptionUp[pos] |= MaskUp;
	}
    }
  pos = this->GenerateStates(nbrBosons, lzMaxUp - 1, lzMaxDown, nbrNUp, nbrNDown, pos);
  return pos;
};


// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU2SpinAllLz::GenerateStates(int nbrBosons, int lzMax, 
						   int currentFermionicPositionUp, int currentFermionicPositionDown, long pos)
{
  if (nbrBosons < 0)
    return pos;
  if (nbrBosons == 0)
    {
      this->StateDescriptionUp[pos] = 0x0ul;
      this->StateDescriptionDown[pos] = 0x0ul;
      return (pos + 1l);
    }
  if (lzMax < 0)
    return pos;

  for (int i = nbrBosons; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul)  << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, lzMax - 1, 
					     currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, pos); 
	  unsigned long MaskDown = ((0x1ul << j) - 0x1ul) << (currentFermionicPositionDown - j - 1);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescriptionUp[pos] |= MaskUp;
	      this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
    }
  return pos;
};


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long BosonOnSphereWithSU2SpinAllLz::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, 
									 int nbrNUp, int nbrNDown)
{
  if ((nbrBosons < 0) || (nbrNUp < 0) || (nbrNDown < 0))
    return 0l;
  if (nbrBosons == 0)
    return 1l;
  if (lzMax < 0)
    return 0l;
  if (nbrBosons == 1)
    {
      return (long) (lzMax + 1);
    }
  long Tmp = 0l;
  for (int i = nbrNUp; i >= 0; --i)
    {
      for (int j = nbrNDown; j >= 0; --j)
	{	  
 	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j), lzMax - 1, 
 							    nbrNUp - i, nbrNDown - j);
	}
    }
  return  Tmp;
}

// evaluate Hilbert space dimension without the Sz constraint
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// return value = Hilbert space dimension

long BosonOnSphereWithSU2SpinAllLz::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax)
{
  if (nbrBosons < 0)
    return 0l;
  if (nbrBosons == 0)
    return 1l;
  if (lzMax < 0)
    return 0l;
  if (nbrBosons == 1)
    {
      return (long) (2 * (lzMax + 1));
    }
  long Tmp = 0l;
  for (int i = nbrBosons; i >= 0; --i)
    Tmp += (i + 1) * this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - i, lzMax - 1);
  return Tmp;
}

