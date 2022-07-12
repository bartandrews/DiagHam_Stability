////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//           class for bosons on sphere with SU(2) spin and both the          //
//                       Lz<->-Lz and Sz<->-Sz symmetries                     //
//                                                                            //
//                        last modification : 25/09/2016                      //
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
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
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

BosonOnSphereWithSU2SpinLzSzSymmetry::BosonOnSphereWithSU2SpinLzSzSymmetry ()
{
  this->SzParitySign = 1.0;
}

// basic constructor without any constraint on Sz
// 
// nbrBosons = number of bosons
// lzMax = twice the maximum Lz value reached by a boson
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity

BosonOnSphereWithSU2SpinLzSzSymmetry::BosonOnSphereWithSU2SpinLzSzSymmetry (int nbrBosons, int lzMax, bool minusSzParity, bool minusLzParity)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  if (minusSzParity == false)
    this->SzParitySign = 1.0;
  else
    this->SzParitySign = -1.0;
  if (minusLzParity == false)
    this->LzParitySign = 1.0;
  else
    this->LzParitySign = -1.0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxOrbitSize = 4;
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
  this->LzSymmetryMaxSwapPosition = (this->LzMax - 1) >> 1;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1);
  
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 
						       this->LzMax + this->NbrBosons + 1, this->LzMax + this->NbrBosons + 1, 0l);

  for (long i = 0; i < TmpHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescriptionUp[i];
      int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[i]);
      this->StateDescriptionUp[i] >>= this->NbrBosons - TmpNbrUps; 
      this->StateDescriptionDown[i] >>= TmpNbrUps; 
    }
  SortDoubleElementArrayDownOrdering<unsigned long>(this->StateDescriptionUp, this->StateDescriptionDown, TmpHilbertSpaceDimension);

  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2SpinLzSzSymmetry!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->GenerateStatesWithDiscreteSymmetry();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(10000000);
      
#ifdef __DEBUG__
      long UsedMemory = 0l;
      UsedMemory += this->LargeHilbertSpaceDimension * ((2 * sizeof(unsigned long)) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024l)
	if (UsedMemory >= 1048576l)
	  cout << (UsedMemory >> 20) << "Mb" << endl;
	else
	  cout << (UsedMemory >> 10) << "kb" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// basic constructor
// 
// nbrBosons = number of bosons
// lzMax = twice the maximum Lz value reached by a boson
// totalSpin = twice the total spin value (not taken into account)
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU2SpinLzSzSymmetry::BosonOnSphereWithSU2SpinLzSzSymmetry (int nbrBosons, int lzMax, int totalSpin, bool minusSzParity, bool minusLzParity, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  if (minusSzParity == false)
    this->SzParitySign = 1.0;
  else
    this->SzParitySign = -1.0;
  if (minusLzParity == false)
    this->LzParitySign = 1.0;
  else
    this->LzParitySign = -1.0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->LzSymmetryMaxSwapPosition = (this->LzMax - 1) >> 1;
  this->MaxOrbitSize = 4;
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
      this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, this->NbrBosonsUp, this->NbrBosonsDown);
    }
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 
						       this->NbrBosonsUp, this->NbrBosonsDown, 0l);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2SpinLzSzSymmetry!" << endl;
      exit(1);
    }

  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->GenerateStatesWithDiscreteSymmetry();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;  
  this->TargetSpace = this;

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0l;
      UsedMemory += this->LargeHilbertSpaceDimension * (2l * sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024l)
	if (UsedMemory >= 1048576l)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// constructor from a binary file that describes the Hilbert space
// 
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU2SpinLzSzSymmetry::BosonOnSphereWithSU2SpinLzSzSymmetry (char* fileName, unsigned long memory)
{
  this->MaxOrbitSize = 4;
  this->ReadHilbertSpace(fileName);

  this->IncNbrBosons = this->NbrBosons + 1;
  this->NbrLzValue = this->LzMax + 1;
  this->LzSymmetryMaxSwapPosition = (this->LzMax - 1) >> 1;
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
    }
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;

  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;  

  this->TargetSpace = this;

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      int UsedMemory = 0;
      UsedMemory += this->HilbertSpaceDimension * (4 * sizeof(unsigned long));
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

BosonOnSphereWithSU2SpinLzSzSymmetry::BosonOnSphereWithSU2SpinLzSzSymmetry(const BosonOnSphereWithSU2SpinLzSzSymmetry& bosons)
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
  this->SzParitySign = bosons.SzParitySign;
  this->LzParitySign = bosons.LzParitySign;
  this->LzSymmetryMaxSwapPosition = bosons.LzSymmetryMaxSwapPosition;
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
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;
  this->RescalingFactors = bosons.RescalingFactors;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  this->MaxOrbitSize = bosons.MaxOrbitSize;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnSphereWithSU2SpinLzSzSymmetry::~BosonOnSphereWithSU2SpinLzSzSymmetry ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU2SpinLzSzSymmetry& BosonOnSphereWithSU2SpinLzSzSymmetry::operator = (const BosonOnSphereWithSU2SpinLzSzSymmetry& bosons)
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
  this->SzParitySign = bosons.SzParitySign;
  this->LzParitySign = bosons.LzParitySign;
  this->LzSymmetryMaxSwapPosition = bosons.LzSymmetryMaxSwapPosition;
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
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;
  this->RescalingFactors = bosons.RescalingFactors;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  this->MaxOrbitSize = bosons.MaxOrbitSize;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereWithSU2SpinLzSzSymmetry::Clone()
{
  return new BosonOnSphereWithSU2SpinLzSzSymmetry(*this);
}

// generate the Hilbert space with the discrete symmetry constraint
//

void BosonOnSphereWithSU2SpinLzSzSymmetry::GenerateStatesWithDiscreteSymmetry()

{
  long TmpLargeHilbertSpaceDimension = 0l;
  unsigned long TmpCanonicalStateDescriptionUp;
  unsigned long TmpCanonicalStateDescriptionDown;
  int NbrSzSymmetry;
  int NbrLzSymmetry;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->FindCanonicalForm(this->StateDescriptionUp[i], this->StateDescriptionDown[i], 
			      TmpCanonicalStateDescriptionUp, TmpCanonicalStateDescriptionDown, NbrLzSymmetry, NbrSzSymmetry);
      if ((this->StateDescriptionUp[i] > TmpCanonicalStateDescriptionUp) ||
	  ((this->StateDescriptionUp[i] == TmpCanonicalStateDescriptionUp) && (this->StateDescriptionDown[i] >= TmpCanonicalStateDescriptionDown)))
	{
	  if (this->TestSymmetryConstraint(this->StateDescriptionUp[i], this->StateDescriptionDown[i]) == true)
	    {
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescriptionUp[i] = 0x0ul;
	      this->StateDescriptionDown[i] = 0x0ul;
	    }
	}
      else
	{
	  this->StateDescriptionUp[i] = 0x0ul;
	  this->StateDescriptionDown[i] = 0x0ul;
	}
    }
  if (TmpLargeHilbertSpaceDimension > 0l)
    {
      unsigned long* TmpStateDescriptionUp = new unsigned long [TmpLargeHilbertSpaceDimension];
      unsigned long* TmpStateDescriptionDown = new unsigned long [TmpLargeHilbertSpaceDimension];
      this->NbrStateInOrbit = new int[TmpLargeHilbertSpaceDimension];
      TmpLargeHilbertSpaceDimension = 0;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if ((this->StateDescriptionUp[i] != 0x0ul) || (this->StateDescriptionDown[i] != 0x0ul))
	    {
	      TmpStateDescriptionUp[TmpLargeHilbertSpaceDimension] = this->StateDescriptionUp[i];
	      TmpStateDescriptionDown[TmpLargeHilbertSpaceDimension] = this->StateDescriptionDown[i];
	      this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] =  this->FindOrbitSize(this->StateDescriptionUp[i], this->StateDescriptionDown[i]);
	      ++TmpLargeHilbertSpaceDimension;
	    }
	}
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      this->StateDescriptionUp = TmpStateDescriptionUp;
      this->StateDescriptionDown = TmpStateDescriptionDown;
      this->RescalingFactors = new double*[this->MaxOrbitSize + 1];
      for (int i = 1; i <= this->MaxOrbitSize; ++i)
	{
	  this->RescalingFactors[i] = new double[this->MaxOrbitSize + 1];
	  for (int j = 1; j <= this->MaxOrbitSize; ++j)
	    {
	      this->RescalingFactors[i][j] = sqrt(((double) i) / ((double) j));
	    }
	}
    }
  this->LargeHilbertSpaceDimension = TmpLargeHilbertSpaceDimension;
}

