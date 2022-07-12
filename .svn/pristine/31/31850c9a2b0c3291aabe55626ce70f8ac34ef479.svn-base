////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//     class for bosons on sphere with SU(2) spin and Sz<->-Sz symmetry       //
//                                                                            //
//                        last modification : 23/09/2016                      //
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
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
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
#include "GeneralTools/Endian.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementMatrixOperation.h"

#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ios;


// default constructor
// 

BosonOnSphereWithSU2SpinSzSymmetry::BosonOnSphereWithSU2SpinSzSymmetry ()
{
  this->SzParitySign = 1.0;
}

// basic constructor without any constraint on Sz
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity

BosonOnSphereWithSU2SpinSzSymmetry::BosonOnSphereWithSU2SpinSzSymmetry (int nbrBosons, int totalLz, int lzMax, bool minusSzParity)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  if (minusSzParity == false)
    this->SzParitySign = 1.0;
  else
    this->SzParitySign = -1.0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->MaxOrbitSize = 2;
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
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2SpinSzSymmetry!" << endl;
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
      UsedMemory += this->LargeHilbertSpaceDimension * (2l * sizeof(unsigned long) + sizeof(int));
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

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// totalSpin = twice the total spin value (not taken into account)
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU2SpinSzSymmetry::BosonOnSphereWithSU2SpinSzSymmetry (int nbrBosons, int totalLz, int lzMax, int totalSpin, bool minusSzParity, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  if (minusSzParity == false)
    this->SzParitySign = 1.0;
  else
    this->SzParitySign = -1.0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxOrbitSize = 2;
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
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2SpinSzSymmetry!" << endl;
      exit(1);
    }

  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->GenerateStatesWithDiscreteSymmetry();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

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

BosonOnSphereWithSU2SpinSzSymmetry::BosonOnSphereWithSU2SpinSzSymmetry (char* fileName, unsigned long memory)
{
  this->MaxOrbitSize = 2;
  this->ReadHilbertSpace(fileName);

  this->IncNbrBosons = this->NbrBosons + 1;
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

BosonOnSphereWithSU2SpinSzSymmetry::BosonOnSphereWithSU2SpinSzSymmetry(const BosonOnSphereWithSU2SpinSzSymmetry& bosons)
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

BosonOnSphereWithSU2SpinSzSymmetry::~BosonOnSphereWithSU2SpinSzSymmetry ()
{
  if ((this->LargeHilbertSpaceDimension != 0l) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->NbrStateInOrbit;
      delete[] this->RescalingFactors;	  
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU2SpinSzSymmetry& BosonOnSphereWithSU2SpinSzSymmetry::operator = (const BosonOnSphereWithSU2SpinSzSymmetry& bosons)
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

AbstractHilbertSpace* BosonOnSphereWithSU2SpinSzSymmetry::Clone()
{
  return new BosonOnSphereWithSU2SpinSzSymmetry(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereWithSU2SpinSzSymmetry::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (BosonOnSphereWithSU2SpinSzSymmetry*) targetSpace;
}

// generate the Hilbert space with the discrete symmetry constraint
//

void BosonOnSphereWithSU2SpinSzSymmetry::GenerateStatesWithDiscreteSymmetry()

{
  long TmpHilbertSpaceDimension = 0l;
  if (this->SzParitySign > 0.0)
    {
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (this->StateDescriptionUp[i] >= this->StateDescriptionDown[i])
	    {
	      ++TmpHilbertSpaceDimension;
	    }
	}
    }
  else
    {
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (this->StateDescriptionUp[i] > this->StateDescriptionDown[i])
	    {
	      ++TmpHilbertSpaceDimension;
	    }
	}
    }
  if (TmpHilbertSpaceDimension > 0l)
    {
      unsigned long* TmpStateDescriptionUp = new unsigned long [TmpHilbertSpaceDimension];
      unsigned long* TmpStateDescriptionDown = new unsigned long [TmpHilbertSpaceDimension];
      this->NbrStateInOrbit = new int[TmpHilbertSpaceDimension];
      TmpHilbertSpaceDimension = 0;
      if (this->SzParitySign > 0.0)
	{
	  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      if (this->StateDescriptionUp[i] >= this->StateDescriptionDown[i])
		{
		  TmpStateDescriptionUp[TmpHilbertSpaceDimension] = this->StateDescriptionUp[i];
		  TmpStateDescriptionDown[TmpHilbertSpaceDimension] = this->StateDescriptionDown[i];
		  if (this->StateDescriptionUp[i] == this->StateDescriptionDown[i])
		    {
		      this->NbrStateInOrbit[TmpHilbertSpaceDimension] = 1;
		    }
		  else
		    {
		      this->NbrStateInOrbit[TmpHilbertSpaceDimension] = 2;
		    }
		  ++TmpHilbertSpaceDimension;
		}
	    }
	}
      else
	{
	  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      if (this->StateDescriptionUp[i] > this->StateDescriptionDown[i])
		{
		  this->NbrStateInOrbit[TmpHilbertSpaceDimension] = 2;
		  TmpStateDescriptionUp[TmpHilbertSpaceDimension] = this->StateDescriptionUp[i];
		  TmpStateDescriptionDown[TmpHilbertSpaceDimension] = this->StateDescriptionDown[i];
		  ++TmpHilbertSpaceDimension;
		}
	    }	 
	}	 
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      this->StateDescriptionUp = TmpStateDescriptionUp;
      this->StateDescriptionDown = TmpStateDescriptionDown;
      this->ComputeRescalingFactors();
    }
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
}

// compute the rescaling factors
//

void BosonOnSphereWithSU2SpinSzSymmetry::ComputeRescalingFactors()
{
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


// find state index
//
// stateDescriptionUp = unsigned integer describing the fermionic state for type up particles
// stateDescriptionDown = unsigned integer describing the fermionic state for type down particles
// return value = corresponding index

int BosonOnSphereWithSU2SpinSzSymmetry::FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  int PosMin = 0;
  int PosMax = this->NbrUniqueStateDescriptionUp - 1;
  if ((stateDescriptionUp > this->UniqueStateDescriptionUp[PosMin]) || (stateDescriptionUp < this->UniqueStateDescriptionUp[PosMax]))
    {
      return this->HilbertSpaceDimension;
    }
  int PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->UniqueStateDescriptionUp[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionUp))
    {
       if (CurrentState > stateDescriptionUp)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->UniqueStateDescriptionUp[PosMid];
    }
  if (CurrentState != stateDescriptionUp)
    {
      PosMid = PosMax;
      if (this->UniqueStateDescriptionUp[PosMid] != stateDescriptionUp)
	return this->HilbertSpaceDimension;
     }

  PosMin = this->FirstIndexUniqueStateDescriptionUp[PosMid];
  PosMax = PosMin + this->UniqueStateDescriptionSubArraySizeUp[PosMid] - 1;
  PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescriptionDown[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionDown))
    {
       if (CurrentState > stateDescriptionDown)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->StateDescriptionDown[PosMid];
    }
  if (CurrentState == stateDescriptionDown)
    return PosMid;
  if ( this->StateDescriptionDown[PosMax] != stateDescriptionDown)
    return this->HilbertSpaceDimension;
  else
    return PosMax;
}


// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::AduAu (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  if (this->TemporaryStateUp[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  coefficient = (double) this->TemporaryStateUp[n];
  --this->TemporaryStateUp[n];
  ++this->TemporaryStateUp[m];
  coefficient *= (double) this->TemporaryStateUp[m];
  coefficient = sqrt(coefficient);  
  unsigned long TmpUp = this->BosonToFermion(this->TemporaryStateUp);
  return this->SymmetrizeAdAdResult(TmpUp, this->StateDescriptionDown[index], coefficient);  
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::AduAd (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  if (this->TemporaryStateDown[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  coefficient = (double) this->TemporaryStateDown[n];
  --this->TemporaryStateDown[n];
  ++this->TemporaryStateUp[m];
  coefficient *= (double) this->TemporaryStateUp[m];
  coefficient = sqrt(coefficient); 
  unsigned long TmpUp = this->BosonToFermion(this->TemporaryStateUp);
  unsigned long TmpDown = this->BosonToFermion(this->TemporaryStateDown);
  return this->SymmetrizeAdAdResult(TmpUp, TmpDown, coefficient);  
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::AddAu (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  if (this->TemporaryStateUp[n] == 0)
    { 
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;      
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  coefficient = (double) this->TemporaryStateUp[n];
  --this->TemporaryStateUp[n];
  ++this->TemporaryStateDown[m];
  coefficient *= (double) this->TemporaryStateDown[m];
  coefficient = sqrt(coefficient);  
  unsigned long TmpUp = this->BosonToFermion(this->TemporaryStateUp);
  unsigned long TmpDown = this->BosonToFermion(this->TemporaryStateDown);
  return this->SymmetrizeAdAdResult(TmpUp, TmpDown, coefficient);  
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::AddAd (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  if (this->TemporaryStateDown[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  coefficient = (double) this->TemporaryStateDown[n];
  --this->TemporaryStateDown[n];
  ++this->TemporaryStateDown[m];
  coefficient *= (double) this->TemporaryStateDown[m];
  coefficient = sqrt(coefficient);  
  unsigned long TmpDown = this->BosonToFermion(this->TemporaryStateDown);
  return this->SymmetrizeAdAdResult(this->StateDescriptionUp[index], TmpDown, coefficient);  
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2SpinSzSymmetry::AuAu (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateUp[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUp[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = this->ProdATemporaryStateUp[n2];
  --this->ProdATemporaryStateUp[n2];
  Coefficient *= this->ProdATemporaryStateUp[n1];
  --this->ProdATemporaryStateUp[n1];
  return sqrt(Coefficient);
}

// apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2SpinSzSymmetry::AuAd (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = this->ProdATemporaryStateDown[n2];
  --this->ProdATemporaryStateDown[n2];
  Coefficient *= this->ProdATemporaryStateUp[n1];
  --this->ProdATemporaryStateUp[n1];
  return sqrt(Coefficient);
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2SpinSzSymmetry::AdAd (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateDown[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDown[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = this->ProdATemporaryStateDown[n2];
  --this->ProdATemporaryStateDown[n2];
  Coefficient *= this->ProdATemporaryStateDown[n1];
  --this->ProdATemporaryStateDown[n1];
  return sqrt(Coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereWithSU2SpinSzSymmetry::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  double Coefficient = 1.0;
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	{
	  unsigned long& Tmp = this->ProdATemporaryStateDown[n[i]];
	  if (Tmp == 0x0ul)
	    {
	      return 0.0;
	    }
	  Coefficient *= Tmp;
	  --Tmp;
	}
      else
	{
	  unsigned long& Tmp = this->ProdATemporaryStateUp[n[i]];
	  if (Tmp == 0x0ul)
	    {
	      return 0.0;
	    }
	  Coefficient *= Tmp;
	  --Tmp;
	}
    }
  return sqrt(Coefficient);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryStateUp[i] = this->ProdATemporaryStateUp[i];
      this->TemporaryStateDown[i] = this->ProdATemporaryStateDown[i];
    }
  coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	{
	  unsigned long& Tmp = this->TemporaryStateDown[m[i]];
	  ++Tmp;
	  coefficient *= Tmp;
	}
      else
	{
	  unsigned long& Tmp = this->TemporaryStateUp[m[i]];
	  ++Tmp;
	  coefficient *= Tmp;
	}
    }
  coefficient = sqrt(coefficient);
  return this->SymmetrizeAdAdResult(this->TemporaryStateUp, this->TemporaryStateDown, coefficient);  
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool BosonOnSphereWithSU2SpinSzSymmetry::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrBosons);
  WriteLittleEndian(File, this->LzMax);
  WriteLittleEndian(File, this->TotalLz);
  WriteLittleEndian(File, this->TotalSpin);
  WriteLittleEndian(File, this->SzParitySign);
  WriteBlockLittleEndian(File, this->StateDescriptionUp, this->LargeHilbertSpaceDimension);
  WriteBlockLittleEndian(File, this->StateDescriptionDown, this->LargeHilbertSpaceDimension);
  WriteBlockLittleEndian(File, this->NbrStateInOrbit, this->LargeHilbertSpaceDimension);
  File.close();
  return true;
}

// read Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description is stored
// return value = true if no error occured

bool BosonOnSphereWithSU2SpinSzSymmetry::ReadHilbertSpace (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->LargeHilbertSpaceDimension = 0l;
      this->HilbertSpaceDimension;
      return false;
    }
  ReadLittleEndian(File, this->LargeHilbertSpaceDimension);
  ReadLittleEndian(File, this->HilbertSpaceDimension);
  ReadLittleEndian(File, this->NbrBosons);
  ReadLittleEndian(File, this->LzMax);
  ReadLittleEndian(File, this->TotalLz);
  ReadLittleEndian(File, this->TotalSpin);
  ReadLittleEndian(File, this->SzParitySign);
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int[this->LargeHilbertSpaceDimension];
  ReadBlockLittleEndian(File, this->StateDescriptionUp, this->LargeHilbertSpaceDimension);
  ReadBlockLittleEndian(File, this->StateDescriptionDown, this->LargeHilbertSpaceDimension);
  ReadBlockLittleEndian(File, this->NbrStateInOrbit, this->LargeHilbertSpaceDimension);
  this->ComputeRescalingFactors();
  File.close();
  return true;
}
  
// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealMatrix BosonOnSphereWithSU2SpinSzSymmetry::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrParticleSector, int lzSector, int szSector, RealVector& groundState)
{
  int nbrOrbitalA = subsytemSize;
  int nbrOrbitalB = this->LzMax + 1 - nbrOrbitalA;  
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;


  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrBosons, this->LzMax);  
  int LzSectorDisk = ConvertLzFromSphereToDisk(lzSector, nbrParticleSector, nbrOrbitalA - 1);
  int ComplementaryLzSectorDisk = (TotalLzDisk - LzSectorDisk) - (ComplementaryNbrParticles * nbrOrbitalA);
  int ComplementaryLzSector = ConvertLzFromDiskToSphere(ComplementaryLzSectorDisk, ComplementaryNbrParticles, nbrOrbitalB - 1);

  if ((abs(ComplementarySzSector) > ComplementaryNbrParticles) || (abs(ComplementaryLzSector) > ((nbrOrbitalB - 1) * ComplementaryNbrParticles)))
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;  
    }
  BosonOnSphereWithSU2Spin SubsytemSpace(nbrParticleSector, lzSector, nbrOrbitalA - 1, szSector);
  BosonOnSphereWithSU2Spin ComplementarySubsytemSpace(ComplementaryNbrParticles, ComplementaryLzSector, nbrOrbitalB - 1, ComplementarySzSector);  
  RealMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySubsytemSpace.GetHilbertSpaceDimension(), true);
  long TmpNbrNonZeroElements = 0l;
  unsigned long** TmpSubsytemSpaceOccupationNumbersUp = new unsigned long* [SubsytemSpace.HilbertSpaceDimension];
  unsigned long** TmpSubsytemSpaceOccupationNumbersDown = new unsigned long* [SubsytemSpace.HilbertSpaceDimension];
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  for (int i = 0; i < SubsytemSpace.HilbertSpaceDimension; ++i)
    {
      TmpSubsytemSpaceOccupationNumbersUp[i] = new unsigned long [SubsytemSpace.NbrLzValue];
      TmpSubsytemSpaceOccupationNumbersDown[i] = new unsigned long [SubsytemSpace.NbrLzValue];
      SubsytemSpace.FermionToBoson(SubsytemSpace.StateDescriptionUp[i], SubsytemSpace.StateDescriptionDown[i], 
				   TmpSubsytemSpaceOccupationNumbersUp[i], TmpSubsytemSpaceOccupationNumbersDown[i]);
    }

  for (int MinIndex = 0; MinIndex < ComplementarySubsytemSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      ComplementarySubsytemSpace.FermionToBoson(ComplementarySubsytemSpace.StateDescriptionUp[MinIndex], ComplementarySubsytemSpace.StateDescriptionDown[MinIndex],  
						ComplementarySubsytemSpace.TemporaryStateUp, ComplementarySubsytemSpace.TemporaryStateDown);
      for (int j = 0; j < SubsytemSpace.HilbertSpaceDimension; ++j)
	{
	  for (int i = 0; i <= SubsytemSpace.LzMax; ++i)
	    {
	      this->TemporaryStateUp[i] = TmpSubsytemSpaceOccupationNumbersUp[j][i];
	      this->TemporaryStateDown[i] = TmpSubsytemSpaceOccupationNumbersDown[j][i];
	    }
	  for (int i = 0; i <= ComplementarySubsytemSpace.LzMax; ++i)
	    {
	      this->TemporaryStateUp[i + nbrOrbitalA] = ComplementarySubsytemSpace.TemporaryStateUp[i];
	      this->TemporaryStateDown[i + nbrOrbitalA] = ComplementarySubsytemSpace.TemporaryStateDown[i];
	    }
	  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpStateUp, TmpStateDown);
	  double TmpCoefficient = 1.0;
	  this->ProdATemporaryNbrStateInOrbit = 1;
	  int TmpPos =  this->SymmetrizeAdAdResult(TmpStateUp, TmpStateDown, TmpCoefficient);  
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      double Tmp = TmpCoefficient * groundState[TmpPos];
	      TmpEntanglementMatrix.SetMatrixElement(j, MinIndex, Tmp);
	      ++TmpNbrNonZeroElements;
	    }
	}
    }

  for (int i = 0; i < SubsytemSpace.HilbertSpaceDimension; ++i)
    {
      delete[] TmpSubsytemSpaceOccupationNumbersUp[i];
      delete[] TmpSubsytemSpaceOccupationNumbersDown[i];
    }
  delete[] TmpSubsytemSpaceOccupationNumbersUp;
  delete[] TmpSubsytemSpaceOccupationNumbersDown;
  if (TmpNbrNonZeroElements > 0l)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix BosonOnSphereWithSU2SpinSzSymmetry::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, RealVector& groundState, bool removeBinomialCoefficient, AbstractArchitecture* architecture)
{
  int nbrOrbitalA = this->LzMax + 1;
  int nbrOrbitalB = this->LzMax + 1;  
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int ComplementaryLzSector = this->TotalLz - lzSector;
  if ((abs(ComplementarySzSector) > ComplementaryNbrParticles) || (abs(ComplementaryLzSector) > (this->LzMax * ComplementaryNbrParticles)))
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;  
    }

  BosonOnSphereWithSU2Spin SubsytemSpace(nbrParticleSector, lzSector, this->LzMax, szSector);
  BosonOnSphereWithSU2Spin ComplementarySubsytemSpace(ComplementaryNbrParticles, ComplementaryLzSector, this->LzMax, ComplementarySzSector);  
  if ((SubsytemSpace.GetHilbertSpaceDimension() == 0) || (ComplementarySubsytemSpace.GetHilbertSpaceDimension() == 0))
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
  RealMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySubsytemSpace.GetHilbertSpaceDimension(), true);

  long TmpNbrNonZeroElements = 0l;

  if (architecture != 0)
    {
      FQHESphereParticleEntanglementMatrixOperation TmpOperation (this, &SubsytemSpace, &ComplementarySubsytemSpace, groundState, TmpEntanglementMatrix, removeBinomialCoefficient);
      TmpOperation.ApplyOperation(architecture);
      TmpNbrNonZeroElements = TmpOperation.GetNbrNonZeroMatrixElements();
    }
  else
    {
      TmpNbrNonZeroElements = this->EvaluatePartialEntanglementMatrixParticlePartitionCore (0, ComplementarySubsytemSpace.GetHilbertSpaceDimension(), &ComplementarySubsytemSpace, 
											    &SubsytemSpace, groundState, &TmpEntanglementMatrix, removeBinomialCoefficient);
    }

  if (TmpNbrNonZeroElements > 0l)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}
   
// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given Lz,Sz=0, Sz parity sectors.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSector = Sz sector in which the density matrix has to be evaluated. It should be equal to zero
// szParitySector = parity sector for the discrete symmetry Sz<->-Sz. Can be either -1 or +1
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix BosonOnSphereWithSU2SpinSzSymmetry::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, int szParity, 
												   RealVector& groundState, 
												   bool removeBinomialCoefficient, AbstractArchitecture* architecture)
{
  int nbrOrbitalA = this->LzMax + 1;
  int nbrOrbitalB = this->LzMax + 1;  
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int ComplementaryLzSector = this->TotalLz - lzSector;
  if ((abs(ComplementarySzSector) > ComplementaryNbrParticles) || (abs(ComplementaryLzSector) > (this->LzMax * ComplementaryNbrParticles)))
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;  
    }
  int ComplementarySzParity = szParity;
  if (this->SzParitySign < 0.0)
    {
      ComplementarySzParity = -szParity;
    }
  BosonOnSphereWithSU2SpinSzSymmetry SubsytemSpace(nbrParticleSector, lzSector, this->LzMax, szSector, (szParity < 0));
  BosonOnSphereWithSU2SpinSzSymmetry ComplementarySubsytemSpace(ComplementaryNbrParticles, ComplementaryLzSector, this->LzMax, ComplementarySzSector, (ComplementarySzParity < 0));  
  if ((SubsytemSpace.GetHilbertSpaceDimension() == 0) || (ComplementarySubsytemSpace.GetHilbertSpaceDimension() == 0))
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
  RealMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySubsytemSpace.GetHilbertSpaceDimension(), true);

  long TmpNbrNonZeroElements = 0l;

  if (architecture != 0)
    {
      FQHESphereParticleEntanglementMatrixOperation TmpOperation (this, &SubsytemSpace, &ComplementarySubsytemSpace, groundState, TmpEntanglementMatrix, removeBinomialCoefficient);
      TmpOperation.ApplyOperation(architecture);
      TmpNbrNonZeroElements = TmpOperation.GetNbrNonZeroMatrixElements();
    }
  else
    {
      TmpNbrNonZeroElements = this->EvaluatePartialEntanglementMatrixParticlePartitionCore (0, ComplementarySubsytemSpace.GetHilbertSpaceDimension(), &ComplementarySubsytemSpace, 
											    &SubsytemSpace, groundState, &TmpEntanglementMatrix, removeBinomialCoefficient);
    }

  if (TmpNbrNonZeroElements > 0l)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}
   
// core part of the entanglement matrix evaluation for the particle partition
// 
// minIndex = first index to consider in the complementary Hilbert space
// nbrIndex = number of indices to consider in the complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// entanglementMatrix = pointer to entanglement matrix
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = number of components that have been added to the entanglement matrix

long BosonOnSphereWithSU2SpinSzSymmetry::EvaluatePartialEntanglementMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace, 
												 ParticleOnSphere* destinationHilbertSpace, RealVector& groundState, 
												 RealMatrix* entanglementMatrix, bool removeBinomialCoefficient)
{
  long TmpNbrNonZeroElements = 0l;
  BosonOnSphereWithSU2Spin* TmpSubsystemSpace = (BosonOnSphereWithSU2Spin*)destinationHilbertSpace;
  bool TmpSzParityFlag = true;
  BosonOnSphereWithSU2Spin* SubsystemSpaceFull = 0;

  if ((this->TotalSpin != 0) || (TmpSubsystemSpace->TotalSpin != 0))
    {
      TmpSzParityFlag = false;
    }
  else
    {
      SubsystemSpaceFull = new BosonOnSphereWithSU2Spin (TmpSubsystemSpace->NbrBosons, TmpSubsystemSpace->TotalLz, TmpSubsystemSpace->LzMax, TmpSubsystemSpace->TotalSpin);
      if (SubsystemSpaceFull->GetHilbertSpaceDimension() == TmpSubsystemSpace->GetHilbertSpaceDimension())
	{
	  delete SubsystemSpaceFull;
	  TmpSzParityFlag = false;
	}
    }

  if (TmpSzParityFlag == false)
    {
      BosonOnSphereWithSU2Spin* SubsystemSpace = (BosonOnSphereWithSU2Spin*)destinationHilbertSpace;
      BosonOnSphereWithSU2Spin* ComplementarySubsystemSpace = (BosonOnSphereWithSU2Spin*) complementaryHilbertSpace;
      
      int nbrOrbitalA = this->LzMax + 1;
      int nbrOrbitalB = this->LzMax + 1;  
      
      double* LogFactorials = new double[this->NbrBosons + 1];
      LogFactorials[0] = 0.0;
      LogFactorials[1] = 0.0;
      for (int i = 2 ; i <= this->NbrBosons; ++i)
	LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
      double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementarySubsystemSpace->NbrBosons] - LogFactorials[SubsystemSpace->NbrBosons];
      if (removeBinomialCoefficient == true)
	TmpLogBinomial = 0.0;
      
      unsigned long** TmpSubsystemSpaceOccupationNumbersUp = new unsigned long* [SubsystemSpace->GetHilbertSpaceDimension()];
      unsigned long** TmpSubsystemSpaceOccupationNumbersDown = new unsigned long* [SubsystemSpace->GetHilbertSpaceDimension()];
      double* TmpSubsystemLogFactorials = new double [SubsystemSpace->GetHilbertSpaceDimension()];
      unsigned long* TmpOccupationUp = new unsigned long [this->LzMax + 1];
      unsigned long* TmpOccupationDown = new unsigned long [this->LzMax + 1];
      int TmpShift = this->LzMax + 1 - nbrOrbitalA;
      int LastIndex = nbrIndex + minIndex;

      for (int i = 0; i < SubsystemSpace->GetHilbertSpaceDimension(); ++i)
	{
	  TmpSubsystemSpaceOccupationNumbersUp[i] = new unsigned long [this->NbrLzValue];
	  TmpSubsystemSpaceOccupationNumbersDown[i] = new unsigned long [this->NbrLzValue];
	  SubsystemSpace->FermionToBoson(SubsystemSpace->StateDescriptionUp[i], SubsystemSpace->StateDescriptionDown[i], 
				   TmpSubsystemSpaceOccupationNumbersUp[i], TmpSubsystemSpaceOccupationNumbersDown[i]);
	  unsigned long* TmpOccupationNumberUp = TmpSubsystemSpaceOccupationNumbersUp[i];
	  unsigned long* TmpOccupationNumberDown = TmpSubsystemSpaceOccupationNumbersDown[i];
	  double TmpFactor = 0.0;
	  for (int k = 0; k <= SubsystemSpace->LzMax; ++k)
	    {
	      TmpFactor += LogFactorials[TmpOccupationNumberUp[k]];
	      TmpFactor += LogFactorials[TmpOccupationNumberDown[k]];
	    }
	  TmpSubsystemLogFactorials[i] = TmpFactor;      
	}
      for (int MinIndex = minIndex; MinIndex < LastIndex; ++MinIndex)    
	{
	  ComplementarySubsystemSpace->FermionToBoson(ComplementarySubsystemSpace->StateDescriptionUp[MinIndex], ComplementarySubsystemSpace->StateDescriptionDown[MinIndex],  
						      ComplementarySubsystemSpace->TemporaryStateUp, ComplementarySubsystemSpace->TemporaryStateDown);
	  double ComplementarySubsystemSpaceFactorial = 0.0;
	  for (int k = 0; k <= ComplementarySubsystemSpace->LzMax; ++k)
	    {
	      ComplementarySubsystemSpaceFactorial += LogFactorials[ComplementarySubsystemSpace->TemporaryStateUp[k]];
	      ComplementarySubsystemSpaceFactorial += LogFactorials[ComplementarySubsystemSpace->TemporaryStateDown[k]];
	    }
	  for (int j = 0; j < SubsystemSpace->GetHilbertSpaceDimension(); ++j)
	    {
	      
	      unsigned long TmpStateUp;
	      unsigned long TmpStateDown;
	      unsigned long* TmpOccupationArrayUp = TmpSubsystemSpaceOccupationNumbersUp[j];
	      unsigned long* TmpOccupationArrayDown = TmpSubsystemSpaceOccupationNumbersDown[j];
	      for (int k = 0; k < nbrOrbitalB; ++k)
		{
		  this->TemporaryStateUp[k] = ComplementarySubsystemSpace->TemporaryStateUp[k];
		  this->TemporaryStateDown[k] = ComplementarySubsystemSpace->TemporaryStateDown[k];
		}
	      for (int k = TmpShift; k < nbrOrbitalB; ++k)
		{
		  this->TemporaryStateUp[k] += TmpOccupationArrayUp[k - TmpShift];
		  this->TemporaryStateDown[k] += TmpOccupationArrayDown[k - TmpShift];
		}
	      for (int k = nbrOrbitalB; k <= this->LzMax; ++k)
		{
		  this->TemporaryStateUp[k] = TmpOccupationArrayUp[k - TmpShift];
		  this->TemporaryStateDown[k] = TmpOccupationArrayDown[k - TmpShift];
		}
	      
	      double TmpCoefficient = 1.0;
	      this->ProdATemporaryNbrStateInOrbit = 1;
	      this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpStateUp, TmpStateDown);
	      int TmpPos =  this->SymmetrizeAdAdResult(TmpStateUp, TmpStateDown, TmpCoefficient);  
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  double TmpFactorial = 0.0;	      
		  for (int k = 0; k <= this->LzMax; ++k)
		    {
		      TmpFactorial += LogFactorials[this->TemporaryStateUp[k]];
		      TmpFactorial += LogFactorials[this->TemporaryStateDown[k]];
		    }
		  TmpFactorial -= ComplementarySubsystemSpaceFactorial + TmpSubsystemLogFactorials[j] + TmpLogBinomial;
		  TmpFactorial *= 0.5; 	      
		  ++TmpNbrNonZeroElements;
		  double Tmp = TmpCoefficient * exp(TmpFactorial) * groundState[TmpPos];
		  entanglementMatrix->SetMatrixElement(j, MinIndex, Tmp);
		}
	    }
	}
      
      for (int i = 0; i < SubsystemSpace->GetHilbertSpaceDimension(); ++i)
	{
	  delete[] TmpSubsystemSpaceOccupationNumbersUp[i];
	  delete[] TmpSubsystemSpaceOccupationNumbersDown[i];
	}
      delete[] TmpSubsystemSpaceOccupationNumbersUp;
      delete[] TmpSubsystemSpaceOccupationNumbersDown;
      delete[] LogFactorials;
    }
  else
    {
      BosonOnSphereWithSU2SpinSzSymmetry* SubsystemSpace = (BosonOnSphereWithSU2SpinSzSymmetry*)destinationHilbertSpace;
      BosonOnSphereWithSU2SpinSzSymmetry* ComplementarySubsystemSpace = (BosonOnSphereWithSU2SpinSzSymmetry*) complementaryHilbertSpace;
      
      int nbrOrbitalA = this->LzMax + 1;
      int nbrOrbitalB = this->LzMax + 1;  
      
      double* LogFactorials = new double[this->NbrBosons + 1];
      LogFactorials[0] = 0.0;
      LogFactorials[1] = 0.0;
      for (int i = 2 ; i <= this->NbrBosons; ++i)
	LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
      double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementarySubsystemSpace->NbrBosons] - LogFactorials[SubsystemSpaceFull->NbrBosons];
      if (removeBinomialCoefficient == true)
	TmpLogBinomial = 0.0;
      
      unsigned long** TmpSubsystemSpaceOccupationNumbersUp = new unsigned long* [SubsystemSpaceFull->GetHilbertSpaceDimension()];
      unsigned long** TmpSubsystemSpaceOccupationNumbersDown = new unsigned long* [SubsystemSpaceFull->GetHilbertSpaceDimension()];
      double* TmpSubsystemLogFactorials = new double [SubsystemSpaceFull->GetHilbertSpaceDimension()];
      unsigned long* TmpOccupationUp = new unsigned long [this->LzMax + 1];
      unsigned long* TmpOccupationDown = new unsigned long [this->LzMax + 1];
      int TmpShift = this->LzMax + 1 - nbrOrbitalA;
      int LastIndex = nbrIndex + minIndex;
      for (int i = 0; i < SubsystemSpaceFull->GetHilbertSpaceDimension(); ++i)
	{
	  TmpSubsystemSpaceOccupationNumbersUp[i] = new unsigned long [this->NbrLzValue];
	  TmpSubsystemSpaceOccupationNumbersDown[i] = new unsigned long [this->NbrLzValue];
	  SubsystemSpaceFull->FermionToBoson(SubsystemSpaceFull->StateDescriptionUp[i], SubsystemSpaceFull->StateDescriptionDown[i], 
					     TmpSubsystemSpaceOccupationNumbersUp[i], TmpSubsystemSpaceOccupationNumbersDown[i]);
	  unsigned long* TmpOccupationNumberUp = TmpSubsystemSpaceOccupationNumbersUp[i];
	  unsigned long* TmpOccupationNumberDown = TmpSubsystemSpaceOccupationNumbersDown[i];
	  double TmpFactor = 0.0;
	  for (int k = 0; k <= SubsystemSpaceFull->LzMax; ++k)
	    {
	      TmpFactor += LogFactorials[TmpOccupationNumberUp[k]];
	      TmpFactor += LogFactorials[TmpOccupationNumberDown[k]];
	    }
	  TmpSubsystemLogFactorials[i] = TmpFactor;      
	}
      for (int MinIndex = minIndex; MinIndex < LastIndex; ++MinIndex)    
	{
	  double TmpRescalingFactor = sqrt((double)  ComplementarySubsystemSpace->NbrStateInOrbit[MinIndex]);
	  ComplementarySubsystemSpace->FermionToBoson(ComplementarySubsystemSpace->StateDescriptionUp[MinIndex], ComplementarySubsystemSpace->StateDescriptionDown[MinIndex],  
						      ComplementarySubsystemSpace->TemporaryStateUp, ComplementarySubsystemSpace->TemporaryStateDown);
	  double ComplementarySubsystemSpaceFactorial = 0.0;
	  for (int k = 0; k <= ComplementarySubsystemSpace->LzMax; ++k)
	    {
	      ComplementarySubsystemSpaceFactorial += LogFactorials[ComplementarySubsystemSpace->TemporaryStateUp[k]];
	      ComplementarySubsystemSpaceFactorial += LogFactorials[ComplementarySubsystemSpace->TemporaryStateDown[k]];
	    }
	  for (int j = 0; j < SubsystemSpaceFull->GetHilbertSpaceDimension(); ++j)
	    {
	      unsigned long TmpStateUp = SubsystemSpaceFull->StateDescriptionUp[j];
	      unsigned long TmpStateDown = SubsystemSpaceFull->StateDescriptionDown[j];
	      SubsystemSpace->ProdATemporaryNbrStateInOrbit = 1;      
	      double TmpDestinationCoefficient = TmpRescalingFactor;
	      int RealDestinationIndex = SubsystemSpace->SymmetrizeAdAdResult(TmpStateUp, TmpStateDown, TmpDestinationCoefficient);
	      if (RealDestinationIndex < SubsystemSpace->GetHilbertSpaceDimension())
		{		  
		  unsigned long* TmpOccupationArrayUp = TmpSubsystemSpaceOccupationNumbersUp[j];
		  unsigned long* TmpOccupationArrayDown = TmpSubsystemSpaceOccupationNumbersDown[j];
		  for (int k = 0; k < nbrOrbitalB; ++k)
		    {
		      this->TemporaryStateUp[k] = ComplementarySubsystemSpace->TemporaryStateUp[k];
		      this->TemporaryStateDown[k] = ComplementarySubsystemSpace->TemporaryStateDown[k];
		    }
		  for (int k = TmpShift; k < nbrOrbitalB; ++k)
		    {
		      this->TemporaryStateUp[k] += TmpOccupationArrayUp[k - TmpShift];
		      this->TemporaryStateDown[k] += TmpOccupationArrayDown[k - TmpShift];
		    }
		  for (int k = nbrOrbitalB; k <= this->LzMax; ++k)
		    {
		      this->TemporaryStateUp[k] = TmpOccupationArrayUp[k - TmpShift];
		      this->TemporaryStateDown[k] = TmpOccupationArrayDown[k - TmpShift];
		    }
		  
		  double TmpCoefficient = 1.0;
		  this->ProdATemporaryNbrStateInOrbit = 1;
		  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpStateUp, TmpStateDown);
		  int TmpPos =  this->SymmetrizeAdAdResult(TmpStateUp, TmpStateDown, TmpCoefficient);  
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      double TmpFactorial = 0.0;	      
		      for (int k = 0; k <= this->LzMax; ++k)
			{
			  TmpFactorial += LogFactorials[this->TemporaryStateUp[k]];
			  TmpFactorial += LogFactorials[this->TemporaryStateDown[k]];
			}
		      TmpFactorial -= ComplementarySubsystemSpaceFactorial + TmpSubsystemLogFactorials[j] + TmpLogBinomial;
		      TmpFactorial *= 0.5; 	      
		      ++TmpNbrNonZeroElements;
		      double Tmp = TmpCoefficient * TmpDestinationCoefficient *  exp(TmpFactorial) * groundState[TmpPos];
		      entanglementMatrix->AddToMatrixElement(RealDestinationIndex, MinIndex, Tmp);
		    }
		}
	    }
	}
      
      for (int i = 0; i < SubsystemSpaceFull->GetHilbertSpaceDimension(); ++i)
	{
	  delete[] TmpSubsystemSpaceOccupationNumbersUp[i];
	  delete[] TmpSubsystemSpaceOccupationNumbersDown[i];
	}
      delete[] TmpSubsystemSpaceOccupationNumbersUp;
      delete[] TmpSubsystemSpaceOccupationNumbersDown;
      delete[] LogFactorials;
      delete SubsystemSpaceFull;
    }
  return TmpNbrNonZeroElements;
}
   
// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrParticleSector = number of particles that belong to the subsystem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& BosonOnSphereWithSU2SpinSzSymmetry::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector,
																   int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
																   int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, RealMatrix& entanglementMatrix)
{
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrBosons, this->LzMax);
  int LzADisk = ConvertLzFromSphereToDisk(lzSector, nbrParticleSector, nbrOrbitalA - 1);
  if ((LzADisk < 0) || (LzADisk > ((nbrOrbitalA - 1) * nbrParticleSector)))
    {
      return entanglementMatrix;	  
    }
  int LzBDisk = (TotalLzDisk - LzADisk) - ComplementaryNbrParticles * (this->LzMax + 1 - nbrOrbitalB);
  if ((LzBDisk < 0) || (LzBDisk > ((nbrOrbitalB - 1) * ComplementaryNbrParticles)))
    {
      return entanglementMatrix;	  
    }
  int ComplementaryLzSector = ConvertLzFromDiskToSphere(LzBDisk, ComplementaryNbrParticles, nbrOrbitalB - 1);
  BosonOnSphereWithSU2Spin SubsystemSpace(nbrParticleSector, lzSector, nbrOrbitalA - 1, szSector);
  BosonOnSphereWithSU2Spin ComplementarySubsystemSpace(ComplementaryNbrParticles, ComplementaryLzSector, nbrOrbitalB - 1, ComplementarySzSector);  

  cout << "subsystem Hilbert space dimension = " << SubsystemSpace.GetHilbertSpaceDimension() << endl;
  unsigned long* TmpMonomialUp1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialDown1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialUp3 = new unsigned long [nbrParticleSector];
  unsigned long* TmpMonomialDown3 = new unsigned long [nbrParticleSector];
  
  for (int i = 0; i < SubsystemSpace.GetHilbertSpaceDimension(); ++i)
    {
      SubsystemSpace.ConvertToMonomial(SubsystemSpace.StateDescriptionUp[i], SubsystemSpace.StateDescriptionDown[i], TmpMonomialUp3, TmpMonomialDown3);
      double Tmp = 1.0;
      for (int j = 0; j < SubsystemSpace.NbrBosonsUp; j++)
	{
	  Tmp *= weightOrbitalAUp[TmpMonomialUp3[j]];
	}
      for (int j = 0; j < SubsystemSpace.NbrBosonsDown; j++)
	{
	  Tmp *= weightOrbitalADown[TmpMonomialDown3[j]];
	}
      for (int j = 0; j < ComplementarySubsystemSpace.GetHilbertSpaceDimension(); ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < ComplementarySubsystemSpace.GetHilbertSpaceDimension(); ++MinIndex)    
    {
      ComplementarySubsystemSpace.ConvertToMonomial(ComplementarySubsystemSpace.StateDescriptionUp[MinIndex], ComplementarySubsystemSpace.StateDescriptionDown[MinIndex], TmpMonomialUp1, TmpMonomialDown1);
      double FormFactor = 1.0;
      for (int i = 0; i < ComplementarySubsystemSpace.NbrBosonsUp; i++)
	FormFactor *= weightOrbitalBUp[TmpMonomialUp1[i]];
      for (int i = 0; i < ComplementarySubsystemSpace.NbrBosonsDown; i++)
	FormFactor *= weightOrbitalBDown[TmpMonomialDown1[i]];
      for (int j = 0; j < SubsystemSpace.GetHilbertSpaceDimension(); ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
  
  delete[] TmpMonomialUp1;
  delete[] TmpMonomialDown1;
  delete[] TmpMonomialUp3;
  delete[] TmpMonomialDown3;
  
  return entanglementMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrParticleSector = number of particles that belong to the subsystem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// szParitySector = parity sector for the discrete symmetry Sz<->-Sz. Can be either -1 or +1
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& BosonOnSphereWithSU2SpinSzSymmetry::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector, 
																   int szParity,
																   int nbrOrbitalA, double* weightOrbitalAUp, 
																   double* weightOrbitalADown, 
																   int nbrOrbitalB, double* weightOrbitalBUp, 
																   double* weightOrbitalBDown, 
																   RealMatrix& entanglementMatrix)
{
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int ComplementarySzParity = szParity;
  if (this->SzParitySign < 0.0)
    {
      ComplementarySzParity = -szParity;
    }
  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrBosons, this->LzMax);
  int LzADisk = ConvertLzFromSphereToDisk(lzSector, nbrParticleSector, nbrOrbitalA - 1);
  if ((LzADisk < 0) || (LzADisk > ((nbrOrbitalA - 1) * nbrParticleSector)))
    {
      return entanglementMatrix;	  
    }
  int LzBDisk = (TotalLzDisk - LzADisk) - ComplementaryNbrParticles * (this->LzMax + 1 - nbrOrbitalB);
  if ((LzBDisk < 0) || (LzBDisk > ((nbrOrbitalB - 1) * ComplementaryNbrParticles)))
    {
      return entanglementMatrix;	  
    }
  int ComplementaryLzSector = ConvertLzFromDiskToSphere(LzBDisk, ComplementaryNbrParticles, nbrOrbitalB - 1);
  BosonOnSphereWithSU2SpinSzSymmetry SubsystemSpace(nbrParticleSector, lzSector, nbrOrbitalA - 1, szSector, (szParity < 0));
  BosonOnSphereWithSU2SpinSzSymmetry ComplementarySubsystemSpace(ComplementaryNbrParticles, ComplementaryLzSector, nbrOrbitalB - 1, ComplementarySzSector, (ComplementarySzParity < 0));  

  cout << "subsystem Hilbert space dimension = " << SubsystemSpace.GetHilbertSpaceDimension() << endl;
  unsigned long* TmpMonomialUp1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialDown1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialUp3 = new unsigned long [nbrParticleSector];
  unsigned long* TmpMonomialDown3 = new unsigned long [nbrParticleSector];
  
  for (int i = 0; i < SubsystemSpace.GetHilbertSpaceDimension(); ++i)
    {
      SubsystemSpace.ConvertToMonomial(SubsystemSpace.StateDescriptionUp[i], SubsystemSpace.StateDescriptionDown[i], TmpMonomialUp3, TmpMonomialDown3);
      double Tmp = 1.0;
      for (int j = 0; j < SubsystemSpace.NbrBosonsUp; j++)
	{
	  Tmp *= weightOrbitalAUp[TmpMonomialUp3[j]];
	}
      for (int j = 0; j < SubsystemSpace.NbrBosonsDown; j++)
	{
	  Tmp *= weightOrbitalADown[TmpMonomialDown3[j]];
	}
      for (int j = 0; j < ComplementarySubsystemSpace.GetHilbertSpaceDimension(); ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < ComplementarySubsystemSpace.GetHilbertSpaceDimension(); ++MinIndex)    
    {
      ComplementarySubsystemSpace.ConvertToMonomial(ComplementarySubsystemSpace.StateDescriptionUp[MinIndex], ComplementarySubsystemSpace.StateDescriptionDown[MinIndex], TmpMonomialUp1, TmpMonomialDown1);
      double FormFactor = 1.0;
      for (int i = 0; i < ComplementarySubsystemSpace.NbrBosonsUp; i++)
	FormFactor *= weightOrbitalBUp[TmpMonomialUp1[i]];
      for (int i = 0; i < ComplementarySubsystemSpace.NbrBosonsDown; i++)
	FormFactor *= weightOrbitalBDown[TmpMonomialDown1[i]];
      for (int j = 0; j < SubsystemSpace.GetHilbertSpaceDimension(); ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
  
  delete[] TmpMonomialUp1;
  delete[] TmpMonomialDown1;
  delete[] TmpMonomialUp3;
  delete[] TmpMonomialDown3;
  
  return entanglementMatrix;
}

// Compute the product of a spinful Slater determinant with a Van der Monde determinant
//
// slaterUp = monomial representation of the Slater spin up part
// slaterDown = monomial representation of the Slater spin up part
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

void BosonOnSphereWithSU2SpinSzSymmetry::VanDerMondeTimesSlater (unsigned long* slaterUp, unsigned long* slaterDown, RealVector& finalState, 
								 double** threeOrbitalOverlaps)
{
  unsigned long TmpNbrStates = 0;
  long CoefUp = 1;
  long CoefDown = 1;
  unsigned long StateUp [this->NbrBosonsUp];
  unsigned long StateDown [this->NbrBosonsDown];
  unsigned long TmpFinalStateUp;
  unsigned long TmpFinalStateDown;
  unsigned long TmpState = 0ul;
  double Sign = 1.0;
  unsigned long Mask = 0ul;
  unsigned long VanDerMonde [this->NbrBosons];
  unsigned long TmpHeapArray [this->NbrBosons];
  unsigned long TmpDim = ((unsigned long) this->NbrBosons);
  for (unsigned long i = 0ul; i < TmpDim; ++i)
    {
      VanDerMonde[i] = i;
      TmpHeapArray[i] = 0ul;
    }
  finalState.ClearVector();


  double TmpFactorUp = 0.0;
  bool ChangedUp = false;
  for (int i = 0; i < this->NbrBosonsUp; ++i)
    {
      StateUp[i] = slaterUp[i] + VanDerMonde[i];
      TmpFactorUp += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
    }
  double TmpFactorDown = 0.0;
  bool ChangedDown = false;
  for (int i = 0; i < this->NbrBosonsDown; ++i)
    {
      StateDown[i] = slaterDown[i] + VanDerMonde[this->NbrBosonsUp + i];
      TmpFactorDown += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
    }
  
  double Coefficient = 1.0;
  this->ProdATemporaryNbrStateInOrbit = 1;
  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
//  int TmpPos =  this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
  int TmpPos =  this->SymmetrizeAdAdResult(TmpFinalStateUp, TmpFinalStateDown, Coefficient);
  if (TmpPos != this->HilbertSpaceDimension)
    {
      finalState[TmpPos] += Sign * Coefficient * exp(TmpFactorUp + TmpFactorDown);
    }
  ++TmpNbrStates;
  unsigned long Tmp = 0;
  while (Tmp < TmpDim)
    {
      if (TmpHeapArray[Tmp] < Tmp)
	{
	  if ((Tmp & 0x1ul) == 0x0ul)
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[0];
	      VanDerMonde[0] = Tmp2;
	      if (Tmp >= this->NbrBosonsUp)
		ChangedDown = true;
	      ChangedUp = true;
	      
	    }
	  else
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[TmpHeapArray[Tmp]];
	      VanDerMonde[TmpHeapArray[Tmp]] = Tmp2;
	      if (TmpHeapArray[Tmp] >= this->NbrBosonsUp)
		{
		  ChangedDown = true;
		}
	      else
		{
		  ChangedUp = true;		  
		  if (Tmp >= this->NbrBosonsUp)
		    {
		      ChangedDown = true;
		    }
		}
	    }
	  if (ChangedUp == true)
	    {
	      TmpFactorUp = 0.0;
	      for (int i = 0; i < this->NbrBosonsUp; ++i)
		{
		  StateUp[i] = slaterUp[i] + VanDerMonde[i];
		  TmpFactorUp += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
		}
	      this->ConvertFromMonomial(StateUp, this->TemporaryStateUp, this->NbrBosonsUp);
	      TmpFinalStateUp = this->BosonToFermion(this->TemporaryStateUp);
	      ChangedUp = false;
	    }
	  if (ChangedDown == true)
	    {
	      TmpFactorDown = 0.0;
	      for (int i = 0; i < this->NbrBosonsDown; ++i)
		{
		  StateDown[i] = slaterDown[i] + VanDerMonde[this->NbrBosonsUp + i];
		  TmpFactorDown += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		}
	      this->ConvertFromMonomial(StateDown, this->TemporaryStateDown, this->NbrBosonsDown);
	      TmpFinalStateDown = this->BosonToFermion(this->TemporaryStateDown);
	      ChangedDown = false;	      
	    }
	  Sign *= -1.0;
	  Coefficient = 1.0;
	  int TmpPos = this->SymmetrizeAdAdResult(TmpFinalStateUp, TmpFinalStateDown, Coefficient);
//	  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      finalState[TmpPos] += Sign* Coefficient * exp(TmpFactorUp + TmpFactorDown);
	    }
	  ++TmpHeapArray[Tmp];
	  Tmp = 0ul;
	}
      else
	{
	  TmpHeapArray[Tmp]= 0ul;
	  ++Tmp;
	}
    }
}

// Compute the product of a spinful Slater determinant with a Van der Monde determinant, assuming a reverse flux attachment
//
// slaterUp = monomial representation of the Slater spin up part
// slaterDown = monomial representation of the Slater spin up part
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

void BosonOnSphereWithSU2SpinSzSymmetry::ReverseVanDerMondeTimesSlater (unsigned long* slaterUp, unsigned long* slaterDown, RealVector& finalState, 
									double** threeOrbitalOverlaps)
{
  unsigned long TmpNbrStates = 0;
  long CoefUp = 1;
  long CoefDown = 1;
  unsigned long StateUp [this->NbrBosonsUp];
  unsigned long StateDown [this->NbrBosonsDown];
  unsigned long TmpFinalStateUp;
  unsigned long TmpFinalStateDown;
  unsigned long TmpState = 0ul;
  double Sign = 1.0;
  unsigned long Mask = 0ul;
  unsigned long VanDerMonde [this->NbrBosons];
  unsigned long TmpHeapArray [this->NbrBosons];
  unsigned long TmpDim = ((unsigned long) this->NbrBosons);
  for (unsigned long i = 0ul; i < TmpDim; ++i)
    {
      VanDerMonde[i] = i;
      TmpHeapArray[i] = 0ul;
    }
  finalState.ClearVector();


  double TmpFactor = 0.0;
  unsigned long Tmp = 0;
  double Coefficient = 1.0;
  this->ProdATemporaryNbrStateInOrbit = 1;
  bool DiscardFlag = false;
  for (int i = 0; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
    {
      StateUp[i] = VanDerMonde[i] - slaterUp[i];
      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
	{
	  TmpFactor += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
	}
      else
	{
	  DiscardFlag = true;
	}
    }
  if (DiscardFlag == false)
    {
      for (int i = 0; (i < this->NbrBosonsDown) && (DiscardFlag == false); ++i)
	{
	  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterDown[i];
	  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
	    {
	      TmpFactor += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
	    }
	  else
	    {
	      DiscardFlag = true;
	    }
	}
      if (DiscardFlag == false)
	{
	  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
	  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
	  Coefficient = 1.0;
//		  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
	  int TmpPos =  this->SymmetrizeAdAdResult(TmpFinalStateUp, TmpFinalStateDown, Coefficient);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      finalState[TmpPos] += Sign * Coefficient * exp(TmpFactor);
	    }
	}
    }
  while (Tmp < TmpDim)
    {
      if (TmpHeapArray[Tmp] < Tmp)
	{
	  if ((Tmp & 0x1ul) == 0x0ul)
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[0];
	      VanDerMonde[0] = Tmp2;
	      
	    }
	  else
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[TmpHeapArray[Tmp]];
	      VanDerMonde[TmpHeapArray[Tmp]] = Tmp2;
	    }
	  Sign *= -1.0;
	  DiscardFlag = false;
	  TmpFactor = 0.0;
	  for (int i = 0; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
	    {
	      StateUp[i] = VanDerMonde[i] - slaterUp[i];
	      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
		{
		  TmpFactor += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
		}
	      else
		{
		  DiscardFlag = true;
		}
	    }
	  if (DiscardFlag == false)
	    {
	      for (int i = 0; (i < this->NbrBosonsDown) && (DiscardFlag == false); ++i)
		{
		  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterDown[i];
		  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
		    {
		      TmpFactor += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		    }
		  else
		    {
		      DiscardFlag = true;
		    }
		}
	      if (DiscardFlag == false)
		{
		  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
		  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
		  Coefficient = 1.0;
		  int TmpPos =  this->SymmetrizeAdAdResult(TmpFinalStateUp, TmpFinalStateDown, Coefficient);
//		  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      finalState[TmpPos] += Sign* Coefficient * exp(TmpFactor);
		    }
		}
	    }
	  ++TmpHeapArray[Tmp];
	  Tmp = 0ul;
	}
      else
	{
	  TmpHeapArray[Tmp]= 0ul;
	  ++Tmp;
	}
    }
}

// convert a given state from a generic basis from the current Sz subspace basis
//
// state = reference on the vector to convert
// space = reference on the basis associated to state
// return value = converted vector

RealVector BosonOnSphereWithSU2SpinSzSymmetry::ConvertToNbodyBasis(RealVector& state, ParticleOnSphereWithSpin* space)
{
  BosonOnSphereWithSU2Spin* TmpSpace = (BosonOnSphereWithSU2Spin*) space;
  RealVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      int Pos = TmpSpace->FindStateIndex(this->StateDescriptionUp[i], this->StateDescriptionDown[i]);
      if (Pos < TmpSpace->HilbertSpaceDimension)
	{
	  TmpVector[i] =  state[Pos] * sqrt((double) this->NbrStateInOrbit[i]);
	}
    }
  return TmpVector;
}
  
// convert a given state from a generic basis to the current Sz subspace basis
//
// state = reference on the vector to convert
// space = reference on the basis associated to state
// return value = converted vector

RealVector BosonOnSphereWithSU2SpinSzSymmetry::ConvertFromNbodyBasis(RealVector& state, ParticleOnSphereWithSpin* space)
{
  BosonOnSphereWithSU2Spin* TmpSpace = (BosonOnSphereWithSU2Spin*) space;
  RealVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      double Coefficient = 1.0;
      this->ProdATemporaryNbrStateInOrbit = 1;
      int TmpPos =  this->SymmetrizeAdAdResult(TmpSpace->StateDescriptionUp[i], TmpSpace->StateDescriptionDown[i], Coefficient);
      if (TmpPos != this->HilbertSpaceDimension)
	{
	  TmpVector[i] = this->SzParitySign * Coefficient * state[TmpPos];
	}
      
    }
  return TmpVector;
}
  
