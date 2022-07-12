////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on the 4D space  torus x sphere              //
//            with magnetic translations and for system size such that        //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 22/02/2017                      //
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
#include "HilbertSpace/BosonOnT2xS2WithMagneticTranslationsShort.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "Architecture/ArchitectureOperation/FQHETorusParticleEntanglementSpectrumOperation.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "MathTools/IntegerAlgebraTools.h"

#include <math.h>
#include <algorithm>
#include <set>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

BosonOnT2xS2WithMagneticTranslationsShort::BosonOnT2xS2WithMagneticTranslationsShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// nbrFluxQuantumTorus = number of flux quanta piercing the torus
// kxMomentum = momentum in the x direction for the torus
// kyMomentum = momentum in the y direction for the torus
// nbrFluxQuantumSphere = number of flux quanta piercing the sphere
// totalLz = projection of the total angular momentum along the z axis for the sphere

BosonOnT2xS2WithMagneticTranslationsShort::BosonOnT2xS2WithMagneticTranslationsShort (int nbrBosons, int nbrFluxQuantumTorus, int kxMomentum, int kyMomentum,
										      int nbrFluxQuantumSphere, int totalLz)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->NbrFluxQuantumTorus = nbrFluxQuantumTorus;
  this->NbrFluxQuantumSphere = nbrFluxQuantumSphere;
  this->TotalLz = totalLz;
  this->ShiftedTotalLz = (this->TotalLz + (this->NbrFluxQuantumSphere * this->NbrBosons)) / 2;
  this->NbrLzValues =  this->NbrFluxQuantumSphere + 1;
  this->MaxMomentum = (this->NbrFluxQuantumSphere + 1) * this->NbrFluxQuantumTorus;
  this->FermionicMaxMomentum = this->MaxMomentum + this->NbrBosons - 1;
  this->NbrKyValue = this->MaxMomentum + 1;
 
  this->MomentumModulo = FindGCD(this->NbrBosons, this->NbrFluxQuantumTorus);
  this->KxMomentum = kxMomentum % this->MomentumModulo;
  this->KyMomentum = kyMomentum % this->NbrFluxQuantumTorus;

  this->StateShift = this->MaxMomentum / this->MomentumModulo;

  this->LastMomentumMask = 0x1ul << (this->MaxMomentum + this->NbrBosons - 1);

  this->TemporaryState = new unsigned long [this->NbrKyValue];
  this->ProdATemporaryState = new unsigned long [this->NbrKyValue];

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrFluxQuantumTorus - 1, this->NbrFluxQuantumSphere, 0, 0);
  if (this->LargeHilbertSpaceDimension != 0l)
    {
      long TmpLargeHilbertSpaceDimension = this->GenerateStates();
      if (TmpLargeHilbertSpaceDimension != 0l)
	{
	  this->LargeHilbertSpaceDimension = TmpLargeHilbertSpaceDimension;
	  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
	  this->Flag.Initialize();
	  this->GenerateLookUpTable(0);
#ifdef __DEBUG__
	  unsigned long UsedMemory = 0;
	  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
	  UsedMemory += this->NbrKyValue * sizeof(int);
	  UsedMemory += this->NbrKyValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
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
      else
	{
	  this->HilbertSpaceDimension = 0;
	  this->LargeHilbertSpaceDimension = 0l;
	}
    }
  else
    {
      this->HilbertSpaceDimension = 0;
    }
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnT2xS2WithMagneticTranslationsShort::BosonOnT2xS2WithMagneticTranslationsShort(const BosonOnT2xS2WithMagneticTranslationsShort& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrKyValue = bosons.NbrKyValue;
  this->MomentumModulo = bosons.MomentumModulo;
  this->FermionicMaxMomentum = bosons.FermionicMaxMomentum;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->NbrFluxQuantumTorus = this->NbrFluxQuantumTorus;
  this->NbrFluxQuantumSphere = this->NbrFluxQuantumSphere;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons.ShiftedTotalLz;
  this->NbrLzValues = bosons.NbrLzValues;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMask = bosons.LastMomentumMask;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->Flag = bosons.Flag;
  this->TemporaryState = new unsigned long [this->NbrKyValue];
  this->ProdATemporaryState = new unsigned long [this->NbrKyValue];
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

}

// destructor
//

BosonOnT2xS2WithMagneticTranslationsShort::~BosonOnT2xS2WithMagneticTranslationsShort ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnT2xS2WithMagneticTranslationsShort& BosonOnT2xS2WithMagneticTranslationsShort::operator = (const BosonOnT2xS2WithMagneticTranslationsShort& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
    }
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->FermionicMaxMomentum = bosons.FermionicMaxMomentum;
  this->NbrKyValue = bosons.NbrKyValue;
  this->MomentumModulo = bosons.MomentumModulo;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->NbrFluxQuantumTorus = this->NbrFluxQuantumTorus;
  this->NbrFluxQuantumSphere = this->NbrFluxQuantumSphere;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons.ShiftedTotalLz;
  this->NbrLzValues = bosons.NbrLzValues;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMask = bosons.LastMomentumMask;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->TemporaryState = new unsigned long [this->NbrKyValue];
  this->ProdATemporaryState = new unsigned long [this->NbrKyValue];
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->Flag = bosons.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnT2xS2WithMagneticTranslationsShort::Clone()
{
  return new BosonOnT2xS2WithMagneticTranslationsShort(*this);
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnT2xS2WithMagneticTranslationsShort::PrintStateMonomial (ostream& Str, long state)
{
  int TmpLz;
  int TmpKz;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  cout << this->StateDescription[state] << endl;
  this->ConvertToMonomial(this->StateDescription[state], this->MaxMomentum + this->NbrBosons - 1, TmpMonomial);
  Str << "[";
  Str << "(" << TmpLz << "," << TmpKz << ")";
  for (int i = 1; i < this->NbrBosons; ++i)
    {
      this->GetLinearizedIndex(TmpMonomial[i], TmpLz, TmpKz);
      Str << "," << "(" << TmpLz << "," << TmpKz << ")";
    }
  Str << "]";
  delete[] TmpMonomial;
  return Str;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnT2xS2WithMagneticTranslationsShort::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescription[state], this->MaxMomentum + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  int TmpLz;
  int TmpKy;
  Str << "[";
  for (int i = 0; i <= this->TemporaryStateKyMax; ++i)
    {
      if (this->TemporaryState[i] > 0)
	{
	  this->GetLinearizedIndex(i, TmpKy, TmpLz);
	  for (int j = 0; j < this->TemporaryState[i]; ++j)
	    Str << "(" << TmpKy << "," << (2 * TmpLz - this->NbrFluxQuantumSphere) << ")";
	}
    }
  Str << "]";
 return Str;
}


// generate all states with both the kx and ky constraint
// 
// return value = new dimension of the Hilbert space

long BosonOnT2xS2WithMagneticTranslationsShort::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  long TmpLargeHilbertSpaceDimension = this->RawGenerateStates(this->NbrBosons, this->NbrFluxQuantumTorus - 1, this->NbrFluxQuantumSphere, 0l, 0, 0);
  if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
    {
      cout << "error : dimension mismatch while generating the Hilbert space ( is" << TmpLargeHilbertSpaceDimension << ", should be " 
	   << this->LargeHilbertSpaceDimension << ")" << endl;
      return 0l;
    }

  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      int NbrTranslation = 0;
      unsigned long TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(this->StateDescription[i], NbrTranslation);
      if (NbrTranslation == 0)
	{
	  ++TmpLargeHilbertSpaceDimension;
	}
      else
	{
	  this->StateDescription[i] = 0x0ul;
	}
    }

  unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	  unsigned long TmpState = this->StateDescription[i];
	  unsigned long TmpReferenceState = TmpState;
	  int TmpOrbitSize = 1;
	  this->ApplySingleTranslation(TmpState);
	  while (TmpState != TmpReferenceState)
	    {
	      this->ApplySingleTranslation(TmpState);
	      ++TmpOrbitSize;
	    }
	  this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = TmpOrbitSize;
	  ++TmpLargeHilbertSpaceDimension;
	}	
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  return TmpLargeHilbertSpaceDimension;
}

// generate all states corresponding to the ky/Lz constraints without taking care of the kx constraint
// 
// nbrBosons = number of bosons
// currentKyMax = torus momentum maximum value for bosons that are still to be placed
// currentLz = current maximum angular momentum for bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentKy = current value of the momentum along y for the torus 
// currentTotalLz = current total angular momentum along the z axis for the sphere  
// return value = position from which new states have to be stored

long BosonOnT2xS2WithMagneticTranslationsShort::RawGenerateStates(int nbrBosons, int currentKyMax, int currentLz, long pos, int currentKy, int currentTotalLz)
{
  if (nbrBosons == 0)
    {
      if (((currentKy % this->NbrFluxQuantumTorus) == this->KyMomentum) && (currentTotalLz == this->ShiftedTotalLz))
	{
	  this->StateDescription[pos] = 0x0ul;
	  return pos + 1l;
	}
      else
	{
	  return pos;
	}
    }

  if (currentLz < 0)
    {
      currentLz = this->NbrFluxQuantumSphere;
      --currentKyMax;
    }
  if (currentKyMax < 0)
    {
      return pos;   
    }

  long TmpPos = pos;
  for (int i = nbrBosons; i > 0; --i)
    {
      TmpPos = this->RawGenerateStates(nbrBosons - i, currentKyMax, currentLz - 1, pos, currentKy + (i * currentKyMax), currentTotalLz + (i * currentLz));
      unsigned long Mask = ((0x1ul << i) - 0x1ul) << ((currentKyMax * this->NbrLzValues) + currentLz + nbrBosons - i);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return this->RawGenerateStates(nbrBosons, currentKyMax, currentLz - 1, pos, currentKy, currentTotalLz);
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKyMax = torus momentum maximum value for bosons that are still to be placed
// currentLz = current maximum angular momentum for bosons that are still to be placed
// currentKy = current value of the momentum along y for the torus 
// currentTotalLz = current total angular momentum along the z axis for the sphere  
// return value = Hilbert space dimension

long BosonOnT2xS2WithMagneticTranslationsShort::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKyMax, int currentLz, int currentKy, int currentTotalLz)
{
  if (nbrBosons == 0)
    {
      if (((currentKy % this->NbrFluxQuantumTorus) == this->KyMomentum) && (currentTotalLz == this->ShiftedTotalLz))
	{
	  return 1l;
	}
      else
	{
	  return 0l;
	}
    }

  if (currentLz < 0)
    {
      currentLz = this->NbrFluxQuantumSphere;
      --currentKyMax;
    }
  if (currentKyMax < 0)
    {
      return 0l;   
    }

  long TmpNbrStates = 0l;
  for (int i = nbrBosons; i >= 0; --i)
    {
      TmpNbrStates += this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKyMax, currentLz - 1, currentKy + (i * currentKyMax), currentTotalLz + (i * currentLz));
    }

  return TmpNbrStates;
}

