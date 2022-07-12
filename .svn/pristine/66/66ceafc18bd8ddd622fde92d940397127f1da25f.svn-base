////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Cecile Repellin                 //
//                                                                            //
//                                                                            //
//               class of Hilbert space for bosons on CP2                     //
//                     including Tz <-> -Tz symmetry                          //
//                                                                            //
//                        last modification : 18/01/2013                      //
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
#include "HilbertSpace/BosonOnCP2TzSymmetry.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

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

BosonOnCP2TzSymmetry::BosonOnCP2TzSymmetry ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// p = number of flux quanta (determines an irreducible representation of SO(5), along with q=0 (LLL))
// totalJz = total value of jz
// totalKz = total value of kz
// memory = amount of memory granted for precalculations

BosonOnCP2TzSymmetry::BosonOnCP2TzSymmetry (int nbrBosons, int nbrFluxQuanta, int totalTz, int totalY, 
									    bool minusTzParity, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->TotalTz = totalTz;
  this->TotalY = totalY;
  this->TotalR = (this->TotalY + 3*this->TotalTz + 2*this->NbrBosons*this->NbrFluxQuanta)/6;
  this->TotalS = (this->TotalY - 3*this->TotalTz + 2*this->NbrBosons*this->NbrFluxQuanta)/6;
  this->TzParitySign = 1.0;
  if (minusTzParity == true)
    this->TzParitySign = -1.0;
  this->NbrLzValue = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2;
  this->LzMax = NbrLzValue - 1;  
  this->Minors = 0;
  this->KeptCoordinates = 0;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->TzStateBosonic = new unsigned long[this->NbrLzValue];
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrFluxQuanta, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0);
  this->quantumNumberTz = new int [this->NbrLzValue];
  this->quantumNumberY = new int [this->NbrLzValue];
  this->quantumNumberR = new int [this->NbrLzValue];
  this->quantumNumberS = new int [this->NbrLzValue];
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(TmpStateDescription, this->NbrBosons, this->NbrFluxQuanta, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0, this->NbrLzValue + this->NbrBosons, 0l);
      this->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberTz, this->quantumNumberY, this->quantumNumberR, this->quantumNumberS);
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space : get " << TmpLargeHilbertSpaceDimension << " , should be " << this->LargeHilbertSpaceDimension << endl;
	}
      this->FermionBasis = new FermionOnSphere(this->NbrBosons, 0, this->LzMax + this->NbrBosons, TmpStateDescription, this->LargeHilbertSpaceDimension);
      TmpLargeHilbertSpaceDimension = 0;
//       cout << "Unsymmetrized Basis : " << endl;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  int TmpTzSymmetry;
// 	  this->PrintState(cout, i);
// 	  cout << endl;
	  if (this->GetCanonicalState(i, TmpTzSymmetry) != this->FermionBasis->StateDescription[i])
	    this->FermionBasis->StateDescription[i] = 0x0ul;
	  else
	  {
	    if (TmpTzSymmetry == 1)
	      {
		if (this->TzParitySign > 0)
		  ++TmpLargeHilbertSpaceDimension;
		else
		  this->FermionBasis->StateDescription[i] = 0x0ul;
	      }
	    else
	      ++TmpLargeHilbertSpaceDimension;
	  }
    }
    TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];
    TmpLargeHilbertSpaceDimension = 0;
//     cout << "Symmetrized basis : " << endl;
    for (int i = 0; i < this->HilbertSpaceDimension; ++i)
      if (this->FermionBasis->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->FermionBasis->StateDescription[i];
// 	  this->PrintState(cout,i);
// 	  cout << endl;
	  ++TmpLargeHilbertSpaceDimension;
	}
    delete[] this->FermionBasis->StateDescription;
    this->FermionBasis->StateDescription = TmpStateDescription;
    this->HilbertSpaceDimension = TmpLargeHilbertSpaceDimension;
    this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
    if (this->HilbertSpaceDimension>0)
    {
      this->FermionBasis = new FermionOnSphere(this->NbrBosons, 0, this->LzMax + this->NbrBosons, TmpStateDescription, this->LargeHilbertSpaceDimension);
      cout << "dim = " << this->LargeHilbertSpaceDimension << endl;
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
      UsedMemory += this->NbrLzValue * this->FermionBasis-> LookUpTableMemorySize * sizeof(int);
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

BosonOnCP2TzSymmetry::BosonOnCP2TzSymmetry(const BosonOnCP2TzSymmetry& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->TotalR = bosons.TotalR;
  this->TotalS = bosons.TotalS;
  this->TzParitySign = bosons.TzParitySign;
  this->ProdASignature = bosons.ProdASignature;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->NbrLzValue = bosons.NbrLzValue;
  this->LzMax = bosons.LzMax;
  this->Minors = 0;
  this->KeptCoordinates = 0;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->TzStateBosonic = new unsigned long [this->NbrLzValue];
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->quantumNumberTz = bosons.quantumNumberTz;
  this->quantumNumberY = bosons.quantumNumberY;  
  this->quantumNumberR = bosons.quantumNumberR;
  this->quantumNumberS = bosons.quantumNumberS;
}

// destructor
//

BosonOnCP2TzSymmetry::~BosonOnCP2TzSymmetry ()
{
  delete[] this->TzStateBosonic;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnCP2TzSymmetry& BosonOnCP2TzSymmetry::operator = (const BosonOnCP2TzSymmetry& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->NbrLzValue = bosons.NbrLzValue;
  this->Minors = 0;
  this->KeptCoordinates = 0;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->TotalR = bosons.TotalR;
  this->TotalS = bosons.TotalS;
  this->TzParitySign = bosons.TzParitySign;
  this->ProdASignature = bosons.ProdASignature;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->TzStateBosonic = new unsigned long [this->NbrLzValue];
  this->quantumNumberTz = bosons.quantumNumberTz;
  this->quantumNumberY = bosons.quantumNumberY;  
  this->quantumNumberR = bosons.quantumNumberR;
  this->quantumNumberS = bosons.quantumNumberS;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnCP2TzSymmetry::Clone()
{
  return new BosonOnCP2TzSymmetry(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured
bool BosonOnCP2TzSymmetry::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrBosons);
  WriteLittleEndian(File, this->TotalTz);
  WriteLittleEndian(File, this->TotalY);
  if (this->HilbertSpaceDimension != 0)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateDescription[i]);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateLzMax[i]);
    }
  else
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateDescription[i]);
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateLzMax[i]);
    }
  File.close();
  return true;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnCP2TzSymmetry::AA (int index, int n1, int n2)
{
  unsigned long TmpState;
  TmpState = this->GetCanonicalState(index, this->ProdASignature);
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  if ((n1 > this->ProdATemporaryStateLzMax) || (n2 > this->ProdATemporaryStateLzMax) || 
      (this->ProdATemporaryState[n1] == 0) || (this->ProdATemporaryState[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState[n1] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n2];
  --this->ProdATemporaryState[n2];
  Coefficient *= this->ProdATemporaryState[n1];
  --this->ProdATemporaryState[n1];
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0ul;
  return sqrt(Coefficient);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnCP2TzSymmetry::AdAd (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  coefficient = this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
//   cout << this->TemporaryStateLzMax << endl;
//   cout << this->ProdASignature << endl;
  return SymmetrizeAdAdResult(coefficient);
}

// convert a given state from symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector  

RealVector BosonOnCP2TzSymmetry::ConvertToNbodyBasis(RealVector& state, BosonOnCP2& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long TmpState2;
  int Signature;  
  int NewLzMax;
  int TmpState2LzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      TmpState2 = nbodyBasis.FermionBasis->StateDescription[i];
      TmpState = this->GetCanonicalStateFromFermionicPartition(TmpState2, Signature);
      NewLzMax = this->LzMax + this->NbrBosons - 1;
      while (TmpState >> NewLzMax == 0x0ul)
     	--NewLzMax;	
      if (Signature != 0)
	{
	  if (this->TzParitySign > 0.0)
	    {
	      TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)];
	    }
	}
      else
	{
	  if (TmpState2 != TmpState)	    
	    TmpVector[i] = this->TzParitySign * M_SQRT1_2 * state[this->FindStateIndex(TmpState, NewLzMax)];
	  else
	    TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2;
	}
    }
  return TmpVector;  
}