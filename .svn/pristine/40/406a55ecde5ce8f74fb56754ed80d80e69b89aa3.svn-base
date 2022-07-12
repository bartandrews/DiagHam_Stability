////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Cecile Repellin                 //
//                                                                            //
//                                                                            //
//               class of Hilbert space for bosons on CP2                     //
//              including Tz <-> -Tz and Z3 symmetry                          //
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
#include "HilbertSpace/BosonOnCP2TzZ3Symmetry.h"
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

BosonOnCP2TzZ3Symmetry::BosonOnCP2TzZ3Symmetry ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// p = number of flux quanta (determines an irreducible representation of SO(5), along with q=0 (LLL))
// totalJz = total value of jz
// totalKz = total value of kz
// memory = amount of memory granted for precalculations

BosonOnCP2TzZ3Symmetry::BosonOnCP2TzZ3Symmetry (int nbrBosons, int nbrFluxQuanta, int totalTz, int totalY, 
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
  this->TzRot1StateBosonic = new unsigned long[this->NbrLzValue];
  this->Rot2StateBosonic = new unsigned long[this->NbrLzValue];
  this->Rot1StateBosonic = new unsigned long[this->NbrLzValue];
  this->TzRot2StateBosonic = new unsigned long[this->NbrLzValue];
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
	    if ((TmpTzSymmetry & 1) != 0)
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

BosonOnCP2TzZ3Symmetry::BosonOnCP2TzZ3Symmetry(const BosonOnCP2TzZ3Symmetry& bosons)
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
  this->TzStateBosonic = new unsigned long[this->NbrLzValue];
  this->Rot1StateBosonic = new unsigned long[this->NbrLzValue];
  this->TzRot1StateBosonic = new unsigned long[this->NbrLzValue];
  this->Rot2StateBosonic = new unsigned long[this->NbrLzValue];
  this->TzRot2StateBosonic = new unsigned long[this->NbrLzValue];
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

BosonOnCP2TzZ3Symmetry::~BosonOnCP2TzZ3Symmetry ()
{
  delete[] this->Rot1StateBosonic;
  delete[] this->TzRot1StateBosonic;
  delete[] this->Rot2StateBosonic;
  delete[] this->TzRot2StateBosonic;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnCP2TzZ3Symmetry& BosonOnCP2TzZ3Symmetry::operator = (const BosonOnCP2TzZ3Symmetry& bosons)
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
  this->TzStateBosonic = new unsigned long[this->NbrLzValue];
  this->Rot1StateBosonic = new unsigned long[this->NbrLzValue];
  this->TzRot1StateBosonic = new unsigned long[this->NbrLzValue];
  this->Rot2StateBosonic = new unsigned long[this->NbrLzValue];
  this->TzRot2StateBosonic = new unsigned long[this->NbrLzValue];
  this->quantumNumberTz = bosons.quantumNumberTz;
  this->quantumNumberY = bosons.quantumNumberY;  
  this->quantumNumberR = bosons.quantumNumberR;
  this->quantumNumberS = bosons.quantumNumberS;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnCP2TzZ3Symmetry::Clone()
{
  return new BosonOnCP2TzZ3Symmetry(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured
bool BosonOnCP2TzZ3Symmetry::WriteHilbertSpace (char* fileName)
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

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& BosonOnCP2TzZ3Symmetry::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0;
  if (reference >= 0l)
    {
      Factor /= state[reference];
    }
  else
    {
      reference = 0l;
    }
  this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
  int dimOrbitalReference = this->GetSymmetryDimension(reference);
  double* SqrtCoefficients = new double [this->NbrLzValue];
  double* InvSqrtCoefficients = new double [this->NbrLzValue];
  FactorialCoefficient Coef;
  for (int k = 0; k <= this->LzMax; ++k)
    {
      int r = quantumNumberR[k];
      int s = quantumNumberS[k];
      int t = this->NbrFluxQuanta - r - s;
      Coef.SetToOne();
      Coef.FactorialDivide(r);
      Coef.FactorialDivide(s);
      Coef.FactorialDivide(t);
      Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
      SqrtCoefficients[k] = sqrt(Coef.GetNumericalValue());
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialDivide(this->TemporaryState[k]);
  
  int dimOrbital;
    for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      dimOrbital = this->GetSymmetryDimension(i);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	{
	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	    if (Index1 < this->NbrBosons)
	    {
	      while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	      {
		Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		++Index2;
	      }	  
	    }
	}
      while (Index1 < this->NbrBosons)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosons)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	  ++Index2;
	}
      if (symmetryFactor == true)
	{
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialMultiply(this->TemporaryState[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      state[i] *= Coefficient*sqrt(dimOrbitalReference)/sqrt(dimOrbital);
    }
  return state;
}
