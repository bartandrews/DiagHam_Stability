////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                          Class author Cecile Repellin                      //
//                                                                            //
//               class of Hilbert space for bosons on CP2                     //
//                                                                            //
//                                                                            //
//                        last modification : 08/01/2013                      //
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
#include "HilbertSpace/BosonOnCP2.h"
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
#include "Architecture/ArchitectureOperation/FQHESquareLatticeSymmetrizeU1U1StateOperation.h"

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

BosonOnCP2::BosonOnCP2 ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// p = number of flux quanta (determines an irreducible representation of SU(3), along with q=0 (LLL))
// totalJz = total value of jz
// totalKz = total value of kz
// memory = amount of memory granted for precalculations

BosonOnCP2::BosonOnCP2 (int nbrBosons, int nbrFluxQuanta, int totalTz, int totalY, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->TotalTz = totalTz;
  this->TotalY = totalY;
  this->TotalR = (this->TotalY + 3*this->TotalTz + 2*this->NbrBosons*this->NbrFluxQuanta)/6;
  this->TotalS = (this->TotalY - 3*this->TotalTz + 2*this->NbrBosons*this->NbrFluxQuanta)/6;
  this->NbrLzValue = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2;
  this->LzMax = NbrLzValue - 1;  
  this->Minors = 0;
  this->KeptCoordinates = 0;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrFluxQuanta, this->NbrFluxQuanta, this->NbrFluxQuanta, 0, 0);
  this->quantumNumberTz = new int [this->NbrLzValue];
  this->quantumNumberY = new int [this->NbrLzValue];
  this->quantumNumberR = new int [this->NbrLzValue];
  this->quantumNumberS = new int [this->NbrLzValue];
  cout << "dim = " << this->LargeHilbertSpaceDimension << endl;
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
// 	if (this->KPauli !=0)
// 	{
// 	 this->EvaluateHilbertSpaceDimensionExclusionPrinciple(); 
// 	}
  }
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnCP2::BosonOnCP2(const BosonOnCP2& bosons)
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
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->NbrLzValue = bosons.NbrLzValue;
  this->LzMax = bosons.LzMax;
  this->Minors = 0;
  this->KeptCoordinates = 0;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  
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

BosonOnCP2::~BosonOnCP2 ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnCP2& BosonOnCP2::operator = (const BosonOnCP2& bosons)
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
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->quantumNumberTz = bosons.quantumNumberTz;
  this->quantumNumberY = bosons.quantumNumberY;  
  this->quantumNumberR = bosons.quantumNumberR;
  this->quantumNumberS = bosons.quantumNumberS;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnCP2::Clone()
{
  return new BosonOnCP2(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured
bool BosonOnCP2::WriteHilbertSpace (char* fileName)
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

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnCP2::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], this->TemporaryState, this->TemporaryStateLzMax);
   //cout << TemporaryStateLzMax << endl;
//   Str << hex << this->FermionBasis->StateDescription[state] << dec << " " <<  this->FermionBasis->StateLzMax[state] << "   ";
//   cout << this->TotalTz << " " << this->TotalY << " ; " << this->TotalR << " " << this->TotalS << " ";
  Str << "[";
  for (int index = 0; index <= this->TemporaryStateLzMax; ++index)
  {
   if (this->TemporaryState[index] > 0)
    {
	for (int i = 0; i < this->TemporaryState[index]; ++i)
	{
	  Str << "(" <<  this->quantumNumberTz[index]   << "," << this->quantumNumberY[index] << ")";
	  }
	}
  }
  Str << "]";
  
//    Str << "  "  << this->FermionBasis->FindStateIndex(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state]) << " = " << state << endl;
  return Str;
}

// generate all states corresponding to the constraints
// 
// stateDescription = array that gives each state description
// nbrBosons = number of bosons
// currentJ = current value of j for a single particle
// currentJz = current value of jz for a single particle
// currentKz = current value of kz for a single particle
// currentTotalTz = current total value of Jz
// currentTotalY = current total value of Kz
// currentFermionicPosition = current fermionic position within the state description
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
  
long BosonOnCP2::GenerateStates(unsigned long* stateDescription, int nbrBosons, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY, int currentFermionicPosition, long pos)
{

  if (nbrBosons < 0)
    return pos;
  if (currentTotalTz + currentTzMax*nbrBosons < this->TotalTz)
    return pos;
  if (currentTotalY + currentY*nbrBosons < this->TotalY)
    return pos;
  
  if (currentTz < -currentTzMax)
   {
     --currentTzMax;
     currentTz = currentTzMax;
     currentY = currentY - 3;
   }
    
  
  if (nbrBosons == 0)
    {
      if ((currentTotalTz == this->TotalTz) && (currentTotalY == this->TotalY))
	{
	  stateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }

  if (currentY < -2*this->NbrFluxQuanta)
    return pos;

  for (int k = nbrBosons; k > 0; --k)
    {
      long TmpPos = this->GenerateStates(stateDescription, nbrBosons - k, currentTz - 2, currentTzMax, currentY, currentTotalTz + k * currentTz, currentTotalY + k * currentY, currentFermionicPosition - k - 1, pos);
      unsigned long Mask = ((0x1ul << k) - 0x1ul) << (currentFermionicPosition - k - 1);
      for (; pos < TmpPos; ++pos)
	stateDescription[pos] |= Mask;
    }
  return this->GenerateStates(stateDescription, nbrBosons, currentTz - 2, currentTzMax, currentY, currentTotalTz, currentTotalY, currentFermionicPosition - 1, pos);
};


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentJz = current value of jz for a single particle
// currentKz = current value of kz for a single particle
// currentTotalTz = current total value of Jz
// currentTotalY = current total value of Kz
// return value = Hilbert space dimension

long BosonOnCP2::EvaluateHilbertSpaceDimension(int nbrBosons, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY)
{
  if (nbrBosons < 0)
    return 0l;
  if (currentTotalTz + currentTzMax*nbrBosons < this->TotalTz)
    return 0l;
  if (currentTotalY + currentY*nbrBosons < this->TotalY)
    return 0l;
  
  if (currentTz < -currentTzMax)
   {
     --currentTzMax;
     currentTz = currentTzMax;
     currentY = currentY - 3;
   }
    
  if (nbrBosons == 0)
    {
      if ((currentTotalTz == this->TotalTz) && (currentTotalY == this->TotalY))
      {
	return 1l;
      }
      else	
	return 0l;
    }
    
  if (currentY < -2*this->NbrFluxQuanta)
    return 0l;
  
  long Count = 0;
  for (int k = nbrBosons; k >= 0; --k)
    Count += this->EvaluateHilbertSpaceDimension(nbrBosons - k, currentTz - 2, currentTzMax, currentY, currentTotalTz + k * currentTz, currentTotalY + k * currentY);
  return Count;
}

// request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
bool BosonOnCP2::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  int* rootPartition = new int[2*this->NbrBosons];
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  int Flag = 0;
  int particleIndex = 0;
  for (int ind = 0; ind <= this->TemporaryStateLzMax; ++ind)
    {
      if (this->TemporaryState[ind] > pauliK)
	return false;
      if (this->TemporaryState[ind] > 0)
      {
	for (int i = 0; i < this->TemporaryState[ind]; ++i)
	  {
	    rootPartition[particleIndex] = this->quantumNumberTz[ind];
	    rootPartition[particleIndex + 1] = this->quantumNumberY[ind];
// 	    rootPartition[particleIndex] = this->quantumNumberR[ind];
// 	    rootPartition[particleIndex + 1] = this->quantumNumberS[ind];
// 	      cout << particleIndex << " " << rootPartition[particleIndex] << " " << rootPartition[particleIndex + 1] << endl;
	    particleIndex += 2;
	  }
	}
    }
  int i = 0;
  while (i < this->NbrBosons - 1)
    {
      int j = i + 1;
      while  (j < this->NbrBosons)
      {
       int distance = 3*abs(rootPartition[i << 1] - rootPartition[j << 1]) + 2*abs(rootPartition[(i << 1) + 1] - rootPartition[(j << 1) + 1]);
// 	  int distance = abs(rootPartition[i << 1] - rootPartition[j << 1]) + abs(rootPartition[(i << 1) + 1] - rootPartition[(j << 1) + 1]) + abs (rootPartition[i << 1] - rootPartition[j << 1] + rootPartition[(i << 1) + 1] - rootPartition[(j << 1) + 1]);
//        cout << (i<<1) << " , " << (j<<1) << " : " << rootPartition[(i << 1)] << " " << rootPartition[j << 1] << " " << rootPartition[(i << 1) + 1] << " " << rootPartition[(j << 1) + 1] << " : " << distance << endl;
       if (distance < 6*pauliR)
	  return false;
       else
	 j += 1;
      }
     i += 1;
    }
  return true; 
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

// RealVector& BosonOnCP2::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
// {
//   unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
//   unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
//   double Factor = 1.0;
//   if (reference >= 0l)
//     {
//       Factor /= state[reference];
//     }
//   else
//     {
//       reference = 0l;
//     }
//   this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
//   double* SqrtCoefficients = new double [this->NbrLzValue];
//   double* InvSqrtCoefficients = new double [this->NbrLzValue];
//   FactorialCoefficient Coef;
//   for (int k = 0; k <= this->LzMax; ++k)
//     {
//       int r = quantumNumberR[k];
//       int s = quantumNumberS[k];
//       int t = this->NbrFluxQuanta - r - s;
//       Coef.SetToOne();
//       Coef.FactorialDivide(r);
//       Coef.FactorialDivide(s);
//       Coef.FactorialDivide(t);
//       Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
//       SqrtCoefficients[k] = sqrt(Coef.GetNumericalValue());
//       InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
//     }
//   FactorialCoefficient ReferenceFactorial;
//   FactorialCoefficient Factorial;
//   this->FermionToBoson(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], 
// 		       this->TemporaryState, this->TemporaryStateLzMax);
//   for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
//     if (this->TemporaryState[k] > 1)
//       ReferenceFactorial.FactorialDivide(this->TemporaryState[k]);
//   for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     {
//       this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
//       
//       int Index1 = 0;
//       int Index2 = 0;
//       double Coefficient = Factor;
//       while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
// 	{
// 	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
// 	    {
// 	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
// 	      ++Index1;
// 	    }
// 	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
// 	    {
// 	      ++Index1;
// 	      ++Index2;
// 	    }
// 	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
// 	    {
// 	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
// 	      ++Index2;
// 	    }	  
// 	}
//       while (Index1 < this->NbrBosons)
// 	{
// 	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
// 	  ++Index1;
// 	}
//       while (Index2 < this->NbrBosons)
// 	{
// 	  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
// 	  ++Index2;
// 	}
//       if (symmetryFactor == true)
// 	{
// 	  Factorial = ReferenceFactorial;
// 	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
// 			       this->TemporaryState, this->TemporaryStateLzMax);
// 	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
// 	    if (this->TemporaryState[k] > 1)
// 	      Factorial.FactorialMultiply(this->TemporaryState[k]);
// 	  Coefficient *= sqrt(Factorial.GetNumericalValue());
// 	}
//       state[i] *= Coefficient;
//     }
//   return state;
// }

RealVector& BosonOnCP2::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
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
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
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
      state[i] *= Coefficient;
    }
  return state;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given (Jz,Kz) sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// jzSector = Jz sector in which the density matrix has to be evaluated 
// kzSector = Kz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnCP2::EvaluatePartialDensityMatrixParticlePartition(int nbrBosonSector, int tzSector, int ySector, RealVector& groundState, AbstractArchitecture* architecture)
{  
  if (nbrBosonSector == 0)
    {
      if ((tzSector == 0) and (ySector ==0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  if (nbrBosonSector == this->NbrBosons)
    {
      if ((tzSector == this->TotalTz) and (ySector == this->TotalY))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;
  int ComplementaryRSector = (this->TotalY - ySector + 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrBosonSector*this->NbrFluxQuanta);
  int ComplementarySSector = (this->TotalY - ySector - 3*(this->TotalTz - tzSector) + 2*ComplementaryNbrBosonSector*this->NbrFluxQuanta);
  if ((ComplementaryRSector < 0) || (ComplementarySSector < 0) || ((ComplementaryRSector % 6) != 0) || ((ComplementarySSector % 6) != 0) || (ComplementaryRSector + ComplementarySSector > 6*this->NbrFluxQuanta*ComplementaryNbrBosonSector))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
  cout << "nbr boson = " << nbrBosonSector << ", tz = " << tzSector << ", y = " << ySector << endl;
  BosonOnCP2 TmpDestinationHilbertSpace(nbrBosonSector, this->NbrFluxQuanta, tzSector, ySector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  BosonOnCP2 TmpHilbertSpace(ComplementaryNbrBosonSector, this->NbrFluxQuanta, this->TotalTz - tzSector, this->TotalY - ySector);

  
  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &TmpDestinationHilbertSpace, &TmpHilbertSpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);

  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to be stored
// return value = number of components that have been added to the density matrix

long BosonOnCP2::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace, RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
   BosonOnCP2* TmpHilbertSpace =  (BosonOnCP2*) complementaryHilbertSpace;
   BosonOnCP2* TmpDestinationHilbertSpace =  (BosonOnCP2*) destinationHilbertSpace;
   int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
   int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
   unsigned long* TmpMonomial2 = new unsigned long [NbrBosonSector];
   unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
   unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
   int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int MaxIndex = minIndex + nbrIndex;
   long TmpNbrNonZeroElements = 0l;

   double* LogFactorials = new double[this->NbrBosons + 1];
   LogFactorials[0] = 0.0;
   LogFactorials[1] = 0.0;
   for (int i = 2 ; i <= this->NbrBosons; ++i)
     LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
   double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];

   for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
	TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
   for (; minIndex < MaxIndex; ++minIndex)    
     {
      int Pos = 0;
      TmpHilbertSpace->ConvertToMonomial(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpMonomial1);
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->FermionBasis->StateDescription[j], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[j], TmpMonomial2);
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < NbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < NbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < NbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 
	      
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      TmpStateCoefficient[Pos] = exp(TmpFactorial);
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		    ++TmpNbrNonZeroElements;
		  }
	    }
	}
     }
   delete[] TmpMonomial2;
   delete[] TmpMonomial1;
   delete[] TmpMonomial3;
   delete[] TmpStatePosition;
   delete[] TmpStatePosition2;
   delete[] TmpStateCoefficient;
   delete[] TmpDestinationLogFactorials;
   return TmpNbrNonZeroElements;
}
