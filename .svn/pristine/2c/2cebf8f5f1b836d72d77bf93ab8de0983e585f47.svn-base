////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of fermion on a torus taking into account magnetic translations    //
//                                                                            //
//                        last modification : 10/09/2003                      //
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
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "Architecture/ArchitectureOperation/FQHETorusParticleEntanglementSpectrumOperation.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "Matrix/Matrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"

#include <math.h>
#include <limits.h>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


// binding to the LAPACK zgetrf routine for LU decomposition and back-substitution
//
#ifdef __LAPACK__
extern "C" void FORTRAN_NAME(zgetrf)(const int* dimensionM, const int* dimensionN, const doublecomplex* matrixA,
				     const int* leadingDimensionA, const int *ipiv, const int *info);
#endif

// default constructor
// 

FermionOnTorusWithMagneticTranslations::FermionOnTorusWithMagneticTranslations ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)

FermionOnTorusWithMagneticTranslations::FermionOnTorusWithMagneticTranslations (int nbrFermions, int maxMomentum, 
										int xMomentum, int yMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;  

  this->NbrMomentum = this->MaxMomentum + 1;
  this->MomentumModulo = FindGCD(this->NbrFermions, this->MaxMomentum);
  cout << "MomentumModulo=" << MomentumModulo<<endl;
  this->XMomentum = xMomentum % this->MomentumModulo;
  this->YMomentum = yMomentum % this->MaxMomentum;

  this->StateShift = this->MaxMomentum / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->MaxMomentum - this->StateShift;
  this->MomentumMask = 0x1ul;
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= 0x1ul;
    }
#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->MaxMomentum + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->MaxMomentum + 1) >> 1);
#endif
  if ((this->MaxMomentum & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;

  this->MaximumSignLookUp = 16;
  this->GenerateSignLookUpTable();
  long TmpDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum);
  cout << "Max dimension: " << TmpDimension << endl;
  if (TmpDimension>INT_MAX)
    cout << "Max-Dimension surpasses integer representation..."<<endl;
  this->LargeHilbertSpaceDimension = this->GenerateStates(TmpDimension);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Actual dimension: " << this->LargeHilbertSpaceDimension << endl;

  this->Flag.Initialize();
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrMomentum * sizeof(int);
//  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnTorusWithMagneticTranslations::FermionOnTorusWithMagneticTranslations(const FermionOnTorusWithMagneticTranslations& fermions)
{
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;
}

// destructor
//

FermionOnTorusWithMagneticTranslations::~FermionOnTorusWithMagneticTranslations ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorusWithMagneticTranslations& FermionOnTorusWithMagneticTranslations::operator = (const FermionOnTorusWithMagneticTranslations& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorusWithMagneticTranslations::Clone()
{
  return new FermionOnTorusWithMagneticTranslations(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnTorusWithMagneticTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->XMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
  List<AbstractQuantumNumber*> TmpList2;
  TmpList2 += new VectorQuantumNumber (TmpList);
  return TmpList2;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnTorusWithMagneticTranslations::GetQuantumNumber (int index)
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->XMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
  return new VectorQuantumNumber (TmpList);
}

// get momemtum value in the y direction of a given state
//
// index = state index
// return value = state momentum in the y direction

int FermionOnTorusWithMagneticTranslations::GetYMomentumValue(int index)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned long State = this->StateDescription[index];
  int Momentum = 0;
  for (int i = 0; i <= StateMaxMomentum; ++i)
    {
      Momentum += ((State >> i ) & 0x1ul) * i;
    }
  return (Momentum % this->MomentumModulo);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnTorusWithMagneticTranslations::ExtractSubspace (AbstractQuantumNumber& q, 
									       SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != (AbstractQuantumNumber::Vector | AbstractQuantumNumber::PeriodicMomentum))
    return 0;
  if (((VectorQuantumNumber&) q).GetQuantumNumbers().GetNbrElement() != 2)
    return 0;
  if ((this->XMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[0]))->GetMomentum()) &&
      (this->YMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[1]))->GetMomentum()))
    return this;
  else 
    return 0;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithMagneticTranslations::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned long State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (0x1ul << n1)) == 0) || 
      ((State & (0x1ul << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(0x1ul << n2);
  if (NewMaxMomentum == n2)
    while ((TmpState >> NewMaxMomentum) == 0)
      --NewMaxMomentum;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(0x1ul << n1);
  if (TmpState == 0x0ul)
    {
      NewMaxMomentum = 0;
    }
  else
    {
      if (NewMaxMomentum == n1)
	while ((TmpState >> NewMaxMomentum) == 0)
	  --NewMaxMomentum;
    }
  if ((TmpState & (0x1ul << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      NewMaxMomentum = m2;
    }
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
  if ((TmpState & (0x1ul << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      NewMaxMomentum = m1;
    }
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
  TmpState = this->FindCanonicalForm(TmpState, NewMaxMomentum, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewMaxMomentum) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewMaxMomentum);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift;
  return TmpIndex;
}

// apply a^+_m a_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation operator
// return value =  resulting multiplicative factor 

double FermionOnTorusWithMagneticTranslations::AdA (int index, int m)
{
  if (this->StateMaxMomentum[index] < m)
    return 0.0;
  if ((this->StateDescription[index] & (0x1ul << m)) == 0)
    return 0.0;
  else
    return 1.0;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnTorusWithMagneticTranslations::AA (int index, int n1, int n2)
{
  this->ProdATemporaryStateMaxMomentum = this->StateMaxMomentum[index];
  this->ProdATemporaryState = this->StateDescription[index];
  if ((n1 >  this->ProdATemporaryStateMaxMomentum) || (n2 >  this->ProdATemporaryStateMaxMomentum) || ((this->ProdATemporaryState & (0x1ul << n1)) == 0x0ul) || 
      ((this->ProdATemporaryState & (0x1ul << n2)) == 0x0ul) || (n1 == n2))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  if (this->ProdATemporaryStateMaxMomentum == n2)
    while ((this->ProdATemporaryState >> this->ProdATemporaryStateMaxMomentum) == 0x0ul)
      --this->ProdATemporaryStateMaxMomentum;
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdATemporaryStateMaxMomentum = 0;
    }
  else
    {
      if (this->ProdATemporaryStateMaxMomentum == n1)
	while ((this->ProdATemporaryState >> this->ProdATemporaryStateMaxMomentum) == 0x0ul)
	  --this->ProdATemporaryStateMaxMomentum;
    }
  return Coefficient;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnTorusWithMagneticTranslations::ProdA (int index, int* n, int nbrIndices)
{
  this->ProdATemporaryStateMaxMomentum = this->StateMaxMomentum[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = n[i];
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
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdATemporaryStateMaxMomentum = 0;
      return Coefficient;      
    }
  while (((this->ProdATemporaryState >> this->ProdATemporaryStateMaxMomentum) == 0) && (this->ProdATemporaryStateMaxMomentum > 0))
    --this->ProdATemporaryStateMaxMomentum;

  return Coefficient;
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithMagneticTranslations::AdAd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if (((TmpState & (0x1ul << m2)) != 0x0ul) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  int NewMaxMomentum =  this->ProdATemporaryStateMaxMomentum;
  if (m2 > NewMaxMomentum)
    {
      NewMaxMomentum = m2;
    }
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
  if ((TmpState & (0x1ul << m1)) != 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      NewMaxMomentum = m1;
    }
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
  TmpState = this->FindCanonicalForm(TmpState, NewMaxMomentum, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewMaxMomentum) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewMaxMomentum);
  coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift;
  return TmpIndex;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithMagneticTranslations::ProdAd (int* m, int nbrIndices, double& coefficient, int& nbrTranslation)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewMaxMomentum = this->ProdATemporaryStateMaxMomentum;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = m[i];
      if ((TmpState & (0x1l << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewMaxMomentum)
	{
	  NewMaxMomentum = Index;
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
  TmpState = this->FindCanonicalForm(TmpState, NewMaxMomentum, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewMaxMomentum) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewMaxMomentum);
  coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift;
  return TmpIndex;
}

// apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator
// return value =  multiplicative factor 

double FermionOnTorusWithMagneticTranslations::A (int index, int n)
{
  this->ProdATemporaryStateMaxMomentum = this->StateMaxMomentum[index];
  this->ProdATemporaryState = this->StateDescription[index];
  if ((n >  this->ProdATemporaryStateMaxMomentum) || ((this->ProdATemporaryState & (0x1ul << n)) == 0x0ul))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n) & this->SignLookUpTableMask[n]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n);
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdATemporaryStateMaxMomentum = 0;
    }
  else
    {
      if (this->ProdATemporaryStateMaxMomentum == n)
	while ((this->ProdATemporaryState >> this->ProdATemporaryStateMaxMomentum) == 0x0ul)
	  --this->ProdATemporaryStateMaxMomentum;
    }
  return Coefficient;
}

// safe version to find canonical form of a state description (with additional checks)
//
// stateDescription = unsigned integer describing the state
// maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description

unsigned long FermionOnTorusWithMagneticTranslations::SafeFindCanonicalForm(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;
  int index = 1;
  // testing
  if (bitcount(stateDescription)!=this->NbrFermions)
    {
      cout << "attention: wrong number of particles in FindCanonicalForm"<<endl;
    }

  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
  if (bitcount(stateDescription)!=this->NbrFermions)
    {
      cout << "attention: wrong number of particles in FindCanonicalForm after shift"<<endl;
    }

  while (stateDescriptionReference != stateDescription)
    {
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = index;
	}
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      ++index;
    }
  if (nbrTranslation != 0)
    {
      maxMomentum = this->MaxMomentum;
      stateDescription = 0x1ul << this->MaxMomentum;
      while ((CanonicalState & stateDescription) ==0)      
	{
	  --maxMomentum;
	  stateDescription >>= 1;
	}
      nbrTranslation = index - nbrTranslation;
    }
  return CanonicalState;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnTorusWithMagneticTranslations::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->MaxMomentum; ++i)
    Str << ((TmpState >> i) & 0x1ul) << " ";
//   Str << "  (" << hex << this->ReorderingSign[state] << dec << ")";
//   Str << " " << this->FindStateIndex(this->StateDescription[state], this->StateMaxMomentum[state]);
//   if (this->FindStateIndex(this->StateDescription[state], this->StateMaxMomentum[state]) != state)
//     {
//       Str << "  error";
//     }
  return Str;
}

#include <string>
using std::string;

string convBase(unsigned long v, long base)
{
  string digits = "0123456789abcdef";
  string result;
  if((base < 2) || (base > 16)) {
    result = "Error: base out of range.";
  }
  else {
    do {
      result = digits[v % base] + result;
      v /= base;
    }
    while(v);
  }
  return result;
}

// generate all states corresponding to the constraints
// tmpDimension = max dimension of Hilbert space (to be reduced by symmetries)
// return value = hilbert space dimension

long FermionOnTorusWithMagneticTranslations::GenerateStates(long tmpDimension)
{
  this->StateDescription = new unsigned long [tmpDimension];
  this->StateMaxMomentum = new int [tmpDimension];
  this->LargeHilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, this->MaxMomentum - 1, 0l, 0);
  long* TmpNbrStateDescription = new long [this->MaxMomentum + 1];  
  for (int i = 0; i <= this->MaxMomentum; ++i)
    {
      TmpNbrStateDescription[i] = 0l;
    }
  int NbrTranslation;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->StateDescription[i] = this->FindCanonicalForm(this->StateDescription[i], this->StateMaxMomentum[i], NbrTranslation);
      ++TmpNbrStateDescription[this->StateMaxMomentum[i]];
    }
  unsigned long** TmpStateDescription = new unsigned long* [this->MaxMomentum + 1];  
  bool** TmpStateDescriptionFlag = new bool* [this->MaxMomentum + 1];  
  for (int i = 0; i <= this->MaxMomentum; ++i)
    {
      if (TmpNbrStateDescription[i] != 0l)
	{
	  TmpStateDescription[i] = new unsigned long [TmpNbrStateDescription[i]];
	  TmpStateDescriptionFlag[i] = new bool [TmpNbrStateDescription[i]];
	}
      else
	{
	  TmpStateDescription[i] = 0;
	  TmpStateDescriptionFlag[i] = 0;
	}
      TmpNbrStateDescription[i] = 0l;
    }
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      TmpStateDescription[this->StateMaxMomentum[i]][TmpNbrStateDescription[this->StateMaxMomentum[i]]++] = this->StateDescription[i];      
    }

  long TmpLargeHilbertSpaceDimension = 0l;
  for (int i = 0; i <= this->MaxMomentum; ++i) 
    {
      int CurrentNbrState = TmpNbrStateDescription[i];
      if (CurrentNbrState > 0)
	{
	  unsigned long* TmpStateArray = TmpStateDescription[i];
	  bool* TmpStateArrayFlag = TmpStateDescriptionFlag[i];
	  SortArrayUpOrdering(TmpStateArray, CurrentNbrState);
	  long TmpNbr = CurrentNbrState;
	  if (this->TestXMomentumConstraint(TmpStateArray[0], i) == true)
	    {
	      TmpStateArrayFlag[0] = true;
	    }
	  else
	    {
	      TmpStateArrayFlag[0] = false;
	      --TmpNbr;	      
	    }
	  for (long j = 1l; j < CurrentNbrState; ++j)
	    {
	      while ((j < CurrentNbrState) && (TmpStateArray[j - 1] == TmpStateArray[j]))
		{
		  TmpStateArrayFlag[j] = false;
		  --TmpNbr;
		  ++j;
		}
	      if (j < CurrentNbrState)
		{
		  if (this->TestXMomentumConstraint(TmpStateArray[j], i) == true)
		    {
		      TmpStateArrayFlag[j] = true;
		    }
		  else
		    {
		      TmpStateArrayFlag[j] = false;
		      --TmpNbr;	      
		    }
		}
	    }
	  TmpLargeHilbertSpaceDimension += TmpNbr;
	}
    }
  delete[] this->StateMaxMomentum;
  delete[] this->StateDescription;
  this->StateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];
  this->StateMaxMomentum = new int [TmpLargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  this->ReorderingSign = new unsigned long [TmpLargeHilbertSpaceDimension];
  long Pos = 0l;
  unsigned long TmpState;
  unsigned long TmpState2;
  int TmpSignature;
  unsigned long TmpReorderingSign;
  int TmpNbrParticle;
  for (int i = 0; i <= this->MaxMomentum; ++i) 
    {
      long CurrentNbrState = TmpNbrStateDescription[i];
      if (CurrentNbrState > 0l)
	{
	  unsigned long* TmpStateArray = TmpStateDescription[i];
	  bool* TmpStateArrayFlag = TmpStateDescriptionFlag[i];
	  for (long j = 0l; j < CurrentNbrState; ++j)
	    {
	      if (TmpStateArrayFlag[j] == true)
		{
		  this->StateDescription[Pos] = TmpStateArray[j];
		  this->StateMaxMomentum[Pos] = i;
		  this->NbrStateInOrbit[Pos] = this->FindNumberXTranslation(this->StateDescription[Pos]);
		  TmpState = this->StateDescription[Pos];
		  TmpSignature = 0;
		  TmpReorderingSign = (unsigned long) 0;
		  if ((this->NbrFermions & 1) == 0)
		    {
		      for (int k = 1; k <= this->NbrStateInOrbit[Pos]; ++k)
			{
			  TmpState2 = TmpState & this->MomentumMask;
			  TmpState =  (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
			  TmpNbrParticle = 0;
			  for (int l = 0; l < this->StateShift; ++l)
			    {
			      if (TmpState2 & 1)
				++TmpNbrParticle;
			      TmpState2 >>= 1;
			    }
			  if (TmpNbrParticle & 1)
			    {
			      TmpReorderingSign |= (((TmpReorderingSign << 1) & (0x1ul << k))) ^ (0x1ul << k);

			      ++TmpSignature;
			    }
			  else
			    {
			      TmpReorderingSign |= ((TmpReorderingSign << 1) & (0x1ul << k));
			    }			  
			}
		    }
		  this->ReorderingSign[Pos] = TmpReorderingSign;
		  ++Pos;
		}
	    }
	  delete[] TmpStateArray;
	  delete[] TmpStateArrayFlag;
	}
    }
  delete[] TmpStateDescription;
  delete[] TmpStateDescriptionFlag;
  delete[] TmpNbrStateDescription;

  return TmpLargeHilbertSpaceDimension;
}

// generate all states corresponding to the constraints (without taking into the canonical form) 
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentYMomentum = current value of the momentum in the y direction
// return value = position from which new states have to be stored

long FermionOnTorusWithMagneticTranslations::RawGenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, long pos, int currentYMomentum)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentMaxMomentum + 1)))
    return pos;
  if (nbrFermions == 1)
    {
      int i = this->YMomentum - (currentYMomentum % this->MaxMomentum);
      if (i < 0)
	i += this->MaxMomentum;
      for (; i <= currentMaxMomentum; i += this->MaxMomentum)
	{
	  this->StateDescription[pos] = 0x1ul << i;
	  this->StateMaxMomentum[pos] = maxMomentum;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  long TmpPos = this->RawGenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos, currentYMomentum + currentMaxMomentum);
  unsigned long Mask = 0x1ul << currentMaxMomentum;
  for (long i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->RawGenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentYMomentum);
  else
    return this->RawGenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentYMomentum);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnTorusWithMagneticTranslations::GenerateLookUpTable(int memory)
{
  // evaluate look-up table size
  memory /= (4 * this->NbrMomentum);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrMomentum)
    this->MaximumLookUpShift = this->NbrMomentum;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrMomentum];
  this->LookUpTableShift = new int [this->NbrMomentum];
  for (int i = 0; i < this->NbrMomentum; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentMaxMomentum = this->StateMaxMomentum[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
  if (CurrentMaxMomentum < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentMaxMomentum] = 0;
  else
    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
  long CurrentLookUpTableValue = 0;
  long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue < TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      ++CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
	  ++CurrentLookUpTableValue;
	  while (CurrentLookUpTableValue <= this->LookUpTableMemorySize)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      ++CurrentLookUpTableValue;
	    }
// 	  --CurrentMaxMomentum;
//  	  while (CurrentMaxMomentum > this->StateMaxMomentum[i])
//  	    {
//  	      this->LookUpTableShift[CurrentMaxMomentum] = -1;
//  	      --CurrentMaxMomentum;
//  	    }
 	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  if (CurrentMaxMomentum < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentMaxMomentum] = 0;
	  else
	    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
	  TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
	  CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = 0;
	  while (CurrentLookUpTableValue < TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      ++CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      ++CurrentLookUpTableValue;
	      while (CurrentLookUpTableValue < TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  ++CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  ++CurrentLookUpTableValue;
  while (CurrentLookUpTableValue <= this->LookUpTableMemorySize)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      ++CurrentLookUpTableValue;
    }
  
  this->RescalingFactors = new double* [this->NbrMomentum];
  for (int i = 1; i <= this->MaxMomentum; ++i)
    {
      this->RescalingFactors[i] = new double [this->NbrMomentum];
      for (int j = 1; j <= this->MaxMomentum; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}

// generate look-up table associated to sign calculations
// 

void FermionOnTorusWithMagneticTranslations::GenerateSignLookUpTable()
{
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  this->NbrParticleLookUpTable = new int [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	{
	  this->SignLookUpTable[j] = -1.0;
	  this->NbrParticleLookUpTable[j] = 1;
	}
      else
	{
	  this->SignLookUpTable[j] = 1.0;
	  this->NbrParticleLookUpTable[j] = 0;
	}
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// return value = Hilbert space dimension

long FermionOnTorusWithMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum)
{
  FactorialCoefficient Dimension; 
  Dimension.PartialFactorialMultiply(maxMomentum - nbrFermions + 1, maxMomentum); 
  Dimension.FactorialDivide(nbrFermions);
  // approximate factor for reduction by symmetries
  if (Dimension.GetIntegerValue()>100000)
    Dimension/=maxMomentum-(maxMomentum>3?3:(maxMomentum>2?2:1));
  return (long)(Dimension.GetNumericalValue());
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnTorusWithMagneticTranslations::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((kxSector == this->XMomentum) && (kySector == this->YMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementaryKxMomentum = (this->XMomentum - kxSector) % this->MomentumModulo;
  int ComplementaryKyMomentum = (this->YMomentum - kySector) % this->MaxMomentum;
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->MomentumModulo;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->MaxMomentum;
  cout << "kx = " << this->XMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
  cout << "ky = " << this->YMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
  FermionOnTorusWithMagneticTranslations SubsytemSpace (nbrParticleSector, this->MaxMomentum, kxSector, kySector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnTorusWithMagneticTranslations ComplementarySpace (ComplementaryNbrParticles, this->MaxMomentum, ComplementaryKxMomentum, ComplementaryKyMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHETorusParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
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
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnTorusWithMagneticTranslations::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithMagneticTranslations* complementaryHilbertSpace,  ParticleOnTorusWithMagneticTranslations* destinationHilbertSpace,
												ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
  FermionOnTorusWithMagneticTranslations* TmpHilbertSpace =  (FermionOnTorusWithMagneticTranslations*) complementaryHilbertSpace;
  FermionOnTorusWithMagneticTranslations* TmpDestinationHilbertSpace =  (FermionOnTorusWithMagneticTranslations*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  BinomialCoefficients TmpBinomial (this->NbrFermions);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
 	      int TmpLzMax = this->MaxMomentum;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  double Coefficient = TmpInvBinomial;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = TmpDestinationHilbertSpace->MaxMomentum;
		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
		    {
		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			--Pos2;
		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState3 ^= TmpState3 >> 32;
#endif	
		      TmpState3 ^= TmpState3 >> 16;
		      TmpState3 ^= TmpState3 >> 8;
		      TmpState3 ^= TmpState3 >> 4;
		      TmpState3 ^= TmpState3 >> 2;
		      TmpState3 ^= TmpState3 >> 1;
		      Sign ^= TmpState3;
		      TmpState2 &= ~(0x1ul << Pos2);
		      --Pos2;
		    }
 		  if ((Sign & 0x1ul) == 0x0ul)		  
 		    Coefficient *= 1.0;
 		  else
 		    Coefficient *= -1.0;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = Coefficient;
		  ++Pos;
		}
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
}

// convert a state defined in the Ky basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector FermionOnTorusWithMagneticTranslations::ConvertToKxKyBasis(ComplexVector& state, ParticleOnTorus* space)  
{
  FermionOnTorus* TmpSpace = (FermionOnTorus*) space;
  ComplexVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      int TmpMaxMomentum = this->StateMaxMomentum[i];
      int Pos = TmpSpace->FindStateIndex(TmpState, TmpMaxMomentum);
      if (Pos < TmpSpace->HilbertSpaceDimension)
	{
	  TmpVector[i] =  state[Pos] * sqrt((double) this->NbrStateInOrbit[i]);
	}
    }
  return TmpVector;
}

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector FermionOnTorusWithMagneticTranslations::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnTorus* space)
{
  FermionOnTorus* TmpSpace = (FermionOnTorus*) space;
  ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  Complex* FourrierCoefficients = new Complex [this->MomentumModulo];
  for (int i = 0; i < this->MomentumModulo; ++i)
    FourrierCoefficients[i] = Phase (-2.0 * M_PI * ((double) (i * this->XMomentum)) / ((double) this->MomentumModulo));
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = TmpSpace->StateDescription[i];
      int NbrTranslation = 0;
      int TmpMaxMomentum = TmpSpace->StateKyMax[i];
      TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(TmpState, TmpMaxMomentum, NbrTranslation);
      if (NbrTranslation >= 0)
	{
	  int Pos = this->FindStateIndex(TmpState, TmpMaxMomentum);
	  if (Pos < this->HilbertSpaceDimension)
	    {
	      TmpVector[i] =  (state[Pos] * (1.0 - (2.0 * ((double) ((this->ReorderingSign[Pos] >> NbrTranslation) & 0x1ul)))) * 
			       FourrierCoefficients[NbrTranslation] / sqrt((double) this->NbrStateInOrbit[Pos]));
	    }
	}
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis, and change its kx quantum number
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// oldKx = vnew value of the relative quantum number kx
// return value = state in the (Kx,Ky) basis

ComplexVector FermionOnTorusWithMagneticTranslations::ConvertFromKxKyBasisAndModifyKx(ComplexVector& state, ParticleOnTorus* space, int oldKx)
{
  FermionOnTorus* TmpSpace = (FermionOnTorus*) space;
  FermionOnTorusWithMagneticTranslations* PreviousSpace = new FermionOnTorusWithMagneticTranslations (this->NbrFermions, this->MaxMomentum, oldKx, this->YMomentum);
  ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  Complex* FourrierCoefficients = new Complex [this->MomentumModulo];
  for (int i = 0; i < this->MomentumModulo; ++i)
    FourrierCoefficients[i] = Phase (-2.0 * M_PI * ((double) (i * this->XMomentum)) / ((double) this->MomentumModulo));
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = TmpSpace->StateDescription[i];
      int NbrTranslation = 0;
      int TmpMaxMomentum = TmpSpace->StateKyMax[i];
      TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(TmpState, TmpMaxMomentum, NbrTranslation);
      if (NbrTranslation >= 0)
	{
	  int TmpNbrTranslation = 0;
	  TmpState = PreviousSpace->FindCanonicalFormAndTestXMomentumConstraint(TmpState, TmpMaxMomentum, TmpNbrTranslation);
	  if (TmpNbrTranslation >= 0)
	  {
	    int Pos = this->FindStateIndex(TmpState, TmpMaxMomentum);
	    if (Pos < this->HilbertSpaceDimension)
	      {
		TmpVector[i] =  (state[Pos] * (1.0 - (2.0 * ((double) ((this->ReorderingSign[Pos] >> NbrTranslation) & 0x1ul)))) * 
			       FourrierCoefficients[NbrTranslation] / sqrt((double) this->NbrStateInOrbit[Pos]));
	      }
	  }
	}
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}

// core part of the C4 rotation
//
// inputState = reference on the state that has to be rotated
// inputSpace = Hilbert space associated to the input state
// outputState = reference on the rotated state
// minIndex = minimum index that has to be computed
// nbrIndices = number of indices that have to be computed
// clockwise = the rotation is done clockwise
// return value = reference on the rotated state

ComplexVector& FermionOnTorusWithMagneticTranslations::CoreC4Rotation (ComplexVector& inputState, ParticleOnTorusWithMagneticTranslations* inputSpace, ComplexVector& outputState, int minIndex, int nbrIndices, 
					       bool clockwise)
{
  FermionOnTorusWithMagneticTranslations* TmpInputSpace = (FermionOnTorusWithMagneticTranslations*) inputSpace;
  unsigned long* TmpInputMonomial = new unsigned long [this->NbrFermions];
  unsigned long* TmpInputMonomial2 = new unsigned long [this->NbrFermions];
  unsigned long* TmpOutputMonomial = new unsigned long [this->NbrFermions];
  int LastIndex = minIndex + nbrIndices;
  Complex Tmp = 0.0;
  Complex Tmp2 = 0.0;
  double TmpCoefficient = pow((double) this->MaxMomentum, -0.5 * ((double) this->NbrFermions));
  double PhaseFactor = 2.0 * M_PI / ((double) this->MaxMomentum);
  if (clockwise == true)
    PhaseFactor *= -1.0;
#ifdef __LAPACK__
  int* Permutation = new int[this->NbrFermions];
  doublecomplex* DeterminantMatrix = new doublecomplex [this->NbrFermions * this->NbrFermions];
#else
  ComplexMatrix DeterminantMatrix (this->NbrFermions, this->NbrFermions);
#endif
  ComplexMatrix PhaseMatrix (this->MaxMomentum, this->MaxMomentum);
  for (int k = 0; k < this->MaxMomentum; ++k)
    for (int l = 0; l < this->MaxMomentum; ++l)
      PhaseMatrix[k][l] = Phase(PhaseFactor * ((double) (k * l)));
 
  for (int i = minIndex ; i < LastIndex; ++i)
    {
      this->ConvertToMonomial(this->StateDescription[i], TmpOutputMonomial);
      Tmp = 0.0;
      for (int j = 0; j < TmpInputSpace->HilbertSpaceDimension; ++j)
	{
	  TmpInputSpace->ConvertToMonomial(TmpInputSpace->StateDescription[j], TmpInputMonomial);
	  unsigned long TmpPhase = 0ul;
#ifdef __LAPACK__
	  for (int k = 0; k < this->NbrFermions; ++k)
	    {
	      ComplexVector& TmpColumn = PhaseMatrix[TmpInputMonomial[k]];
	      for (int l = 0; l < this->NbrFermions; ++l)
		{
		  Complex& Tmp = TmpColumn[(int) TmpOutputMonomial[l]];
		  DeterminantMatrix[k + l * this->NbrFermions].r = Tmp.Re;
		  DeterminantMatrix[k + l * this->NbrFermions].i = Tmp.Im;
		}
	    }
	  int Information = 0;
	  int TmpDimension = this->NbrFermions;
	  FORTRAN_NAME(zgetrf)(&TmpDimension, &TmpDimension, DeterminantMatrix, &TmpDimension, Permutation, &Information);
	  int Sign = 0;
	  Complex Determinant (1.0,0.0);
	  for (int k = 0; k < TmpDimension; ++k)
	    {
	      if (Permutation[k] != k + 1)
		Sign ^= 1;
	      Determinant *= Complex(DeterminantMatrix[k + TmpDimension * k].r, DeterminantMatrix[k + TmpDimension * k].i);
	    }
	  if (Sign & 1)
	    Determinant *= -1.0;
	  Tmp += inputState[j] * TmpCoefficient * Determinant * sqrt((double) TmpInputSpace->NbrStateInOrbit[j]);
#else
	  for (int k = 0; k < this->NbrFermions; ++k)
	    for (int l = 0; l < this->NbrFermions; ++l)
	      DeterminantMatrix[k][l] = PhaseMatrix[TmpInputMonomial[k]][(long)TmpOutputMonomial[l]];
	  Tmp += inputState[j] * TmpCoefficient * DeterminantMatrix.Determinant() * sqrt((double) TmpInputSpace->NbrStateInOrbit[j]);
#endif
	}
      outputState[i] = Tmp * sqrt((double) (this->NbrStateInOrbit[i]));      
    }
#ifdef __LAPACK__
  delete[] DeterminantMatrix;
  delete[] Permutation;
#endif
  delete[] TmpInputMonomial;
  delete[] TmpInputMonomial2;
  delete[] TmpOutputMonomial;
  return outputState;
}

// request whether state with given index satisfies a general Pauli exclusion principle
//
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals
// return value = true if teh state satisfies the general Pauli exclusion principle

bool FermionOnTorusWithMagneticTranslations::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  unsigned long TmpState = this->StateDescription[index];
  unsigned long RMask = (0x1ul << pauliR) - 0x1ul;
  TmpState |= (TmpState & RMask) << this->MaxMomentum;
  unsigned long UnsignedPauliK = (unsigned long) pauliK;
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      unsigned long TmpOccupation = 0l;
      unsigned long TmpState2 = TmpState & RMask;
      while (TmpState2 != 0x0ul)
	{
	  TmpOccupation += (TmpState2 & 0x1ul);
	  TmpState2 >>= 1;
	}
      if (TmpOccupation > UnsignedPauliK)
	return false;
      TmpState >>= 1;
    }
  return true;
}
