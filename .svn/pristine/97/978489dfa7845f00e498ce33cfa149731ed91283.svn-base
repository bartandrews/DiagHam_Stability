////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                      class author: Cecile Repellin                         //
//                                                                            //
//                                                                            //
//                 class of fermion with spin on a torus with time            //
//                           reversal symmetry, taking                        //
//                      into account magnetic translations                    //
//                                                                            //
//                        last modification : 20/02/2015                      //
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
#include "HilbertSpace/FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "Matrix/Matrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ArrayTools.h"

#include <cmath>
#include <cstdlib>


using std::cout;
using std::endl;



// basic constructor
// 
// nbrFermions= number of fermions
// totalSpin = twice the total spin value
// maxMomentum = momentum maximum value for a fermion
// xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)

FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations::FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations (int nbrFermions, int totalSpin, int maxMomentum, 
											      int xMomentum, int yMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;  
  this->MaxMomentum = maxMomentum;  

  this->NbrMomentum = this->MaxMomentum + 1;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  if (this->TotalSpin == 0)
    this->MomentumModulo = this->MaxMomentum;
  else
    this->MomentumModulo = FindGCD(abs(this->TotalSpin), this->MaxMomentum);
  cout << "MomentumModulo=" << MomentumModulo<<endl;
  this->XMomentum = xMomentum % this->MomentumModulo;
  this->YMomentum = yMomentum % this->MaxMomentum;

  this->StateShift = 2 * this->MaxMomentum / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = 0x1ul;
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= 0x1ul;
    }

  
  this->MaximumSignLookUp = 16;
  this->GenerateSignLookUpTable();
  this->HilbertSpaceDimension = ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum - 1, 0, this->NbrFermionsUp);
  cout << "intermediate Hilbert space dimension = " << HilbertSpaceDimension << endl;
  this->HilbertSpaceDimension = this->GenerateStates();
  
  
  cout << "Hilbert space dimension = "<< HilbertSpaceDimension << endl;
  this->Flag.Initialize();

  if (this->HilbertSpaceDimension !=0)
    this->GenerateLookUpTable(1000000);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
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

FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations::FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations(const FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations& fermions)
{
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->IncNbrFermions = fermions.IncNbrFermions;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations::~FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations& FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations::operator = (const FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;

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
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

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

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations::Clone()
{
  return new FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations(*this);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalMomentum = momentum total value
// totalSpinUp = number of particles with spin up
// return value = Hilbert space dimension

long FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalMomentum, int totalSpinUp)
{
  if ((nbrFermions < 0) || (totalSpinUp < 0) || (totalSpinUp > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpinUp) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpinUp)) )
    return 0l;

  if (nbrFermions == 1)
    {
      long Tmp = 0;
      if (totalSpinUp == 1)
      {
	for (int k = 0; k <= lzMax; ++k)
	  if (((k+totalMomentum) % MaxMomentum) == YMomentum)
	  {
	    ++Tmp;
	  }
      }
      else
      {
	for (int k = 0; k <= lzMax; ++k)
	  if (((-k+totalMomentum + MaxMomentum) % MaxMomentum) == YMomentum)
	  {
	    ++Tmp;
	  }
      }
      return Tmp;
    }
  

  if ((lzMax == 0)  && ( (totalMomentum % MaxMomentum) != YMomentum))
    return 0l;

  unsigned long Tmp = 0l;  
  if (nbrFermions > 2)
    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalMomentum, totalSpinUp - 1);
  else
    if ((totalSpinUp == 1) && (((totalMomentum) % MaxMomentum) == YMomentum) )
      ++Tmp;
    return  (Tmp + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalMomentum + lzMax, totalSpinUp - 1)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, ((totalMomentum - lzMax + MaxMomentum) % MaxMomentum), totalSpinUp)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalMomentum, totalSpinUp));
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalMomentum = momentum total value
// totalSpinUp = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = Hilbert space dimension

long FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations::RawGenerateStates(int nbrFermions, int lzMax, int totalMomentum, int totalSpinUp, long pos)
{  
  if ((nbrFermions < 0) || (totalSpinUp < 0) || (totalSpinUp > nbrFermions))
    return pos;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpinUp) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpinUp)) )
    return pos;
  if ((nbrFermions == 0) && (lzMax == 0) && (totalSpinUp == 0) && ((totalMomentum % MaxMomentum)== YMomentum))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  
  if (nbrFermions == 1)
    {
      for (int k = 0; k <= lzMax; ++k)
      {
	if (totalSpinUp == 1)
	{
	  if (((k+totalMomentum + MaxMomentum) % MaxMomentum) == YMomentum)
	    {	    
	      this->StateDescription[pos++] = 0x1ul << ((k << 1) + totalSpinUp);
	    }
	}
	else
	{
	  if (((-k+totalMomentum) % MaxMomentum) == YMomentum)
	    {	    
	      this->StateDescription[pos++] = 0x1ul << (k << 1);
	    }
	}
      }
      return (pos);
    }
  
  if ((lzMax == 0)  && ( (totalMomentum % MaxMomentum) != YMomentum))
    return pos;

  // put last two particles:
  if ((lzMax==0) && (nbrFermions == 2) && (totalSpinUp == 1))
    {
      this->StateDescription[pos] = 0x3ul;
      return (pos + 1l);
    }

  // enter recursion, here:
  // put two particles
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, lzMax - 1, totalMomentum, totalSpinUp - 1,  pos);
  unsigned long Mask = 0x3ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  // put one particle with spin up
  TmpPos = this->RawGenerateStates(nbrFermions - 1, lzMax - 1, totalMomentum + lzMax, totalSpinUp - 1,  pos);
  Mask = 0x2ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  // put one particle with spin down
  TmpPos = this->RawGenerateStates(nbrFermions - 1, lzMax - 1, totalMomentum - lzMax, totalSpinUp,  pos);
  Mask = 0x1ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  // do not put any particles
  return this->RawGenerateStates(nbrFermions, lzMax - 1, totalMomentum, totalSpinUp, pos);  
  
}