////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of spin 1/2 chain without any Sz contraint            //
//                  and 2d translations plus inversion symmetry               //
//                                                                            //
//                        last modification : 20/12/2015                      //
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


#include "HilbertSpace/Spin1_2ChainFullInversionAnd2DTranslation.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"

#include <math.h>
#include <iostream>


using std::cout;
using std::dec;
using std::hex;
using std::endl;



// default constructor
//

Spin1_2ChainFullInversionAnd2DTranslation::Spin1_2ChainFullInversionAnd2DTranslation ()
{
  this->InversionSector = 1.0;
  this->InversionShift = 0;
  this->InversionUnshift = 0;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// inversionSector = inversion sector (can be either +1 or -1)
// xMomentum = momentum along the x direction
// maxXMomentum = number of sites in the x direction
// yMomentum = momentum along the y direction
// maxYMomentum = number of sites in the y direction
// memory = amount of memory granted for precalculations

Spin1_2ChainFullInversionAnd2DTranslation::Spin1_2ChainFullInversionAnd2DTranslation (int inversionSector, int xMomentum, int maxXMomentum, 
										      int yMomentum, int maxYMomentum, unsigned long memory) 
{
  this->Flag.Initialize();
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->NbrSite = this->MaxXMomentum * this->MaxYMomentum;
  this->ChainLength = this->NbrSite;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  if (inversionSector == 1)
    {
      this->InversionSector = 1.0;
    }
  else
    {
      this->InversionSector = -1.0;
    }
  this->Sz = 0;
  this->FixedQuantumNumberFlag = false;
  
  this->StateXShift = this->NbrSite / this->MaxXMomentum;
  this->ComplementaryStateXShift = this->NbrSite - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->NbrYMomentumBlocks = this->NbrSite / this->StateXShift;
  this->StateYShift = (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 

#ifdef __64_BITS__
  this->InversionShift = 32 - ((this->YMomentumBlockSize - 1) >> 1) -1;
#else
  this->InversionUnshift = 16 - ((this->YMomentumBlockSize - 1) >> 1) -1;
#endif
  if ((this->YMomentumBlockSize & 1) == 0)
    this->InversionUnshift = this->InversionShift - 1;
  else
    this->InversionUnshift = this->InversionShift;
  
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0)
    {
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory +=  this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->NbrSite * sizeof(int);
      UsedMemory += this->NbrSite * this->LookUpTableMemorySize * sizeof(int);
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

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainFullInversionAnd2DTranslation::Spin1_2ChainFullInversionAnd2DTranslation (const Spin1_2ChainFullInversionAnd2DTranslation& chain)
{
  this->Flag = chain.Flag;
  if (chain.NbrSite != 0)
    {
      this->NbrSite = chain.NbrSite;
      this->ChainLength = chain.NbrSite;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->XMomentum = chain.XMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;
      this->MaxYMomentum = chain.MaxYMomentum;
      this->YMomentum = chain.YMomentum;
      this->NbrYMomentumBlocks = chain.NbrYMomentumBlocks;
      this->StateYShift = chain.StateYShift;
      this->YMomentumBlockSize = chain.YMomentumBlockSize;
      this->ComplementaryStateYShift = chain.ComplementaryStateYShift;
      this->YMomentumMask = chain.YMomentumMask;
      this->YMomentumBlockMask = chain.YMomentumBlockMask;  
      this->YMomentumFullMask = chain.YMomentumFullMask;
      this->ComplementaryYMomentumFullMask = chain.ComplementaryYMomentumFullMask; 
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;

      this->Sz = chain.Sz;
      this->InversionSector = chain.InversionSector;

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
    }
  else
    {
      this->NbrSite = 0;
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->LookUpTableSize = 0;
      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;
    }
}

// destructor
//

Spin1_2ChainFullInversionAnd2DTranslation::~Spin1_2ChainFullInversionAnd2DTranslation () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainFullInversionAnd2DTranslation& Spin1_2ChainFullInversionAnd2DTranslation::operator = (const Spin1_2ChainFullInversionAnd2DTranslation& chain)
{
  if ((this->Flag.Used() == true) && (this->NbrSite != 0) && (this->Flag.Shared() == false))
    {
      delete[] this->StateDescription;
    }
  this->Flag = chain.Flag;
  if (chain.NbrSite != 0)
    {
      this->NbrSite = chain.NbrSite;
      this->ChainLength = chain.NbrSite;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->XMomentum = chain.XMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;
      this->MaxYMomentum = chain.MaxYMomentum;
      this->YMomentum = chain.YMomentum;
      this->NbrYMomentumBlocks = chain.NbrYMomentumBlocks;
      this->StateYShift = chain.StateYShift;
      this->YMomentumBlockSize = chain.YMomentumBlockSize;
      this->ComplementaryStateYShift = chain.ComplementaryStateYShift;
      this->YMomentumMask = chain.YMomentumMask;
      this->YMomentumBlockMask = chain.YMomentumBlockMask;  
      this->YMomentumFullMask = chain.YMomentumFullMask;
      this->ComplementaryYMomentumFullMask = chain.ComplementaryYMomentumFullMask; 
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;

      this->Sz = chain.Sz;
      this->InversionSector = chain.InversionSector;

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
    }
  else
    {
      this->NbrSite = 0;
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;
   }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainFullInversionAnd2DTranslation::Clone()
{
  return new Spin1_2ChainFullInversionAnd2DTranslation (*this);
}

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long Spin1_2ChainFullInversionAnd2DTranslation::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates();
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  double TmpSign;
  int InversionFlip;
#ifdef  __64_BITS__
  unsigned long Discard = 0xfffffffffffffffful;
#else
  unsigned long Discard = 0xfffffffful;
#endif
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
// 	  unsigned long Tmp = this->StateDescription[i];
// 	  cout << "| ";
// 	  for (int j = 0; j < this->MaxXMomentum; ++j)
// 	    {
// 	      for (int k = 0; k < this->MaxYMomentum; ++k)
// 		{
// 		  if (((Tmp >> (k + (j * this->MaxYMomentum))) & 0x1ul) == 0x0ul)
// 		    cout << "+ ";
// 		  else
// 		    cout << "- ";
// 		}
// 	      cout << "| ";
// 	    }
// 	  cout << " -> ";
// 	  cout << "| ";
// 	  Tmp = this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY, TmpSign);
// 	  for (int j = 0; j < this->MaxXMomentum; ++j)
// 	    {
// 	      for (int k = 0; k < this->MaxYMomentum; ++k)
// 		{
// 		  if (((Tmp >> (k + (j * this->MaxYMomentum))) & 0x1ul) == 0x0ul)
// 		    cout << "+ ";
// 		  else
// 		    cout << "- ";
// 		}
// 	      cout << "| ";
// 	    }
// 	  cout << endl;
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY, TmpSign) == this->StateDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
//	      cout << "check" << endl;
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = Discard;
	    }
	}
      else
	{
	  this->StateDescription[i] = Discard;
	}
    }
  if (TmpLargeHilbertSpaceDimension > 0)
    {
      unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
      this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
      TmpLargeHilbertSpaceDimension = 0l;
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (this->StateDescription[i] != Discard)
	    {
	      TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	      this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);

// 	      cout << "| ";
// 	      for (int j = 0; j < this->MaxXMomentum; ++j)
// 		{
// 		  for (int k = 0; k < this->MaxYMomentum; ++k)
// 		    {
// 		      if (((TmpStateDescription[TmpLargeHilbertSpaceDimension] >> (k + (j * this->MaxYMomentum))) & 0x1ul) == 0x0ul)
// 			cout << "+ ";
// 		      else
// 			cout << "- ";
// 		    }
// 		  cout << "| ";
// 		}
// 	      cout << " orb=" << this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] << endl;
	      ++TmpLargeHilbertSpaceDimension;
	    }
	}
      delete[] this->StateDescription;
      this->StateDescription = TmpStateDescription;
    }
  else
    {
      delete[] this->StateDescription;
    }
  return TmpLargeHilbertSpaceDimension;
}

// convert a state defined in the (Kx,Ky) basis into a state in the real space basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector Spin1_2ChainFullInversionAnd2DTranslation::ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  Spin1_2ChainFull* TmpSpace = (Spin1_2ChainFull*) space;
  ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  Complex** FourrierCoefficients = new Complex* [this->MaxXMomentum];
  
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficients[i][j] = Phase (-2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
	}
    }
  
  int NbrTranslationX;
  int NbrTranslationY;
  double TmpSign;
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = TmpSpace->StateDescription[i];
      TmpState = this->FindCanonicalForm(TmpState, NbrTranslationX, NbrTranslationY, TmpSign);
      NbrTranslationX = (this->MaxXMomentum - NbrTranslationX) % this->MaxXMomentum;
      NbrTranslationY = (this->MaxYMomentum - NbrTranslationY) % this->MaxYMomentum;
      int TmpMaxMomentum = this->NbrSite;
      while (((TmpState >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
	--TmpMaxMomentum;
      
      int Pos = this->FindStateIndex(TmpState, TmpMaxMomentum);
      if (Pos < this->HilbertSpaceDimension)
	{
	  TmpVector[i] =  (state[Pos] * FourrierCoefficients[NbrTranslationX][NbrTranslationY]) * (TmpSign  / sqrt((double) this->NbrStateInOrbit[Pos]));
	}
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}
  
// compute the rescaling factors
//

void Spin1_2ChainFullInversionAnd2DTranslation::ComputeRescalingFactors()
{
  int Tmp = 2 * this->NbrSite;
  this->RescalingFactors = new double* [Tmp + 1];
  for (int i = 1; i <= Tmp; ++i)
    {
      this->RescalingFactors[i] = new double [Tmp + 1];
      for (int j = 1; j <= Tmp; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}

