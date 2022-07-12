////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a fixed Sz value                //
//         and a fixed parity under the parity symmetry (aka Sz<->-Sz)        //
//                                                                            //
//                        last modification : 21/07/2015                      //
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


#include "HilbertSpace/Spin1_2ChainParitySymmetry.h"
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
using std::endl;
using std::hex;
using std::dec;



// default constructor
//

Spin1_2ChainParitySymmetry::Spin1_2ChainParitySymmetry ()
{
  this->Parity = 0;
  this->ParitySign = 1.0;
  this->ParityMask = 0x0ul;
}

// constructor for complete Hilbert space with a given total spin projection Sz
//
// chainLength = number of spin 1/2
// parity = parity of the total (Sz + 1/2)
// memorySize = memory size in bytes allowed for look-up table

Spin1_2ChainParitySymmetry::Spin1_2ChainParitySymmetry (int chainLength, int parity, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
#ifdef __64_BITS__
  if (this->ChainLength  > 64)
#else
  if (this->ChainLength  > 32)
#endif
    {
      this->ChainLength = 1;
    }
  this->FixedQuantumNumberFlag = true;
  this->Sz = 0;
  this->Parity = parity;
  if (this->Parity == 0)
    this->ParitySign = 1.0;
  else
    this->ParitySign = -1.0;
  this->ParityMask = (0x1ul << this->ChainLength) - 1ul;;
  unsigned long TmpLargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);
  this->StateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = this->GenerateStates ((this->Sz + this->ChainLength) >> 1, this->ChainLength - 1, 0l);
  this->LargeHilbertSpaceDimension = 0ul;
  if (this->Parity == 0)
    {
      for (long i = 0l; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->ApplyParitySymmetry(this->StateDescription[i]);
	  if (TmpState >= this->StateDescription[i])
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	}
      unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->LargeHilbertSpaceDimension = 0ul;
      for (long i = 0l; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->ApplyParitySymmetry(this->StateDescription[i]);
	  if (TmpState >= this->StateDescription[i])
	    {
	      if (TmpState != this->StateDescription[i])
		{
		  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->StateDescription[i] | SPIN1_2CHAIN_PARITYSYMMETRIC_BIT;
		}
	      else
		{
		  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->StateDescription[i];
		}
	      ++this->LargeHilbertSpaceDimension;	  
	    }
	}
      delete[] this->StateDescription;
      this->StateDescription = TmpStateDescription;
    }
  else
    {
      for (long i = 0l; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->ApplyParitySymmetry(this->StateDescription[i]);
	  if (TmpState > this->StateDescription[i])
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	}
      if (this->LargeHilbertSpaceDimension > 0)
	{
	  unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
	  this->LargeHilbertSpaceDimension = 0ul;
	  for (long i = 0l; i < TmpLargeHilbertSpaceDimension; ++i)
	    {
	      unsigned long TmpState = this->ApplyParitySymmetry(this->StateDescription[i]);
	      if (TmpState > this->StateDescription[i])
		{
		  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->StateDescription[i] | SPIN1_2CHAIN_PARITYSYMMETRIC_BIT;
		  ++this->LargeHilbertSpaceDimension;	  
		}
	    }
	  delete[] this->StateDescription;
	  this->StateDescription = TmpStateDescription;
	}
      else
	{
	  this->StateDescription = 0;
	}
    }
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->MaxStateDescription = this->StateDescription[0l] & SPIN1_2CHAIN_FULLSYMMETRY_MASK;
      this->MinStateDescription = this->StateDescription[this->LargeHilbertSpaceDimension - 1l] & SPIN1_2CHAIN_FULLSYMMETRY_MASK;
      this->GenerateLookUpTable(memorySize, SPIN1_2CHAIN_FULLSYMMETRY_MASK);
    }
  this->Flag.Initialize();
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainParitySymmetry::Spin1_2ChainParitySymmetry (const Spin1_2ChainParitySymmetry& chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
      this->Sz = chain.Sz;
      this->StateDescription = chain.StateDescription;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
      this->Parity = chain.Parity;
      this->ParitySign = chain.ParitySign;
      this->ParityMask = chain.ParityMask;
      this->MaxStateDescription = chain.MaxStateDescription;
      this->MinStateDescription = chain.MinStateDescription;
      this->Flag = chain.Flag;
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->LookUpTableMemorySize = 0;
      this->MaximumLookUpShift = 0;
      this->Parity = 0;
      this->ParitySign = 1.0;
      this->ParityMask = 0x0ul;
      this->MaxStateDescription = 0x0ul;
      this->MinStateDescription = 0x0ul;
    }
}

// destructor
//

Spin1_2ChainParitySymmetry::~Spin1_2ChainParitySymmetry () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainParitySymmetry& Spin1_2ChainParitySymmetry::operator = (const Spin1_2ChainParitySymmetry& chain)
{
  if ((this->Flag.Used() == true) && (this->ChainLength != 0) && (this->Flag.Shared() == false))
    {
      delete[] this->StateDescription;
    }
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
      this->Sz = chain.Sz;
      this->StateDescription = chain.StateDescription;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
      this->Parity = chain.Parity;
      this->ParitySign = chain.ParitySign;
      this->ParityMask = chain.ParityMask;
      this->MaxStateDescription = chain.MaxStateDescription;
      this->MinStateDescription = chain.MinStateDescription;
      this->Flag = chain.Flag;
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->LookUpTableMemorySize = 0;
      this->MaximumLookUpShift = 0;
      this->Parity = 0;
      this->ParitySign = 1.0;
      this->ParityMask = 0x0ul;
      this->MaxStateDescription = 0x0ul;
      this->MinStateDescription = 0x0ul;
   }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainParitySymmetry::Clone()
{
  return new Spin1_2ChainParitySymmetry (*this);
}



