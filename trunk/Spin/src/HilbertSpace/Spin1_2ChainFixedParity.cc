////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of spin 1/2 chain without any Sz contraint           //
//                                                                            //
//                        last modification : 11/12/2013                      //
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


#include "HilbertSpace/Spin1_2ChainFixedParity.h"
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


#define M_SQRT3 1.73205080756888
#define M_SQRT6 2.44948974278318


// default constructor
//

Spin1_2ChainFixedParity::Spin1_2ChainFixedParity ()
{
  this->Flag.Initialize();
  this->ChainLength = 0;
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->Sz = 0;
  this->StateDescription = 0;
  this->LookUpTable = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
  this->SzParity = 0;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1/2
// parity = parity of the total (Sz + 1/2) (can be 0 or 1)

Spin1_2ChainFixedParity::Spin1_2ChainFixedParity (int chainLength, int parity) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  if (this->ChainLength  > 63)
    {
      this->ChainLength = 1;
    }
  this->FixedQuantumNumberFlag = false;
  this->SzParity = (unsigned long) parity;

  this->LargeHilbertSpaceDimension = 1l << (this->ChainLength - 1);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];

  unsigned long TmpMax = 1ul << this->ChainLength;
  for (unsigned long i = 0l ; i < TmpMax; ++i)
    {
      unsigned long TmpState = i;
#ifdef __64_BITS__
      TmpState ^= TmpState >> 32;
#endif
      TmpState ^= TmpState >> 16;
      TmpState ^= TmpState >> 8;
      TmpState ^= TmpState >> 4;
      TmpState ^= TmpState >> 2;
      TmpState ^= TmpState >> 1;
      if (((int) (TmpState & 0x1ul)) == this->SzParity)
	this->StateDescription[i >> 1] = i;
    }
  this->LookUpTable = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainFixedParity::Spin1_2ChainFixedParity (const Spin1_2ChainFixedParity& chain)
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
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
      this->SzParity = chain.SzParity;
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->LookUpTableSize = 0;
      this->SzParity = 0;
    }
}

// destructor
//

Spin1_2ChainFixedParity::~Spin1_2ChainFixedParity () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainFixedParity& Spin1_2ChainFixedParity::operator = (const Spin1_2ChainFixedParity& chain)
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
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainFixedParity::Clone()
{
  return new Spin1_2ChainFixedParity (*this);
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainFixedParity::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex((State | (0x1ul << j)) & ~(0x1ul << i));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x0ul)
    {
      coefficient = -0.25;
      return state;
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainFixedParity::SpiSpj (int i, int j, int state, double& coefficient)
{  
  if (i != j)
    {
      unsigned long TmpState = this->StateDescription[state];
      unsigned long Mask = ((0x1ul << i) | (0x1ul << j));
      if ((TmpState & Mask) == 0x0ul)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex(TmpState | Mask);
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainFixedParity::SmiSmj (int i, int j, int state, double& coefficient)
{  
  if (i != j)
    {
      unsigned long TmpState = this->StateDescription[state];
      unsigned long Mask = ((0x1ul << i) | (0x1ul << j));
      if ((TmpState & Mask) == Mask)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex(TmpState & ~Mask);
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  return this->HilbertSpaceDimension;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix Spin1_2ChainFixedParity::EvaluatePartialDensityMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture)
{ 
  if (nbrSites == 0)
    {
      if (szSector == 0)
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
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
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
  Spin1_2ChainFixedParity TmpDestinationHilbertSpace(nbrSites, szSector);
  int ComplementaryQSector = this->Sz - szSector;
  if (ComplementaryQSector < 0)
    ComplementaryQSector += 2;
   if (ComplementaryQSector >= 2)
    ComplementaryQSector -= 2;
  Spin1_2ChainFixedParity TmpHilbertSpace(this->ChainLength - nbrSites, ComplementaryQSector);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  long TmpNbrNonZeroElements = 0l;
  unsigned long Mask = ((0x1ul << shift) - 0x1ul);
  unsigned long ComplementaryMask = ((0x1ul << this->ChainLength) - 0x1ul) ^ Mask;
  int ComplementaryShift = this->ChainLength - shift;
  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
	  TmpState2 = ((TmpState2 << shift) & ComplementaryMask) | ((TmpState2 >> ComplementaryShift) & Mask);
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]];
	      for (int k = 0; k < Pos; ++k)
		{
		  if (TmpStatePosition2[k] >= Pos2)
		    {
		      TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
		    }
		}
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  return TmpDensityMatrix;
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix Spin1_2ChainFixedParity::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
      
    }
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}      
    }
  Spin1_2ChainFixedParity TmpDestinationHilbertSpace(nbrSites, szSector);
  int ComplementaryQSector = this->Sz - szSector;
  if (ComplementaryQSector < 0)
    ComplementaryQSector += 2;
  if (ComplementaryQSector >= 2)
    ComplementaryQSector -= 2;
  Spin1_2ChainFixedParity TmpHilbertSpace(this->ChainLength - nbrSites, ComplementaryQSector);
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
  int Shift = nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  long TmpNbrNonZeroElements = 0l;
  unsigned long Mask = ((0x1ul << shift) - 0x1ul);
  unsigned long ComplementaryMask = ((0x1ul << this->ChainLength) - 0x1ul) ^ Mask;
  int ComplementaryShift = this->ChainLength - shift;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
	  TmpState2 = ((TmpState2 << shift) & ComplementaryMask) | ((TmpState2 >> ComplementaryShift) & Mask);
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpEntanglementMatrix.AddToMatrixElement(j, MinIndex, groundState[TmpPos]);
	    }
	}
    }
  return TmpEntanglementMatrix;
}
