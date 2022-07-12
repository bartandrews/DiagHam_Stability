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


#include "HilbertSpace/Spin1_2ChainFull.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"
#include "Matrix/SparseComplexMatrix.h"

#include <math.h>
#include <iostream>


using std::cout;
using std::endl;


#define M_SQRT3 1.73205080756888
#define M_SQRT6 2.44948974278318


// default constructor
//

Spin1_2ChainFull::Spin1_2ChainFull ()
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
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1/2

Spin1_2ChainFull::Spin1_2ChainFull (int chainLength) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  if (this->ChainLength  > 32)
    {
      this->ChainLength = 1;
    }
  this->FixedQuantumNumberFlag = false;

  this->LargeHilbertSpaceDimension = 1l << this->ChainLength;
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];

  for (unsigned long i = 0l ; i < this->LargeHilbertSpaceDimension; ++i)
    this->StateDescription[i] = i;
  this->LookUpTable = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainFull::Spin1_2ChainFull (const Spin1_2ChainFull& chain)
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
    }
}

// destructor
//

Spin1_2ChainFull::~Spin1_2ChainFull () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainFull& Spin1_2ChainFull::operator = (const Spin1_2ChainFull& chain)
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

AbstractHilbertSpace* Spin1_2ChainFull::Clone()
{
  return new Spin1_2ChainFull (*this);
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainFull::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = state;
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
	  return (int) ((State | (0x1ul << j)) & ~(0x1ul << i));
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

int Spin1_2ChainFull::SpiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = state;
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
      if (tmpState2 == 0x3ul)
	{
	  coefficient = 1.0;
	  return (int) (State & ~((0x1ul << j) | (0x1ul << i)));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
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

int Spin1_2ChainFull::SmiSmj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = state;
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
      if (tmpState2 == 0x0ul)
	{
	  coefficient = 1.0;
	  return (int) (State | ((0x1ul << j) | (0x1ul << i)));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  return this->HilbertSpaceDimension;
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix Spin1_2ChainFull::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      RealMatrix TmpEntanglementMatrix(1, 1);
      TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
      return TmpEntanglementMatrix;
    }
  if (nbrSites == this->ChainLength)
    {
      RealMatrix TmpEntanglementMatrix(1, 1);
      TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
      return TmpEntanglementMatrix;
    }
  Spin1_2ChainFull TmpDestinationHilbertSpace(nbrSites);
  Spin1_2ChainFull TmpHilbertSpace(this->ChainLength - nbrSites);

  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int Shift = nbrSites;
  unsigned long Mask = (0x1ul << Shift) - 0x1ul;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift) & 0xfffffffful;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask);
	  TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[(int) TmpState2]);
	}
    }

   return TmpEntanglementMatrix;
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix Spin1_2ChainFull::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      ComplexMatrix TmpEntanglementMatrix(1, 1);
      Complex Tmp(1.0, 0.0);
      TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
      return TmpEntanglementMatrix;
    }
  if (nbrSites == this->ChainLength)
    {
      ComplexMatrix TmpEntanglementMatrix(1, 1);
      Complex Tmp(1.0, 0.0);
      TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
      return TmpEntanglementMatrix;
    }
  Spin1_2ChainFull TmpDestinationHilbertSpace(nbrSites);
  Spin1_2ChainFull TmpHilbertSpace(this->ChainLength - nbrSites);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int Shift = nbrSites;
  unsigned long Mask = (0x1ul << Shift) - 0x1ul;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift) & 0xfffffffful;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask);
	  TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[(int) TmpState2]);
	}
    }
   return TmpEntanglementMatrix;
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
// 
// sites = list of sites that define the A subsystem 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix Spin1_2ChainFull::EvaluatePartialEntanglementMatrix (int* sites, int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      RealMatrix TmpEntanglementMatrix(1, 1);
      TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
      return TmpEntanglementMatrix;
    }
  if (nbrSites == this->ChainLength)
    {
      RealMatrix TmpEntanglementMatrix(1, 1);
      TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
      return TmpEntanglementMatrix;
    }
  Spin1_2ChainFull TmpDestinationHilbertSpace(nbrSites);
  Spin1_2ChainFull TmpHilbertSpace(this->ChainLength - nbrSites);

  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int NbrBSites = this->ChainLength - nbrSites;
  int* BSites = new int [NbrBSites];
  NbrBSites = 0;
  for (int i = 0; i < this->ChainLength; ++i)
    {
      int j = 0;
      while ((j < nbrSites) && (sites[j] != i))
	++j;
      if (j == nbrSites)
	{
	  BSites[NbrBSites] = i;
	  ++NbrBSites;
	}
    }

  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState3 = TmpHilbertSpace.StateDescription[MinIndex];
      unsigned long TmpState = 0x0ul;
      for (int i = 0; i < NbrBSites; ++i)
	TmpState |= ((TmpState3 >> i) & 0x1ul) << BSites[i];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState4 = TmpDestinationHilbertSpace.StateDescription[j];
	  unsigned long TmpState2 = 0x0ul;
	  for (int i = 0; i < nbrSites; ++i)
	    TmpState2 |= ((TmpState4 >> i) & 0x1ul) << sites[i];	  
	  TmpState2 |= TmpState;
	  TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[(int) TmpState2]);
	}
    }
  delete[] BSites;
  return TmpEntanglementMatrix;
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
// 
// sites = list of sites that define the A subsystem 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix Spin1_2ChainFull::EvaluatePartialEntanglementMatrix (int* sites, int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      ComplexMatrix TmpEntanglementMatrix(1, 1);
      Complex Tmp(1.0, 0.0);
      TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
      return TmpEntanglementMatrix;
    }
  if (nbrSites == this->ChainLength)
    {
      ComplexMatrix TmpEntanglementMatrix(1, 1);
      Complex Tmp(1.0, 0.0);
      TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
      return TmpEntanglementMatrix;
    }
  Spin1_2ChainFull TmpDestinationHilbertSpace(nbrSites);
  Spin1_2ChainFull TmpHilbertSpace(this->ChainLength - nbrSites);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int NbrBSites = this->ChainLength - nbrSites;
  int* BSites = new int [NbrBSites];
  NbrBSites = 0;
  for (int i = 0; i < this->ChainLength; ++i)
    {
      int j = 0;
      while ((j < nbrSites) && (sites[j] != i))
	++j;
      if (j == nbrSites)
	{
	  BSites[NbrBSites] = i;
	  ++NbrBSites;
	}
    }

  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState3 = TmpHilbertSpace.StateDescription[MinIndex];
      unsigned long TmpState = 0x0ul;
      for (int i = 0; i < NbrBSites; ++i)
	TmpState |= ((TmpState3 >> i) & 0x1ul) << BSites[i];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState4 = TmpDestinationHilbertSpace.StateDescription[j];
	  unsigned long TmpState2 = 0x0ul;
	  for (int i = 0; i < nbrSites; ++i)
	    TmpState2 |= ((TmpState4 >> i) & 0x1ul) << sites[i];	  
	  TmpState2 |= TmpState;
	  TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[(int) TmpState2]);
	}
    }
  return TmpEntanglementMatrix;
}

// create a state from its MPS description
//
// bMatrices = array that gives the B matrices 
// state = reference to vector that will contain the state description
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void Spin1_2ChainFull::CreateStateFromMPSDescription (SparseComplexMatrix* bMatrices, ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
						      long initialIndex, long nbrComponents)
{

  ComplexVector TmpVector1 (bMatrices->GetNbrColumn());
  ComplexVector TmpVector2 (bMatrices->GetNbrColumn());
  unsigned long MinIndex = (unsigned long) initialIndex;
  unsigned long MaxIndex = (unsigned long) this->LargeHilbertSpaceDimension;
  if (nbrComponents > 0)
    {
      MaxIndex = MinIndex + ((unsigned long) nbrComponents);
    }
  if (mPSRowIndex >= 0)
    {
      for (; MinIndex < MaxIndex; ++MinIndex)    
	{
	  state[(long) MinIndex] = 0.0;
	  TmpVector1.ClearVector();
	  TmpVector1[mPSColumnIndex] = 1.0;
	  for (int i =  this->ChainLength - 1 ; i >= 0; --i)
	    {
	      if ((MinIndex & (0x1ul << i)) == 0x0ul)
		{
		  bMatrices[0].RightMultiply(TmpVector1, TmpVector2);
		}
	      else
		{
		  bMatrices[1].RightMultiply(TmpVector1, TmpVector2);
		}
	      ComplexVector TmpVector3 = TmpVector2;
	      TmpVector2 = TmpVector1;
	      TmpVector1 = TmpVector3;
	    }
	  state[(long) MinIndex] = TmpVector1[mPSRowIndex];
	}
    }
  else
    {
      for (; MinIndex < MaxIndex; ++MinIndex)    
	{
	  state[(long) MinIndex] = 0.0;
	  for (int j = 0; j < bMatrices->GetNbrColumn(); ++j)
	    {
	      TmpVector1.ClearVector();
	      TmpVector1[j] = 1.0;
	      for (int i =  this->ChainLength - 1 ; i >= 0; --i)
		{
		  if ((MinIndex & (0x1ul << i)) == 0x0ul)
		    {
		      bMatrices[0].RightMultiply(TmpVector1, TmpVector2);
		    }
		  else
		    {
		      bMatrices[1].RightMultiply(TmpVector1, TmpVector2);
		    }
		  ComplexVector TmpVector3 = TmpVector2;
		  TmpVector2 = TmpVector1;
		  TmpVector1 = TmpVector3;
		}
	      state[(long) MinIndex] += TmpVector1[j];
	    }
	}
    }
}

// create a state from its MPS description
//
// bMatrices = array that gives the B matrices 
// state = reference to vector that will contain the state description
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void Spin1_2ChainFull::CreateStateFromMPSDescription (ComplexMatrix* bMatrices, ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
						      long initialIndex, long nbrComponents)
{

  ComplexVector TmpVector1 (bMatrices->GetNbrColumn());
  ComplexVector TmpVector2 (bMatrices->GetNbrColumn());
  unsigned long MinIndex = (unsigned long) initialIndex;
  unsigned long MaxIndex = (unsigned long) this->LargeHilbertSpaceDimension;
  if (nbrComponents > 0)
    {
      MaxIndex = MinIndex + ((unsigned long) nbrComponents);
    }
  if (mPSRowIndex >= 0)
    {
      for (; MinIndex < MaxIndex; ++MinIndex)    
	{
	  state[(long) MinIndex] = 0.0;
	  TmpVector1.ClearVector();
	  TmpVector1[mPSColumnIndex] = 1.0;
	  for (int i = 0; i < this->ChainLength; ++i)
	    {
	      if ((MinIndex & (0x1ul << i)) == 0x0ul)
		{
		  TmpVector2.Multiply(bMatrices[0], TmpVector1);
		}
	      else
		{
		  TmpVector2.Multiply(bMatrices[1], TmpVector1);
		}
	      ComplexVector TmpVector3 = TmpVector2;
	      TmpVector2 = TmpVector1;
	      TmpVector1 = TmpVector3;
	    }
	  state[(long) MinIndex] = TmpVector1[mPSRowIndex];
	}
    }
  else
    {
      for (; MinIndex < MaxIndex; ++MinIndex)    
	{
	  state[(long) MinIndex] = 0.0;
	  for (int j = 0; j < bMatrices->GetNbrColumn(); ++j)
	    {
	      TmpVector1.ClearVector();
	      TmpVector1[j] = 1.0;
	      for (int i = 0; i < this->ChainLength; ++i)
		{
		  if ((MinIndex & (0x1ul << i)) == 0x0ul)
		    {
		      TmpVector2.Multiply(bMatrices[0], TmpVector1);
		    }
		  else
		    {
		      TmpVector2.Multiply(bMatrices[1], TmpVector1);
		    }
		  ComplexVector TmpVector3 = TmpVector2;
		  TmpVector2 = TmpVector1;
		  TmpVector1 = TmpVector3;
		}
	      state[(long) MinIndex] += TmpVector1[j];
	    }
	}
    }
}

