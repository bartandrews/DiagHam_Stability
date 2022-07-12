////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            /
//             class of doubled spin 0 +1/2 chain with translations           //
//                                                                            //
//                        last modification : 21/01/2016                      //
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

#include "HilbertSpace/VirtualSpacePEPSWithTranslations.h"
#include "HilbertSpace/VirtualSpaceTransferMatrixWithTranslations.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/ArrayTools.h"
#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


// default constructor
//

VirtualSpaceTransferMatrixWithTranslations::VirtualSpaceTransferMatrixWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->BondDimension = 0;
  this->FixedSpinProjectionFlag = false;
  this->CompatibilityWithMomentum = 0;
  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->PowerD = 0;
  this->TranslationPhase=0;
}

// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

VirtualSpaceTransferMatrixWithTranslations::VirtualSpaceTransferMatrixWithTranslations (int chainLength, int bondDimension, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->BondDimension = bondDimension;
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;
  this->TranslationPhase=0;

  this->PowerD = new int [this->ChainLength];
  this->PowerD[0]=1;
  this->PowerD[1]=  this->BondDimension* this->BondDimension;
  for(int i = 2 ; i <this->ChainLength;i++)
    this->PowerD[i]=this->PowerD[i-1]*this->PowerD[1];
  
  this->LargeHilbertSpaceDimension =  this->PowerD[this->ChainLength-1]*this->PowerD[1];
  cout << "Hilbert Space Dimension = " <<this->LargeHilbertSpaceDimension<<endl;
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, 0l);

  SortArrayDownOrdering(this->ChainDescription ,TmpHilbertSpaceDimension);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in DoubledSpin0_1_2_ChainWithTranslations!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
    }
  this->RescalingFactors = 0;
  int MemoryCost = sizeof(unsigned long) * this->LargeHilbertSpaceDimension + this->ChainLength *  sizeof(int);
  cout <<" Memory cost of the Hilbert Space " <<MemoryCost <<endl;
}
 


// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

VirtualSpaceTransferMatrixWithTranslations::VirtualSpaceTransferMatrixWithTranslations (int chainLength,  int bondDimension, int momentum,  int translationStep,  int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->BondDimension = bondDimension;

  this->MaxXMomentum = this->ChainLength/ translationStep;
  this->Momentum = momentum %  this->MaxXMomentum;
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;


  this->PowerD = new int [this->ChainLength];
  this->PowerD[0]=1;
  this->PowerD[1]=  this->BondDimension* this->BondDimension;
  for(int i = 2 ; i <this->ChainLength;i++)
    this->PowerD[i]=this->PowerD[i-1]*this->PowerD[1] ;

  
  this->LargeHilbertSpaceDimension =  this->PowerD[this->ChainLength-1]*this->PowerD[1];
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  cout << "Hilbert Space Dimension = " <<this->LargeHilbertSpaceDimension<<endl;
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1,0l);

  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2Spin!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long TmpCanonicalState;

  int NbrTranslation;
  int CurrentNbrStateInOrbit;
  unsigned long DicardFlag = ~0x0ul;
  unsigned long TmpState;
  this->CreatePrecalculationTable();
  SortArrayDownOrdering(this->ChainDescription ,TmpHilbertSpaceDimension);
  
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpState = this->ChainDescription[i];
      this->FindCanonicalForm(TmpState,TmpCanonicalState, NbrTranslation);

      if ((TmpState  == TmpCanonicalState))
	{
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpCanonicalState);
	  if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->ChainDescription[i] = DicardFlag;
	    }
	}
      else
	{
	  this->ChainDescription[i] = DicardFlag;
	}
    }
  
  unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];

  this->NbrStateInOrbit = new int [this->LargeHilbertSpaceDimension];

  this->LargeHilbertSpaceDimension = 0l;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      if ( this->ChainDescription[i] != DicardFlag)
	{
	  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->ChainDescription[i];
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(this->ChainDescription[i]);
	  this->NbrStateInOrbit[this->LargeHilbertSpaceDimension] = CurrentNbrStateInOrbit;
	  ++this->LargeHilbertSpaceDimension;
	}
    }
  
  delete[]  this->ChainDescription;
  this->ChainDescription = TmpStateDescription;
  
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  
  this->LookUpTable =0;
  
  
  if (this->HilbertSpaceDimension > 0)
    this->GenerateLookUpTable(memorySize);  
  this->EvaluateExponentialFactors();

  unsigned long MemoryCost = (sizeof(unsigned long) + sizeof(int)) * this->LargeHilbertSpaceDimension + this->ChainLength *  sizeof(int) + this->MaxXMomentum*sizeof(Complex);
  cout <<" Memory cost of the Hilbert Space " <<MemoryCost <<endl;
}


// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

VirtualSpaceTransferMatrixWithTranslations::VirtualSpaceTransferMatrixWithTranslations (const VirtualSpaceTransferMatrixWithTranslations & chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescription = chain.ChainDescription;
      this->BondDimension = chain.BondDimension;
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->UniqueStateDescriptionBra = chain.UniqueStateDescriptionBra;
      this->UniqueStateDescriptionSubArraySizeBra = chain.UniqueStateDescriptionSubArraySizeBra;
      this->NbrUniqueStateDescriptionBra = chain.NbrUniqueStateDescriptionBra;
      this->FirstIndexUniqueStateDescriptionBra = chain.FirstIndexUniqueStateDescriptionBra;
      this->PowerD = chain.PowerD;
      this->TranslationPhase=chain.TranslationPhase;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescriptionBra = 0;
      this->ChainDescriptionKet = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->BondDimension = 0;
      this->ChainDescription = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->PowerD = 0;
      this->TranslationPhase=0;
      this->MaxXMomentum = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

VirtualSpaceTransferMatrixWithTranslations::~VirtualSpaceTransferMatrixWithTranslations () 
{
  delete [] this->PowerD;
  delete [] this->ChainDescription;
  delete [] TranslationPhase;
  this->LargeHilbertSpaceDimension = 0;
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

VirtualSpaceTransferMatrixWithTranslations & VirtualSpaceTransferMatrixWithTranslations::operator = (const VirtualSpaceTransferMatrixWithTranslations & chain)
{
  AbstractDoubledSpinChainWithTranslations::operator =(chain);
  this->PowerD = chain.PowerD;
  this->BondDimension = chain.BondDimension;
  this->ChainDescription = chain.ChainDescription;
  this->TranslationPhase = chain.TranslationPhase;
  this->MaxXMomentum = chain.MaxXMomentum;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* VirtualSpaceTransferMatrixWithTranslations::Clone()
{
  return new VirtualSpaceTransferMatrixWithTranslations (*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& VirtualSpaceTransferMatrixWithTranslations::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned int tmpBra,tmpKet;
  int BraNumber, KetNumber;
  unsigned long StateDescription = this->ChainDescription[state];  

  Str << this->FindStateIndex(StateDescription) << " : " << StateDescription<< " : ";
  for (int j = this->ChainLength; j >0; j--)
    {
      this->GetBraAndKetIndicesFromCommonIndex(tmpBra,tmpKet, StateDescription%this->PowerD[1]);
      StateDescription/=this->PowerD[1];
      Str << "("<< tmpBra <<" " << tmpKet<< ") ";
    }
  return Str;
}


// generate all states corresponding to the constraints
// 
// lengthBra = length of the chain to be decided for bra spins
// lengthBra = length of the chain to be decided for ket spins
// diffSz = difference of spin projection between bra and ket chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
long VirtualSpaceTransferMatrixWithTranslations::GenerateStates(int length, long pos)
{
  if (length == 0 )
    {
      for (int i =0; i < this->PowerD[1]; i++)
	{
	  this->ChainDescription[pos] = i;
	  pos++;
	}
      return pos;
    }
  
  if(length > 0)
    { 
      long TmpPos;
      for (int i = 0; i < this->PowerD[1]; i++)
	{
	  TmpPos = this->GenerateStates(length-1, pos); 
	  for (; pos < TmpPos; ++pos)
	    {
	      this->ChainDescription[pos] += i*this->PowerD[length];
	    }
	}
      return pos;
    }
  return pos;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long VirtualSpaceTransferMatrixWithTranslations::ShiftedEvaluateHilbertSpaceDimension(int length)
{
  long Tmp =1;
  for(int i =0; i < length; i++)
    {
      Tmp *= this->PowerD[1];
    }
  return Tmp;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void VirtualSpaceTransferMatrixWithTranslations::GenerateLookUpTable(unsigned long memory)
{  
  this->LookUpTable = new long [this->PowerD[2]+1];
  int CurrentBeginning = this->ChainDescription[0] / this->PowerD[this->ChainLength-2];
  for (int i =CurrentBeginning+1; i <this->PowerD[2]+1;i++)
    {
      this->LookUpTable[i]=0;
    }
  for(int i =0; i < this->HilbertSpaceDimension; i++)
    {
      if( this->ChainDescription[i]/this->PowerD[this->ChainLength-2] >= CurrentBeginning)
	{
	  
	}
      else
	{
	  for(int p = this->ChainDescription[i] / this->PowerD[this->ChainLength-2] + 1 ;p <= CurrentBeginning;p++)
	    this->LookUpTable[p] = i;
	  CurrentBeginning = this->ChainDescription[i] / this->PowerD[this->ChainLength-2];
	}
    }
  for(int p = 0 ;p <= CurrentBeginning;p++)
  this->LookUpTable[p] = this->HilbertSpaceDimension;
}
 

int VirtualSpaceTransferMatrixWithTranslations::FindStateIndex(unsigned long stateDescription)
{
  int PosMax = this->LookUpTable[stateDescription/this->PowerD[this->ChainLength-2]];
  int PosMin = this->LookUpTable[stateDescription/this->PowerD[this->ChainLength-2]+1];
  int PosMid = (PosMin + PosMax) >>1;
  while( (PosMax - PosMin) >1 ) 
    {
      if (this->ChainDescription[PosMid]  >  stateDescription)
	{
	  PosMin = PosMid;
	}
      else
	{
	  if (this->ChainDescription[PosMid]  <  stateDescription)
	    {
	      PosMax = PosMid;
	    }
	  else
	    return PosMid;
	}
      
      PosMid = (PosMin + PosMax) >>1;
    }
  
  if (this->ChainDescription[PosMid] == stateDescription)
    return PosMid;
  else
    {
      if ((PosMax == this->HilbertSpaceDimension) ||(this->ChainDescription[PosMax] != stateDescription))
	return this->HilbertSpaceDimension;
      else
	return PosMax;
    }
}
 
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix VirtualSpaceTransferMatrixWithTranslations::EvaluatePartialDensityMatrix (RealVector& groundState)
{
  VirtualSpacePEPSWithTranslations TmpDestinationHilbertSpace (this->ChainLength,this->BondDimension, 10000,10000); 
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  unsigned long TmpBra,TmpKet,TmpState;

  ComplexMatrix HRep(TmpDestinationHilbertSpace.HilbertSpaceDimension,TmpDestinationHilbertSpace.HilbertSpaceDimension,true);
  for(int i=0;i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      for(int j=0;j < TmpDestinationHilbertSpace.HilbertSpaceDimension;j++)
	{
	  TmpState=0;
	  TmpBra = TmpDestinationHilbertSpace.ChainDescription[i];
	  TmpKet = TmpDestinationHilbertSpace.ChainDescription[j];
	  for (int p = 0;p <this->ChainLength;p++)
	    {
	      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra %this->BondDimension , TmpKet %this->BondDimension);
	      TmpBra/=this->BondDimension;
	      TmpKet/=this->BondDimension;
	    }
	  int Index = this->FindStateIndex (TmpState);
	  if (Index < this->HilbertSpaceDimension ) 
	    HRep.SetMatrixElement(i,j,groundState[Index]);
	}
    }
  
  ComplexMatrix SquareRho( TmpDestinationHilbertSpace.HilbertSpaceDimension,  TmpDestinationHilbertSpace.HilbertSpaceDimension,true);
  SquareRho = HRep*HRep;
  
  double Trace = 0.0;
  double Tmp = 0.0;
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      SquareRho.GetMatrixElement(i,i,Tmp);
      Trace+=Tmp;
    }
  cout <<"Trace "<< Trace<<endl;

  
  double Tmp1;
  double Tmp2;
  cout << "check hermiticity" << endl;

  double Error = 5e-8;
  
  for (int i = 0; i <  TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    for (int j = i; j <  TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
      {
	HRep.GetMatrixElement(i, j, Tmp1);
	HRep.GetMatrixElement(j, i, Tmp2);
	if (Norm(Tmp1 - Tmp2) > Error )
	  {
	    cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << " " << Norm(Tmp1 - Tmp2) << " (should be lower than " << (Error ) << ")" << endl;
	  }
	HRep.GetMatrixElement(i,j,Tmp);
	TmpDensityMatrix.SetMatrixElement(i,j,Tmp);
      }  
  return TmpDensityMatrix;
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix VirtualSpaceTransferMatrixWithTranslations::EvaluatePartialDensityMatrix (ComplexVector& groundState)
{
  VirtualSpacePEPSWithTranslations TmpDestinationHilbertSpace (this->ChainLength,this->BondDimension, 10000,10000);
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  unsigned long TmpBra,TmpKet,TmpState;

  ComplexMatrix HRep(TmpDestinationHilbertSpace.HilbertSpaceDimension,TmpDestinationHilbertSpace.HilbertSpaceDimension,true);
  for(int i=0;i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      for(int j=0;j < TmpDestinationHilbertSpace.HilbertSpaceDimension;j++)
	{
	  TmpState=0;
	  TmpBra = TmpDestinationHilbertSpace.ChainDescription[i];
	  TmpKet = TmpDestinationHilbertSpace.ChainDescription[j];
	  for (int p = 0;p <this->ChainLength;p++)
	    {
	      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra %this->BondDimension , TmpKet %this->BondDimension);
	      TmpBra/=this->BondDimension;
	      TmpKet/=this->BondDimension;
	    }
	  int Index = this->FindStateIndex (TmpState);
	  if (Index < this->HilbertSpaceDimension ) 
	    HRep.SetMatrixElement(i,j,groundState[Index]);
	}
    }

  ComplexMatrix SquareRho( TmpDestinationHilbertSpace.HilbertSpaceDimension,  TmpDestinationHilbertSpace.HilbertSpaceDimension,true);
  SquareRho = HRep*HRep;
  
  Complex Trace = 0.0;
  Complex Tmp = 0.0;
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      SquareRho.GetMatrixElement(i,i,Tmp);
      Trace+=Tmp;
    }
  cout <<"Trace "<< Trace<<endl;
  
  Complex Tmp1;
  Complex Tmp2;
  cout << "check hermiticity" << endl;

  double Error = 5e-8;
  
  for (int i = 0; i <  TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    for (int j = i; j <  TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
      {
	HRep.GetMatrixElement(i, j, Tmp1);
	HRep.GetMatrixElement(j, i, Tmp2);
	if (Norm(Tmp1 - Conj(Tmp2)) > Error )
	  {
	    cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << " " << Norm(Tmp1 - Conj(Tmp2)) << " (should be lower than " << (Error ) << ")" << endl;
	  }
	HRep.GetMatrixElement(i,j,Tmp);
	TmpDensityMatrix.SetMatrixElement(i,j,Tmp);
      }  
  return TmpDensityMatrix;
}



void VirtualSpaceTransferMatrixWithTranslations::ConvertToGeneralSpace(ComplexVector vSource,ComplexVector & vDestination)
{
  for(int i =0; i <this->HilbertSpaceDimension; i++)
    {
      vDestination[(long) this->ChainDescription[i]] = vSource[i];
    }
}

void VirtualSpaceTransferMatrixWithTranslations::AddConvertFromGeneralSpace(ComplexVector vSource,ComplexVector & vDestination)
{
  for(int i =0; i <this->HilbertSpaceDimension; i++)
    {
      vDestination[i] = vSource[(long) this->ChainDescription[i]];
    }
}


void VirtualSpaceTransferMatrixWithTranslations::ConvertToGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination)
{
  for(int i =0; i <this->HilbertSpaceDimension; i++)
    {
      vDestination[(long) this->ChainDescription[i]] = vSource[i]* sqrt ( ((double) this->NbrStateInOrbit[i]));
    }
}

void VirtualSpaceTransferMatrixWithTranslations::AddConvertFromGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination)
{
  unsigned long TmpState;
  for(int i = 0; i <this->HilbertSpaceDimension; i++)
    {
      TmpState = (unsigned long) this->ChainDescription[i];
      for(int p = 0; p < this->NbrStateInOrbit[i]; p++)
	{
	  vDestination[i] += this->TranslationPhase[p] * vSource[(int ) TmpState] / sqrt ( ((double) this->NbrStateInOrbit[i]));
 	  this->ApplySingleXTranslation(TmpState);
	}
    }
}

// evaluate all exponential factors
//   

void VirtualSpaceTransferMatrixWithTranslations::EvaluateExponentialFactors()
{
  this->TranslationPhase = new Complex[this->MaxXMomentum];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      this->TranslationPhase[i] = Phase(2.0 * M_PI * ((this->Momentum * ((double) i) / ((double) this->MaxXMomentum))));
    }
}

HermitianMatrix VirtualSpaceTransferMatrixWithTranslations::EvaluatePartialDensityMatrix ( int momentumSector, ComplexVector& groundState)
{
  cout <<"creating ket Hilbert space"<<endl;
  VirtualSpacePEPSWithTranslations TmpDestinationHilbertSpace (this->ChainLength,this->BondDimension, momentumSector,  this->ChainLength / this->MaxXMomentum, 10000,10000);
  int ComplementaryKSector = (this->Momentum - momentumSector) %  this->MaxXMomentum;
  if (ComplementaryKSector < 0)
    ComplementaryKSector +=  this->MaxXMomentum;
  

  cout <<"creating bra Hilbert space"<<endl;
  VirtualSpacePEPSWithTranslations TmpHilbertSpace (this->ChainLength,this->BondDimension, ComplementaryKSector,  this->ChainLength / this->MaxXMomentum, 10000,10000);
  
  int MaxDimension = TmpDestinationHilbertSpace.HilbertSpaceDimension;
  
  if ( MaxDimension < TmpHilbertSpace.HilbertSpaceDimension )
    {
      MaxDimension = TmpHilbertSpace.HilbertSpaceDimension;
    }
  
  HermitianMatrix TmpDensityMatrix(MaxDimension, true);
  unsigned long TmpState,TmpBra,TmpKet,TmpCanonicalState,ReferenceBra,ReferenceKet; 
  int NbrTranslation;
  cout<<" filling the matrix up"<<endl;
  cout <<TmpDestinationHilbertSpace.HilbertSpaceDimension<<" " << TmpHilbertSpace.HilbertSpaceDimension<<endl;
  ComplexMatrix HRep( MaxDimension, MaxDimension,true); 
  for(int i = 0 ; i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      for(int j = 0 ; j < TmpHilbertSpace.HilbertSpaceDimension;j++)
	{
	  ReferenceBra = TmpDestinationHilbertSpace.ChainDescription[i];
	  ReferenceKet = TmpHilbertSpace.ChainDescription[j];
	  
	  for(int t = 0; t < TmpDestinationHilbertSpace.NbrStateInOrbit[i]; t++)
	    {
	      for(int k = 0; k < TmpHilbertSpace.NbrStateInOrbit[j]; k++)
		{
		  TmpState=0;
		  TmpBra =  ReferenceBra;
		  TmpKet =  ReferenceKet;
		  for (int p = 0;p <this->ChainLength;p++)
		    {
		      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra %this->BondDimension , TmpKet %this->BondDimension);
		      TmpBra/=this->BondDimension;
		      TmpKet/=this->BondDimension;
		    }
		  
		  this->FindCanonicalForm(TmpState,TmpCanonicalState, NbrTranslation);
		  int Index = this->FindStateIndex (TmpCanonicalState);
		  if (Index < this->HilbertSpaceDimension ) 
		    {
		      double TmpFactor =sqrt( (double) (TmpDestinationHilbertSpace.NbrStateInOrbit[i] * TmpHilbertSpace.NbrStateInOrbit[j] * ( double) this->NbrStateInOrbit[Index])); 
		      double Argument =  2.0 * M_PI * (NbrTranslation * this->Momentum  - t * momentumSector  - k * ComplementaryKSector) /  this->MaxXMomentum ;
		      HRep.AddToMatrixElement(i,j,groundState[Index]/TmpFactor*Phase(Argument));
		    }
		  TmpHilbertSpace.ApplySingleXTranslation(ReferenceKet);
		}
	      TmpDestinationHilbertSpace.ApplySingleXTranslation(ReferenceBra);
	    }
	}
    }
  
  ComplexMatrix SquareRho( MaxDimension, MaxDimension,true);
  SquareRho = HRep*HRep;

  Complex Trace = 0.0;
  Complex Tmp = 0.0;
  for (int i = 0; i < MaxDimension;i++)
    {
      SquareRho.GetMatrixElement(i,i,Tmp);
      Trace+=Tmp;
    }
  cout <<"Trace "<< Trace<<endl;
  HRep/= Phase(0.5*Arg(Trace));
  Complex Tmp1;
  Complex Tmp2;
  cout << "check hermiticity" << endl;
  
  double Error = 5e-8;
  
  for (int i = 0; i <   MaxDimension; ++i)
    for (int j = i; j <  MaxDimension; ++j)
      {
	HRep.GetMatrixElement(i, j, Tmp1);
	HRep.GetMatrixElement(j, i, Tmp2);
	if (Norm(Tmp1 - Conj(Tmp2)) > Error )
	  {
	    cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << " " << Norm(Tmp1 - Conj(Tmp2)) << " (should be lower than " << (Error ) << ")" << endl;
	  }
	HRep.GetMatrixElement(i,j,Tmp);
	TmpDensityMatrix.SetMatrixElement(i,j,Tmp);
      }
  return TmpDensityMatrix;
}

/*
void  DoubledSpin0_1_2_ChainWithTranslations::ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector)
{
//  cout <<"using void  DoubledSpin0_1_2_ChainWithTranslations::ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector)"<<endl;
  RealDiagonalMatrix ZMatrix(9,true);
  double  Tmp[3];
  Tmp[0] = -1.0;Tmp[1] = -1.0;Tmp[2] = 1.0;
  for(int i =0;i <9; i++)
    {
      ZMatrix.SetMatrixElement(i,i,Tmp[i%3] *Tmp[i/3] ); 
    }
//  cout <<ZMatrix<<endl;
  double TmpFactor;
  unsigned long OldState;
  int TmpIndice;
  unsigned long  TmpStateDescription;


  int tmp = this->FindNextInversionSymetricIndice(0);

  if((tmp ==this->HilbertSpaceDimension) )
    {
      tmp=0;
      while(Norm(sourceVector[tmp]) < 1e-8 )
	tmp++;
    }
  
  while(Norm(sourceVector[tmp]) < 1e-8 )
    {
      tmp=this->FindNextInversionSymetricIndice(tmp);
    }
  
  cout <<tmp<<" "<< Norm(sourceVector[tmp])<<endl;
//  sourceVector*=Conj(sourceVector[tmp])/Norm(sourceVector[tmp]);
  sourceVector/=Phase(Arg(sourceVector[tmp]));

//  this->NormalizeDensityMatrix(sourceVector);
  for (int i =0; i <this->HilbertSpaceDimension; i++)
    {
      cout <<" i = " <<i<< " : "<<this->ChainDescription[i]<<endl;;
      TmpFactor = 1.0;
      OldState = 0ul;
      TmpStateDescription =  this->ChainDescription[i]; 
      for (int p = 0 ; p <this->ChainLength ; p++)
	{
	  TmpIndice = (TmpStateDescription/this->PowerD[p]) % this->PowerD[1];
	  cout <<TmpIndice<<" ";
	  if  (p%2 == 0 ) 
	    TmpFactor*=ZMatrix(TmpIndice,TmpIndice);
	  OldState+= TmpIndice*this->PowerD[this->ChainLength-p-1];
	}

      cout <<endl<<" "<<OldState<< " "<<this->FindStateIndex(OldState)<<" "<<sourceVector[this->FindStateIndex(OldState)]<<endl;
      destinationVector[i] = Conj(TmpFactor*sourceVector[this->FindStateIndex(OldState)]) ;
    }
}
*/

int VirtualSpaceTransferMatrixWithTranslations::FindNextInversionSymetricIndice(int previousOne)
{
  unsigned long OldState;
  int TmpIndice;
  unsigned long  TmpStateDescription;
  for (int i =previousOne+1; i <this->HilbertSpaceDimension; i++)
    {
      OldState = 0ul;
      TmpStateDescription =  this->ChainDescription[i]; 
      for (int p = 0 ; p <this->ChainLength ; p++)
	{
	  TmpIndice = (TmpStateDescription/this->PowerD[p]) % this->PowerD[1];
	  OldState+= TmpIndice*this->PowerD[this->ChainLength-p-1];
	} 
      if (OldState ==  this->ChainDescription[i]) 
	return i;
    } 
  return this->HilbertSpaceDimension;
}

/*
void  DoubledSpin0_1_2_ChainWithTranslations::ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector, bool translationFlag)
{
  if (translationFlag == false ) 
    {
      this->ApplyInversionSymmetry(sourceVector, destinationVector) ;
    }
  else
    {
      //  cout <<"using void  DoubledSpin0_1_2_ChainWithTranslations::ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector)"<<endl;
      RealDiagonalMatrix ZMatrix(9,true);
      double  Tmp[3];
      Tmp[0] = -1.0;Tmp[1] = -1.0;Tmp[2] = 1.0;
      for(int i =0;i <9; i++)
	{
	  ZMatrix.SetMatrixElement(i,i,Tmp[i%3] *Tmp[i/3] ); 
	}
      //  cout <<ZMatrix<<endl;
      double TmpFactor;
      unsigned long OldState;
      int TmpIndice;
      unsigned long  TmpStateDescription;
      unsigned long  TmpCanonicalState;
      int NbrTranslation;

      int tmp = this->FindNextInversionSymetricIndice(0);
      
      if((tmp ==this->HilbertSpaceDimension) )
	{
	  tmp=0;
	  while(Norm(sourceVector[tmp]) < 1e-8 )
	    tmp++;
	}
      
      while(Norm(sourceVector[tmp]) < 1e-8 )
	{
	  tmp=this->FindNextInversionSymetricIndice(tmp);
	}
      

      sourceVector*=Conj(sourceVector[tmp])/Norm(sourceVector[tmp]);
      for (int i =0; i <this->HilbertSpaceDimension; i++)
	{
	  cout <<" i = " <<i<< " : "<<this->ChainDescription[i]<<endl;;
	  TmpFactor = 1.0;
	  OldState = 0ul;
	  TmpStateDescription =  this->ChainDescription[i]; 
	  for (int p = 0 ; p <this->ChainLength ; p++)
	    {
	      TmpIndice = (TmpStateDescription/this->PowerD[p]) % this->PowerD[1];
	      cout <<TmpIndice<<" ";
	      if  (p%2 == 0 ) 
		TmpFactor*=ZMatrix(TmpIndice,TmpIndice);
	      OldState+= TmpIndice*this->PowerD[this->ChainLength-p-1];
	    }

	  this->FindCanonicalForm(OldState,TmpCanonicalState, NbrTranslation);	  
//	  cout <<endl<<" "<<OldState<< " "<<" "<<this->FindStateIndex(TmpCanonicalState)<<" "<< NbrTranslation <<" "<<sourceVector[this->FindStateIndex(TmpCanonicalState)]<<endl;


	  destinationVector[i] = Conj(TmpFactor*sourceVector[this->FindStateIndex(TmpCanonicalState)]*this->TranslationPhase[NbrTranslation]) ;
	}
    }
}
*/

void VirtualSpaceTransferMatrixWithTranslations::NormalizeDensityMatrix(ComplexVector & sourceVector)
{
  unsigned long SourceState,TmpState;
  unsigned int TmpBra,TmpKet;
  Complex Trace = 0.0;
  for(int i=0;i < this->HilbertSpaceDimension;i++)
    {
      if ( Norm(sourceVector[i]) > 1e-8 ) 
	{
	  TmpState=0;
	  SourceState = this->ChainDescription[i];
	  for (int p = 0;p <this->ChainLength;p++)
	    {
	      this->GetBraAndKetIndicesFromCommonIndex(TmpBra,TmpKet, SourceState%this->PowerD[1]);
	      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpKet, TmpBra );
	      SourceState/=this->PowerD[1];
	    }
	  int Index = this->FindStateIndex (TmpState);
	  if (Index < this->HilbertSpaceDimension ) 
	    {
	      sourceVector/=Phase(0.5*(Arg(sourceVector[i]) + Arg(sourceVector[Index])));
	      return;
	    }
	}
    }
}
