////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
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


#include "HilbertSpace/DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers.h"
#include "HilbertSpace/Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers.h"
// #include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry.cc"

class DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry;

#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


// default constructor
//

DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers::DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers () 
{
}

 

// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers::DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers (int chainLength,  int diffSz, int subLatticeDifferenceKet, int subLatticeDifferenceBra, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->SubLatticeDifferenceKet = subLatticeDifferenceKet;
  this->SubLatticeDifferenceBra = subLatticeDifferenceBra;

  
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
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
  for(int i = 1 ; i <this->ChainLength;i++)
    this->PowerD[i]=this->PowerD[i-1]*9;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, this->DiffSz);
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);

  SortArrayDownOrdering(this->ChainDescription,TmpHilbertSpaceDimension);
  
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in DoubledSpin0_1_2_ChainWithTranslations!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  

  this->LargeHilbertSpaceDimension = 0l;
  unsigned long DicardFlag = ~0x0ul;
  int BraNumber, KetNumber;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      this->ComputeDiffereenceSubLatticeNumberZero(this->ChainDescription[i], BraNumber, KetNumber);
      if (( BraNumber == this->SubLatticeDifferenceBra ) && ( KetNumber == this->SubLatticeDifferenceKet ) )
	{
	  ++this->LargeHilbertSpaceDimension;
	}
      else
	{
	  this->ChainDescription[i] = DicardFlag;
	}
    }
  
  unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  
  this->LargeHilbertSpaceDimension = 0l;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      if ( this->ChainDescription[i] != DicardFlag)
	{
	  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->ChainDescription[i];
	  ++this->LargeHilbertSpaceDimension;
	}
    }
  
  delete[]  this->ChainDescription;
  
  this->ChainDescription = TmpStateDescription;
  
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  
  this->LookUpTable =0;
  
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
    }
  this->RescalingFactors = 0;
  cout <<"Hilbert Space dimension = "<< this->GetHilbertSpaceDimension()<<endl;  
}



// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers::DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers (int chainLength, int momentum, int sz, int subLatticeDifferenceKet, int subLatticeDifferenceBra, int memorySize, int memorySlice)
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = sz;
  this->FixedSpinProjectionFlag = true;
  this->Momentum = momentum;
  this->SubLatticeDifferenceKet = subLatticeDifferenceKet;
  this->SubLatticeDifferenceBra = subLatticeDifferenceBra;
  
  this->MaxXMomentum = this->ChainLength / 2;  
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);

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
  for(int i = 1 ; i <this->ChainLength;i++)
    this->PowerD[i]=this->PowerD[i-1]*9;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, this->DiffSz);
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);
  SortArrayDownOrdering(this->ChainDescription,TmpHilbertSpaceDimension);
  
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
  int BraNumber,KetNumber;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpState = this->ChainDescription[i];
      this->FindCanonicalForm(TmpState,TmpCanonicalState, NbrTranslation);
      
      if (TmpState  == TmpCanonicalState)
	{
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpCanonicalState);
	  if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)
	    {
	      this->ComputeDiffereenceSubLatticeNumberZero(this->ChainDescription[i], BraNumber, KetNumber);
	      if ((  BraNumber == this->SubLatticeDifferenceBra ) && ( KetNumber == this->SubLatticeDifferenceKet  ) )
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
  
  cout <<"Hilbert Space dimension = "<< this->GetHilbertSpaceDimension()<<endl;  
}


// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers::DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers (const DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers & chain) :    DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers(chain)
{
}

// destructor
//

DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers::~DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers & DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers::operator = (const DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers & chain)
{
  DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers::operator =(chain);
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace*  DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers::Clone()
{
  return new  DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers (*this);
}


HermitianMatrix DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers::EvaluatePartialDensityMatrix (int szSector, int momentumSector, int sublatticeQuantumNumberSector,  int complementarySublatticeQuantumNumberSector, ComplexVector& groundState)
{

  Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers TmpDestinationHilbertSpace(this->ChainLength, momentumSector ,2, szSector,sublatticeQuantumNumberSector, 10000,10000);
 
  int ComplementaryKSector = (this->Momentum - momentumSector) %  this->MaxXMomentum;


  if (ComplementaryKSector < 0)
    ComplementaryKSector +=  this->MaxXMomentum;
  
  Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers TmpHilbertSpace(this->ChainLength,ComplementaryKSector,2 ,szSector, complementarySublatticeQuantumNumberSector, 10000,10000);
  cout << TmpDestinationHilbertSpace.HilbertSpaceDimension << " "<<TmpHilbertSpace.HilbertSpaceDimension<<endl;
  
  int MaxDimension = TmpDestinationHilbertSpace.HilbertSpaceDimension;
  
  if ( MaxDimension < TmpHilbertSpace.HilbertSpaceDimension )
    {
      MaxDimension = TmpHilbertSpace.HilbertSpaceDimension;
    }
  
  HermitianMatrix TmpDensityMatrix(MaxDimension, true);
  unsigned long TmpState,TmpBra,TmpKet,TmpCanonicalState,ReferenceBra,ReferenceKet; 
  int NbrTranslation;

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
		  TmpState = 0;
		  TmpBra = ReferenceBra;
		  TmpKet = ReferenceKet;
		  for (int p = 0 ; p < this->ChainLength ; p++)
		    {
		      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra &0x3ul, TmpKet &0x3ul);
		      TmpBra>>=2;
		      TmpKet>>=2;
		    }
		  
		  this->FindCanonicalForm(TmpState,TmpCanonicalState, NbrTranslation);
		  
		  int Index = this->FindStateIndex (TmpCanonicalState);
		  if (Index < this->HilbertSpaceDimension ) 
		    {
		      double TmpFactor = sqrt( (double) (TmpDestinationHilbertSpace.NbrStateInOrbit[i] * TmpHilbertSpace.NbrStateInOrbit[j] * ( double) this->NbrStateInOrbit[Index])); 
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

