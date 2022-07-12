////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of bosons on disk                         //
//                                                                            //
//                        last modification : 04/06/2002                      //
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
#include "HilbertSpace/BosonOnDisk.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"

#include <math.h>
#include <iostream>


using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum angular momentum that a single particle can reach (negative if it has to be deduced from nbrBosons and totalLz)

BosonOnDisk::BosonOnDisk (int nbrBosons, int totalLz, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  if (lzMax < 0)
    this->LzMax = totalLz;
  else
    this->LzMax = lzMax;
  this->ShiftedTotalLz = this->TotalLz;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, this->TotalLz);
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  int TmpLzMax = this->LzMax;
  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->GenerateStates(this->NbrBosons, TmpLzMax, TmpLzMax, this->ShiftedTotalLz, 0);
  this->KeyMultiplicationTable = new int [this->LzMax + 1];
  this->GenerateLookUpTable(0);
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateLzMax[i] + 1) * sizeof(int) + sizeof(int*);
  UsedMemory += (this->TotalLz + 1) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += ((this->TotalLz + 1) * this->IncNbrBosons) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
/*  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;*/
#endif
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnDisk::BosonOnDisk(const BosonOnDisk& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->Flag = bosons.Flag;
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->LzMaxPosition = bosons.LzMaxPosition;
  this->KeyInvertSectorSize = bosons.KeyInvertSectorSize;
  this->KeyInvertTable = bosons.KeyInvertTable;
  this->KeyInvertTableNbrIndices = bosons.KeyInvertTableNbrIndices;
  this->KeyInvertIndices = bosons.KeyInvertIndices;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
}

// destructor
//

BosonOnDisk::~BosonOnDisk ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnDisk& BosonOnDisk::operator = (const BosonOnDisk& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Keys;
      delete[] this->KeyMultiplicationTable;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      int Size = (this->LzMax + 2) * this->IncNbrBosons;
      for (int i = 0; i < Size; ++i)
	{
	  if (this->KeyInvertSectorSize[i] > 0)
	    {
	      for (int j= 0; j < this->KeyInvertSectorSize[i]; ++j)
		delete[] this->KeyInvertIndices[i][j];
	      delete[] this->KeyInvertTable[i];
	      delete[] this->KeyInvertTableNbrIndices[i];
	      delete[] this->KeyInvertIndices[i];
	    }
	}      
      delete[] this->KeyInvertSectorSize;
      delete[] this->KeyInvertTable;
      delete[] this->KeyInvertTableNbrIndices;
      delete[] this->KeyInvertIndices;
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->Flag = bosons.Flag;
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->LzMaxPosition = bosons.LzMaxPosition;
  this->KeyInvertSectorSize = bosons.KeyInvertSectorSize;
  this->KeyInvertTable = bosons.KeyInvertTable;
  this->KeyInvertTableNbrIndices = bosons.KeyInvertTableNbrIndices;
  this->KeyInvertIndices = bosons.KeyInvertIndices;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  return *this;
}


// forge an eigenstate from a description given by a file
//
// filename = name of the file that contains the state description
// state = reference on the vector where the state has to be stored
// return value = true if no error occured

bool BosonOnDisk::ForgeEigenstate(char* filename, RealVector& state)
{
  int NbrComponents;
  Complex* Coefficients;
  int** ComponentDescription;
  if (this->LoadEigenstateDescrition(filename, this->TotalLz + 1, true, NbrComponents, Coefficients, ComponentDescription) == false)    
    return false;

  state.ResizeAndClean(this->HilbertSpaceDimension);
  
  double TmpNorm = 0.0;
  for (int i = 0; i < NbrComponents; ++i)
    {
      TmpNorm += Coefficients[i].Re * Coefficients[i].Re;
    }
  TmpNorm = 1.0 / sqrt(TmpNorm);
  //  int Pos = 0;
  bool SuccessfullParsing = true;
  for (int i = 0; i < NbrComponents; ++i)
    {
      int TmpNbrParticles = 0;
      int TmpMomentum = 0;
      bool ErrorFlag = false;
      int* TmpDescription = ComponentDescription[i];
      int TmpLzMax = 0;
      for (int j = 0; j <= this->TotalLz; ++j)
	{
	  if (TmpDescription[j] < 0)
	    {
	      ErrorFlag = true;
	    }
	  else
	    {
	      this->TemporaryState[j] = TmpDescription[j];
	      if (TmpDescription[j] > 0)
		{
		  TmpLzMax = j;
		  TmpMomentum += (TmpDescription[j] * j);
		  TmpNbrParticles += TmpDescription[j];
		}
	    }
	}
      if (TmpMomentum != this->TotalLz)
	{
	  cout << "wrong total momentum in state description ";
	  for (int j = 0; j <= this->TotalLz; ++j)
	    cout << TmpDescription[j] << " ";
	  cout << endl;
	  SuccessfullParsing = false;
	}
      else
	if (TmpNbrParticles != this->NbrBosons)
	  {
	    cout << "wrong number of particles in state description ";
	    for (int j = 0; j <= this->TotalLz; ++j)
	      cout << TmpDescription[j] << " ";
	    cout << endl;
	    SuccessfullParsing = false;
	  }
	else
	  {
	    int TmpIndex = this->FindStateIndex(this->TemporaryState, TmpLzMax);
	    if (TmpIndex < this->HilbertSpaceDimension)
	      {
		state[TmpIndex] += Coefficients[i].Re * TmpNorm;
	      }
	  }
      delete[] TmpDescription;
   }
  state /= state.Norm();
  delete[] Coefficients;
  delete[] ComponentDescription;
  return SuccessfullParsing;
}

