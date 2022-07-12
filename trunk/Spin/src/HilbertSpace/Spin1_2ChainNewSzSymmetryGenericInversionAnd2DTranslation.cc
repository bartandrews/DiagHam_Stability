////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 1/2 chain with Sz contraint                //
//                and 2d translations plus the Sz<->-Sz symmetry              //
//                       and any generic inversion symmetry                   //
//                                                                            //
//                        last modification : 14/06/2017                      //
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


#include "HilbertSpace/Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <iostream>


using std::cout;
using std::dec;
using std::hex;
using std::endl;
using std::ios;



// default constructor
//

Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation::Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation ()
{
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// nbrSite = total number or spins
// sz = value of the total magnetization
// szSymmetrySector = Sz<->-Sz symmetry sector (can be either +1 or -1)
// inversionSymmetrySectorSector = inversion symmetry sector (can be either +1 or -1)
// inversionList = array that give the mapped indices under the inversion symmetry
// xMomentum = momentum along the x direction
// maxXMomentum = number of sites in the x direction
// yMomentum = momentum along the y direction
// maxYMomentum = number of sites in the y direction
// memory = amount of memory granted for precalculations

Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation::Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation (int nbrSite, int sz, int szSymmetrySector, 
														      int inversionSymmetrySector, int* inversionList, 
														      int xMomentum, int maxXMomentum, 
														      int yMomentum, int maxYMomentum, unsigned long memory) 
{
  this->Flag.Initialize();
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->NbrSite = nbrSite;
  this->ChainLength = nbrSite;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  if (szSymmetrySector == 1)
    {
      this->SzSymmetrySector = 1.0;
    }
  else
    {
      this->SzSymmetrySector = -1.0;
    }
  this->SzSymmetryMask = (0x1ul << this->ChainLength) - 0x1ul;
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

  if (inversionSymmetrySector == 1)
    {
      this->InversionSymmetrySector = 1.0;
    }
  else
    {
      this->InversionSymmetrySector = -1.0;
    }
  this->InversionSymmetryTable = new int[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    {
      this->InversionSymmetryTable[i] = inversionList[i];
    }  

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->NbrSite);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
//   cout << (this->LargeHilbertSpaceDimension) << endl;
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
//     for (int i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//       this->PrintState(cout, i) << endl;
}



// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation::Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation (char* fileName, unsigned long memory)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->HilbertSpaceDimension = 0;
      return;
    }
  File.seekg (0l, ios::end);
  unsigned long FileSize = File.tellg ();
  File.close();

  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->HilbertSpaceDimension = 0;
      return;
    }
 
  ReadLittleEndian(File, this->HilbertSpaceDimension);
  ReadLittleEndian(File, this->LargeHilbertSpaceDimension);  
  ReadLittleEndian(File, this->NbrSite);
  ReadLittleEndian(File, this->ChainLength);
  ReadLittleEndian(File, this->MaxXMomentum);
  ReadLittleEndian(File, this->MaxYMomentum);
  ReadLittleEndian(File, this->Sz);
  ReadLittleEndian(File, this->XMomentum);
  ReadLittleEndian(File, this->YMomentum);
  ReadLittleEndian(File, this->SzSymmetrySector);
  ReadLittleEndian(File, this->SzSymmetryMask);
  ReadLittleEndian(File, this->FixedQuantumNumberFlag);
  ReadLittleEndian(File, this->StateXShift);
  ReadLittleEndian(File, this->ComplementaryStateXShift);
  ReadLittleEndian(File, this->XMomentumMask);
  ReadLittleEndian(File, this->NbrYMomentumBlocks);
  ReadLittleEndian(File, this->StateYShift);
  ReadLittleEndian(File, this->YMomentumBlockSize);
  ReadLittleEndian(File, this->ComplementaryStateYShift);
  ReadLittleEndian(File, this->YMomentumMask);
  ReadLittleEndian(File, this->YMomentumBlockMask);
  ReadLittleEndian(File, this->YMomentumFullMask);
  ReadLittleEndian(File, this->ComplementaryYMomentumFullMask);
  ReadLittleEndian(File, this->InversionSymmetrySector);
  this->InversionSymmetryTable = new int[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    ReadLittleEndian(File, this->InversionSymmetryTable[i]);
  
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [this->LargeHilbertSpaceDimension];
  if (this->HilbertSpaceDimension != 0)
  {
    for (int i = 0; i < this->HilbertSpaceDimension; ++i)
      ReadLittleEndian(File, this->StateDescription[i]);
    for (int i = 0; i < this->HilbertSpaceDimension; ++i)
      ReadLittleEndian(File, this->NbrStateInOrbit[i]);
  }
  else
  {
    for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
      ReadLittleEndian(File, this->StateDescription[i]);
    for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
      ReadLittleEndian(File, this->NbrStateInOrbit[i]);
  }
//   FileSize -= this->LargeHilbertSpaceDimension * sizeof (unsigned long);
//   if (FileSize == 0ul)
//     {
// 
//       int TmpLzMax = this->LzMax;
//       this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
//       for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
// 	{
// 	  unsigned long TmpState = this->StateDescription[i];
// 	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
// 	    --TmpLzMax;
// 	  this->StateLzMax[i] = TmpLzMax;
// 	}
//    }
//   else
//     {
//       this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
//       for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
// 	ReadLittleEndian(File, this->StateLzMax[i]);
//     }
  File.close();


  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
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

Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation::Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation (const Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation& chain)
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

      this->Sz = chain.Sz;
      this->SzSymmetrySector = chain.SzSymmetrySector;
      this->SzSymmetryMask = chain.SzSymmetryMask;

      this->InversionSymmetrySector = chain.InversionSymmetrySector;
      this->InversionSymmetryTable = new int [this->NbrSite];
      for (int i = 0; i < this->NbrSite; ++i)
	{
	  this->InversionSymmetryTable[i] = chain.InversionSymmetryTable[i];
	}

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->LookUpTableShift = chain.LookUpTableShift;
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
      this->LookUpTableShift = 0;
      this->SzSymmetrySector = 0.0;
      this->SzSymmetryMask = 0x0ul;
      this->InversionSymmetrySector = 0.0;
      this->InversionSymmetryTable = 0;
    }
}

// destructor
//

Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation::~Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation& Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation::operator = (const Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation& chain)
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

      this->Sz = chain.Sz;
      this->SzSymmetrySector = chain.SzSymmetrySector;
      this->SzSymmetryMask = chain.SzSymmetryMask;

      this->InversionSymmetrySector = chain.InversionSymmetrySector;
      this->InversionSymmetryTable = new int [this->NbrSite];
      for (int i = 0; i < this->NbrSite; ++i)
	{
	  this->InversionSymmetryTable[i] = chain.InversionSymmetryTable[i];
	}

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->LookUpTableShift = chain.LookUpTableShift;
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
      this->LookUpTableShift = 0;
      this->SzSymmetrySector = 0.0;
      this->SzSymmetryMask = 0x0ul;
      this->InversionSymmetrySector = 0.0;
      this->InversionSymmetryTable = 0;
   }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation::Clone()
{
  return new Spin1_2ChainNewSzSymmetryGenericInversionAnd2DTranslation (*this);
}

