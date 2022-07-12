////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                    class author: Cecile Repellin                           //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a Sz = 0                        //
//        a pseudospin 1/2 (not a conserved quantity), spin inversion         //
//                     symmetry and 2D translations                           //
//                                                                            //
//                        last modification : 11/12/2016                      //
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


#include "HilbertSpace/Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainNewSzSymmetryAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainWithPseudospin.h"
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

Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation ()
{
  this->SzSymmetrySector = 1.0;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// inversionSector = inversion sector (can be either +1 or -1)
// xMomentum = momentum along the x direction
// maxXMomentum = number of sites in the x direction
// yMomentum = momentum along the y direction
// maxYMomentum = number of sites in the y direction
// memory = amount of memory granted for precalculations

Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation (int nbrSite, int sz, int szSymmetrySector, int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, unsigned long memory) 
{
  this->Flag.Initialize();
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
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
    
  this->SzSymmetryMask = 0x0ul;
  for (int i = 0; i < this->ChainLength; ++i)
    this->SzSymmetryMask |= (0x1ul << (2*i + 1));
  unsigned long TmpMask= (0x1ul << (2*this->ChainLength)) - 0x1ul;
  this->SzSymmetryComplementaryMask = (~this->SzSymmetryMask) & TmpMask;
  
  
//   cout << (this->SzSymmetryMask) << " " << (this->SzSymmetryComplementaryMask) << endl;
  
  this->Sz = 0;
  this->FixedQuantumNumberFlag = false;
  
  this->StateXShift = 2 * this->ChainLength / this->MaxXMomentum;
  this->ComplementaryStateXShift = 2 * this->ChainLength- this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->NbrYMomentumBlocks = (2 * this->ChainLength) / this->StateXShift;
  this->StateYShift = ((2 * this->ChainLength) / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 
  
  this->FillMirrorSymmetryTable();

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << (this->LargeHilbertSpaceDimension) << endl;
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
      UsedMemory = this->ChainLength * sizeof(int);
      UsedMemory += this->ChainLength * this->LookUpTableMemorySize * sizeof(int);
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

Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation (char* fileName, unsigned long memory)
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
  ReadLittleEndian(File, this->ChainLength);
  ReadLittleEndian(File, this->ChainLength);
  ReadLittleEndian(File, this->MaxXMomentum);
  ReadLittleEndian(File, this->MaxYMomentum);
  ReadLittleEndian(File, this->Sz);
  ReadLittleEndian(File, this->XMomentum);
  ReadLittleEndian(File, this->YMomentum);
  ReadLittleEndian(File, this->SzSymmetrySector);
  ReadLittleEndian(File, this->SzSymmetryMask);
  ReadLittleEndian(File, this->SzSymmetryComplementaryMask);
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
  ReadLittleEndian(File, this->MirrorTransformationTable);
  
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
      UsedMemory = this->ChainLength * sizeof(int);
      UsedMemory += this->ChainLength * this->LookUpTableMemorySize * sizeof(int);
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

Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation (const Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation& chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
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
      this->SzSymmetryComplementaryMask = chain.SzSymmetryComplementaryMask;
      
      this->MirrorTransformationTable = chain.MirrorTransformationTable;

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
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->LookUpTableSize = 0;
      this->SzSymmetrySector = 0.0;
      this->SzSymmetryMask = 0.0;
      this->MirrorTransformationTable = 0;
    }
}

// destructor
//

Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::~Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation () 
{
  delete[] this->MirrorTransformationTable;
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation& Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::operator = (const Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation& chain)
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
      this->SzSymmetryComplementaryMask = chain.SzSymmetryComplementaryMask;
      
      this->MirrorTransformationTable = chain.MirrorTransformationTable;

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
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->SzSymmetrySector = 0.0;
      this->SzSymmetryMask = 0.0;
      this->MirrorTransformationTable = 0;
   }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::Clone()
{
  return new Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation (*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->ChainLength);
  WriteLittleEndian(File, this->MaxXMomentum);
  WriteLittleEndian(File, this->MaxYMomentum);
  WriteLittleEndian(File, this->Sz);
  WriteLittleEndian(File, this->XMomentum);
  WriteLittleEndian(File, this->YMomentum);
  WriteLittleEndian(File, this->SzSymmetrySector);
  WriteLittleEndian(File, this->SzSymmetryMask);
  WriteLittleEndian(File, this->SzSymmetryComplementaryMask);
  WriteLittleEndian(File, this->FixedQuantumNumberFlag);
  WriteLittleEndian(File, this->StateXShift);
  WriteLittleEndian(File, this->ComplementaryStateXShift);
  WriteLittleEndian(File, this->XMomentumMask);
  WriteLittleEndian(File, this->NbrYMomentumBlocks);
  WriteLittleEndian(File, this->StateYShift);
  WriteLittleEndian(File, this->YMomentumBlockSize);
  WriteLittleEndian(File, this->ComplementaryStateYShift);
  WriteLittleEndian(File, this->YMomentumMask);
  WriteLittleEndian(File, this->YMomentumBlockMask);
  WriteLittleEndian(File, this->YMomentumFullMask);
  WriteLittleEndian(File, this->ComplementaryYMomentumFullMask);
  WriteLittleEndian(File, this->MirrorTransformationTable);
    
  if (this->HilbertSpaceDimension != 0)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->NbrStateInOrbit[i]);
    }
  else
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->NbrStateInOrbit[i]);
    }
  File.close();
  return true;
}



// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::GenerateStates()
{
  cout << "Intermediary Hilbert space dimension = " << (this->LargeHilbertSpaceDimension) << endl;
  this->HilbertSpaceDimension = this->LargeHilbertSpaceDimension;
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates((this->Sz + this->ChainLength) >> 1, this->ChainLength - 1, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  double TmpSign;
#ifdef  __64_BITS__
  unsigned long Discard = 0xfffffffffffffffful;
#else
  unsigned long Discard = 0xfffffffful;
#endif
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY, TmpSign) == this->StateDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
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

ComplexVector Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  Spin1_2ChainWithPseudospin* TmpSpace = (Spin1_2ChainWithPseudospin*) space;
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
      
      int Pos = this->FindStateIndex(TmpState);
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

void Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::ComputeRescalingFactors()
{
  int Tmp = 2 * this->ChainLength;
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


// convert a state from a SU(2) basis to another one, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, RealMatrix oneBodyBasis, AbstractSpinChain* space, 
							 long firstComponent, long nbrComponents)
{
  Spin1_2ChainNewSzSymmetryAnd2DTranslation* TmpSpace = (Spin1_2ChainNewSzSymmetryAnd2DTranslation*) space;
  int* TmpSpinIndices = new int [this->ChainLength];
  int* TmpPseudospinIndices = new int [this->ChainLength];
  int* TmpPseudospinIndices2 = new int [this->ChainLength];
  targetState.ClearVector();
  int NbrStateInOrbit;
  
  Complex** ExponentialFactors = new Complex*[this->MaxXMomentum];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      ExponentialFactors[i] = new Complex[this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{ 
	  ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->MaxXMomentum))
							       + (this->YMomentum * ((double) j) / ((double) this->MaxYMomentum))));
	}
    }

  for (long i = 0; i < (TmpSpace->HilbertSpaceDimension); ++i)
    {
      unsigned long TmpState = TmpSpace->StateDescription[i];
      unsigned long Tmp;
      unsigned long Tmp1;
      int TmpIndex = 0;
      bool flag = true;
      int j;
      while ((TmpIndex < this->ChainLength) && (flag = true))
	{
	  j = this->ChainLength - TmpIndex - 1;
	  Tmp = (TmpState >> (3 * TmpIndex)) & 0x7ul;
	  Tmp1 = (Tmp & 0x1ul) + ((Tmp >> 1) & 0x1ul) + ((Tmp >> 2) & 0x1ul);
// 	  cout << TmpState << " " << (3*j) << " " << (TmpState >> (3 * j)) << " " << Tmp << " " << Tmp1 << endl;
	  if ((Tmp == 0x7ul) || (Tmp == 0x0ul))
	  {
	    flag = false;
	    TmpIndex = this->ChainLength;
	  }
	  else
	  {
	    if (Tmp1 == 0x1ul)
	      TmpSpinIndices[TmpIndex] = 0;
	    if (Tmp1 == 0x2ul)
	      TmpSpinIndices[TmpIndex] = 1;
	    
	    if ((Tmp & 0x1ul) == ((Tmp >>1) & 0x1ul))
	      TmpPseudospinIndices[TmpIndex] = 0;
	    if ((Tmp & 0x1ul) == ((Tmp >>2) & 0x1ul))
	      TmpPseudospinIndices[TmpIndex] = 1;
	    if (((Tmp >>1) & 0x1ul) == ((Tmp >>2) & 0x1ul))
	      TmpPseudospinIndices[TmpIndex] = 2;
	    
	    ++TmpIndex;
	  }
	  
	  
// 	  cout << TmpState << " " << j << " " << TmpIndex << " " << flag << endl;
	}

      if (flag == true)
      {
// 	for (int k = 0; k < this->ChainLength; ++k)
// 	  cout <<  TmpPseudospinIndices[k] << endl;
	NbrStateInOrbit = TmpSpace->GetNbrStateinOrbit(i);
	this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpSpinIndices, TmpPseudospinIndices, TmpPseudospinIndices2, oneBodyBasis, NbrStateInOrbit, ExponentialFactors);
      }
//     cout << endl;
    }
  delete[] TmpSpinIndices;
  delete[] TmpPseudospinIndices;
  delete[] TmpPseudospinIndices2;
  for (int i = 0; i < this->MaxXMomentum; ++i)
    delete[] ExponentialFactors[i];
  delete[] ExponentialFactors;
}

// recursive part of the convertion from a state from a SU(2) basis to another one, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSpinIndices = array that gives the spin dressing the initial n-body state
// currentSpinIndices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector

void Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
								int position, int* spinIndices, int* initialPseudospinIndices, int* currentPseudospinIndices, RealMatrix oneBodyBasis, int nbrStateInOrbit, Complex** ExponentialFactors) 
{
//   cout << position << " : " << endl;
//   for (int i = 0; i < position; ++i)
//     cout << currentSpinIndices[i] << " ";
//   cout << endl;
  if (position == this->ChainLength)
    {
      unsigned long TmpState = 0x0ul;
      unsigned long Mask = 0x0ul;
      for (int i = 0; i < this->ChainLength; ++i)
	{
	  Mask =  (((currentPseudospinIndices[i]) << (2*i)) | ((spinIndices[i]) << (2*i + 1))); // Mask = 00...0100...0 : one fermion state in the second quantized basis
// 	  cout << i << " mask = " << ((currentPseudospinIndices[i]) << (2*i)) << " " << ((spinIndices[i]) << (2*i + 1)) << " " << Mask << endl;
	  if ((TmpState & Mask) != 0x0ul)
	    return;
	  TmpState |= Mask; //set bit corresponding to the current fermion state to 1 in TmpState
	}
      int nbrTranslationX;
      int nbrTranslationY;
      double TmpCoefficient = 1.0;
      int Index = this->SymmetrizeResult(TmpState, nbrStateInOrbit, TmpCoefficient, nbrTranslationX, nbrTranslationY);
//       cout << TmpCoefficient << endl;
      if (Index < this->HilbertSpaceDimension)
	{
	  targetState[Index] += coefficient * TmpCoefficient * ExponentialFactors[nbrTranslationX][nbrTranslationY];
	}
      return;      
    }
  else
    {
      currentPseudospinIndices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[0][initialPseudospinIndices[position]]), position + 1, spinIndices, initialPseudospinIndices, currentPseudospinIndices, oneBodyBasis, nbrStateInOrbit, ExponentialFactors);
      currentPseudospinIndices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[1][initialPseudospinIndices[position]]), position + 1, spinIndices, initialPseudospinIndices, currentPseudospinIndices, oneBodyBasis, nbrStateInOrbit, ExponentialFactors);
    }
}
