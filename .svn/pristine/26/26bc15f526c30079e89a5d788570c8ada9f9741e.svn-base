////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 1/2 chain with Sz contraint                //
//                  and 2d translations plus Sz<->-Sz symmetry                //
//                                                                            //
//                        last modification : 24/07/2018                      //
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


#include "HilbertSpace/Spin3_2ChainSzSymmetryAnd2DTranslation.h"
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

Spin3_2ChainSzSymmetryAnd2DTranslation::Spin3_2ChainSzSymmetryAnd2DTranslation ()
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

Spin3_2ChainSzSymmetryAnd2DTranslation::Spin3_2ChainSzSymmetryAnd2DTranslation (int nbrSite, int sz, int szSymmetrySector, int xMomentum, int maxXMomentum, 
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
  this->SzSymmetryMask = (0x1ul << (2 * this->ChainLength)) - 0x1ul;
  this->Sz = 0;
  this->FixedQuantumNumberFlag = false;
  
  this->StateXShift = 2 * (this->NbrSite / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this->NbrSite) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->NbrYMomentumBlocks = (2 * this->NbrSite) / this->StateXShift;
  this->StateYShift = 2 * (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
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

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->NbrSite);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
//   cout << (this->LargeHilbertSpaceDimension) << endl;
  if (this->LargeHilbertSpaceDimension > 0)
    {
      this->GenerateLookUpTable(memory);
      this->ComputeRescalingFactors();
      
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

Spin3_2ChainSzSymmetryAnd2DTranslation::Spin3_2ChainSzSymmetryAnd2DTranslation (char* fileName, unsigned long memory)
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
  File.close();
  

  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0)
    {
      this->GenerateLookUpTable(memory);
      this->ComputeRescalingFactors();
      
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

Spin3_2ChainSzSymmetryAnd2DTranslation::Spin3_2ChainSzSymmetryAnd2DTranslation (const Spin3_2ChainSzSymmetryAnd2DTranslation& chain)
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
      this->SzSymmetryMask = 0.0;
    }
}

// destructor
//

Spin3_2ChainSzSymmetryAnd2DTranslation::~Spin3_2ChainSzSymmetryAnd2DTranslation () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin3_2ChainSzSymmetryAnd2DTranslation& Spin3_2ChainSzSymmetryAnd2DTranslation::operator = (const Spin3_2ChainSzSymmetryAnd2DTranslation& chain)
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
      this->SzSymmetryMask = 0.0;
   }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin3_2ChainSzSymmetryAnd2DTranslation::Clone()
{
  return new Spin3_2ChainSzSymmetryAnd2DTranslation (*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool Spin3_2ChainSzSymmetryAnd2DTranslation::WriteHilbertSpace (char* fileName)
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
  WriteLittleEndian(File, this->NbrSite);
  WriteLittleEndian(File, this->ChainLength);
  WriteLittleEndian(File, this->MaxXMomentum);
  WriteLittleEndian(File, this->MaxYMomentum);
  WriteLittleEndian(File, this->Sz);
  WriteLittleEndian(File, this->XMomentum);
  WriteLittleEndian(File, this->YMomentum);
  WriteLittleEndian(File, this->SzSymmetrySector);
  WriteLittleEndian(File, this->SzSymmetryMask);
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

long Spin3_2ChainSzSymmetryAnd2DTranslation::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->NbrSite - 1, 0);
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

ComplexVector Spin3_2ChainSzSymmetryAnd2DTranslation::ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  Spin3_2Chain* TmpSpace = (Spin3_2Chain*) space;
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

void Spin3_2ChainSzSymmetryAnd2DTranslation::ComputeRescalingFactors()
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

