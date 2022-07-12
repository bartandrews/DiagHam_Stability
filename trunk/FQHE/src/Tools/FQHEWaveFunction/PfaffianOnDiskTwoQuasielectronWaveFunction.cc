////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of Pfaffian wave function with two quasiholes on disk          //
//                                                                            //
//                        last modification : 23/10/2008                      //
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
#include "Tools/FQHEWaveFunction/PfaffianOnDiskTwoQuasielectronWaveFunction.h"
#include "Vector/RealVector.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "GeneralTools/Endian.h"
#include "MathTools/BinomialCoefficients.h"

#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// group maximum size in bits
#ifdef __64_BITS__
#define GROUP_MAXSIZE 12
#else
#define GROUP_MAXSIZE 6
#endif


// constructor
//
// nbrParticles = number of particles
// zElectron1 = position of the first quasielectron
// zElectron2 = position of the second quasielectron (spherical coordinates, theta angle)
// fermions = flag indicating whether to calculate bosonic or fermionic pfaffian

PfaffianOnDiskTwoQuasielectronWaveFunction::PfaffianOnDiskTwoQuasielectronWaveFunction(int nbrParticles,  Complex zElectron1, Complex zElectron2, bool fermions)
{
  this->NbrParticles = nbrParticles;
  this->ZElectron1 = zElectron1;
  this->ZElectron2 = zElectron2;
  this->GaussianWeight = exp (-0.125 * (SqrNorm(this->ZElectron1) + SqrNorm(this->ZElectron2)));
  this->InvScale = 1.0 / 10.0;

  this->TmpPfaffian = new Complex* [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpPfaffian[i] = new Complex [this->NbrParticles];
  this->TmpSqrPfaffian = new Complex* [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpSqrPfaffian[i] = new Complex [this->NbrParticles];
  this->TmpIndexArray = new int [this->NbrParticles];
  this->TmpWeights1 = new Complex [this->NbrParticles];
  this->TmpWeights2 = new Complex [this->NbrParticles];
  this->NextCoordinate = -1;
  this->FermionFlag = fermions;

  this->EvaluatePermutations();
  this->CurrentWaveFunctionPart1 = new Complex [this->NbrPermutations];
  this->CurrentWaveFunctionPart2 = new Complex [this->NbrPermutations];
  this->PreviousWaveFunctionPart1 = new Complex [this->NbrPermutations];
  this->PreviousWaveFunctionPart2 = new Complex [this->NbrPermutations];
  this->Flag.Initialize();
}


// constructor from data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// zElectron1 = position of the first quasielectron
// zElectron2 = position of the second quasielectron (spherical coordinates, theta angle)
// fermions = flag indicating whether to calculate bosonic or fermionic pfaffian

PfaffianOnDiskTwoQuasielectronWaveFunction::PfaffianOnDiskTwoQuasielectronWaveFunction(char* filename, Complex zElectron1, Complex zElectron2, bool fermions)
{
  this->NbrParticles = 0;
  this->ReadPermutations(filename);
  this->ZElectron1 = zElectron1;
  this->ZElectron2 = zElectron2;
  this->GaussianWeight = exp (-0.125 * (SqrNorm(this->ZElectron1) + SqrNorm(this->ZElectron2)));
  this->InvScale = 1.0 / 10.0;

  this->TmpPfaffian = new Complex* [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpPfaffian[i] = new Complex [this->NbrParticles];
  this->TmpSqrPfaffian = new Complex* [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpSqrPfaffian[i] = new Complex [this->NbrParticles];
  this->TmpIndexArray = new int [this->NbrParticles];
  this->TmpWeights1 = new Complex [this->NbrParticles];
  this->TmpWeights2 = new Complex [this->NbrParticles];
  this->NextCoordinate = -1;
  this->CurrentWaveFunctionPart1 = new Complex [this->NbrPermutations];
  this->CurrentWaveFunctionPart2 = new Complex [this->NbrPermutations];
  this->PreviousWaveFunctionPart1 = new Complex [this->NbrPermutations];
  this->PreviousWaveFunctionPart2 = new Complex [this->NbrPermutations];

  this->FermionFlag = fermions;
  this->Flag.Initialize();
}


// copy constructor
//
// function = reference on the wave function to copy

PfaffianOnDiskTwoQuasielectronWaveFunction::PfaffianOnDiskTwoQuasielectronWaveFunction(const PfaffianOnDiskTwoQuasielectronWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->ZElectron1 = function.ZElectron1;
  this->ZElectron2 = function.ZElectron2;
  this->GaussianWeight = exp (-0.125 * (SqrNorm(this->ZElectron1) + SqrNorm(this->ZElectron2)));
  this->InvScale = function.InvScale;

  this->TmpPfaffian = new Complex* [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpPfaffian[i] = new Complex [this->NbrParticles];
  this->TmpSqrPfaffian = new Complex* [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpSqrPfaffian[i] = new Complex [this->NbrParticles];
  this->TmpIndexArray = new int [this->NbrParticles];
  this->TmpWeights1 = new Complex [this->NbrParticles];
  this->TmpWeights2 = new Complex [this->NbrParticles];
  this->NextCoordinate = -1;

  this->FermionFlag=function.FermionFlag;

  this->NbrPermutations = function.NbrPermutations;
  this->Permutations1= function.Permutations1;
  this->Permutations2= function.Permutations2;
  this->Flag = function.Flag;

  this->CurrentWaveFunctionPart1 = new Complex [this->NbrPermutations];
  this->CurrentWaveFunctionPart2 = new Complex [this->NbrPermutations];
  this->PreviousWaveFunctionPart1 = new Complex [this->NbrPermutations];
  this->PreviousWaveFunctionPart2 = new Complex [this->NbrPermutations];
}

// destructor
//

PfaffianOnDiskTwoQuasielectronWaveFunction::~PfaffianOnDiskTwoQuasielectronWaveFunction()
{
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      delete[] this->TmpPfaffian[i];
      delete[] this->TmpSqrPfaffian[i];
    }
  delete[] this->TmpPfaffian;
  delete[] this->TmpSqrPfaffian;
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete[] this->Permutations1;
      delete[] this->Permutations2;
    }
  delete[] this->TmpIndexArray;
  delete[] this->TmpWeights1;
  delete[] this->TmpWeights2;
  delete[] this->CurrentWaveFunctionPart1;
  delete[] this->CurrentWaveFunctionPart2;
  delete[] this->PreviousWaveFunctionPart1;
  delete[] this->PreviousWaveFunctionPart2;  
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianOnDiskTwoQuasielectronWaveFunction::Clone ()
{
  return new PfaffianOnDiskTwoQuasielectronWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianOnDiskTwoQuasielectronWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp, T1, T2;
  Complex WaveFunction(1.0);
  Complex TmpZ;
  Complex TmpZElectron1;
  Complex TmpZElectron2;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpZ.Re = x[i << 1];
      TmpZ.Im = x[1 + (i << 1)];
      this->TmpWeights1[i] = exp (0.25 * (Conj(this->ZElectron1) * TmpZ));
      this->TmpWeights2[i] = exp (0.25 * (Conj(this->ZElectron2) * TmpZ));
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = TmpZ.Re - x[j << 1];
	  Tmp.Im = TmpZ.Im - x[1 + (j << 1)];
	  Tmp *= this->InvScale;
	  this->TmpPfaffian[i][j] = Tmp;
	  this->TmpPfaffian[j][i] = -Tmp;
	  this->TmpSqrPfaffian[i][j] = Tmp * Tmp;
	  this->TmpSqrPfaffian[j][i] = this->TmpSqrPfaffian[i][j];
	  WaveFunction *= Tmp;
	}
    }

  if (this->FermionFlag == false)
    WaveFunction = 1.0;

  int NbrParticlesPerColor = this->NbrParticles >> 1;
  int Shift = ((NbrParticlesPerColor - 1) * 5);
  Complex WaveFunctionSum = 0.0;
  if (this->NextCoordinate == -1)
    for (unsigned long i = 0ul; i < this->NbrPermutations; ++i)
      {
	unsigned long TmpPerm = this->Permutations1[i];
	unsigned long TmpPerm2 = this->Permutations2[i];
	Complex WaveFunctionPart1 = 0.0;
	for (int j = 0; j < NbrParticlesPerColor; ++ j)	
	  {
	    Complex WaveFunctionPart12 = 0.0;	  
	    Complex WaveFunctionPart12b = 1.0;	  
	    for (int k = 0; k < NbrParticlesPerColor; ++k)
	      TmpIndexArray[k] = (TmpPerm >> (k * 5)) & 0x1ful;
	    TmpPerm = (TmpPerm  >> 5) | ((TmpPerm & 0x1ful) << Shift);
	    
	    Complex* TmpArray = this->TmpPfaffian[TmpIndexArray[0]];
	    for (int k = 1; k < NbrParticlesPerColor; ++k)
	      {
		Tmp = 1.0; 
		for (int l = 1; l < k; ++l)
		  Tmp *= TmpArray[TmpIndexArray[l]];
		Tmp *= this->TmpWeights1[TmpIndexArray[k]];
		Complex* TmpArraySqr = this->TmpSqrPfaffian[TmpIndexArray[k]];
		for (int l = k + 1; l < NbrParticlesPerColor; ++l)
		  {
		    Tmp *= TmpArray[TmpIndexArray[l]];
		    WaveFunctionPart12b *= TmpArraySqr[TmpIndexArray[l]];
		  }
		WaveFunctionPart12 += Tmp;
	      }
	    WaveFunctionPart1 += WaveFunctionPart12 * WaveFunctionPart12b;
	  }
	this->CurrentWaveFunctionPart1[i] = WaveFunctionPart1;
	
	Complex WaveFunctionPart2 = 0.0;
	for (int j = 0; j < NbrParticlesPerColor; ++ j)
	  {
	    Complex WaveFunctionPart22 = 0.0;	  
	    Complex WaveFunctionPart22b = 1.0;	  
	    for (int k = 0; k < NbrParticlesPerColor; ++k)
	      TmpIndexArray[k] = (TmpPerm2 >> (k * 5)) & 0x1ful;
	    TmpPerm2 = (TmpPerm2  >> 5) | ((TmpPerm2 & 0x1ful) << Shift);
	    Complex* TmpArray = this->TmpPfaffian[TmpIndexArray[0]];
	    for (int k = 1; k < NbrParticlesPerColor; ++k)
	      {
		Tmp = 1.0; 
		for (int l = 1; l < k; ++l)
		  Tmp *= TmpArray[TmpIndexArray[l]];	      
		Tmp *= this->TmpWeights2[TmpIndexArray[k]];
		Complex* TmpArraySqr = this->TmpSqrPfaffian[TmpIndexArray[k]];
		for (int l = k + 1; l < NbrParticlesPerColor; ++l)
		  {
		    Tmp *= TmpArray[TmpIndexArray[l]];	      
		    WaveFunctionPart22b *=  TmpArraySqr[TmpIndexArray[l]];
		  }
		WaveFunctionPart22 += Tmp;
	      }
	    WaveFunctionPart2 += WaveFunctionPart22 * WaveFunctionPart22b;
	  }      
	this->CurrentWaveFunctionPart2[i] = WaveFunctionPart2;
	
	WaveFunctionSum += WaveFunctionPart1 * WaveFunctionPart2;
      }
  else
    for (unsigned long i = 0ul; i < this->NbrPermutations; ++i)
      {
	unsigned long TmpPerm = this->Permutations1[i];
	unsigned long TmpPerm2 = this->Permutations2[i];
	Complex WaveFunctionPart1 = 0.0;

	for (int j = 0; j < NbrParticlesPerColor; ++j)	
	  if (((TmpPerm >> (5 * j)) & 0x1ful)  == ((unsigned long) this->NextCoordinate))
	    TmpPerm = 0ul;
	if (TmpPerm != 0ul)
	  for (int j = 0; j < NbrParticlesPerColor; ++j)	
	    {
	      Complex WaveFunctionPart12 = 0.0;	  
	      Complex WaveFunctionPart12b = 1.0;	  
	      for (int k = 0; k < NbrParticlesPerColor; ++k)
		TmpIndexArray[k] = (TmpPerm >> (k * 5)) & 0x1ful;
	      TmpPerm = (TmpPerm  >> 5) | ((TmpPerm & 0x1ful) << Shift);
	      
	      Complex* TmpArray = this->TmpPfaffian[TmpIndexArray[0]];
	      for (int k = 1; k < NbrParticlesPerColor; ++k)
		{
		  Tmp = 1.0; 
		  for (int l = 1; l < k; ++l)
		    Tmp *= TmpArray[TmpIndexArray[l]];
		  Tmp *= this->TmpWeights1[TmpIndexArray[k]];
		  Complex* TmpArraySqr = this->TmpSqrPfaffian[TmpIndexArray[k]];
		  for (int l = k + 1; l < NbrParticlesPerColor; ++l)
		    {
		      Tmp *= TmpArray[TmpIndexArray[l]];
		      WaveFunctionPart12b *= TmpArraySqr[TmpIndexArray[l]];
		    }
		  WaveFunctionPart12 += Tmp;
		}
	      WaveFunctionPart1 += WaveFunctionPart12 * WaveFunctionPart12b;
	    }
	else
	  WaveFunctionPart1 = this->PreviousWaveFunctionPart1[i];
	this->CurrentWaveFunctionPart1[i] = WaveFunctionPart1;

	Complex WaveFunctionPart2 = 0.0;
	if (TmpPerm == 0ul)
	  for (int j = 0; j < NbrParticlesPerColor; ++ j)
	    {
	      Complex WaveFunctionPart22 = 0.0;	  
	      Complex WaveFunctionPart22b = 1.0;	  
	      for (int k = 0; k < NbrParticlesPerColor; ++k)
		TmpIndexArray[k] = (TmpPerm2 >> (k * 5)) & 0x1ful;
	      TmpPerm2 = (TmpPerm2  >> 5) | ((TmpPerm2 & 0x1ful) << Shift);
	      Complex* TmpArray = this->TmpPfaffian[TmpIndexArray[0]];
	      for (int k = 1; k < NbrParticlesPerColor; ++k)
		{
		  Tmp = 1.0; 
		  for (int l = 1; l < k; ++l)
		    Tmp *= TmpArray[TmpIndexArray[l]];	      
		  Tmp *= this->TmpWeights2[TmpIndexArray[k]];
		  Complex* TmpArraySqr = this->TmpSqrPfaffian[TmpIndexArray[k]];
		  for (int l = k + 1; l < NbrParticlesPerColor; ++l)
		    {
		      Tmp *= TmpArray[TmpIndexArray[l]];	      
		      WaveFunctionPart22b *=  TmpArraySqr[TmpIndexArray[l]];
		    }
		  WaveFunctionPart22 += Tmp;
		}
	      WaveFunctionPart2 += WaveFunctionPart22 * WaveFunctionPart22b;
	    }      
	else
	  WaveFunctionPart2 = this->PreviousWaveFunctionPart2[i];
	this->CurrentWaveFunctionPart2[i] = WaveFunctionPart2;
	
	WaveFunctionSum += WaveFunctionPart1 * WaveFunctionPart2;
      }
  WaveFunction *= WaveFunctionSum;
  WaveFunction *= this->GaussianWeight;
  return WaveFunction;
}

// indicate that only a given coordinate will be changed during the next wave function evaluation
//
// coordinate = coordinate index (-1 to reset the wave function memory)

void PfaffianOnDiskTwoQuasielectronWaveFunction::SetNextCoordinate(int coordinate)
{
  this->NextCoordinate = coordinate;
  Complex* Tmp = this->CurrentWaveFunctionPart1;
  this->CurrentWaveFunctionPart1 = this->PreviousWaveFunctionPart1;
  this->PreviousWaveFunctionPart1 = Tmp;
  Tmp = this->CurrentWaveFunctionPart2;
  this->CurrentWaveFunctionPart2 = this->PreviousWaveFunctionPart2;
  this->PreviousWaveFunctionPart2 = Tmp;
}

// cancel data that have been modified during the last wave function evaluation
//

void PfaffianOnDiskTwoQuasielectronWaveFunction::RestorePreviousData()
{
  Complex* Tmp = this->CurrentWaveFunctionPart1;
  this->CurrentWaveFunctionPart1 = this->PreviousWaveFunctionPart1;
  this->PreviousWaveFunctionPart1 = Tmp;
  Tmp = this->CurrentWaveFunctionPart2;
  this->CurrentWaveFunctionPart2 = this->PreviousWaveFunctionPart2;
  this->PreviousWaveFunctionPart2 = Tmp;
}

// evaluate all permutations requested to symmetrize the state
//

void PfaffianOnDiskTwoQuasielectronWaveFunction::EvaluatePermutations()
{
  BinomialCoefficients Binomial(this->NbrParticles);
  int NbrParticlesPerColor = this->NbrParticles >> 1;
  this->NbrPermutations = Binomial(this->NbrParticles, NbrParticlesPerColor);
  this->Permutations1 = new unsigned long[this->NbrPermutations];
  this->Permutations2 = new unsigned long[this->NbrPermutations];
  unsigned long MinValue = (0x1ul << NbrParticlesPerColor) - 0x1ul;
  unsigned long MaxValue = MinValue << (this->NbrParticles - NbrParticlesPerColor);
  unsigned long* TmpArrayPerm = new unsigned long [this->NbrParticles];
  this->NbrPermutations = 0;
  for (; MinValue <= MaxValue; ++MinValue)
    {
      int Count = 0;
      int Pos = 0;
      while ((Pos < this->NbrParticles) && (Count <= NbrParticlesPerColor))
	{
	  if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
	    ++Count;
	  ++Pos;
	}
      if (Count == NbrParticlesPerColor)
	{
	  int Pos1 = 0;
	  int Pos2 = NbrParticlesPerColor;
	  for (Pos = 0; Pos < this->NbrParticles; ++Pos)
	    if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
	      {
		TmpArrayPerm[Pos1] = (unsigned long) Pos;
		++Pos1;
	      }
	    else
	      {
		TmpArrayPerm[Pos2] = (unsigned long) Pos;
		++Pos2;
	      }
	  unsigned long TmpPerm2 = 0ul;
	  unsigned long TmpPerm3 = 0ul;
	  for (int i = 0; i < NbrParticlesPerColor; ++i)
	    {
	      TmpPerm2 |= TmpArrayPerm[i] << (i * 5);
	      TmpPerm3 |= TmpArrayPerm[i + NbrParticlesPerColor] << (i *5);
	    }
	  this->Permutations1[this->NbrPermutations] = TmpPerm2;
	  this->Permutations2[this->NbrPermutations] = TmpPerm3;	      
	  ++this->NbrPermutations;
	}
    }
  delete[] TmpArrayPerm;
  return;
}


// get all permutations requested to symmetrize the SU(K) state from data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool PfaffianOnDiskTwoQuasielectronWaveFunction::ReadPermutations(char* filename)
{
  ifstream File;
  File.open(filename, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << filename << endl;
      return false;
    }
  ReadLittleEndian(File, this->NbrParticles);
  ReadLittleEndian(File, this->NbrPermutations);
  this->Permutations1 = new unsigned long[this->NbrPermutations];
  this->Permutations2 = new unsigned long[this->NbrPermutations];
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    ReadLittleEndian(File, this->Permutations1[i]);
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    ReadLittleEndian(File, this->Permutations2[i]);
  File.close();
  return true;
}

// write all permutations requested to symmetrize the SU(K) state to data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool PfaffianOnDiskTwoQuasielectronWaveFunction::WritePermutations(char* filename)
{
  ofstream File;
  File.open(filename, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << filename << endl;
      return false;
    }
  WriteLittleEndian(File, this->NbrParticles);
  WriteLittleEndian(File, this->NbrPermutations);
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    WriteLittleEndian(File, this->Permutations1[i]);
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    WriteLittleEndian(File, this->Permutations2[i]);
  File.close();
  return true;
}

