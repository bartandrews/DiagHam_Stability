////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of Pfaffian wave function with two quasielectrons on sphere       //
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
#include "Tools/FQHEWaveFunction/PfaffianOnSphereTwoQuasielectronWaveFunction.h"
#include "Vector/RealVector.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
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
// theta1 = position of the first quasielectron (spherical coordinates, theta angle)
// phi1 = position of the first quasielectron (spherical coordinates, phi angle)
// theta2 = position of the second quasielectron (spherical coordinates, theta angle)
// phi2 = position of the second quasielectron (spherical coordinates, phi angle)
// fermions = flag indicating whether to calculate bosonic or fermionic pfaffian

PfaffianOnSphereTwoQuasielectronWaveFunction::PfaffianOnSphereTwoQuasielectronWaveFunction(int nbrParticles, double theta1, double phi1, 
											   double theta2, double phi2, bool fermions)
{
  this->NbrParticles = nbrParticles;

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
  this->TmpPreviousPfaffian = new Complex [this->NbrParticles];
  this->TmpPreviousSqrPfaffian = new Complex [this->NbrParticles];
  this->FermionFlag = fermions;

  this->EvaluatePermutations();
  this->Flag.Initialize();

  this->UElectron1.Re = cos(0.5*phi1);
  this->UElectron1.Im= -sin(0.5*phi1);
  this->UElectron1 *= cos(0.5*theta1);

  this->VElectron1.Re = cos(0.5*phi1);
  this->VElectron1.Im = sin(0.5*phi1);
  this->VElectron1 *= sin(0.5*theta1);
  
  this->UElectron2.Re = cos(0.5*phi2);
  this->UElectron2.Im = -sin(0.5*phi2);
  this->UElectron2 *= cos(0.5*theta2);

  this->VElectron2.Re = cos(0.5*phi2);
  this->VElectron2.Im = sin(0.5*phi2);
  this->VElectron2 *= sin(0.5*theta2);

  this->ConjUElectron1 = Conj(this->UElectron1);
  this->ConjVElectron1 = Conj(this->VElectron1);
  this->ConjUElectron2 = Conj(this->UElectron2);
  this->ConjVElectron2 = Conj(this->VElectron2);

  this->FermionFlag = fermions;

}

// constructor using permutation description stored in a file
//
// filename = pointer to the file name that described the symmetrization procedure
// theta1 = position of the first quasielectron (spherical coordinates, theta angle)
// phi1 = position of the first quasielectron (spherical coordinates, phi angle)
// theta2 = position of the second quasielectron (spherical coordinates, theta angle)
// phi2 = position of the second quasielectron (spherical coordinates, phi angle)
// fermions = flag indicating whether to calculate bosonic or fermionic pfaffian

PfaffianOnSphereTwoQuasielectronWaveFunction::PfaffianOnSphereTwoQuasielectronWaveFunction(char* filename, 
											   double theta1, double phi1, 
											   double theta2, double phi2, bool fermions)
{
  this->NbrParticles = 0;

  this->ReadPermutations(filename);
  this->Flag.Initialize();

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
  this->TmpPreviousPfaffian = new Complex [this->NbrParticles];
  this->TmpPreviousSqrPfaffian = new Complex [this->NbrParticles];
  this->FermionFlag = fermions;

  this->UElectron1.Re = cos(0.5*phi1);
  this->UElectron1.Im= -sin(0.5*phi1);
  this->UElectron1 *= cos(0.5*theta1);

  this->VElectron1.Re = cos(0.5*phi1);
  this->VElectron1.Im = sin(0.5*phi1);
  this->VElectron1 *= sin(0.5*theta1);
  
  this->UElectron2.Re = cos(0.5*phi2);
  this->UElectron2.Im = -sin(0.5*phi2);
  this->UElectron2 *= cos(0.5*theta2);

  this->VElectron2.Re = cos(0.5*phi2);
  this->VElectron2.Im = sin(0.5*phi2);
  this->VElectron2 *= sin(0.5*theta2);

  this->ConjUElectron1 = Conj(this->UElectron1);
  this->ConjVElectron1 = Conj(this->VElectron1);
  this->ConjUElectron2 = Conj(this->UElectron2);
  this->ConjVElectron2 = Conj(this->VElectron2);

  this->FermionFlag = fermions;
}



// copy constructor
//
// function = reference on the wave function to copy

PfaffianOnSphereTwoQuasielectronWaveFunction::PfaffianOnSphereTwoQuasielectronWaveFunction(const PfaffianOnSphereTwoQuasielectronWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;

  this->UElectron1 = function.UElectron1;
  this->VElectron1 = function.VElectron1;
  this->UElectron2 = function.UElectron2;
  this->VElectron2 = function.VElectron2;

  this->ConjUElectron1 = Conj(this->UElectron1);
  this->ConjVElectron1 = Conj(this->VElectron1);
  this->ConjUElectron2 = Conj(this->UElectron2);
  this->ConjVElectron2 = Conj(this->VElectron2);

  this->FermionFlag = function.FermionFlag;
  this->NbrPermutations = function.NbrPermutations;
  this->Permutations1= function.Permutations1;
  this->Permutations2= function.Permutations2;
  this->Flag = function.Flag;

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
  this->TmpPreviousPfaffian = new Complex [this->NbrParticles];
  this->TmpPreviousSqrPfaffian = new Complex [this->NbrParticles];
}

// destructor
//

PfaffianOnSphereTwoQuasielectronWaveFunction::~PfaffianOnSphereTwoQuasielectronWaveFunction()
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
  delete[] this->TmpPreviousPfaffian;
  delete[] this->TmpPreviousSqrPfaffian;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianOnSphereTwoQuasielectronWaveFunction::Clone ()
{
  return new PfaffianOnSphereTwoQuasielectronWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianOnSphereTwoQuasielectronWaveFunction::operator ()(RealVector& x)
{
  Complex WaveFunction(1.0);
  return WaveFunction;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex PfaffianOnSphereTwoQuasielectronWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex Tmp;
  Complex TmpU1;
  Complex TmpV1;
  Complex TmpU2;
  Complex TmpV2;
  Complex WaveFunction(1.0);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpU1 = uv[i << 1];
      TmpV1 = uv[1 + (i << 1)];
      this->TmpWeights1[i] = ((this->ConjUElectron1 * TmpU1) + (this->ConjVElectron1 * TmpV1));
      this->TmpWeights2[i] =  ((this->ConjUElectron2 * TmpU1) + (this->ConjVElectron2 * TmpV1));
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  TmpU2 = uv[j << 1];
	  TmpV2 = uv[1 + (j << 1)];
	  Tmp = (TmpU1 * TmpV2) - (TmpU2 * TmpV1);
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
  unsigned long Mask = (0x1ul << (NbrParticlesPerColor * 5)) - 0x1ul;
  Complex WaveFunctionSum = 0.0;
  for (unsigned long i = 0ul; i < this->NbrPermutations; ++i)
    {
      unsigned long TmpPerm = this->Permutations1[i];
      unsigned long TmpPerm2 = this->Permutations2[i];
      TmpPerm &= Mask;
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

      WaveFunctionSum += WaveFunctionPart1 * WaveFunctionPart2;
    }
  WaveFunction *= WaveFunctionSum;
  return WaveFunction;
}

// evaluate all permutations requested to symmetrize the state
//

void PfaffianOnSphereTwoQuasielectronWaveFunction::EvaluatePermutations()
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


// get all permutations requested to symmetrize the state from data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool PfaffianOnSphereTwoQuasielectronWaveFunction::ReadPermutations(char* filename)
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

// write all permutations requested to symmetrize the state to data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool PfaffianOnSphereTwoQuasielectronWaveFunction::WritePermutations(char* filename)
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

