////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of U(1) wave function obtained form symmetrization of a   //
//                        SU(2) wave function on sphere                       //
//                                                                            //
//                        last modification : 20/10/2008                      //
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
#include "MathTools/BinomialCoefficients.h"
#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSU2ToU1WaveFunction.h"
#include "Tools/FQHEWaveFunction/HalperinOnSphereWaveFunction.h"
#include "Vector/RealVector.h"
#include "GeneralTools/Endian.h"

#include <iostream>
#include <math.h>


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

// default constructor
//

FQHESphereSymmetrizedSU2ToU1WaveFunction::FQHESphereSymmetrizedSU2ToU1WaveFunction()
{
}

// constructor
//
// nbrParticles = number of particles
// nbrParticlesUp = number of particles with spin up
// sU2Wavefunction = pointer to the base SU(2) wave function
// fermionFlag = true if the final state should be a fermionic state
 
FQHESphereSymmetrizedSU2ToU1WaveFunction::FQHESphereSymmetrizedSU2ToU1WaveFunction(int nbrParticles, int nbrParticlesUp, Abstract1DComplexFunctionOnSphere* sU2Wavefunction, bool fermionFlag)
{
  this->NbrParticles = nbrParticles;
  this->NbrParticlesUp = nbrParticlesUp;
  this->NbrParticlesDown = this->NbrParticles - this->NbrParticlesUp;
  this->SU2WaveFunction = (Abstract1DComplexFunctionOnSphere*) (sU2Wavefunction->Clone());
  this->FermionFlag = fermionFlag;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  this->FullySymmetrize = false;  
  this->EvaluatePermutations();
  this->Flag.Initialize();
}

// constructor from data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// sU2Wavefunction = pointer to the base SU(2) wave function
// fermionFlag = true if the final state should be a fermionic state
 
FQHESphereSymmetrizedSU2ToU1WaveFunction::FQHESphereSymmetrizedSU2ToU1WaveFunction(char* filename, Abstract1DComplexFunctionOnSphere* sU2Wavefunction, bool fermionFlag)
{
  this->NbrParticles = 0;
  this->NbrParticlesUp = 0;
  this->NbrParticlesDown = 0;
  this->ReadPermutations(filename);
  this->SU2WaveFunction = (Abstract1DComplexFunctionOnSphere*) (sU2Wavefunction->Clone());
  this->FermionFlag = fermionFlag;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  this->FullySymmetrize = false;  
  this->Flag.Initialize();
}

// copy constructor
//
// function = reference on the wave function to copy

FQHESphereSymmetrizedSU2ToU1WaveFunction::FQHESphereSymmetrizedSU2ToU1WaveFunction(const FQHESphereSymmetrizedSU2ToU1WaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrParticlesUp = function.NbrParticlesUp;
  this->NbrParticlesDown = function.NbrParticlesDown;
  this->Permutations = function.Permutations;
  this->NbrPermutations = function.NbrPermutations;
  this->SU2WaveFunction = (Abstract1DComplexFunctionOnSphere*) (function.SU2WaveFunction->Clone());
  this->FermionFlag = function.FermionFlag;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  this->FullySymmetrize = false;  
  this->Flag = function.Flag;
}

// destructor
//

FQHESphereSymmetrizedSU2ToU1WaveFunction::~FQHESphereSymmetrizedSU2ToU1WaveFunction()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete[] this->Permutations;
    }
  if (this->SU2WaveFunction != 0)
    delete this->SU2WaveFunction;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHESphereSymmetrizedSU2ToU1WaveFunction::Clone ()
{
  return new FQHESphereSymmetrizedSU2ToU1WaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex FQHESphereSymmetrizedSU2ToU1WaveFunction::operator ()(RealVector& x)
{  
  ComplexVector TmpUV (this->NbrParticles * 2);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpUV[2 * i].Re = cos(0.5 * x[i << 1]);
      TmpUV[2 * i].Im = TmpUV[2 * i].Re;
      TmpUV[2 * i].Re *= cos(0.5 * x[1 + (i << 1)]);
      TmpUV[2 * i].Im *= sin(0.5 * x[1 + (i << 1)]);
      TmpUV[(2 * i) + 1].Re = sin(0.5 * x[i << 1]);
      TmpUV[(2 * i) + 1].Im = TmpUV[(2 * i) + 1].Re;
      TmpUV[(2 * i) + 1].Re *= cos(0.5 * x[1 + (i << 1)]);
      TmpUV[(2 * i) + 1].Im *= -sin(0.5 * x[1 + (i << 1)]);
    }
  return this->CalculateFromSpinorVariables(TmpUV);
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex FQHESphereSymmetrizedSU2ToU1WaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{  
  Complex TotalValue = 0.0;
  for (unsigned long i = 0ul; i < this->NbrPermutations; ++i)
    {
      unsigned long Tmp = this->Permutations[i];
      int TmpIndex2 = 0;
      for (int k = 0; k < this->NbrParticles; ++k)
	{
	  unsigned long TmpIndex = ((Tmp >> (k << 2)) & 0xful) << 1;
	  this->TemporaryUV.Re(TmpIndex2) = uv.Re((int) TmpIndex);
	  this->TemporaryUV.Im(TmpIndex2) = uv.Im((int) TmpIndex);
	  ++TmpIndex2;
	  ++TmpIndex;
	  this->TemporaryUV.Re(TmpIndex2) = uv.Re((int) TmpIndex);
	  this->TemporaryUV.Im(TmpIndex2) = uv.Im((int) TmpIndex);
	  ++TmpIndex2;
	}
      TotalValue += this->LocalCalculateFromSpinorVariables(this->TemporaryUV);
    }
  if (this->FermionFlag == true)
    {
      Complex TmpU;
      Complex TmpV;
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->NbrParticles; ++j)
	    WaveFunction *=  ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      TotalValue *= WaveFunction;
    }
  return TotalValue;
}  


// evaluate all permutations requested to symmetrize the SU(2) state
//

void FQHESphereSymmetrizedSU2ToU1WaveFunction::EvaluatePermutations()
{
  unsigned long Fact = 2;
  for (unsigned long i = 3; i <= ((unsigned long) this->NbrParticles); ++i)
    Fact *= i;
  if (this->FullySymmetrize == false)
    {
      unsigned long TmpNbrPermutation = Fact;  
      for (unsigned long i = 2; i <= ((unsigned long) this->NbrParticlesUp); ++i)
	TmpNbrPermutation /= i;
      for (unsigned long i = 2; i <= ((unsigned long) this->NbrParticlesDown); ++i)
	TmpNbrPermutation /= i;
      
      this->Permutations = new unsigned long[TmpNbrPermutation];
      
      unsigned long TmpPerm =  0x0ul;
      unsigned long* TmpArrayPerm = new unsigned long [this->NbrParticles];
      for (int k = 0; k < this->NbrParticles; ++k) 
	{
	  TmpPerm |= ((unsigned long) k) << (k << 2);
	  TmpArrayPerm[k] = (unsigned long) k;
	}
      
      this->Permutations[0] = TmpPerm;
      TmpNbrPermutation = 1ul;
      for (unsigned long j = 1; j < Fact; ++j)
	{
	  int Pos1 = this->NbrParticles - 1;
	  while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	    --Pos1;
	  --Pos1;
	  int Pos2 = this->NbrParticles - 1;      
	  while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	    --Pos2;
	  unsigned long TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  Pos2 = this->NbrParticles - 1;   
	  Pos1++;
	  while (Pos1 < Pos2)
	    {
	      TmpIndex = TmpArrayPerm[Pos1];
	      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	      TmpArrayPerm[Pos2] = TmpIndex;
	      ++Pos1;
	      --Pos2;
	    }
	  bool Flag = true;
	  int TmpIndex2 = 0;
	  int k = 1;
	  while ((k < this->NbrParticlesUp) && (Flag == true))
	    {
	      if (TmpArrayPerm[TmpIndex2] > TmpArrayPerm[TmpIndex2 + 1])
		Flag = false;
	      ++TmpIndex2;
	      ++k;
	    }
	  ++TmpIndex2;
	  k = 1;
	  while ((k < this->NbrParticlesDown) && (Flag == true))
	    {
	      if (TmpArrayPerm[TmpIndex2] > TmpArrayPerm[TmpIndex2 + 1])
		Flag = false;
	      ++TmpIndex2;
	      ++k;
	    } 
	  if (Flag == true)
	    {
	      TmpPerm =  0x0ul;
	      for (int k = 0; k < this->NbrParticles; ++k) 
		TmpPerm |= TmpArrayPerm[k] << (k << 2);
	      this->Permutations[TmpNbrPermutation] = TmpPerm;
	      ++TmpNbrPermutation;
	    }
	}
      this->NbrPermutations = TmpNbrPermutation;
    }
  else
    {
      this->Permutations = new unsigned long[Fact];
      this->NbrPermutations = Fact;
      unsigned long TmpPerm =  0x0ul;
      unsigned long* TmpArrayPerm = new unsigned long [this->NbrParticles];
      for (int k = 0; k < this->NbrParticles; ++k) 
	{
	  TmpPerm |= ((unsigned long) k) << (k << 2);
	  TmpArrayPerm[k] = (unsigned long) k;
	}
      this->Permutations[0] = TmpPerm;
      for (unsigned long j = 1; j < Fact; ++j)
	{
	  int Pos1 = this->NbrParticles - 1;
	  while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	    --Pos1;
	  --Pos1;
	  int Pos2 = this->NbrParticles - 1;      
	  while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	    --Pos2;
	  unsigned long TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  Pos2 = this->NbrParticles - 1;   
	  Pos1++;
	  while (Pos1 < Pos2)
	    {
	      TmpIndex = TmpArrayPerm[Pos1];
	      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	      TmpArrayPerm[Pos2] = TmpIndex;
	      ++Pos1;
	      --Pos2;
	    }
	  TmpPerm =  0x0ul;
	  for (int k = 0; k < this->NbrParticles; ++k) 
	    TmpPerm |= TmpArrayPerm[k] << (k << 2);
	  this->Permutations[j] = TmpPerm;
	}
    }
  return;
}


// get all permutations requested to symmetrize the SU(2) state from data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool FQHESphereSymmetrizedSU2ToU1WaveFunction::ReadPermutations(char* filename)
{
  ifstream File;
  File.open(filename, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << filename << endl;
      return false;
    }
  ReadLittleEndian(File, this->NbrParticles);
  ReadLittleEndian(File, this->NbrParticlesUp);
  ReadLittleEndian(File, this->NbrParticlesDown);
  ReadLittleEndian(File, this->NbrPermutations);
  this->Permutations = new unsigned long[this->NbrPermutations];
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    ReadLittleEndian(File, this->Permutations[i]);
  File.close();
  return true;
}

// write all permutations requested to symmetrize the SU(2) state to data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool FQHESphereSymmetrizedSU2ToU1WaveFunction::WritePermutations(char* filename)
{
  ofstream File;
  File.open(filename, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << filename << endl;
      return false;
    }
  WriteLittleEndian(File, this->NbrParticles);
  WriteLittleEndian(File, this->NbrParticlesUp);
  WriteLittleEndian(File, this->NbrParticlesDown);
  WriteLittleEndian(File, this->NbrPermutations);
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    WriteLittleEndian(File, this->Permutations[i]);
  File.close();
  return true;
}

// evaluate function at a given point(the first 2*N1 coordinates correspond to the position of the type 1 particles, 
//                                     the following 2*N2 coordinates correspond to the position of the type 2 particles,
//                                     last the 2*N3 coordinates correspond to the position of the type 3 particles)
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
//      this method is for only for internal class usage
// return value = function value at (uv)

Complex FQHESphereSymmetrizedSU2ToU1WaveFunction::LocalCalculateFromSpinorVariables(ComplexVector& uv)
{
  return this->SU2WaveFunction->CalculateFromSpinorVariables(uv);
}
