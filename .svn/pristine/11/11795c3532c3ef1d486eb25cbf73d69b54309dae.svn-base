////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of U(1) wave function obtained form symmetrization of a   //
//                        SU(K) wave function on sphere                       //
//                                                                            //
//                        last modification : 17/03/2008                      //
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
#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSUKToU1WaveFunction.h"
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

FQHESphereSymmetrizedSUKToU1WaveFunction::FQHESphereSymmetrizedSUKToU1WaveFunction()
{
}

// constructor
//
// nbrParticles = number of particles
// kValue = number of particle per cluster
// sUKWavefunction = pointer to the base SU(K) wave function
// fermionFlag = true if the final state should be a fermionic state
 
FQHESphereSymmetrizedSUKToU1WaveFunction::FQHESphereSymmetrizedSUKToU1WaveFunction(int nbrParticles, int kValue, Abstract1DComplexFunctionOnSphere* sUKWavefunction, bool fermionFlag)
{
  this->NbrParticles = nbrParticles;
  this->KValue = kValue;
  this->NbrParticlesPerColor = this->NbrParticles / this->KValue;
  this->SUKWaveFunction = (Abstract1DComplexFunctionOnSphere*) (sUKWavefunction->Clone());
  this->FermionFlag = fermionFlag;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  this->FullySymmetrize = false;  
  this->EvaluatePermutations();
  this->Flag.Initialize();
}

// constructor from data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// sUKWavefunction = pointer to the base SU(K) wave function
// fermionFlag = true if the final state should be a fermionic state
 
FQHESphereSymmetrizedSUKToU1WaveFunction::FQHESphereSymmetrizedSUKToU1WaveFunction(char* filename, Abstract1DComplexFunctionOnSphere* sUKWavefunction, bool fermionFlag)
{
  this->NbrParticles = 0;
  this->KValue = 0;
  this->ReadPermutations(filename);
  this->NbrParticlesPerColor = this->NbrParticles / this->KValue;
  this->SUKWaveFunction = (Abstract1DComplexFunctionOnSphere*) (sUKWavefunction->Clone());
  this->FermionFlag = fermionFlag;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  this->FullySymmetrize = false;  
  this->Flag.Initialize();
}

// copy constructor
//
// function = reference on the wave function to copy

FQHESphereSymmetrizedSUKToU1WaveFunction::FQHESphereSymmetrizedSUKToU1WaveFunction(const FQHESphereSymmetrizedSUKToU1WaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrParticlesPerColor = function.NbrParticlesPerColor;
  this->KValue = function.KValue;
  this->Permutations = function.Permutations;
  this->NbrPermutations = function.NbrPermutations;
  this->SUKWaveFunction = (Abstract1DComplexFunctionOnSphere*) (function.SUKWaveFunction->Clone());
  this->FermionFlag = function.FermionFlag;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  this->FullySymmetrize = false;  
  this->Flag = function.Flag;
}

// destructor
//

FQHESphereSymmetrizedSUKToU1WaveFunction::~FQHESphereSymmetrizedSUKToU1WaveFunction()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete[] this->Permutations;
    }
  if (this->SUKWaveFunction != 0)
    delete this->SUKWaveFunction;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHESphereSymmetrizedSUKToU1WaveFunction::Clone ()
{
  return new FQHESphereSymmetrizedSUKToU1WaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex FQHESphereSymmetrizedSUKToU1WaveFunction::operator ()(RealVector& x)
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

Complex FQHESphereSymmetrizedSUKToU1WaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
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


// evaluate all permutations requested to symmetrize the SU(K) state
//

void FQHESphereSymmetrizedSUKToU1WaveFunction::EvaluatePermutations()
{
  unsigned long Fact = 2;
  for (unsigned long i = 3; i <= ((unsigned long) this->NbrParticles); ++i)
    Fact *= i;
  if (this->FullySymmetrize == false)
    {
      unsigned long TmpNbrPermutation = Fact;  
      unsigned long FactNbrParticlesPerColor = 1;
      for (unsigned long i = 2; i <= ((unsigned long) this->NbrParticlesPerColor); ++i)
	FactNbrParticlesPerColor *= i;
      for (int i = 0; i < this->KValue; ++i)
	TmpNbrPermutation /= FactNbrParticlesPerColor;
      unsigned long FactKValue = 1;
      for (unsigned long i = 2; i <= ((unsigned long) this->KValue); ++i)
	FactKValue *= i;
      TmpNbrPermutation /= FactKValue;
      
      this->Permutations = new unsigned long[TmpNbrPermutation];
      
      unsigned long TmpPerm =  0x0ul;
      unsigned long* TmpArrayPerm = new unsigned long [this->NbrParticles];
      unsigned long* TmpArrayPerm2 = new unsigned long [this->KValue];
      for (int k = 0; k < this->NbrParticles; ++k) 
	{
	  TmpPerm |= ((unsigned long) k) << (k << 2);
	  TmpArrayPerm[k] = (unsigned long) k;
	}
      
      this->Permutations[0] = TmpPerm;
      TmpNbrPermutation = 1ul;
      Fact /= ((unsigned long) this->NbrParticles);
      Fact *= ((unsigned long) (this->NbrParticles - this->NbrParticlesPerColor + 1));
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
	  for (int i = 0; ((i < this->KValue) && (Flag == true)); ++i)
	    {
	      int k = 1;
	      while ((k < this->NbrParticlesPerColor) && (Flag == true))
		{
		  if (TmpArrayPerm[TmpIndex2] > TmpArrayPerm[TmpIndex2 + 1])
		    Flag = false;
		  ++TmpIndex2;
		  ++k;
		}
	      ++TmpIndex2;
	    } 
	  if (Flag == true)
	    {
	      for (int k = 0; k < this->KValue; ++k)
		{
		  TmpIndex = 0ul;
		  for (int i = 0; i < this->NbrParticlesPerColor; ++i)
		    TmpIndex |= TmpArrayPerm[i + (k * this->NbrParticlesPerColor)] << (i << 2);
		  TmpArrayPerm2[k] = TmpIndex;
		}
	      for (int i = 1; ((i < this->KValue) && (Flag == true)); ++i)
		if (TmpArrayPerm2[i - 1] > TmpArrayPerm2[i])
		  Flag = false;
	      if (Flag == true)
		{
		  TmpPerm =  0x0ul;
		  for (int i = 0; i < this->KValue; ++i)
		    TmpPerm |= (TmpArrayPerm2[i] << ((i * this->NbrParticlesPerColor) << 2));
		  this->Permutations[TmpNbrPermutation] = TmpPerm;	      
		  ++TmpNbrPermutation;
		}
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


// get all permutations requested to symmetrize the SU(K) state from data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool FQHESphereSymmetrizedSUKToU1WaveFunction::ReadPermutations(char* filename)
{
  ifstream File;
  File.open(filename, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << filename << endl;
      return false;
    }
  ReadLittleEndian(File, this->NbrParticles);
  ReadLittleEndian(File, this->KValue);
  ReadLittleEndian(File, this->NbrPermutations);
  this->Permutations = new unsigned long[this->NbrPermutations];
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    ReadLittleEndian(File, this->Permutations[i]);
  File.close();
  return true;
}

// write all permutations requested to symmetrize the SU(K) state to data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool FQHESphereSymmetrizedSUKToU1WaveFunction::WritePermutations(char* filename)
{
  ofstream File;
  File.open(filename, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << filename << endl;
      return false;
    }
  WriteLittleEndian(File, this->NbrParticles);
  WriteLittleEndian(File, this->KValue);
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

Complex FQHESphereSymmetrizedSUKToU1WaveFunction::LocalCalculateFromSpinorVariables(ComplexVector& uv)
{
  return this->SUKWaveFunction->CalculateFromSpinorVariables(uv);
}
