////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on square lattice with spin            //
//                                  in momentum space                         //
//                                                                            //
//                        last modification : 16/02/2011                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceTruncated.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


#ifdef HAVE_FTI
// constructors relying on FTI technology...
// options as above, except:
//
// tightBindingModel = pointer to the relevant tight-binding model
// cutoff = energy cut-off in configuration space, with respect to the groundstate energy
//
FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::FermionOnSquareLatticeWithSpinMomentumSpaceTruncated (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, Abstract2DTightBindingModel *tightBindingModel, double cutoff, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->TightBindingEnergies = new double[2*this->LzMax];
  double MinimumEnergy = 0.0;
  if (tightBindingModel!=NULL)
    {
      MinimumEnergy = tightBindingModel->SingleParticleGroundstateEnergy();
      for (int kx=0; kx<NbrSiteX; ++kx)
	for (int ky=0; ky<NbrSiteY; ++ky)
	  {
	    int InternalIndex = (kx * this->NbrSiteY) + ky;
	    int TightBindingIndex = tightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	    this->TightBindingEnergies[InternalIndex<<1] = tightBindingModel->GetEnergy(0, TightBindingIndex)-MinimumEnergy;
	    this->TightBindingEnergies[(InternalIndex<<1)+1] = tightBindingModel->GetEnergy(1, TightBindingIndex)-MinimumEnergy;
	  }
    }
  else
    {
      for (int i=0; i<LzMax; ++i)
	{
	  this->TightBindingEnergies[i<<1] = 0.0;
	  this->TightBindingEnergies[(i<<1)+1] = 1.0;
	}
    }
  int FilledNbrBands;
  this->EnergyCutoff = tightBindingModel->ComputeGroundstateEnergy(NbrFermions, FilledNbrBands)-NbrFermions*MinimumEnergy+cutoff;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0.0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0.0, 0l);
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      PrintMemorySize(cout,UsedMemory)<<endl;
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      PrintMemorySize(cout,UsedMemory)<<endl;
#endif
    }

}


FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::FermionOnSquareLatticeWithSpinMomentumSpaceTruncated (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, Abstract2DTightBindingModel *tightBindingModel, double cutoff, unsigned long memory)
{
    this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = nbrSpinUp;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->TightBindingEnergies = new double[2*this->LzMax];
  double MinimumEnergy = 0.0;
  if (tightBindingModel!=NULL)
    {
      MinimumEnergy = tightBindingModel->SingleParticleGroundstateEnergy();
      for (int kx=0; kx<NbrSiteX; ++kx)
	for (int ky=0; ky<NbrSiteY; ++ky)
	  {
	    int InternalIndex = (kx * this->NbrSiteY) + ky;
	    int TightBindingIndex = tightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	    this->TightBindingEnergies[InternalIndex<<1] = tightBindingModel->GetEnergy(0, TightBindingIndex)-MinimumEnergy;
	    this->TightBindingEnergies[(InternalIndex<<1)+1] = tightBindingModel->GetEnergy(1, TightBindingIndex)-MinimumEnergy;
	  }
    }
  else
    {
      for (int i=0; i<LzMax; ++i)
	{
	  this->TightBindingEnergies[i<<1] = 0.0;
	  this->TightBindingEnergies[(i<<1)+1] = 1.0;
	}
    }
  int FilledNbrBands;
  this->EnergyCutoff = tightBindingModel->ComputeGroundstateEnergy(NbrFermions, FilledNbrBands) - NbrFermions*MinimumEnergy + cutoff;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrFermionsUp, 0.0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrFermionsUp, 0.0, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      PrintMemorySize(cout,UsedMemory)<<endl;

      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      PrintMemorySize(cout,UsedMemory)<<endl;
#endif
    }
  
}
  
#endif



// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// tightBindingEnergies = pointer to the energies of a tight binding model of the corresponding band structure
// cutoff = energy cut-off in configuration space, in absolute units
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::FermionOnSquareLatticeWithSpinMomentumSpaceTruncated (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, double *tightBindingEnergies, double cutoff, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->TightBindingEnergies = new double[2*this->LzMax];
  double MinimumEnergy = tightBindingEnergies[0];
  for (int i=1; i<this->LzMax<<1; ++i)
    if (tightBindingEnergies[i]<MinimumEnergy)
      MinimumEnergy = tightBindingEnergies[i];
  for (int i=0; i<this->LzMax<<1; ++i)
    this->TightBindingEnergies[i] = tightBindingEnergies[i] - MinimumEnergy;
  this->EnergyCutoff = cutoff - NbrFermions * MinimumEnergy;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0.0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0.0, 0l);
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      PrintMemorySize(cout,UsedMemory)<<endl;
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      PrintMemorySize(cout,UsedMemory)<<endl;
#endif
    }
}

// basic constructor when Sz is preserved
// 
// nbrFermions = number of fermions
// nbrSpinUp = number of particles with spin up
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// tightBindingEnergies = pointer to the energies of a tight binding model of the corresponding band structure
// cutoff = energy cut-off in configuration space, in absolute units
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::FermionOnSquareLatticeWithSpinMomentumSpaceTruncated (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, double *tightBindingEnergies, double cutoff, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = nbrSpinUp;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->TightBindingEnergies = new double[2*this->LzMax];
  double MinimumEnergy = tightBindingEnergies[0];
  for (int i=1; i<this->LzMax<<1; ++i)
    if (tightBindingEnergies[i]<MinimumEnergy)
      MinimumEnergy = tightBindingEnergies[i];
  for (int i=0; i<this->LzMax<<1; ++i)
    this->TightBindingEnergies[i] = tightBindingEnergies[i] - MinimumEnergy;
  this->EnergyCutoff = cutoff - NbrFermions * MinimumEnergy;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrFermionsUp, 0.0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrFermionsUp, 0.0, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      PrintMemorySize(cout,UsedMemory)<<endl;

      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      PrintMemorySize(cout,UsedMemory)<<endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::FermionOnSquareLatticeWithSpinMomentumSpaceTruncated(const FermionOnSquareLatticeWithSpinMomentumSpaceTruncated& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->HighestBit = fermions.HighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->TightBindingEnergies = fermions.TightBindingEnergies;
  this->EnergyCutoff = fermions.EnergyCutoff;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::~FermionOnSquareLatticeWithSpinMomentumSpaceTruncated ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->TightBindingEnergies;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSpinMomentumSpaceTruncated& FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::operator = (const FermionOnSquareLatticeWithSpinMomentumSpaceTruncated& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
      delete[] this->TightBindingEnergies;

    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->SzFlag = fermions.SzFlag;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->TightBindingEnergies = fermions.TightBindingEnergies;
  this->EnergyCutoff = fermions.EnergyCutoff;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::Clone()
{
  return new FermionOnSquareLatticeWithSpinMomentumSpaceTruncated(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i << 1));
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if ((Tmp & 0x2l) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",+)";
      if ((Tmp & 0x1l) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",-)";
    }
  Str << "]";
//   Str << " " << TmpState; 
//   Str << " " << hex << TmpState << dec; 
//   Str << " " << state;
//   int TmpLzMax = (this->LzMax << 1) + 1;
//   while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
//     --TmpLzMax;
//   Str << " " << this->FindStateIndex(TmpState, TmpLzMax);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentEnergy = current energy of particles in configuration
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, double currentEnergy, long pos)
{
  if (currentEnergy > EnergyCutoff)
    return pos;
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
    
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(currentKx, j, 1) < EnergyCutoff))
	    {
	      this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 1);
	      ++pos;
	    }
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(currentKx, j, 0) < EnergyCutoff))
	    {
	      this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 1);
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(i, j, 1) < EnergyCutoff))
		{
		  this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 1);
		  ++pos;
		}
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(i, j, 0) < EnergyCutoff))
		{
 		  this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 1);
 		  ++pos;
		}
	    }
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentEnergy + GetEnergy(currentKx, currentKy, 0) + GetEnergy(currentKx, currentKy, 1), pos);
  unsigned long Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentEnergy + GetEnergy(currentKx, currentKy, 1), pos);
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
   TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentEnergy + GetEnergy(currentKx, currentKy, 0), pos);
   Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
   return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, currentEnergy, pos);
};

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of fermions with spin up
// currentEnergy = current energy of particles in configuration
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrSpinUp, double currentEnergy, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }

  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrFermions) || (currentEnergy > EnergyCutoff))
    return pos;

  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      if (nbrSpinUp == 1)
	{
	  for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(currentKx, j, nbrSpinUp) < EnergyCutoff))
		{
		  this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 1);
		  ++pos;
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(i, j, nbrSpinUp) < EnergyCutoff))
		    {
		      this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 1);
		      ++pos;
			}
		}
	    }
	}
      else
	{
	  for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(currentKx, j, nbrSpinUp) < EnergyCutoff))
		{
		  this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 1);
		  ++pos;
		}
	    }
	  for (int i = currentKx - 1; i >= 0; --i)
	    {
	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(i, j, nbrSpinUp) < EnergyCutoff))
		    {
		      this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 1);
		      ++pos;
		    }
		}
	    }
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrSpinUp - 1, currentEnergy + GetEnergy(currentKx, currentKy, 0) + GetEnergy(currentKx, currentKy, 1), pos);
  unsigned long Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrSpinUp - 1, currentEnergy + GetEnergy(currentKx, currentKy, 1), pos);
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrSpinUp, currentEnergy + GetEnergy(currentKx, currentKy, 0), pos);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrSpinUp, currentEnergy, pos);  
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentEnergy = current energy of particles in configuration
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, double currentEnergy)
{
  if (currentEnergy > EnergyCutoff)
    return 0l;
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(currentKx, j, 0) < EnergyCutoff))
	    Count += 1l;
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(currentKx, j, 1) < EnergyCutoff))
	    Count += 1l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(i, j, 0) < EnergyCutoff))
		Count += 1l;
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(i, j, 1) < EnergyCutoff))
		Count += 1l;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentEnergy + this->GetEnergy(currentKx, currentKy, 0) + this->GetEnergy(currentKx, currentKy, 1));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentEnergy + this->GetEnergy(currentKx, currentKy, 0));
    Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentEnergy + this->GetEnergy(currentKx, currentKy, 1));
    Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, currentEnergy);
  return Count;
}


// evaluate Hilbert space dimension with a fixed number of fermions with spin up
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of fermions with spin up
// currentEnergy = current energy of particles in configuration
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrSpinUp, double currentEnergy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrFermions) || (currentEnergy > EnergyCutoff))
    return 0l;

  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      int BandIndex = nbrSpinUp;

      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(currentKx, j, BandIndex) < EnergyCutoff))
	    Count++;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && (currentEnergy + GetEnergy(i, j, BandIndex) < EnergyCutoff))
		Count++;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrSpinUp - 1, currentEnergy + this->GetEnergy(currentKx, currentKy, 0) + this->GetEnergy(currentKx, currentKy, 1) );
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrSpinUp, currentEnergy + this->GetEnergy(currentKx, currentKy, 0));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrSpinUp - 1, currentEnergy + this->GetEnergy(currentKx, currentKy, 1));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrSpinUp, currentEnergy);
  return Count;
}


// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::FindStateIndex(unsigned long stateDescription, int lzmax)
{
//   cout << "Need to fix search algorithm!"<<endl;
//   long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
//   long PosMin = this->LookUpTable[lzmax][PosMax];
//   PosMax = this->LookUpTable[lzmax][PosMax + 1];

  // skip look-up table for now:
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    {
      return this->HilbertSpaceDimension;
    }
  
  long PosMax = 0;
  long PosMin = this->HilbertSpaceDimension - 1;  
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    if (stateDescription == this->StateDescription[PosMin])
      return PosMin;
    else
      return this->HilbertSpaceDimension;
}
