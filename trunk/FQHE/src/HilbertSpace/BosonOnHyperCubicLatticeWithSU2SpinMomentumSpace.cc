////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of bosons on a 4D hypercubic lattice with SU(2) spin      //
//                                in momentum space                           //
//                                                                            //
//                        last modification : 04/04/2012                      //
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
#include "HilbertSpace/BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace.h"
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
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

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


// basic constructor
// 
// nbrBosons = number of bosons
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// nbrSiteT = number of sites in the t direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// ktMomentum = momentum along the t direction
// memory = amount of memory granted for precalculations

BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace (int nbrBosons, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT,
												    int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrBosonsUp = 0;
  this->NbrBosonsDown = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteT = nbrSiteT;
  this->NbrSiteYZT = this->NbrSiteZ * this->NbrSiteY * this->NbrSiteT;
  this->NbrSiteZT = this->NbrSiteZ * this->NbrSiteT;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->KzMomentum = kzMomentum;
  this->KtMomentum = ktMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ * this->NbrSiteT;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];

  this->NUpLzMax = this->LzMax + this->NbrBosons - 1;
  this->NDownLzMax = this->LzMax + this->NbrBosons - 1;
  this->FermionicLzMax = this->NUpLzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, this->NbrSiteT - 1, 0, 0, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, this->NbrSiteT - 1, 0, 0, 0, 0, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      for (long i = 0; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->StateDescriptionUp[i];
	  unsigned long Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionUp[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescriptionDown[i];
	  Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionDown[i] >>= this->NbrBosons - Tmp; 
	}
      SortDoubleElementArrayDownOrdering<unsigned long>(this->StateDescriptionUp, this->StateDescriptionDown, TmpLargeHilbertSpaceDimension);
      this->GenerateLookUpTable(memory);
//    for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
//       {
//         cout << i << " : ";
//         this->PrintState(cout, i);
//         cout << this->FindStateIndex(this->StateDescriptionUp[i], this->StateDescriptionDown[i]);
//         cout << endl;
//       }
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (4 * sizeof(unsigned long));
      cout << "memory requested for Hilbert space = ";
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

// basic constructor when Sz is conserved
// 
// nbrBosons = number of bosons
// nbrSpinUp = number of particles with spin up
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// nbrSiteT = number of sites in the t direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// ktMomentum = momentum along the t direction
// memory = amount of memory granted for precalculations

BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace (int nbrBosons, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT,
												    int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrBosonsUp = nbrSpinUp;
  this->NbrBosonsDown = this->NbrBosons - this->NbrBosonsUp;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteT = nbrSiteT;
  this->NbrSiteYZT = this->NbrSiteZ * this->NbrSiteY * this->NbrSiteT;
  this->NbrSiteZT = this->NbrSiteZ * this->NbrSiteT;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->KzMomentum = kzMomentum;
  this->KtMomentum = ktMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ * this->NbrSiteT;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];

  this->NUpLzMax = this->LzMax + this->NbrBosonsUp - 1;
  this->NDownLzMax = this->LzMax + this->NbrBosonsDown - 1;
  this->FermionicLzMax = this->NUpLzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, this->NbrSiteT - 1, 0, 0, 0, 0, this->NbrBosonsUp);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, this->NbrSiteT - 1, 0, 0, 0, 0, this->LzMax + this->NbrBosonsUp, this->LzMax + this->NbrBosonsDown, this->NbrBosonsUp, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      SortDoubleElementArrayDownOrdering<unsigned long>(this->StateDescriptionUp, this->StateDescriptionDown, TmpLargeHilbertSpaceDimension);
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (4 * sizeof(unsigned long));
      cout << "memory requested for Hilbert space = ";
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
// bosons = reference on the hilbert space to copy to copy

BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace(const BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace& bosons)
{
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->NbrSiteZ = bosons.NbrSiteZ;
  this->NbrSiteT = bosons.NbrSiteT;
  this->NbrSiteYZT = bosons.NbrSiteYZT;
  this->NbrSiteZT = bosons.NbrSiteZT;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->KzMomentum = bosons.KzMomentum;
  this->KtMomentum = bosons.KtMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->SzFlag = bosons.SzFlag;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
}

// destructor
//

BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::~BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace& BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::operator = (const BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      delete[] this->UniqueStateDescriptionUp;
      delete[] this->UniqueStateDescriptionSubArraySizeUp;
      delete[] this->FirstIndexUniqueStateDescriptionUp;
    }
  delete[] this->TemporaryStateUp;
  delete[] this->TemporaryStateDown;
  delete[] this->ProdATemporaryStateUp;
  delete[] this->ProdATemporaryStateDown;
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->NbrSiteZ = bosons.NbrSiteZ; 
  this->NbrSiteT = bosons.NbrSiteT;
  this->NbrSiteYZT = bosons.NbrSiteYZT;
  this->NbrSiteZT = bosons.NbrSiteZT;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->KzMomentum = bosons.KzMomentum;
  this->KtMomentum = bosons.KtMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->SzFlag = bosons.SzFlag;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::Clone()
{
  return new BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state], TemporaryStateUp, TemporaryStateDown); 

  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int TmpKx = i / this->NbrSiteYZT;
      int TmpKy = i % this->NbrSiteYZT;
      int TmpKz = TmpKy % this->NbrSiteZT;
      TmpKy /= this->NbrSiteZT;
      int TmpKt = TmpKz % this->NbrSiteT;
      TmpKz /= this->NbrSiteT;
      if (this->TemporaryStateUp[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateUp[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << "," << TmpKt << ",up)";
	}
      if (this->TemporaryStateDown[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateDown[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << "," << TmpKt << ",down)";
	}
    }
  Str << "]";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentKz, int currentKt,
								      int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, int currentFermionicPositionUp, int currentFermionicPositionDown, long pos)
{
  if (currentKt < 0)
    {
      currentKt = this->NbrSiteT - 1;
      currentKz--;
      if (currentKz < 0)
	{
	  currentKz = this->NbrSiteZ - 1;
	  currentKy--;
	  if (currentKy < 0)
	    {
	      currentKy = this->NbrSiteY - 1;
	      currentKx--;
	    }
	}
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
	  && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum) && ((currentTotalKt % this->NbrSiteT) == this->KtMomentum))
	{
 	  this->StateDescriptionUp[pos] = 0x0ul;
 	  this->StateDescriptionDown[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  for (int i = nbrBosons; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul) << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + ((i + j) * currentKx), currentTotalKy + ((i + j) * currentKy), currentTotalKz + ((i + j) * currentKz), currentTotalKt + ((i + j) * currentKt), currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, pos);
	  unsigned long MaskDown = ((0x1ul << j) - 0x1ul) << (currentFermionicPositionDown - j - 1);
	  for (; pos < TmpPos; ++pos)
	    {
 		  this->StateDescriptionUp[pos] |= MaskUp;
 		  this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
    }
  return pos;
};

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// nbrSpinUp = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentKz, int currentKt,
								      int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, 
								      int currentFermionicPositionUp, int currentFermionicPositionDown, int nbrSpinUp, long pos)
{
  if (currentKt < 0)
    {
      currentKt = this->NbrSiteT - 1;
      currentKz--;
      if (currentKz < 0)
	{
	  currentKz = this->NbrSiteZ - 1;
	  currentKy--;
	  if (currentKy < 0)
	    {
	      currentKy = this->NbrSiteY - 1;
	      currentKx--;
	    }
	}
    }
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
    return 0l;
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
	  && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum) && ((currentTotalKt % this->NbrSiteT) == this->KtMomentum))
	{
 	  this->StateDescriptionUp[pos] = 0x0ul;
 	  this->StateDescriptionDown[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  for (int i = nbrSpinUp; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul) << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + ((i + j) * currentKx), currentTotalKy + ((i + j) * currentKy), currentTotalKz + ((i + j) * currentKz), currentTotalKt + ((i + j) * currentKt), currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, nbrSpinUp - i, pos);
	  unsigned long MaskDown = ((0x1ul << j) - 0x1ul) << (currentFermionicPositionDown - j - 1);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescriptionUp[pos] |= MaskUp;
	      this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
    }
  return pos;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// return value = Hilbert space dimension

long BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentKz, int currentKt,
										     int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt)
{
  if (currentKt < 0)
    {
      currentKt = this->NbrSiteT - 1;
      currentKz--;
      if (currentKz < 0)
	{
	  currentKz = this->NbrSiteZ - 1;
	  currentKy--;
	  if (currentKy < 0)
	    {
	      currentKy = this->NbrSiteY - 1;
	      currentKx--;
	    }
	}
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
	  && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum)&& ((currentTotalKt % this->NbrSiteT) == this->KtMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  for (int i = nbrBosons; i >= 0; --i)
    Count += (((long) i) + 1l) * this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy), currentTotalKz + (i * currentKz), currentTotalKt + (i * currentKt));
  return Count;
}

// evaluate Hilbert space dimension with a fixed number of bosons with spin up
//
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// nbrSpinUp = number of particles with spin up
// return value = Hilbert space dimension

long BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentKz, int currentKt,
										int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, int nbrSpinUp)
{
  if (currentKt < 0)
    {
      currentKt = this->NbrSiteT - 1;
      currentKz--;
      if (currentKz < 0)
	{
	  currentKz = this->NbrSiteZ - 1;
	  currentKy--;
	  if (currentKy < 0)
	    {
	      currentKy = this->NbrSiteY - 1;
	      currentKx--;
	    }
	}
    }
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
    return 0l;

  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
	  && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum)&& ((currentTotalKt % this->NbrSiteT) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  for (int i = nbrBosons; i >= 0; --i)
    for (int j = i; j >= 0; --j)
      Count += this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy), currentTotalKz + (i * currentKz), currentTotalKt + (i * currentKt), nbrSpinUp -j);
  return Count;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = ky sector in which the density matrix has to be evaluated 
// kzSector = kz sector in which the density matrix has to be evaluated 
// kzSector = kt sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, int kzSector, int ktSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0) && (kzSector == 0) && (ktSector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum) && (kzSector == this->KzMomentum)
	  && (ktSector == this->KxMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementaryKxMomentum = (this->KxMomentum - kxSector) % this->NbrSiteX;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector) % this->NbrSiteY;
  int ComplementaryKzMomentum = (this->KzMomentum - kzSector) % this->NbrSiteZ;
  int ComplementaryKtMomentum = (this->KtMomentum - ktSector) % this->NbrSiteT;
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->NbrSiteX;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrSiteY;
  if (ComplementaryKzMomentum < 0)
    ComplementaryKzMomentum += this->NbrSiteZ;
  if (ComplementaryKtMomentum < 0)
    ComplementaryKtMomentum += this->NbrSiteT;
  cout << "kx = " << this->KxMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
  cout << "ky = " << this->KyMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
  cout << "kz = " << this->KzMomentum << " " << kzSector << " " << ComplementaryKzMomentum << endl;
  cout << "kt = " << this->KtMomentum << " " << ktSector << " " << ComplementaryKtMomentum << endl;
  BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace SubsytemSpace (nbrParticleSector, this->NbrSiteX, this->NbrSiteY, this->NbrSiteZ, this->NbrSiteT, kxSector, kySector, kzSector, ktSector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, this->NbrSiteZ, this->NbrSiteT, ComplementaryKxMomentum, ComplementaryKyMomentum, ComplementaryKzMomentum, ComplementaryKtMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = ky sector in which the density matrix has to be evaluated 
// kzSector = kz sector in which the density matrix has to be evaluated 
// ktSector = kt sector in which the density matrix has to be evaluated 
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, int kzSector, int ktSector, 
													    int nbrGroundStates, ComplexVector* groundStates, double* weights, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0) && (kzSector == 0) && (ktSector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 0.0;
	  for (int i = 0; i < nbrGroundStates; ++i)
	    TmpDensityMatrix(0, 0) += weights[i];
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum) && (kzSector == this->KzMomentum)
	   && (ktSector == this->KtMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 0.0;
	  for (int i = 0; i < nbrGroundStates; ++i)
	    TmpDensityMatrix(0, 0) += weights[i];
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementaryKxMomentum = (this->KxMomentum - kxSector) % this->NbrSiteX;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector) % this->NbrSiteY;
  int ComplementaryKzMomentum = (this->KzMomentum - kzSector) % this->NbrSiteZ;
  int ComplementaryKtMomentum = (this->KtMomentum - ktSector) % this->NbrSiteT;
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->NbrSiteX;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrSiteY;
  if (ComplementaryKzMomentum < 0)
    ComplementaryKzMomentum += this->NbrSiteZ;
  if (ComplementaryKtMomentum < 0)
    ComplementaryKtMomentum += this->NbrSiteT;
  cout << "kx = " << this->KxMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
  cout << "ky = " << this->KyMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
  cout << "kz = " << this->KzMomentum << " " << kzSector << " " << ComplementaryKzMomentum << endl;
  cout << "kt = " << this->KtMomentum << " " << ktSector << " " << ComplementaryKtMomentum << endl;
  BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace SubsytemSpace (nbrParticleSector, this->NbrSiteX, this->NbrSiteY, this->NbrSiteZ, this->NbrSiteT, 
								  kxSector, kySector, kzSector, ktSector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, this->NbrSiteZ, this->NbrSiteT,
								 ComplementaryKxMomentum, ComplementaryKyMomentum, ComplementaryKzMomentum, ComplementaryKtMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, nbrGroundStates, groundStates, weights, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

