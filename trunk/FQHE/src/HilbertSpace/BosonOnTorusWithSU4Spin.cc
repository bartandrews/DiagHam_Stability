////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of bosons on a torus with SU(4) spin                //
//                                                                            //
//                        last modification : 26/06/2012                      //
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
#include "HilbertSpace/BosonOnTorusWithSU4Spin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "HilbertSpace/BosonOnTorusShort.h" 
#include "GeneralTools/ArrayTools.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;


// constructor with a constraint on total spin momentum and total momentum
// 
// nbrBosons = number of bosons
// totalSpin = twice the total spin value
// totalIsospin = twice the total isospin value
// totalEntanglement = twice the total entanglement value
// maxMomentum = momentum maximum value for a boson
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnTorusWithSU4Spin::BosonOnTorusWithSU4Spin (int nbrBosons, int totalSpin, int totalIsospin, int totalEntanglement, int maxMomentum, int kyMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = maxMomentum;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = totalIsospin;
  this->TotalEntanglement = totalEntanglement;

  this->KyMomentum = kyMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = maxMomentum;
  this->Flag.Initialize();
  this->TemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->NbrLzValue];

  int NUpPlus = this->NbrBosons + this->TotalSpin + this->TotalIsospin + this->TotalEntanglement;
  int NUpMinus = this->NbrBosons + this->TotalSpin - this->TotalIsospin - this->TotalEntanglement;
  int NDownPlus = this->NbrBosons - this->TotalSpin + this->TotalIsospin - this->TotalEntanglement;
  int NDownMinus = this->NbrBosons - this->TotalSpin - this->TotalIsospin + this->TotalEntanglement;
  if ((NUpPlus < 0) || (NUpMinus < 0) || (NDownPlus < 0) || (NDownMinus < 0) || ((NUpPlus & 3) != 0) || ((NUpMinus & 3) != 0) || ((NDownPlus & 3) != 0) || ((NDownMinus & 3) != 0))
    this->LargeHilbertSpaceDimension = 0l;
  else
    {
      NUpPlus >>= 2;
      NUpMinus >>= 2;
      NDownPlus  >>= 2;
      NDownMinus  >>= 2;
      this->NUpPlusLzMax = this->LzMax + NUpPlus - 1;
      this->NUpMinusLzMax = this->LzMax + NUpMinus - 1;
      this->NDownPlusLzMax = this->LzMax + NDownPlus - 1;
      this->NDownMinusLzMax = this->LzMax + NDownMinus - 1;
      this->FermionicLzMax = this->NUpPlusLzMax;
      if (this->NUpMinusLzMax > this->FermionicLzMax)
	this->FermionicLzMax = this->NUpMinusLzMax;
      if (this->NDownPlusLzMax > this->FermionicLzMax)
	this->FermionicLzMax = this->NDownPlusLzMax;
      if (this->NDownMinusLzMax > this->FermionicLzMax)
	this->FermionicLzMax = this->NDownMinusLzMax;
      this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, 0, NUpPlus, NUpMinus, NDownPlus, NDownMinus);
    }
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUpPlus = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionUpMinus = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDownPlus = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDownMinus = new unsigned long [this->LargeHilbertSpaceDimension];
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, this->LzMax, this->LzMax, 0, 
								NUpPlus, NUpMinus, NDownPlus, NDownMinus, 0l);
      cout  << "Dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      SortQuadElementArrayDownOrdering<unsigned long>(this->StateDescriptionUpPlus, this->StateDescriptionUpMinus, this->StateDescriptionDownPlus, this->StateDescriptionDownMinus, TmpLargeHilbertSpaceDimension);
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

BosonOnTorusWithSU4Spin::BosonOnTorusWithSU4Spin(const BosonOnTorusWithSU4Spin& bosons)
{
  this->KyMomentum = bosons.KyMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->TotalIsospin = bosons.TotalIsospin;
  this->TotalEntanglement = bosons.TotalEntanglement;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpPlusLzMax = bosons.NUpPlusLzMax;
  this->NUpMinusLzMax = bosons.NUpMinusLzMax;
  this->NDownPlusLzMax = bosons.NDownPlusLzMax;
  this->NDownMinusLzMax = bosons.NDownMinusLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->StateDescriptionUpPlus = bosons.StateDescriptionUpPlus;
  this->StateDescriptionUpMinus = bosons.StateDescriptionUpMinus;
  this->StateDescriptionDownPlus = bosons.StateDescriptionDownPlus;
  this->StateDescriptionDownMinus = bosons.StateDescriptionDownMinus;
  this->NbrUniqueStateDescriptionUpPlus = bosons.NbrUniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionUpPlus = bosons.UniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionSubArraySizeUpPlus = bosons.UniqueStateDescriptionSubArraySizeUpPlus;
  this->NbrUniqueStateDescriptionUpMinus = bosons.NbrUniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionUpMinus = bosons.UniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionSubArraySizeUpMinus = bosons.UniqueStateDescriptionSubArraySizeUpMinus;
  this->FirstIndexUniqueStateDescriptionUpMinus = bosons.FirstIndexUniqueStateDescriptionUpMinus;
  this->NbrUniqueStateDescriptionDownPlus = bosons.NbrUniqueStateDescriptionDownPlus;
  this->UniqueStateDescriptionDownPlus = bosons.UniqueStateDescriptionDownPlus;
  this->UniqueStateDescriptionSubArraySizeDownPlus = bosons.UniqueStateDescriptionSubArraySizeDownPlus;
  this->FirstIndexUniqueStateDescriptionDownPlus = bosons.FirstIndexUniqueStateDescriptionDownPlus;
}

// destructor
//

BosonOnTorusWithSU4Spin::~BosonOnTorusWithSU4Spin ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusWithSU4Spin& BosonOnTorusWithSU4Spin::operator = (const BosonOnTorusWithSU4Spin& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionUpPlus;
      delete[] this->StateDescriptionUpMinus;
      delete[] this->StateDescriptionDownPlus;
      delete[] this->StateDescriptionDownMinus;
      delete[] this->UniqueStateDescriptionUpPlus;
      delete[] this->UniqueStateDescriptionSubArraySizeUpPlus;
      delete[] this->NbrUniqueStateDescriptionUpMinus;
      for (long i = 0l; i < this->NbrUniqueStateDescriptionUpPlus; ++i)
	{
	  delete[] this->UniqueStateDescriptionUpMinus[i];
	  delete[] this->UniqueStateDescriptionSubArraySizeUpMinus[i];
	  delete[] this->FirstIndexUniqueStateDescriptionUpMinus[i];
	}
      delete[] this->UniqueStateDescriptionUpMinus;
      delete[] this->UniqueStateDescriptionSubArraySizeUpMinus;
      delete[] this->FirstIndexUniqueStateDescriptionUpMinus;
    }
  delete[] this->TemporaryStateUpPlus;
  delete[] this->TemporaryStateUpMinus;
  delete[] this->TemporaryStateDownPlus;
  delete[] this->TemporaryStateDownMinus;
  delete[] this->ProdATemporaryStateUpPlus;
  delete[] this->ProdATemporaryStateUpMinus;
  delete[] this->ProdATemporaryStateDownPlus;
  delete[] this->ProdATemporaryStateDownMinus;
  this->KyMomentum = bosons.KyMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->TotalIsospin = bosons.TotalIsospin;
  this->TotalEntanglement = bosons.TotalEntanglement;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpPlusLzMax = bosons.NUpPlusLzMax;
  this->NUpMinusLzMax = bosons.NUpMinusLzMax;
  this->NDownPlusLzMax = bosons.NDownPlusLzMax;
  this->NDownMinusLzMax = bosons.NDownMinusLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->StateDescriptionUpPlus = bosons.StateDescriptionUpPlus;
  this->StateDescriptionUpMinus = bosons.StateDescriptionUpMinus;
  this->StateDescriptionDownPlus = bosons.StateDescriptionDownPlus;
  this->StateDescriptionDownMinus = bosons.StateDescriptionDownMinus;
  this->NbrUniqueStateDescriptionUpPlus = bosons.NbrUniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionUpPlus = bosons.UniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionSubArraySizeUpPlus = bosons.UniqueStateDescriptionSubArraySizeUpPlus;
  this->NbrUniqueStateDescriptionUpMinus = bosons.NbrUniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionUpMinus = bosons.UniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionSubArraySizeUpMinus = bosons.UniqueStateDescriptionSubArraySizeUpMinus;
  this->FirstIndexUniqueStateDescriptionUpMinus = bosons.FirstIndexUniqueStateDescriptionUpMinus;
  this->NbrUniqueStateDescriptionDownPlus = bosons.NbrUniqueStateDescriptionDownPlus;
  this->UniqueStateDescriptionDownPlus = bosons.UniqueStateDescriptionDownPlus;
  this->UniqueStateDescriptionSubArraySizeDownPlus = bosons.UniqueStateDescriptionSubArraySizeDownPlus;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusWithSU4Spin::Clone()
{
  return new BosonOnTorusWithSU4Spin(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithSU4Spin::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[state], this->StateDescriptionUpMinus[state], this->StateDescriptionDownPlus[state], this->StateDescriptionDownMinus[state],
		       this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, this->TemporaryStateDownPlus, this->TemporaryStateDownMinus); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = 0; i <= this->LzMax; ++i)
    {
      Str << "(" << this->TemporaryStateUpPlus[i] << "," << this->TemporaryStateUpMinus[i] << "," << this->TemporaryStateDownPlus[i] << "," << this->TemporaryStateDownMinus[i] << ") | ";
    }
  Str << " : " << this->FindStateIndex(this->StateDescriptionUpPlus[state], this->StateDescriptionUpMinus[state], this->StateDescriptionDownPlus[state], this->StateDescriptionDownMinus[state]) << " = " << state;
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKyUpPlus = current momentum along y for a single particle with up-plus
// currentKyUpMinus = current momentum along y for a single particle with up-minus
// currentKyDownPlus = current momentum along y for a single particle with down-plus
// currentKyDownMinus = current momentum along y for a single particle with down-minus
// currentTotalKy = current total momentum along y
// nbrNUpPlus = number of particles with quantum number up-plus
// nbrNUpMinus = number of particles with quantum number up-minus
// nbrNDownPlus = number of particles with quantum number down-plus
// nbrNDownMinus = number of particles with quantum number down-minus
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSU4Spin::GenerateStates(int nbrBosons, int currentKyUpPlus, int currentKyUpMinus, int currentKyDownPlus, int currentKyDownMinus, int currentTotalKy, 
					     int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus, long pos)
{
  if ((nbrBosons < 0) || (nbrNUpPlus < 0) || (nbrNUpMinus < 0) || (nbrNDownPlus < 0) || (nbrNDownMinus < 0))
    return pos;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	{
	  this->StateDescriptionUpPlus[pos] = 0x0ul;
	  this->StateDescriptionUpMinus[pos] = 0x0ul;
	  this->StateDescriptionDownPlus[pos] = 0x0ul;
	  this->StateDescriptionDownMinus[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else
	return pos;
    }
  if ((currentKyUpPlus < 0) || (currentKyUpMinus < 0) || (currentKyDownPlus < 0) || (currentKyDownMinus < 0))
    return pos;

  long TmpPos;
  for (int i = nbrNUpPlus; i >= 0; --i)    
    {
      unsigned long MaskUpPlus = ((0x1ul << i) - 1ul) << (currentKyUpPlus + nbrNUpPlus - i);
      for (int j = nbrNUpMinus; j >= 0; --j)
	{
    	  unsigned long MaskUpMinus = ((0x1ul << j) - 1ul) << (currentKyUpMinus + nbrNUpMinus - j);	  
	  for (int k = nbrNDownPlus; k >= 0; --k)
	    {
	      unsigned long MaskDownPlus = ((0x1ul << k) - 1ul) << (currentKyDownPlus + nbrNDownPlus - k);	  
	      for (int l = nbrNDownMinus; l >= 0; --l)
		{
		  unsigned long MaskDownMinus = ((0x1ul << l) - 1ul) << (currentKyDownMinus + nbrNDownMinus - l);	  
		  TmpPos = this->GenerateStates(nbrBosons - i - j - k - l, currentKyUpPlus - 1, currentKyUpMinus - 1, currentKyDownPlus - 1, currentKyDownMinus - 1, 
						currentTotalKy + (currentKyUpPlus * i) + (currentKyUpMinus * j) + (currentKyDownPlus * k) + (currentKyDownMinus * l), 
						nbrNUpPlus - i, nbrNUpMinus - j, nbrNDownPlus - k, nbrNDownMinus - l, pos); 
		  for (; pos < TmpPos; ++pos)
		    {
		      this->StateDescriptionUpPlus[pos] |= MaskUpPlus;
		      this->StateDescriptionUpMinus[pos] |= MaskUpMinus;
		      this->StateDescriptionDownPlus[pos] |= MaskDownPlus;
		      this->StateDescriptionDownMinus[pos] |= MaskDownMinus;
		    }
		}
	    }
	}
    }
  return pos;
};


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// nbrNUpPlus = number of particles with quantum number up-plus
// nbrNUpMinus = number of particles with quantum number up-minus
// nbrNDownPlus = number of particles with quantum number down-plus
// nbrNDownMinus = number of particles with quantum number down-minus
// return value = Hilbert space dimension

long BosonOnTorusWithSU4Spin::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus)
{
  if ((nbrBosons < 0) || (nbrNUpPlus < 0) || (nbrNUpMinus < 0) || (nbrNDownPlus < 0) || (nbrNDownMinus < 0))
    return 0l;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long Tmp = 0l;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->NbrLzValue) == this->KyMomentum)
	    ++Tmp;
	}
      return Tmp;
    }
  for (int i = nbrNUpPlus; i >= 0; --i)
    for (int j = nbrNUpMinus; j >= 0; --j)
      for (int k = nbrNDownPlus; k >= 0; --k)
	for (int l = nbrNDownMinus; l >= 0; --l)
	  Tmp += this->EvaluateHilbertSpaceDimension(nbrBosons - (i + j + k + l), currentKy - 1, currentTotalKy + (currentKy * (i + j + k + l)), 
						     nbrNUpPlus - i, nbrNUpMinus - j, nbrNDownPlus - k, nbrNDownMinus - l);
  return  Tmp;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrNUpPlusSector = number of particles with quantum number up-plus that belong to the subsytem 
// nbrNUpMinusSector = number of particles with quantum number up-minus that belong to the subsytem 
// nbrNDownPlusSector = number of particles with quantum number down-plus that belong to the subsytem 
// nbrNDownMinusSector = number of particles with quantum number down-plus that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnTorusWithSU4Spin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
											    int nbrNUpPlusSector, int nbrNUpMinusSector, int nbrNDownPlusSector, int nbrNDownMinusSector, 
											    RealVector& groundState, AbstractArchitecture* architecture)
{
  int NUpPlus = this->NbrBosons + this->TotalSpin + this->TotalIsospin + this->TotalEntanglement;
  int NUpMinus = this->NbrBosons + this->TotalSpin - this->TotalIsospin - this->TotalEntanglement;
  int NDownPlus = this->NbrBosons - this->TotalSpin + this->TotalIsospin - this->TotalEntanglement;
  int NDownMinus = this->NbrBosons - this->TotalSpin - this->TotalIsospin + this->TotalEntanglement;
  NUpPlus >>= 2;
  NUpMinus >>= 2;
  NDownPlus  >>= 2;
  NDownMinus  >>= 2;
 if (nbrParticleSector == 0)
    {
      if ((lzSector == 0) && (nbrNUpPlusSector == 0) && (nbrNUpMinusSector == 0) && (nbrNDownPlusSector == 0) && (nbrNDownMinusSector == 0))
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
       if ((lzSector == this->KyMomentum) && (nbrNUpPlusSector == NUpPlus) && (nbrNUpMinusSector == NUpMinus) && (nbrNDownPlusSector == NDownPlus) && (nbrNDownMinusSector == NDownMinus))
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
  int ComplementaryKyMomentum = (this->KyMomentum - lzSector) % this->NbrLzValue;
  int ComplementaryNUpPlusSector = NUpPlus - nbrNUpPlusSector;
  int ComplementaryNUpMinusSector = NUpMinus - nbrNUpMinusSector;
  int ComplementaryNDownPlusSector = NDownPlus - nbrNDownPlusSector;
  int ComplementaryNDownMinusSector = NDownMinus - nbrNDownMinusSector;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrLzValue;

  int SubsytemTotalSz = (nbrNUpPlusSector + nbrNUpMinusSector) - (nbrNDownPlusSector + nbrNDownMinusSector);
  int SubsytemTotalIz = (nbrNUpPlusSector + nbrNDownPlusSector) - (nbrNUpMinusSector + nbrNDownMinusSector);
  int SubsytemTotalPz = (nbrNUpPlusSector + nbrNDownMinusSector) - (nbrNUpMinusSector + nbrNDownPlusSector);
  int ComplementaryTotalSz = (ComplementaryNUpPlusSector + ComplementaryNUpMinusSector) - (ComplementaryNDownPlusSector + ComplementaryNDownMinusSector);
  int ComplementaryTotalIz = (ComplementaryNUpPlusSector + ComplementaryNDownPlusSector) - (ComplementaryNUpMinusSector + ComplementaryNDownMinusSector);
  int ComplementaryTotalPz = (ComplementaryNUpPlusSector + ComplementaryNDownMinusSector) - (ComplementaryNUpMinusSector + ComplementaryNDownPlusSector);
  cout << "ky = " << this->KyMomentum << " " << lzSector << " " << ComplementaryKyMomentum << endl;
  cout << "nupplus = " << NUpPlus << " " << nbrNUpPlusSector << " " << ComplementaryNUpPlusSector << endl;
  cout << "nupminus = " << NUpMinus << " " << nbrNUpMinusSector << " " << ComplementaryNUpMinusSector << endl;
  cout << "ndownplus = " << NDownPlus << " " << nbrNDownPlusSector << " " << ComplementaryNDownPlusSector << endl;
  cout << "ndownminus = " << NDownPlus << " " << nbrNDownPlusSector << " " << ComplementaryNDownPlusSector << endl;
  BosonOnTorusWithSU4Spin SubsytemSpace (nbrParticleSector, SubsytemTotalSz, SubsytemTotalIz, SubsytemTotalPz, this->NbrLzValue, lzSector);
  RealSymmetricMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnTorusWithSU4Spin ComplementarySpace (ComplementaryNbrParticles, ComplementaryTotalSz, ComplementaryTotalIz, ComplementaryTotalPz, this->NbrLzValue, ComplementaryKyMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnTorusWithSU4Spin::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
										 RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
  BosonOnTorusWithSU4Spin* TmpHilbertSpace =  (BosonOnTorusWithSU4Spin*) complementaryHilbertSpace;
  BosonOnTorusWithSU4Spin* TmpDestinationHilbertSpace =  (BosonOnTorusWithSU4Spin*) destinationHilbertSpace;
  int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
  int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  
  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
  double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];
  for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescriptionUpPlus[i], TmpDestinationHilbertSpace->StateDescriptionUpMinus[i], TmpDestinationHilbertSpace->StateDescriptionDownPlus[i], TmpDestinationHilbertSpace->StateDescriptionDownMinus[i], TmpDestinationHilbertSpace->TemporaryStateUpPlus, TmpDestinationHilbertSpace->TemporaryStateUpMinus, TmpDestinationHilbertSpace->TemporaryStateDownPlus, TmpDestinationHilbertSpace->TemporaryStateDownMinus); 

      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->LzMax; ++k)
	{
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryStateUpPlus[k]];
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryStateUpMinus[k]];
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryStateDownPlus[k]];
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryStateDownMinus[k]];
	}
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }

  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->StateDescriptionUpPlus[minIndex], TmpHilbertSpace->StateDescriptionUpMinus[minIndex], TmpHilbertSpace->StateDescriptionDownPlus[minIndex], TmpHilbertSpace->StateDescriptionDownMinus[minIndex], 
				      TmpHilbertSpace->TemporaryStateUpPlus, TmpHilbertSpace->TemporaryStateUpMinus, TmpHilbertSpace->TemporaryStateDownPlus, TmpHilbertSpace->TemporaryStateDownMinus);
       double TmpHilbertSpaceFactorial = 0.0;
       for (int k = 0; k <= TmpHilbertSpace->LzMax; ++k)
	 {
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryStateUpPlus[k]];
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryStateUpMinus[k]];
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryStateDownPlus[k]];
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryStateDownMinus[k]];
	 }
       for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	 {
	   TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescriptionUpPlus[j], TmpDestinationHilbertSpace->StateDescriptionUpMinus[j], TmpDestinationHilbertSpace->StateDescriptionDownPlus[j], TmpDestinationHilbertSpace->StateDescriptionDownMinus[j], 
						      TmpDestinationHilbertSpace->TemporaryStateUpPlus, TmpDestinationHilbertSpace->TemporaryStateUpMinus, TmpDestinationHilbertSpace->TemporaryStateDownPlus, TmpDestinationHilbertSpace->TemporaryStateDownMinus);
	   for (int k = 0; k <=  TmpDestinationHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryStateUpPlus[k] = TmpDestinationHilbertSpace->TemporaryStateUpPlus[k];
	       this->TemporaryStateUpMinus[k] = TmpDestinationHilbertSpace->TemporaryStateUpMinus[k];
	       this->TemporaryStateDownPlus[k] = TmpDestinationHilbertSpace->TemporaryStateDownPlus[k];
	       this->TemporaryStateDownMinus[k] = TmpDestinationHilbertSpace->TemporaryStateDownMinus[k];
	     }
	   for (int k = TmpDestinationHilbertSpace->LzMax + 1; k <=  this->LzMax; ++k)
	     {
	       this->TemporaryStateUpPlus[k] = 0x0ul;
	       this->TemporaryStateUpMinus[k] = 0x0ul;
	       this->TemporaryStateDownPlus[k] = 0x0ul;
	       this->TemporaryStateDownMinus[k] = 0x0ul;
	     }	   
	   for (int k = 0; k <=  TmpHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryStateUpPlus[k] += TmpHilbertSpace->TemporaryStateUpPlus[k];
	       this->TemporaryStateUpMinus[k] += TmpHilbertSpace->TemporaryStateUpMinus[k];
	       this->TemporaryStateDownPlus[k] += TmpHilbertSpace->TemporaryStateDownPlus[k];
	       this->TemporaryStateDownMinus[k] += TmpHilbertSpace->TemporaryStateDownMinus[k];
	     }

	   int TmpPos = this->FindStateIndex(this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, this->TemporaryStateDownPlus, this->TemporaryStateDownMinus);
	   if (TmpPos != this->HilbertSpaceDimension)
	     {
	       double TmpFactorial = 0.0;	      
	       for (int k = 0; k <= this->LzMax; ++k)
		 {
		   TmpFactorial += LogFactorials[this->TemporaryStateUpPlus[k]];
		   TmpFactorial += LogFactorials[this->TemporaryStateUpMinus[k]];
		   TmpFactorial += LogFactorials[this->TemporaryStateDownPlus[k]];
		   TmpFactorial += LogFactorials[this->TemporaryStateDownMinus[k]];
		 }
	       TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	       TmpFactorial *= 0.5; 
	       
	       TmpStatePosition[Pos] = TmpPos;
	       TmpStatePosition2[Pos] = j;
	       TmpStateCoefficient[Pos] = exp(TmpFactorial);
	       ++Pos;
	     }
	 }
       if (Pos != 0)
 	{
 	  ++TmpNbrNonZeroElements;
 	  for (int j = 0; j < Pos; ++j)
 	    {
 	      int Pos2 = TmpStatePosition2[j];
 	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
 	      for (int k = 0; k < Pos; ++k)
 		if (TmpStatePosition2[k] >= Pos2)
 		  {
 		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
 		  }
 	    }
 	}
     }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
  
}

