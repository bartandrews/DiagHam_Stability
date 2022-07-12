////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//      NbrStates + NbrBosons - 1 < 63 or 31 (64 bit or 32 bit systems)       //
//                                                                            //
//                        last modification : 10/02/2008                      //
//                                                                            //
//                                                                            //
//    this program is free software; you can redistribute it and/or modify    //
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
#include "HilbertSpace/BosonOnLatticeKyLong.h"
#include "HilbertSpace/BosonOnLattice.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "Architecture/ArchitectureOperation/FQHELatticeParticleEntanglementSpectrumOperation.h"

#include <cmath>
#include <cstdlib>

#include <bitset>

using std::bitset;


using std::cout;
using std::endl;

// switch some debugging output
// #define DEBUG_OUTPUT


// default constructor
//

BosonOnLatticeKyLong::BosonOnLatticeKyLong ()
{
  this->HilbertSpaceDimension=0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// basic constructor -> yields a square lattice in Landau gauge
// 
// nbrBosons = number of bosons
// lx = length of simulation cell in x-direction
// ly = length of simulation cell in y-direction
// ky = many-body momentum in y-direction
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// memory = memory that can be allocated for precalculations
BosonOnLatticeKyLong::BosonOnLatticeKyLong (int nbrBosons, int lx, int ly, int ky, int nbrFluxQuanta, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->Lx = lx;
  this->Ly = ly;
  this->NbrSublattices = 1;  
  this->NbrStates = Lx*Ly;
  
  this->NbrCodingBits = NbrStates + NbrBosons - 1;
  
#ifdef __128_BIT_LONGLONG__
  if (NbrCodingBits>128)    
#else
    if (NbrCodingBits>64)    
#endif
      {
	cout<<"BosonOnLatticeKyLong: Cannot represent the "<<NbrCodingBits<<" coding bits requested in a single word"<<endl;
	exit(1);
      }
  
  
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->FluxDensity = ((double)NbrFluxQuanta)/this->NbrStates;
  this->LxTranslationPhase = Polar(1.0, -2.0*M_PI*FluxDensity*this->Lx);
  this->LyTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge...
  
  for (int i=1; i<=Ly; ++i)
    if ((i*NbrFluxQuanta)%Ly == 0)
      {
	this->TranslationCell = i;
	break;
      }
  if (Ly % TranslationCell != 0)
    {
      cout << "System has no translation symmetries along the y-direction, please calculate in real-space!"<<endl;
    }
  // maximum many-body momentum (Ly / TranslationCell)
  this->Kmax = Ly / TranslationCell;
  
  if( ky >= Kmax )
    {
      cout << "Attention, requested total momentum Ky="<<ky<<" is outside of fundamental BZ!"<<endl;
    }
  this->Ky = ky % Kmax;
  this->HilbertSpaceDimension = EvaluateHilbertSpaceDimension(nbrBosons, NbrStates-1, NbrStates-1, 0);
  
  this->StateDescription=new ULONGLONG [this->HilbertSpaceDimension];
  this->StateHighestBit=new int[this->HilbertSpaceDimension];
  
  this->TemporaryState = new unsigned long [this->NbrStates];
  this->ShiftedState = new unsigned long [this->NbrStates];
  this->ProdATemporaryState = new unsigned long [this->NbrStates];
  this->Flag.Initialize();
  
  this->GenerateStates(nbrBosons, NbrStates-1, NbrStates-1, 0, 0);
  
  this->TargetSpace=this;
  this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  
#ifdef DEBUG_OUTPUT
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    {
      PrintState(cout, i);
      cout << endl;
    }
#endif
  
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnLatticeKyLong::BosonOnLatticeKyLong (const BosonOnLatticeKyLong& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->Lx = bosons.Lx;
  this->Ly = bosons.Ly;
  this->NbrSublattices = bosons.NbrSublattices;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->FluxDensity = bosons.FluxDensity;
  this->LxTranslationPhase = bosons.LxTranslationPhase;
  this->LyTranslationPhase = bosons.LyTranslationPhase;
  this->NbrStates = bosons.NbrStates;
  this->NbrCodingBits = bosons.NbrCodingBits;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;

  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrStates];
  this->ShiftedState = new unsigned long [this->NbrStates];
  this->ProdATemporaryState = new unsigned long [this->NbrStates];  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

BosonOnLatticeKyLong::~BosonOnLatticeKyLong ()
{
  if (this->HilbertSpaceDimension != 0) 
    {
      if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete[] this->StateDescription;
	  if (this->StateHighestBit != 0)
	    delete[] this->StateHighestBit;
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < this->NbrCodingBits; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	}
      delete[] this->TemporaryState;
      delete[] this->ShiftedState;
      delete[] this->ProdATemporaryState;
    }
  
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeKyLong& BosonOnLatticeKyLong::operator = (const BosonOnLatticeKyLong& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrCodingBits; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
      delete[] this->TemporaryState;
      delete[] this->ShiftedState;
      delete[] this->ProdATemporaryState;
    }
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->Lx = bosons.Lx;
  this->Ly = bosons.Ly;
  this->NbrSublattices = bosons.NbrSublattices;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->FluxDensity = bosons.FluxDensity;
  this->LxTranslationPhase = bosons.LxTranslationPhase;
  this->LyTranslationPhase = bosons.LyTranslationPhase;
  this->NbrStates = bosons.NbrStates;
  this->NbrCodingBits = bosons.NbrCodingBits;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrStates];
  this->ShiftedState = new unsigned long [this->NbrStates];
  this->ProdATemporaryState = new unsigned long [this->NbrStates];  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeKyLong::Clone()
{
  return new BosonOnLatticeKyLong(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnLatticeKyLong::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(this->NbrBosons);
  TmpList += new PeriodicMomentumQuantumNumber(this->Ky, this->Kmax);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnLatticeKyLong::GetQuantumNumber (int index)
{
  return new NumberParticleQuantumNumber(this->NbrBosons);
}

// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id
int BosonOnLatticeKyLong::GetHilbertSpaceAdditionalSymmetry()
{
  return ParticleOnLattice::YTranslations;
}


// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnLatticeKyLong::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  return 0;
}

// get the number of sites
//
// return value = number of sites
int BosonOnLatticeKyLong::GetNbrSites()
{
  return this->NbrStates;
}

// it is possible to change the flux through the simulation cell
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
void BosonOnLatticeKyLong::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  cout << "Attention, do not use SetNbrFluxQuanta for Bosons on Lattice with translation symmetry!"<<endl;
  exit(-1);
}

// change flux through cell and periodic boundary conditions
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
void BosonOnLatticeKyLong::SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY)
{
  cout << "Attention, do not use SetNbrFluxQuanta for Bosons on Lattice with translation symmetry!"<<endl;
  exit(-1);
}


// request solenoid fluxes
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
//
void BosonOnLatticeKyLong::GetSolenoidFluxes(double &solenoidX, double &solenoidY)
{
  solenoidX=0.0;
  solenoidY=0.0;
}

// obtain the current setting of the flux piercing this lattice
int BosonOnLatticeKyLong::GetNbrFluxQuanta()
{
  return this->NbrFluxQuanta;
}


// get maximum possible momentum for this geometry
// return = maximum value of Ky
int BosonOnLatticeKyLong::GetMaximumKy()
{
  return this->Kmax;
}

// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added

long unsigned int BosonOnLatticeKyLong::Ad (long unsigned int state, int q, double &coefficient)
{
  cout <<"Using wrong function BosonOnLatticeKyLong::Ad"<<endl;
  exit(1);
  return 0x0ul;
}

// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
ULONGLONG BosonOnLatticeKyLong::Ad (ULONGLONG state, int q, double &coefficient)
{
  int StateHighestBit = getHighestBit(state)-1;
  this->FermionToBoson(state, StateHighestBit, this->TemporaryState, this->TemporaryStateHighestBit);
  if (q>this->TemporaryStateHighestBit)
    {
      for (int i = this->TemporaryStateHighestBit + 1; i < q; ++i)
	this->TemporaryState[i] = 0ul;
      this->TemporaryState[q] = 1ul;
      this->TemporaryStateHighestBit = q;
      coefficient = 1.0;
    }
  else
    {
      ++this->TemporaryState[q];
      coefficient = sqrt((double)this->TemporaryState[q]);
    }
  return this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit);
}


// apply a^+_q1 a^+_q2 a_r1 a_r2 operator to a given state (with q1+q2=n1+r2)
//
// index = index of the state on which the operator has to be applied
// q1 = first index for creation operator
// q2 = second index for creation operator
// r1 = first index for annihilation operator
// r2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLatticeKyLong::AdAdAA (int index, int q1, int q2, int r1, int r2, double& coefficient)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  
  if ((r1 > this->TemporaryStateHighestBit) || (r2 > this->TemporaryStateHighestBit) || (this->TemporaryState[r1] == 0) || (this->TemporaryState[r2] == 0) || ((r1 == r2) && (this->TemporaryState[r1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = this->TemporaryState[r2];
  --this->TemporaryState[r2];
  if ((r2==TemporaryStateHighestBit)&&(this->TemporaryState[r2]==0)&&(TemporaryStateHighestBit>0))
    {
      --TemporaryStateHighestBit;
      while ((TemporaryStateHighestBit>0)&&(TemporaryState[TemporaryStateHighestBit] == 0))
	--TemporaryStateHighestBit;
    }
  coefficient *= this->TemporaryState[r1];
  --this->TemporaryState[r1];
  if ((r1==TemporaryStateHighestBit)&&(this->TemporaryState[r1]==0)&&(TemporaryStateHighestBit>0))
    {
      --TemporaryStateHighestBit;
      while ((TemporaryStateHighestBit>0)&&(TemporaryState[TemporaryStateHighestBit] == 0))
	--TemporaryStateHighestBit;
    }
  if (q2 > TemporaryStateHighestBit)
    {
      for (int i = this->TemporaryStateHighestBit + 1; i <= q2; ++i)
	this->TemporaryState[i] = 0ul;
      TemporaryStateHighestBit = q2;      
    }
  ++this->TemporaryState[q2];
  coefficient *= this->TemporaryState[q2];
  if (q1 > TemporaryStateHighestBit)
    {
      for (int i = this->TemporaryStateHighestBit + 1; i <= q1; ++i)
	this->TemporaryState[i] = 0ul;
      TemporaryStateHighestBit = q1;
    }
  ++this->TemporaryState[q1];
  coefficient *= this->TemporaryState[q1];  
  coefficient = sqrt(coefficient);

  //   cout << "Checking correspondance: " <<endl;  
  //   for (int i=0; i<=TemporaryStateHighestBit; ++i) cout << " "<< this->TemporaryState[i];
//   cout <<endl;
//   int tmp = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
//   PrintState(cout,tmp); cout << " prefactor: "<<coefficient<<endl;

  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
}


// apply a_r1 a_r2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// r1 = first index for annihilation operator
// r2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeKyLong::AA (int index, int r1, int r2)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->ProdATemporaryState, this->ProdATemporaryStateHighestBit);
  
  if ((r1 > this->ProdATemporaryStateHighestBit) || (r2 > this->ProdATemporaryStateHighestBit) || (this->ProdATemporaryState[r1] == 0) || (this->ProdATemporaryState[r2] == 0) || ((r1 == r2) && (this->ProdATemporaryState[r1] == 1)))
    {
      return 0.0;
    }
  for (int i = this->ProdATemporaryStateHighestBit + 1; i < this->NbrStates; ++i)
    this->ProdATemporaryState[i] = 0ul;
  
  double Coefficient = this->ProdATemporaryState[r2];
  --this->ProdATemporaryState[r2];
  if ((r2==this->ProdATemporaryStateHighestBit)&&(this->ProdATemporaryState[r2]==0)&&(this->ProdATemporaryStateHighestBit>0))
    {
      --this->ProdATemporaryStateHighestBit;
      while (ProdATemporaryState[this->ProdATemporaryStateHighestBit] == 0)
	--ProdATemporaryStateHighestBit;
    }
  Coefficient *= this->ProdATemporaryState[r1];
  --this->ProdATemporaryState[r1];
  if ((r1==this->ProdATemporaryStateHighestBit)&&(this->ProdATemporaryState[r1]==0)&&(this->ProdATemporaryStateHighestBit>0))
    {
      --this->ProdATemporaryStateHighestBit;
      while (ProdATemporaryState[this->ProdATemporaryStateHighestBit] == 0)
	--this->ProdATemporaryStateHighestBit;
    }
  return sqrt(Coefficient);
}


// apply a^+_q1 a^+_q2 operator to the state produced using AA method (without destroying it)
//
// q1 = first index for creation operator
// q2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLatticeKyLong::AdAd (int q1, int q2, double& coefficient)
{
  for (int i = 0; i < this->NbrStates; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  
  this->TemporaryStateHighestBit = this->ProdATemporaryStateHighestBit;
  ++this->TemporaryState[q2];
  coefficient = this->TemporaryState[q2];
  if (q2 > TemporaryStateHighestBit)
    TemporaryStateHighestBit = q2;
  ++this->TemporaryState[q1];
  coefficient *= this->TemporaryState[q1];
  if (q1 > TemporaryStateHighestBit)
    TemporaryStateHighestBit = q1;
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor   

double BosonOnLatticeKyLong::ProdA (int index, int* n, int nbrIndices)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->ProdATemporaryState, this->ProdATemporaryStateHighestBit);
  
  for (int i = this->ProdATemporaryStateHighestBit + 1; i < this->NbrStates; ++i)
    this->ProdATemporaryState[i] = 0;
  
  unsigned long TmpCoefficient = 1ul;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (n[i] > this->ProdATemporaryStateHighestBit)
	return 0.0;
      unsigned long& Tmp = this->ProdATemporaryState[n[i]];
      if (Tmp == 0)
	return 0.0;
      TmpCoefficient *= Tmp;
      --Tmp;
      if ((n[i] == this->ProdATemporaryStateHighestBit)&&(this->ProdATemporaryState[n[i]]==0)&&(this->ProdATemporaryStateHighestBit>0))
	{
	  --this->ProdATemporaryStateHighestBit;
	  while (this->ProdATemporaryState[this->ProdATemporaryStateHighestBit] == 0)
	    this->ProdATemporaryStateHighestBit--;
	}
    }
  return sqrt((double) TmpCoefficient);
}
  
// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLatticeKyLong::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  this->TemporaryStateHighestBit=this->ProdATemporaryStateHighestBit;
  for (int i = 0; i < this->NbrStates; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  
  int TmpCoefficient = 1;
  for (int i = 0; i < nbrIndices; ++i)
    {
      if(m[i] > this->TemporaryStateHighestBit)
	this->TemporaryStateHighestBit = m[i];
      TmpCoefficient *= ++this->TemporaryState[m[i]];
    }
  coefficient = sqrt((double) TmpCoefficient);
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
}

// apply a^+_q a_q operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_q a_q

double BosonOnLatticeKyLong::AdA (int index, int q)
{ 
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  if (q>this->TemporaryStateHighestBit)
    return 0.0;
  else
    return (double) (this->TemporaryState[q]);  
}

// apply a^+_q a_r operator to a given state 
//
// index = index of the state on which the operator has to be applied
// q = index of the creation operator
// r = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnLatticeKyLong::AdA (int index, int q, int r, double& coefficient)
{
  // cout << "AdA: " << index << ", " << q << ", " << r <<endl;
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
//   cout << "Temporary state before: ";
//   for (int i=0; i<NbrStates; ++i)
//      cout << TemporaryState[i]<<" ";
//   cout << ", h.b.=" <<TemporaryStateHighestBit<< endl;
  // bitset<32>b2=this->StateDescription[index];
  // cout << "converted: "<<b2 << endl;
  if (q==r)
    {
      if (q>this->TemporaryStateHighestBit)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      else
	{
	  coefficient = (double) (this->TemporaryState[q]);
	  return index;
	}
    }
  if ((r > this->TemporaryStateHighestBit) || (this->TemporaryState[r] == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = this->TemporaryState[r];
  --this->TemporaryState[r];
  
  if ((TemporaryStateHighestBit==r)&&(this->TemporaryState[r]==0)&&(TemporaryStateHighestBit>0))
    {
      --TemporaryStateHighestBit;
      while ((TemporaryState[TemporaryStateHighestBit] == 0)&&(TemporaryStateHighestBit>0))
	--TemporaryStateHighestBit;
    }
  if (q > TemporaryStateHighestBit)
    {
      for (int i = this->TemporaryStateHighestBit + 1; i <= q; ++i)
	this->TemporaryState[i] = 0ul;
      TemporaryStateHighestBit = q;
    }
  ++this->TemporaryState[q];
  coefficient *= this->TemporaryState[q];
  coefficient = sqrt(coefficient);

//   cout << "Temporary state after:  ";
//   for (int i=0; i<NbrStates; ++i)
//     cout << TemporaryState[i]<<" ";
//   cout << ", h.b.=" <<TemporaryStateHighestBit<< endl;
//   unsigned long tmpL=this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit);
//   bitset<32>b=tmpL;
//   cout << "converted: "<<b << endl;
  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
}

// apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
// index = index of the state on which the operator has to be applied
// NbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
// QValues = array of quantum numbers where an interaction is present
// InteractionPerQ = coefficient U_q of the interaction
//
double BosonOnLatticeKyLong::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  int BosonsLeft=this->NbrBosons;
  double result=0.0;
  if (nbrInteraction==NbrStates)
    {
      int Q=TemporaryStateHighestBit;
      while ((Q>-1) && (BosonsLeft>0))
	{
	  if (this->TemporaryState[Q]!=0)
	    {
	      result+=interactionPerQ[Q]*TemporaryState[Q]*(TemporaryState[Q]-1);
	      BosonsLeft-=TemporaryState[Q];
	    }
	  --Q;
	}
    }
  else // cannot assume ordered set of qValues
    {
      for (int i=0; (i<nbrInteraction)&&(BosonsLeft>0); ++i)
	{
	  int Q=qValues[i];
	  result+=interactionPerQ[Q]*TemporaryState[Q]*(TemporaryState[Q]-1);
	  BosonsLeft-=TemporaryState[Q];	  
	}
    }
  
//   cout << "state:";
//   for (int i=0; i<NbrStates; ++i)
//     cout << " " << TemporaryState[i];
//   cout << " (hb "<<TemporaryStateHighestBit<<") -> result = " << result<<endl;
  
  return result;
}



// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxQ = maximum value for the quantum number of a boson in the state
// currentMaxQ = current max value for the quantum number of bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

int BosonOnLatticeKyLong::GenerateStates(int nbrBosons, int maxQ, int currentMaxQ, int pos, int currentMomentum)
//, int debugLevel)
{
  if (nbrBosons == 0)
    {
      // create an entry for a state and set it to zero
      if ((currentMomentum % this->Kmax) == this->Ky)
	{	  
	  this->StateHighestBit[pos] = 0;
	  this->StateDescription[pos] =  ((ULONGLONG)0x0l);	  
	  return pos + 1;
	}
      else
	{
	  return pos;
	}
    }
  if (currentMaxQ == 0) // assumes k_y(Q=0)=0
    {
      // create an entry for the current state and fill the lowest bit with the remaining bosons
      if ((currentMomentum % this->Kmax) == this->Ky)
	{
	  TemporaryState[0] = nbrBosons;
	  TemporaryStateHighestBit = 0;
	  this->StateDescription[pos] = this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit, this->StateHighestBit[pos]);
	  return pos + 1;
	}
      else
	{
	  return pos;
	}
    }

  int TmpNbrBosons = 0;
  int ReducedCurrentMaxQ = currentMaxQ - 1;
  int TmpPos = pos;
  int currentQMomentum = this->DecodeKy(currentMaxQ);
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = this->GenerateStates(TmpNbrBosons, maxQ, ReducedCurrentMaxQ, pos, currentMomentum + (nbrBosons - TmpNbrBosons) * currentQMomentum);
      
      for (int i = pos; i < TmpPos; i++)
	{
	  
	  this->FermionToBoson(this->StateDescription[i], this->StateHighestBit[i], this->TemporaryState, this->TemporaryStateHighestBit);
	  
	  for (int bit=TemporaryStateHighestBit+1; bit<currentMaxQ; ++bit)
	    this->TemporaryState[bit] = 0;	  
	  
	  this->TemporaryState[currentMaxQ] = nbrBosons - TmpNbrBosons;
	  this->StateDescription[i] = this->BosonToFermion(this->TemporaryState, currentMaxQ, this->StateHighestBit[i]);  
	}
      
      ++TmpNbrBosons;
      pos = TmpPos;
    }
  
  if (maxQ == currentMaxQ)
    return this->GenerateStates(nbrBosons, ReducedCurrentMaxQ, ReducedCurrentMaxQ, pos, currentMomentum);
  else
    {
      return this->GenerateStates(nbrBosons, maxQ, ReducedCurrentMaxQ, pos, currentMomentum);
    }
}



// recursively evaluate Hilbert space dimension 
//
// nbrBosons = number of bosons
// maxQ = maximum value for the quantum number of a boson in the state
// currentMaxQ = current max value for the quantum number of bosons that are still to be placed
// currentMomentum = current value of the momentum
// return value = Hilbert space dimension

int BosonOnLatticeKyLong::EvaluateHilbertSpaceDimension(int nbrBosons, int maxQ, int currentMaxQ, int currentMomentum)
{
  if ((nbrBosons == 0)||(currentMaxQ == 0))
    {
      if ((currentMomentum % this->Kmax) == this->Ky)
	{
	  return 1;
	}
      else
	return 0;
    }
  
  int TmpNbrBosons = 0;
  int ReducedCurrentMaxQ = currentMaxQ - 1;
  int currentQMomentum = this->DecodeKy(currentMaxQ);
  int rst=0;
  while (TmpNbrBosons < nbrBosons)
    {
      rst += this->EvaluateHilbertSpaceDimension(TmpNbrBosons, maxQ, ReducedCurrentMaxQ, currentMomentum + (nbrBosons - TmpNbrBosons) * currentQMomentum);
      ++TmpNbrBosons;
    }
  if (maxQ == currentMaxQ)
    {
      return rst + this->EvaluateHilbertSpaceDimension(nbrBosons, ReducedCurrentMaxQ, ReducedCurrentMaxQ, currentMomentum);
    }
  else
    {
      return  rst + this->EvaluateHilbertSpaceDimension(nbrBosons, maxQ, ReducedCurrentMaxQ, currentMomentum);
    }
}

// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// ky = codes momentum in y-direction and 'subcell' index
// sublattice = sublattice index
// 

int BosonOnLatticeKyLong::EncodeQuantumNumber(int posx, int ky, int sublattice, Complex &translationPhase)
{
  //cout << "Encoding " << posx<<", "<<posy<<": ";
  int numXTranslations=0, numYTranslations=0;  
  while (posx<0)
    {
      posx+=this->Lx;
      ++numXTranslations;
    }
  while (posx>=this->Lx)
    {
      posx-=this->Lx;
      --numXTranslations;
    }
  while (ky<0)
    {
      ky+=this->Ly;
      ++numYTranslations;      
    }
  while (ky>=this->Ly)
    {
      ky-=this->Ly;
      --numYTranslations;
    }
  int rst = posx + this->Lx*ky;
  rst*=this->NbrSublattices;
  rst+=sublattice;
  // determine phase for shifting site to the simulation cell:
  Complex tmpPhase(1.0,0.0);
  Complex tmpPhase2;
  translationPhase=tmpPhase;
  if (numXTranslations>0)
    tmpPhase2=LxTranslationPhase;
  else
    tmpPhase2=Conj(LxTranslationPhase);
  for (int i=0; i<abs(numXTranslations); ++i)
    tmpPhase*=tmpPhase2;


  for (int y=1;y<=(ky%TranslationCell); ++y)
    translationPhase*=tmpPhase;
  // in Landau gauge chosen here, translations in y-direction are trivial...
//   tmpPhase=1.0;
//   if (numYTranslations>0)
//     tmpPhase2=LyTranslationPhase;
//   else
//     tmpPhase2=Conj(LyTranslationPhase);
//   for (int i=0; i<abs(numYTranslations); ++i)
//     tmpPhase*=tmpPhase2;
//   //cout<<" tmpPhaseY="<<tmpPhase;
//   for (int x=1;x<=posx; ++x)
//     translationPhase*=tmpPhase;
  //cout << "tX="<<numXTranslations<< ", tY="<<numYTranslations<<", translationPhase= " <<translationPhase<<endl;
  return rst;
}

// decode a single encoded quantum number q to the set of quantum numbers posx, posy
// posx = position along x-direction
// posy = codes momentum in y-direction and 'subcell' index

void BosonOnLatticeKyLong::DecodeQuantumNumber(int q, int &posx, int &ky, int &sublattice)
{
  int m = this->NbrSublattices;
  sublattice=q%m;
  posx=(q%(m*Lx))/m;
  m*=Lx;
  ky=(q%(m*Ly))/m;
}

// extract the momentum ky from a quantum number q
// return: momentum ky (in range 0...Kmax-1)

int BosonOnLatticeKyLong::DecodeKy(int q)
{
  int posx, ky, subl, m=this->NbrSublattices;
  subl=q%m;
  posx=(q%(m*Lx))/m;
  m*=Lx;
  ky=(q%(m*Ly))/m;
  return ky/TranslationCell;
}


// ky = true momentum in y-direction
// fluxSubLattice = 'sublattice' index remaining after translation symmetry

int BosonOnLatticeKyLong::EncodeCompositeMomentum(int ky, int fluxSubLattice)
{
  return ky*TranslationCell+fluxSubLattice;
}

// decode composite ky-momentum
// cK = composite momentum (momentum plus flux sublattice)
// ky = true momentum in y-direction
// fluxSubLattice = 'sublattice' index remaining after translation symmetry

void BosonOnLatticeKyLong::DecodeCompositeMomentum(int cK, int &ky, int &fluxSubLattice)
{
  ky = cK/TranslationCell;
  fluxSubLattice = cK%TranslationCell;
}

// check whether HilbertSpace implements ordering of operators
//

bool BosonOnLatticeKyLong::HaveOrder ()
{
  return true;
}


// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value

int BosonOnLatticeKyLong::CheckOrder (int* m, int* n, int nbrIndices)
{
  if (nbrIndices>this->NbrBosons)
    return 0;
  this->TemporaryStateHighestBit=0;
  for (int i=0; i<this->NbrStates; ++i)
    this->TemporaryState[i]=0;
  for (int i=0; i<nbrIndices; ++i)
    {
      ++this->TemporaryState[m[i]];
      if (m[i]>TemporaryStateHighestBit)
	TemporaryStateHighestBit=m[i];
    }
  
  ULONGLONG CreationValue=this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit);
  this->TemporaryStateHighestBit=0;
  for (int i=0; i<this->NbrStates; ++i)
    this->TemporaryState[i] = 0;
  for (int i=0; i<nbrIndices; ++i)
    {
      ++this->TemporaryState[n[i]];
      if (n[i]>TemporaryStateHighestBit)
	TemporaryStateHighestBit=n[i];
    }
  ULONGLONG AnnihilationValue=this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit);
 
  if (CreationValue > AnnihilationValue)
    return 1;
  else if (CreationValue < AnnihilationValue)
    return -1;
  else return 0;
}

// obtain a list of quantum numbers in state
// index = index of many-body state to be considered
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
// normalization = indicating the multiplicity of the state for bosonic spaces

void BosonOnLatticeKyLong::ListQuantumNumbers(int index, int *quantumNumbers, double &normalization)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  int NbrQ = 0;
  normalization = 1.0;
  
  for (int q = 0; q <= TemporaryStateHighestBit; ++q)
    if (this->TemporaryState[q] > 0)
      {		
	for (unsigned int l=this->TemporaryState[q];l>0; --l)
	  {
	    normalization *= l;
	    quantumNumbers[NbrQ]=q;
	    ++NbrQ;
	  }
      }
  normalization=sqrt(normalization);
}

// obtain a list of quantum numbers in state
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles

void BosonOnLatticeKyLong::ListQuantumNumbers(int index, int *quantumNumbers)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  int NbrQ = 0;
  for (int q=0; q<=TemporaryStateHighestBit; ++q)
    if (this->TemporaryState[q]>0)
      {		
	for (unsigned int l=this->TemporaryState[q];l>0; --l)
	  {
	    quantumNumbers[NbrQ]=q;
	    ++NbrQ;
	  }
      }
}



// translate a state by a multiple of the lattice vectors
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// translationPhase = returns phase inccurred by translation
// return value = index of translated state

int BosonOnLatticeKyLong::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  int BosonsLeft=this->NbrBosons;
  int Q = TemporaryStateHighestBit;
  int OldX, OldY, OldSl;
  int NewQ;
  int CountYCoordinates = 0; // total phase is shiftX * sum_i y_i in Landau gauge
  Complex PeriodicPhase;
  Complex CumulatedPhase = 1.0;

//cout << "TS:";
//for (int i=0; i<=TemporaryStateHighestBit; ++i)
//  cout << " "<<TemporaryState[i];
//cout << endl;

  this->ShiftedStateHighestBit=0;
  this->ShiftedState[0]=0;
  while ((Q>-1) && (BosonsLeft>0))
    {
      if (this->TemporaryState[Q] != 0)
	{
	  this->DecodeQuantumNumber(Q,OldX, OldY, OldSl);
	  CountYCoordinates+=this->TemporaryState[Q]*OldY;
	  NewQ=this->EncodeQuantumNumber(OldX+shiftX, OldY+shiftY, OldSl, PeriodicPhase);

//cout << "PeriodicPhase for shift ("<<OldX<<","<<OldY<<")->("<<OldX+shiftX<<","<<OldY+shiftY<<")="<< PeriodicPhase<<endl;
	  
	  for (unsigned i = 0; i<this->TemporaryState[Q]; ++i)
	    CumulatedPhase *= PeriodicPhase;
	  if (NewQ > ShiftedStateHighestBit)
	    {
	      for (int q=ShiftedStateHighestBit+1; q<NewQ;++q)
		this->ShiftedState[q]=0;
	      ShiftedStateHighestBit=NewQ;
	    }
	  this->ShiftedState[NewQ]=this->TemporaryState[Q];
	  BosonsLeft-=TemporaryState[Q];
	}
      --Q;
    }
  
// verify sign of phase!
//cout << "TranslationPhase for shift by ("<<shiftX<<","<<shiftY<<")="<<Polar(1.0, 2.0*M_PI*FluxDensity*shiftX*CountYCoordinates)<<endl;
  
  translationPhase = Polar(1.0, 2.0*M_PI*FluxDensity*shiftX*CountYCoordinates)* Conj(CumulatedPhase);
  
//cout<<"Cumulated Phase="<<Arg(CumulatedPhase)/M_PI<<"pi"<<endl;
//cout << "NS:";
//for (int i=0; i<=ShiftedStateHighestBit; ++i)
//  cout << " "<<ShiftedState[i];
//cout << endl;

  return this->FindStateIndex(this->BosonToFermion(this->ShiftedState, this->ShiftedStateHighestBit), this->ShiftedStateHighestBit + this->NbrBosons - 1);
}

// find whether there is a translation vector from state i to state f
// i = index of initial state
// f = index of final state
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// return value = final state can be reached by translation

bool BosonOnLatticeKyLong::IsTranslation(int i, int f, int &shiftX, int &shiftY)
{
  this->FermionToBoson(this->StateDescription[i], this->StateHighestBit[i], this->TemporaryState, this->TemporaryStateHighestBit);
  this->FermionToBoson(this->StateDescription[f], this->StateHighestBit[f], this->ShiftedState, this->ShiftedStateHighestBit);
  int *count = new int[this->NbrBosons];
  int *count2 = new int[this->NbrBosons];
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      count[j]=0;
      count2[j]=0;
    }  
  for (int j=0; j<=this->TemporaryStateHighestBit; ++j) 
    if (TemporaryState[j]>0) 
      ++count[TemporaryState[j]-1];
  
  for (int j=0; j<=this->ShiftedStateHighestBit; ++j) 
    if (ShiftedState[j]>0) 
      ++count2[ShiftedState[j]-1];
  
  for (int j=0; j<this->NbrBosons; ++j)    
    if (count[j]!=count2[j]) 
      return false;
  
  delete [] count2;
  
  // more tests needed...
  
  delete [] count;
  return true;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnLatticeKyLong::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescription[state], this->StateHighestBit[state], this->TemporaryState, this->TemporaryStateHighestBit);
  int i = 0;
  for (; i <= this->TemporaryStateHighestBit; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i < this->NbrStates; ++i)
    Str << "0 ";
#ifdef DEBUG_OUTPUT
  bitset<32> b = this->StateDescription[state];
  Str << "   " << b << "   highestBit = " << this->StateHighestBit[state];
  cout << "   highestState = "<<this->TemporaryStateHighestBit;
#endif
  return Str;
}


// find state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int BosonOnLatticeKyLong::FindStateIndex(ULONGLONG stateDescription, int highestBit)
{
//   cout << "in FindStateIndex: desc=" << b<<", highest bit=" << highestBit << endl;
//   cout << "LookUpTable = "<<this->LookUpTableShift[highestBit]<<"   "; cout.flush();

  long PosMax = stateDescription >> this->LookUpTableShift[highestBit];  
  
  //   cout << "PosMax = "<<PosMax<<"   "; cout.flush();
//   cout << "highestBit = "<<highestBit<<endl; cout.flush();
  
  long PosMin = this->LookUpTable[highestBit][PosMax];
  PosMax = this->LookUpTable[highestBit][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  ULONGLONG CurrentState = this->StateDescription[PosMid];
  
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
    return PosMin;
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnLatticeKyLong::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrCodingBits);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrCodingBits)
    this->MaximumLookUpShift = this->NbrCodingBits;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;
  
  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrCodingBits];
  this->LookUpTableShift = new int [this->NbrCodingBits];
  for (int i = 0; i < this->NbrCodingBits; ++i)
    {
      this->LookUpTableShift[i]=0;
      this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
      for (unsigned int j=0; j<this->LookUpTableMemorySize + 1; ++j)
	this->LookUpTable[i][j]=0;
    }
  
  int CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
  ULONGLONG CurrentLookUpTableValue = this->LookUpTableMemorySize;
  ULONGLONG TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentHighestBit != this->StateHighestBit[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  	 
/*	  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
	    cout << TmpLookUpTable[j] << " ";
	    cout << endl << "-------------------------------------------" << endl;*/
 	  
	  CurrentHighestBit = this->StateHighestBit[i];
	  TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
	  if (CurrentHighestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentHighestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentHighestBit];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;
  
/*  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
      cout << TmpLookUpTable[j] << " ";
    cout << endl << "-------------------------------------------" << endl;*/
}

// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found

int BosonOnLatticeKyLong::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
{
  cout <<"Calling wrong function BosonOnLatticeKyLong::CarefulFindStateIndex"<<endl;
  exit(1);
 return this->HilbertSpaceDimension;
}

// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found

int BosonOnLatticeKyLong::CarefulFindStateIndex(ULONGLONG stateDescription, int highestBit)
{
  if (bitcount(stateDescription) != this->NbrBosons)
    {
      return this->HilbertSpaceDimension;
    }
  if (highestBit<0)
    {
      highestBit = getHighestBit(stateDescription);
    }
  if (highestBit >= this->NbrStates+this->NbrBosons - 1)
    {
      return this->HilbertSpaceDimension;
    }      
  return this->FindStateIndex(stateDescription, highestBit);
}


// conversion to generic (full) many-body representation in real-space basis
// state: many-body state in Ky-momentum basis
// nbodyBasis: full Hilbert-space in real-space representation
// returns: vector in many-body basis of targetSpace

ComplexVector& BosonOnLatticeKyLong::ConvertToNbodyBasis(ComplexVector& state, ParticleOnLattice &nbodyBasis)
{
  ComplexVector *TmpV=new ComplexVector();
  *TmpV=this->ConvertToNbodyBasis(state, nbodyBasis, 0, this->GetHilbertSpaceDimension());
  cout << "Norm was:" << TmpV->Norm() << endl;
  *TmpV/=TmpV->Norm();
  return *TmpV;
}


// conversion to generic (full) many-body representation in real-space basis
// state: many-body state in Ky-momentum basis
// nbodyBasis: full Hilbert-space in real-space representation
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// returns: vector in many-body basis of targetSpace

ComplexVector& BosonOnLatticeKyLong::ConvertToNbodyBasis(ComplexVector& state, ParticleOnLattice &nbodyBasis, int firstComponent, int nbrComponent)
{
  this->TargetVector.ResizeAndClean(nbodyBasis.GetHilbertSpaceDimension());
  this->FullSpace = ((BosonOnLattice*)&nbodyBasis);
  int *QuantumNumbers = new int[this->NbrBosons];
  double Normalization;
  int LastComponent=firstComponent+nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      if (Norm(state[i])>1e-15)
	{
	  this->ListQuantumNumbers(i,QuantumNumbers,Normalization);
	  this->ExpandBasisState(this->NbrBosons, QuantumNumbers, ((ULONGLONG) 0x0ul), Conj(state[i])/Normalization);
	}
    }
  return TargetVector;
}


// recursive expansion of an operator-product of creation operators in k
// in the full n-body basis
// nbrOperators = number of operators N remaining to be applied
// quantumNumbers = array of quantum numbers q1,...qN of creation operators
// state = state to be acted upon
// prefactor = previous coefficients applied to state
// 
// in last stage of recursion, writes to this->TargetVector using the Hilbert-Space this->FullSpace

void BosonOnLatticeKyLong::ExpandBasisState (int nbrOperators, int *quantumNumbers, ULONGLONG state, Complex prefactor)
{
  int Index;
  ULONGLONG ResultingState;
  int cK, K, N, S, Subl, TargetQ;
  double NewFactor = sqrt((double)1.0/(double)Kmax);
  double AdFactor;
  Complex TranslationPhase;
  this->DecodeQuantumNumber(quantumNumbers[nbrOperators-1], N, cK, Subl);
  
  //cout << "Decoded step 1 QN "<<quantumNumbers[nbrOperators-1]<<" -> N="<<N<<", cK="<<cK<<", Subl="<<Subl<<endl;

  this->DecodeCompositeMomentum(cK, K, S);
  
  //cout << "Decoded step 2 QN "<<quantumNumbers[nbrOperators-1]<<" -> N="<<N<<", K="<<K<<", S="<<S<<", Subl="<<Subl<<endl;
  double ExpFactor = -2.0*M_PI*(double)K/(double)Kmax;
  
  if (nbrOperators>1)
    {
      for (int r=0; r<this->Kmax; ++r)
	{
	  TargetQ=this->FullSpace->EncodeQuantumNumber(N, TranslationCell*r+S, Subl, TranslationPhase);
	  
//cout << "r="<<r<<": TargetQ="<<TargetQ<< " (x="<<N<<" y="<<TranslationCell*r+S<<"), applying onto state " << state << endl;	  

	  ResultingState = this->FullSpace->Ad(state,TargetQ, AdFactor);

//cout << "Recursing with state: "<<ResultingState<<endl;
	  
	  if (ResultingState != ((ULONGLONG) 0x0ul))
	    this->ExpandBasisState(nbrOperators-1, quantumNumbers, ResultingState, prefactor*AdFactor*NewFactor*Polar(1.0,ExpFactor*r));
	}
    }
  else
    {
      for (int r=0; r<this->Kmax; ++r)
	{
	  TargetQ=this->FullSpace->EncodeQuantumNumber(N, TranslationCell*r+S, Subl, TranslationPhase);
	  
//cout << "r="<<r<<": TargetQ="<<TargetQ<< " (x="<<N<<" y="<<TranslationCell*r+S<<"), applying onto state " << state << endl;*

	  ResultingState = this->FullSpace->Ad(state,TargetQ, AdFactor);
	  
//cout << "Finishing with state: "<<ResultingState<<endl;
	  
	  if (ResultingState!= ((ULONGLONG) 0x0ul) )
	    {	      
	      if ((Index=FullSpace->CarefulFindStateIndex(ResultingState,-1))<FullSpace->GetHilbertSpaceDimension())
		{
		  this->TargetVector[Index]+= prefactor*AdFactor*NewFactor*Polar(1.0,ExpFactor*r);

//cout << "Adding "<< prefactor*AdFactor*NewFactor*Polar(1.0,ExpFactor*r) <<" to component "<<Index<<endl;
		}
	    }
	}
    }
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// kySector = Ky sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnLatticeKyLong::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  
  if (nbrParticleSector == this->NbrBosons)
    {
      if (kySector == this->Ky)
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
  BosonOnLatticeKyLong SubsytemSpace (nbrParticleSector, this->Lx, this->Ly , kySector , this->NbrFluxQuanta, 10000000);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  
  int ComplementaryKySector = this->Ky - kySector;
  
  if (ComplementaryKySector < 0)
    ComplementaryKySector += this->Kmax;
  
  if (ComplementaryKySector >= this->Kmax)
    ComplementaryKySector -= this->Kmax;
  
  BosonOnLatticeKyLong ComplementarySpace (ComplementaryNbrParticles, this->Lx, this->Ly, ComplementaryKySector , this->NbrFluxQuanta, 10000000);
  
  FQHELatticeParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  
  Operation.ApplyOperation(architecture);
  
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnLatticeKyLong::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnLattice* complementaryHilbertSpace,  ParticleOnLattice* destinationHilbertSpace, ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
  BosonOnLatticeKyLong* TmpHilbertSpace =  (BosonOnLatticeKyLong*) complementaryHilbertSpace;
  BosonOnLatticeKyLong* TmpDestinationHilbertSpace =  (BosonOnLatticeKyLong*) destinationHilbertSpace;
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
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescription[i], TmpDestinationHilbertSpace->StateHighestBit[i], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateHighestBit);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateHighestBit; ++k)
	TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->StateDescription[minIndex], TmpHilbertSpace->StateHighestBit[minIndex], TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateHighestBit);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateHighestBit; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  for (int i = 0; i < this->NbrStates; i++)
	    this->TemporaryState[i] = 0;
	  
	  TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescription[j], TmpDestinationHilbertSpace->StateHighestBit[j], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateHighestBit);
	  
	  if (TmpHilbertSpace->TemporaryStateHighestBit >TmpDestinationHilbertSpace->TemporaryStateHighestBit)
	    this->TemporaryStateHighestBit = TmpHilbertSpace->TemporaryStateHighestBit;
	  else
	    this->TemporaryStateHighestBit = TmpDestinationHilbertSpace->TemporaryStateHighestBit;
	  
	  for (int i = 0; i <= TmpHilbertSpace->TemporaryStateHighestBit; i++)
	    this->TemporaryState[i] += TmpHilbertSpace->TemporaryState[i];
	  for (int i = 0; i <= TmpDestinationHilbertSpace->TemporaryStateHighestBit; i++)
	    this->TemporaryState[i] += TmpDestinationHilbertSpace->TemporaryState[i];
	  
	  int TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
	  
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateHighestBit; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
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
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		    ++TmpNbrNonZeroElements;
		  }
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  delete[] TmpDestinationLogFactorials;
  return TmpNbrNonZeroElements;
}
