////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//      NbrStates + NbrBosons - 1 < 63 or 31 (64 bit or 32bit systems)        //
//                                                                            //
//                        last modification : 10/02/2008                      //
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
#include "HilbertSpace/BosonOnLatticeGeneric.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <cmath>
#include <cstdlib>

#include <bitset>
using std::bitset;


using std::cout;
using std::endl;


// switch verbosity:
//#define DEBUG_OUTPUT



// default constructor
//

BosonOnLatticeGeneric::BosonOnLatticeGeneric ()
{
  this->HardCoreBasis=0;
}

// basic constructor -> yields a square lattice in Landau gauge
// 
// nbrBosons = number of bosons
// latticeGeometry = geometry of lattice system is living on
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// memory = memory that can be allocated for precalculations
// solenoidX = solenoid flux through lattice in x-direction (in units of pi)
// solenoidY = solenoid flux through lattice in y-direction (in units of pi)
BosonOnLatticeGeneric::BosonOnLatticeGeneric (int nbrBosons, LatticePhases *latticeGeometry, int nbrFluxQuanta, unsigned long memory, double solenoidX, double solenoidY)
{
  this->NbrBosons = nbrBosons;
  this->LatticeGeometry = latticeGeometry;
  this->TrueDimension = LatticeGeometry->GetLatticeDimension();
  this->NbrSublattices = LatticeGeometry->GetNbrSubLattices();
  this->Length = new int[TrueDimension];
  for (int i=0; i<TrueDimension; ++i)
    this->Length[i] = LatticeGeometry->GetLatticeLength(i);
  this->NbrStates = this->LatticeGeometry->GetNbrSites();
  
  TmpTranslations = new int[TrueDimension];
  TmpCoordinates = new int[TrueDimension];
  for (int i=0; i<TrueDimension; ++i)
    {
      TmpTranslations[i]=0;
      TmpCoordinates[i]=0;
    }
  
  this->SetNbrFluxQuanta(nbrFluxQuanta, solenoidX, solenoidY);

  this->HardCoreBasis = new HardCoreBoson(nbrBosons, NbrStates + nbrBosons - 1, memory);
  this->HilbertSpaceDimension = this->HardCoreBasis->GetHilbertSpaceDimension();

  this->TemporaryState = new unsigned long [this->NbrStates];
  this->ShiftedState = new unsigned long [this->NbrStates];
  this->ProdATemporaryState = new unsigned long [this->NbrStates];
  this->Flag.Initialize();

  this->CurrentTranslation = new int[2];
  this->CurrentMappings = new int[NbrStates]; 
  this->CurrentTranslationPhases = new Complex[NbrStates];
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

BosonOnLatticeGeneric::BosonOnLatticeGeneric(const BosonOnLatticeGeneric& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->LatticeGeometry = bosons.LatticeGeometry;
  this->Length = bosons.Length;
  this->NbrSublattices = bosons.NbrSublattices;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->FluxDensity = bosons.FluxDensity;
  this->TrueDimension = bosons.TrueDimension;
  this->SolenoidX = bosons.SolenoidX;
  this->SolenoidY = bosons.SolenoidY;
  this->LxTranslationPhase = bosons.LxTranslationPhase;
  this->LyTranslationPhase = bosons.LyTranslationPhase;
  this->NbrStates = bosons.NbrStates;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->HardCoreBasis = (HardCoreBoson*) bosons.HardCoreBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrStates];
  this->ShiftedState = new unsigned long [this->NbrStates];
  this->ProdATemporaryState = new unsigned long [this->NbrStates];
  this->TmpTranslations = new int[TrueDimension];
  this->TmpCoordinates = new int[TrueDimension];
  for (int i=0; i<TrueDimension; ++i)
    {
      TmpTranslations[i]=0;
      TmpCoordinates[i]=0;
    }
  this->CurrentTranslation = new int[2];
  this->CurrentTranslation[0]=bosons.CurrentTranslation[0];
  this->CurrentTranslation[1]=bosons.CurrentTranslation[1];
  this->CurrentMappings = new int[NbrStates]; 
  this->CurrentTranslationPhases = new Complex[NbrStates];
  for (int i=0; i<NbrStates; ++i)
    {
      this->CurrentMappings[i] = bosons.CurrentMappings[i];
      this->CurrentTranslationPhases[i] = bosons.CurrentTranslationPhases[i];
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

BosonOnLatticeGeneric::~BosonOnLatticeGeneric ()
{
  if (this->HardCoreBasis != 0)
    {
      delete this->HardCoreBasis;
      delete[] this->TemporaryState;
      delete[] this->ShiftedState;
      delete[] this->ProdATemporaryState;
      delete[] this->TmpTranslations;
      delete[] this->TmpCoordinates;
      delete[] this->CurrentTranslation;
      delete[] this->CurrentMappings;
      delete[] this->CurrentTranslationPhases;
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
      delete [] Length;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeGeneric& BosonOnLatticeGeneric::operator = (const BosonOnLatticeGeneric& bosons)
{
  if (this->HardCoreBasis != 0)
    {
      delete this->HardCoreBasis;
      delete[] this->TemporaryState;
      delete[] this->ShiftedState;
      delete[] this->ProdATemporaryState;
      delete[] this->TmpTranslations;
      delete[] this->TmpCoordinates;
      delete[] this->CurrentTranslation;
      delete[] this->CurrentMappings;
      delete[] this->CurrentTranslationPhases;
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
      delete [] Length;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->LatticeGeometry = bosons.LatticeGeometry;
  this->Length = bosons.Length;
  this->NbrSublattices = bosons.NbrSublattices;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->FluxDensity = bosons.FluxDensity;
  this->TrueDimension = bosons.TrueDimension;
  this->SolenoidX = bosons.SolenoidX;
  this->SolenoidY = bosons.SolenoidY;
  this->LxTranslationPhase = bosons.LxTranslationPhase;
  this->LyTranslationPhase = bosons.LyTranslationPhase;
  this->NbrStates = bosons.NbrStates;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->HardCoreBasis = (HardCoreBoson*) bosons.HardCoreBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrStates];
  this->ShiftedState = new unsigned long [this->NbrStates];
  this->ProdATemporaryState = new unsigned long [this->NbrStates];
  this->TmpTranslations = new int[TrueDimension];
  this->TmpCoordinates = new int[TrueDimension];
  for (int i=0; i<TrueDimension; ++i)
    {
      TmpTranslations[i]=0;
      TmpCoordinates[i]=0;
    }
  this->CurrentTranslation = new int[2];
  this->CurrentTranslation[0]=bosons.CurrentTranslation[0];
  this->CurrentTranslation[1]=bosons.CurrentTranslation[1];
  this->CurrentMappings = new int[NbrStates]; 
  this->CurrentTranslationPhases = new Complex[NbrStates];
  for (int i=0; i<NbrStates; ++i)
    {
      this->CurrentMappings[i] = bosons.CurrentMappings[i];
      this->CurrentTranslationPhases[i] = bosons.CurrentTranslationPhases[i];
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeGeneric::Clone()
{
  return new BosonOnLatticeGeneric(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnLatticeGeneric::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(this->NbrBosons);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnLatticeGeneric::GetQuantumNumber (int index)
{
  return new NumberParticleQuantumNumber(this->NbrBosons);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnLatticeGeneric::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// get the number of sites
//
// return value = number of sites
int BosonOnLatticeGeneric::GetNbrSites()
{
  return this->NbrStates;
}

// it is possible to change the flux through the simulation cell
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
void BosonOnLatticeGeneric::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->FluxDensity = ((double)NbrFluxQuanta)/this->NbrStates;  
  cout << "FluxDensity="<<this->FluxDensity<<endl;
  // xxx do something about translation phases!
}

// change flux through cell and periodic boundary conditions
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
void BosonOnLatticeGeneric::SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY)
{
  this->SolenoidX = M_PI*solenoidX;
  this->SolenoidY = M_PI*solenoidY;
  this->SetNbrFluxQuanta(nbrFluxQuanta);
}

// request solenoid fluxes
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
//
void BosonOnLatticeGeneric::GetSolenoidFluxes(double &solenoidX, double &solenoidY)
{
  solenoidX=this->SolenoidX;
  solenoidY=this->SolenoidY;
}

// obtain the current setting of the flux piercing this lattice
int BosonOnLatticeGeneric::GetNbrFluxQuanta()
{
  return this->NbrFluxQuanta;
}


/*void print_array2(int length, long unsigned int*array)
{
  if (length>0)
    {
      cout << array[0];
      for (int i=1; i<length; ++i) cout << " " << array[i];
      cout << " (length "<<length<<")"<<endl;
    }
}
*/


// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
unsigned long BosonOnLatticeGeneric::Ad (unsigned long state, int q, double& coefficient)
{
  int StateHighestBit = getHighestBit(state)-1;
  this->FermionToBoson(state, StateHighestBit, this->TemporaryState, this->TemporaryStateHighestBit);
  if (q>this->TemporaryStateHighestBit)
    {
      for (int i = this->TemporaryStateHighestBit + 1; i < q; ++i)
	this->TemporaryState[i] = 0x0ul;
      this->TemporaryState[q] = 0x1ul;
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

int BosonOnLatticeGeneric::AdAdAA (int index, int q1, int q2, int r1, int r2, double& coefficient)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
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
  return this->HardCoreBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
}

// apply a_r1 a_r2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// r1 = first index for annihilation operator
// r2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeGeneric::AA (int index, int r1, int r2)
{
#ifdef DEBUG_OUTPUT
  cout << "AA ("<<index <<" " <<r1 <<" " <<r2<<") ";
  this->PrintState(cout,index);
#endif
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->ProdATemporaryState, this->ProdATemporaryStateHighestBit);
  if ((r1 > this->ProdATemporaryStateHighestBit) || (r2 > this->ProdATemporaryStateHighestBit) || 
      (this->ProdATemporaryState[r1] == 0) || (this->ProdATemporaryState[r2] == 0) || ((r1 == r2) && (this->ProdATemporaryState[r1] == 1)))
    {
#ifdef DEBUG_OUTPUT
      if ((index==4)||(index==2))
	cout << "Vanishing on "<<index<<": "<<r1<<" "<<r2<<endl;
#endif
      return 0.0;
    }
  for (int i = this->ProdATemporaryStateHighestBit + 1; i < this->NbrStates; ++i)
    this->ProdATemporaryState[i] = 0ul;
  double Coefficient = this->ProdATemporaryState[r2];
  --this->ProdATemporaryState[r2];
  if ((r2==ProdATemporaryStateHighestBit)&&(ProdATemporaryStateHighestBit>0)&&(this->ProdATemporaryState[r2]==0))
    {
      --ProdATemporaryStateHighestBit;
      while ((ProdATemporaryStateHighestBit>0)&&(ProdATemporaryState[ProdATemporaryStateHighestBit] == 0))
	--ProdATemporaryStateHighestBit;
    }
  Coefficient *= this->ProdATemporaryState[r1];
  --this->ProdATemporaryState[r1];
  if ((r1==ProdATemporaryStateHighestBit)&&(ProdATemporaryStateHighestBit>0)&&(this->ProdATemporaryState[r1]==0))
    {
      --ProdATemporaryStateHighestBit;
      while ((ProdATemporaryStateHighestBit>0)&&(ProdATemporaryState[ProdATemporaryStateHighestBit] == 0))
	--ProdATemporaryStateHighestBit;
    }
#ifdef DEBUG_OUTPUT
  if ((index==4)||(index==2))
    cout << "Acting on "<<index<<": "<<r1<<" "<<r2<<endl;
  cout << "ProdATemporaryState="<<this->ProdATemporaryState[0];
  for (int i=1; i<NbrStates; ++i) cout << " " <<this->ProdATemporaryState[i];
  cout << " (h.b. "<<ProdATemporaryStateHighestBit<<")"<<endl;
#endif
  return sqrt(Coefficient);
}


// apply a^+_q1 a^+_q2 operator to the state produced using AA method (without destroying it)
//
// q1 = first index for creation operator
// q2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLatticeGeneric::AdAd (int q1, int q2, double& coefficient)
{
  for (int i = 0; i < this->NbrStates; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  TemporaryStateHighestBit=ProdATemporaryStateHighestBit;
  ++this->TemporaryState[q2];
  coefficient = this->TemporaryState[q2];
  if (q2 > TemporaryStateHighestBit)
    TemporaryStateHighestBit = q2;
  ++this->TemporaryState[q1];
  coefficient *= this->TemporaryState[q1];
  if (q1 > TemporaryStateHighestBit)
    TemporaryStateHighestBit = q1;
  coefficient = sqrt(coefficient);
#ifdef DEBUG_OUTPUT
  cout << "Temporary State="<<this->TemporaryState[0];
  for (int i=1; i<NbrStates; ++i) cout << " " <<this->TemporaryState[i];
  cout << " (h.b. "<<TemporaryStateHighestBit<<")"<<endl;
  unsigned long TargetState=this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit);
  int TargetIndex=this->HardCoreBasis->FindStateIndex(TargetState, this->TemporaryStateHighestBit + this->NbrBosons - 1);
  bitset<32> b=TargetState;
  bitset<32> b2=HardCoreBasis->StateDescription[TargetIndex];
  cout << q1 << " " << q2 <<" Binary Target: "<<b<<" index "<<TargetIndex<<" vs "<<b2<<endl;
  return TargetIndex;
#else
  return this->HardCoreBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
#endif
}


// apply a^+_q a_q operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_q a_q

double BosonOnLatticeGeneric::AdA (int index, int q)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  if (q > this->NbrStates)  
    return 0.0;
  return (double) (this->TemporaryState[q]);  
}

// apply a^+_q a_r operator to a given state 
//
// index = index of the state on which the operator has to be applied
// q = index of the creation operator
// r = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnLatticeGeneric::AdA (int index, int q, int r, double& coefficient)
{
  // cout << "AdA: " << index << ", " << q << ", " << r <<endl;
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  // cout << "Temporary state before: ";
//   for (int i=0; i<NbrStates; ++i)
//      cout << TemporaryState[i]<<" ";
//   cout << ", h.b.=" <<TemporaryStateHighestBit<< endl;
  // bitset<32>b2=this->HardCoreBasis->StateDescription[index];
  // cout << "converted: "<<b2 << endl;

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
  return this->HardCoreBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
}

// apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
// index = index of the state on which the operator has to be applied
// NbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
// QValues = array of quantum numbers where an interaction is present
// InteractionPerQ = coefficient U_q of the interaction
//
double BosonOnLatticeGeneric::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
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



// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// sublattice = sublattice index
// 
int BosonOnLatticeGeneric::EncodeQuantumNumber(int posx, int posy, int sublattice, Complex &translationPhase)
{
  TmpCoordinates[0]=posx;
  TmpCoordinates[1]=posy;
  int rst = LatticeGeometry->GetSiteNumber(TmpCoordinates, sublattice, TmpTranslations);
  double phase = LatticeGeometry->GetTunnellingPhaseFromGauge(rst, rst, TmpTranslations);
  translationPhase = Polar(1.0,2.0*M_PI*this->FluxDensity*phase);
  return rst;
}

// decode a single encoded quantum number q to the set of quantum numbers posx, posy
// posx = position along x-direction
// posy = position along y-direction
void BosonOnLatticeGeneric::DecodeQuantumNumber(int q, int &posx, int &posy, int &sublattice)
{
  LatticeGeometry->GetSiteCoordinates(q, TmpCoordinates, sublattice);
  posx=TmpCoordinates[0];
  posy=TmpCoordinates[1];
}

// check whether HilbertSpace implements ordering of operators
//
bool BosonOnLatticeGeneric::HaveOrder ()
{
  return true;
}


// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
int BosonOnLatticeGeneric::CheckOrder (int* m, int* n, int nbrIndices)
{
  if (nbrIndices>this->NbrBosons) return 0;
  this->TemporaryStateHighestBit=0;
  for (int i=0; i<this->NbrStates; ++i)
    this->TemporaryState[i]=0;
  for (int i=0; i<nbrIndices; ++i)
    {
      ++this->TemporaryState[m[i]];
      if (m[i]>TemporaryStateHighestBit)
	TemporaryStateHighestBit=m[i];
    }
  unsigned long CreationValue=this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit);
  this->TemporaryStateHighestBit=0;
  for (int i=0; i<this->NbrStates; ++i)
    this->TemporaryState[i]=0;
  for (int i=0; i<nbrIndices; ++i)
    {
      ++this->TemporaryState[n[i]];
      if (n[i]>TemporaryStateHighestBit)
	TemporaryStateHighestBit=n[i];
    }
  unsigned long AnnihilationValue=this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit);
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
void BosonOnLatticeGeneric::ListQuantumNumbers(int index, int *quantumNumbers, double &normalization)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  int NbrQ=0;
  normalization=1.0;
  for (int q=0; q<=TemporaryStateHighestBit; ++q)
    if (this->TemporaryState[q]>0)
      {		
	for (unsigned long l=this->TemporaryState[q];l>0; --l)
	  {
	    normalization*=l;
	    quantumNumbers[NbrQ]=q;
	    ++NbrQ;
	  }
      }
  normalization=sqrt(normalization);
}

// obtain a list of quantum numbers in state
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
void BosonOnLatticeGeneric::ListQuantumNumbers(int index, int *quantumNumbers)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  int NbrQ=0;
  for (int q=0; q<=TemporaryStateHighestBit; ++q)
    if (this->TemporaryState[q]!=0)
      {		
	for (unsigned int l=this->TemporaryState[q];l>0; --l)
	  quantumNumbers[NbrQ++]=q;
      }
}


// translate a state by a multiple of the lattice vectors
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// translationPhase = returns phase inccurred by translation
// return value = index of translated state
int BosonOnLatticeGeneric::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase,
					  double *translationX, double *translationY)
{
  cout << "Need to implement translation with non-standard boundary conditions for BosonOnLatticeGeneric::TranslateState!"<<endl;
  return this->TranslateState(index, shiftX, shiftY, translationPhase);
}

// translate a state by a multiple of the lattice vectors
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// translationPhase = returns phase inccurred by translation
// return value = index of translated state
int BosonOnLatticeGeneric::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase)
{
  if ((this->CurrentTranslation[0]!=shiftX)||(this->CurrentTranslation[1]!=shiftY))
    {
      double SolenoidFluxes[2];
      SolenoidFluxes[0]=this->SolenoidX;
      SolenoidFluxes[1]=this->SolenoidY;
      this->CurrentTranslation[0]=shiftX;
      this->CurrentTranslation[1]=shiftY;
      this->LatticeGeometry->GetTranslations(this->CurrentTranslation, this->CurrentMappings,
					     this->CurrentTranslationPhases, SolenoidFluxes);
    }
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  int BosonsLeft=this->NbrBosons;
  int Q=TemporaryStateHighestBit;
  int NewQ;
  Complex PeriodicPhase;
  Complex CumulatedPhase=1.0;
  //cout << "TS:";
  //for (int i=0; i<=TemporaryStateHighestBit; ++i)
  //  cout << " "<<TemporaryState[i];
  //cout << endl;
  this->ShiftedStateHighestBit=0;
  this->ShiftedState[0]=0;
  while ((Q>-1) && (BosonsLeft>0))
    {
      if (this->TemporaryState[Q]>0)
	{
	  NewQ=CurrentMappings[Q];
	  PeriodicPhase=CurrentTranslationPhases[Q];
	  for (unsigned i=0; i<this->TemporaryState[Q]; ++i) CumulatedPhase*=PeriodicPhase;
	  if (NewQ>ShiftedStateHighestBit)
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
  translationPhase = Conj(CumulatedPhase);
  return this->HardCoreBasis->FindStateIndex(this->BosonToFermion(this->ShiftedState, this->ShiftedStateHighestBit), this->ShiftedStateHighestBit + this->NbrBosons - 1);
}


// find whether there is a translation vector from state i to state f
// i = index of initial state
// f = index of final state
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// return value = final state can be reached by translation
bool BosonOnLatticeGeneric::IsTranslation(int i, int f, int &shiftX, int &shiftY)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[i], this->HardCoreBasis->StateHighestBit[i], this->TemporaryState, this->TemporaryStateHighestBit);
  this->FermionToBoson(this->HardCoreBasis->StateDescription[f], this->HardCoreBasis->StateHighestBit[f], this->ShiftedState, this->ShiftedStateHighestBit);
  int *count = new int[this->NbrBosons];
  int *count2 = new int[this->NbrBosons];
  for (int j=0; j<this->NbrBosons; ++j)
    {
      count[j]=0;
      count2[j]=0;
    }  
  for (int j=0; j<=this->TemporaryStateHighestBit; ++j) if (TemporaryState[j]>0) ++count[TemporaryState[j]-1];
  for (int j=0; j<=this->ShiftedStateHighestBit; ++j) if (ShiftedState[j]>0) ++count2[ShiftedState[j]-1];
  for (int j=0; j<this->NbrBosons; ++j)    
    if (count[j]!=count2[j]) return false;
  delete [] count2;
  
  // more tests needed...
  
  delete [] count;
  return true;
}

// apply a gauge transformation
// phases = phases in array ordered according to the quantum number q
// input = vector that has to be transformed according to that gauge
ComplexVector& BosonOnLatticeGeneric::GaugeTransformVector(double *phases, ComplexVector& input)
{
  double CumulatedPhase;
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->HardCoreBasis->StateDescription[i], this->HardCoreBasis->StateHighestBit[i], this->TemporaryState, this->TemporaryStateHighestBit);
      CumulatedPhase=0.0;
      int BosonsLeft=this->NbrBosons;
      int Q=TemporaryStateHighestBit;
      while ((Q>-1) && (BosonsLeft>0))
	{
	  if (this->TemporaryState[Q]>0)
	    CumulatedPhase+=this->TemporaryState[Q]*phases[Q];
	  --Q;
	}
      input[i]*=Polar(1.0, CumulatedPhase);
    }
  return input;
}



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnLatticeGeneric::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[state], this->HardCoreBasis->StateHighestBit[state], this->TemporaryState, this->TemporaryStateHighestBit);
  int i = 0;
  for (; i <= this->TemporaryStateHighestBit; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i < this->NbrStates; ++i)
    Str << "0 ";
#ifdef DEBUG_OUTPUT
  Str << "   highestBit = " << this->TemporaryStateHighestBit;
#endif
  return Str;
}




// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found
int BosonOnLatticeGeneric::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
{
  if (bitcount(stateDescription)!=this->NbrBosons)
    {
      return this->HilbertSpaceDimension;
    }
  if (highestBit<0)
    {
      highestBit = getHighestBit(stateDescription)-1;
    }
  if (highestBit >= this->NbrStates+this->NbrBosons - 1)
    {
      return this->HilbertSpaceDimension;
    }
  int Index = this->HardCoreBasis->FindStateIndex(stateDescription, highestBit);  
  if (this->HardCoreBasis->StateDescription[Index] == stateDescription)
    return Index;
  else
    {
      cout << "Unexpected situation in CarefulFindStateIndex!"<<endl;
      for (int i=0; i<HilbertSpaceDimension; ++i)
	if (this->HardCoreBasis->StateDescription[i] == stateDescription)
	  cout << "Element now found at i="<<i<<", "<<this->HardCoreBasis->StateDescription[i]
	       <<"="<<stateDescription<<"!"<<endl;
      
      return this->HilbertSpaceDimension;
    }
}
