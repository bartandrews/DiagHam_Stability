////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//          class of hard-core boson with a single quantum number             //
//                                                                            //
//                        last modification : 11/02/2008                      //
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
#include "HilbertSpace/HardCoreBosonOnLatticeGeneric.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"

#include <bitset>
#include <iostream>
#include <cstdlib>

using std::bitset;
using std::cout;
using std::endl;

// default constructor
HardCoreBosonOnLatticeGeneric::HardCoreBosonOnLatticeGeneric()
{
  this->NbrBosons=0;
  this->HilbertSpaceDimension=0;
}


// basic constructor -> yields a square lattice in Landau gauge
// 
// nbrBosons = number of bosons
// latticeGeometry = geometry of lattice system is living on
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// memory = memory that can be allocated for precalculations
// solenoidX = solenoid flux through lattice in x-direction (in units of pi)
// solenoidY = solenoid flux through lattice in y-direction (in units of pi)
HardCoreBosonOnLatticeGeneric::HardCoreBosonOnLatticeGeneric (int nbrBosons, LatticePhases *latticeGeometry, int nbrFluxQuanta, unsigned long memory, double solenoidX, double solenoidY)
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

#ifdef __64_BITS__  
  if (this->NbrStates>64)    
#else
  if (this->NbrStates>32)    
#endif
    {
      cout<<"HardCoreBosonOnLatticeGeneric: Cannot represent the "<<NbrStates<<" states requested in a single word"<<endl;
      exit(1);
    }

  this->SetNbrFluxQuanta(nbrFluxQuanta, solenoidX, solenoidY);
  this->Flag.Initialize();

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(nbrBosons,NbrStates);
  
  this->StateDescription=new unsigned long[this->HilbertSpaceDimension];
  this->StateHighestBit=new int[this->HilbertSpaceDimension];
  this->GenerateStates(nbrBosons,NbrStates);
  this->TargetSpace=this;
  this->GenerateLookUpTable(memory);

  this->CurrentTranslation = new int[2];
  this->CurrentTranslation[0] = 0;
  this->CurrentTranslation[1] = 0;
  this->CurrentMappings = new int[NbrStates]; 
  this->CurrentTranslationPhases = new Complex[NbrStates];

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

//   for (int i=0; i<HilbertSpaceDimension; ++i)
//     {
//       this->PrintState(cout,i);
//       cout << endl;
//     }

  this->Flag.Initialize();
#ifdef __DEBUG__
  // for (int i=0; i<HilbertSpaceDimension; ++i)
//     {
//       this->PrintState(cout,i);
//       cout << endl;
//     }
  unsigned long UsedMemory = 0;  
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  // review here!
  //  UsedMemory += this->NbrLzValue * sizeof(int);
  // UsedMemory += this->NbrLzValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
  // UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension<<endl;
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
HardCoreBosonOnLatticeGeneric::HardCoreBosonOnLatticeGeneric(const HardCoreBosonOnLatticeGeneric& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->LatticeGeometry = bosons.LatticeGeometry;
  this->Length = bosons.Length;
  this->TrueDimension = bosons.TrueDimension;
  this->NbrSublattices = bosons.NbrSublattices;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->FluxDensity = bosons.FluxDensity;  
  this->LxTranslationPhase = bosons.LxTranslationPhase;
  this->LyTranslationPhase = bosons.LyTranslationPhase;  
  this->NbrStates = bosons.NbrStates;
  this->SolenoidX = bosons.SolenoidX;
  this->SolenoidY = bosons.SolenoidY;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->Flag = bosons.Flag;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
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


// virtual destructor
//

HardCoreBosonOnLatticeGeneric::~HardCoreBosonOnLatticeGeneric ()
{
  delete[] this->TmpTranslations;
  delete[] this->TmpCoordinates;
  delete[] this->CurrentTranslation;
  delete[] this->CurrentMappings;
  delete[] this->CurrentTranslationPhases;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrStates; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
      delete[] Length;
    }
}

// assignment (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space
HardCoreBosonOnLatticeGeneric& HardCoreBosonOnLatticeGeneric::operator = (const HardCoreBosonOnLatticeGeneric& bosons)
{
  delete[] this->TmpTranslations;
  delete[] this->TmpCoordinates;
  delete[] this->CurrentTranslation;
  delete[] this->CurrentMappings;
  delete[] this->CurrentTranslationPhases;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrStates; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
      delete[] Length;
    }
  if (this->TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->LatticeGeometry = bosons.LatticeGeometry;
  this->Length = bosons.Length;
  this->TrueDimension = bosons.TrueDimension;
  this->NbrSublattices = bosons.NbrSublattices;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->FluxDensity = bosons.FluxDensity;
  this->LxTranslationPhase = bosons.LxTranslationPhase;
  this->LyTranslationPhase = bosons.LyTranslationPhase;  
  this->NbrStates = bosons.NbrStates;
  this->SolenoidX = bosons.SolenoidX;
  this->SolenoidY = bosons.SolenoidY;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->Flag = bosons.Flag;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
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
AbstractHilbertSpace* HardCoreBosonOnLatticeGeneric::Clone()
{
  return new HardCoreBosonOnLatticeGeneric(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void HardCoreBosonOnLatticeGeneric::SetTargetSpace(ParticleOnLattice* targetSpace)
{
  this->TargetSpace=(HardCoreBosonOnLatticeGeneric*)targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int HardCoreBosonOnLatticeGeneric::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->GetHilbertSpaceDimension();
}

// get the particle statistic 
//
// return value = particle statistic
int HardCoreBosonOnLatticeGeneric::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}

// get the quantization axis 
//
// return value = particle statistic
char HardCoreBosonOnLatticeGeneric::GetLandauGaugeAxis()
{
  return 'y';
}


// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int HardCoreBosonOnLatticeGeneric::GetHilbertSpaceAdditionalSymmetry()
{
  return HardCoreBosonOnLatticeGeneric::NoSymmetry;
}

// check whether HilbertSpace implements ordering of operators
//
bool HardCoreBosonOnLatticeGeneric::HaveOrder ()
{
  return true;
}


// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
int HardCoreBosonOnLatticeGeneric::CheckOrder (int* m, int* n, int nbrIndices)
{
  unsigned long CreationValue=0x0ul;
  unsigned long AnnihilationValue=0x0ul;
  for (int i=0; i<nbrIndices; ++i)
    {
      CreationValue |= 0x1ul << m[i];
      AnnihilationValue |= 0x1ul << n[i];
    }
  if (CreationValue > AnnihilationValue)
    return 1;
  else if (CreationValue < AnnihilationValue)
    return -1;
  else return 0;
}


// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number
List<AbstractQuantumNumber*> HardCoreBosonOnLatticeGeneric::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(this->NbrBosons);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number
AbstractQuantumNumber* HardCoreBosonOnLatticeGeneric::GetQuantumNumber (int index)
{
  return new NumberParticleQuantumNumber(this->NbrBosons);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace
AbstractHilbertSpace* HardCoreBosonOnLatticeGeneric::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// get the number of sites
//
// return value = number of sites
int HardCoreBosonOnLatticeGeneric::GetNbrSites()
{
  return this->NbrStates;
}

// it is possible to change the flux through the simulation cell
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
void HardCoreBosonOnLatticeGeneric::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->FluxDensity = ((double)NbrFluxQuanta)/this->NbrStates;
  // cout << "FluxDensity="<<this->FluxDensity<<endl;
  // xxx do something about translation phases!
  // this->LxTranslationPhase = Polar(1.0, -2.0*M_PI*FluxDensity*this->Lx);
  // cout << "LxTranslationPhase= exp(I*"<<2.0*M_PI*FluxDensity*this->Lx<<")="<<LxTranslationPhase<<endl;
  // this->LyTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge ...
}

// change flux through cell and periodic boundary conditions
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
void HardCoreBosonOnLatticeGeneric::SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY)
{
  this->SolenoidX = M_PI*solenoidX;
  this->SolenoidY = M_PI*solenoidY;
  this->SetNbrFluxQuanta(nbrFluxQuanta);
}


// request solenoid fluxes
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
//
void HardCoreBosonOnLatticeGeneric::GetSolenoidFluxes(double &solenoidX, double &solenoidY)
{
  solenoidX=this->SolenoidX;
  solenoidY=this->SolenoidY;
}


// obtain the current setting of the flux piercing this lattice
int HardCoreBosonOnLatticeGeneric::GetNbrFluxQuanta()
{
  return this->NbrFluxQuanta;
}



// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
unsigned long HardCoreBosonOnLatticeGeneric::Ad (unsigned long state, int q, double& coefficient)
{  
  if ((state & (((unsigned long) (0x1)) << q))!= 0)
    {
      coefficient=0.0;
      return 0x0l;
    }
  coefficient=1.0;
  state |= (((unsigned long) (0x1)) << q);
  return state;
}

// apply annihilation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
// coefficient = reference on the double where the multiplicative factor has to be stored
unsigned long HardCoreBosonOnLatticeGeneric::A (unsigned long state, int q, double &coefficient)
{
  if (state & (((unsigned long) (0x1)) << q))
    {
      coefficient=1.0;
      state |= ~(((unsigned long) (0x1)) << q);
      return state;
    }
  coefficient=0.0;
  return 0x0l;
}


// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int HardCoreBosonOnLatticeGeneric::AdAdAA (int index, int m1, int m2, int n1, int n2, double &coefficient)
{
  coefficient = 1.0;
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((State & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  // perform annihilation operators
  TmpState &= ~(((unsigned long) (0x1)) << n2);
  if (NewHighestBit == n2)
    while ((TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  TmpState &= ~(((unsigned long) (0x1)) << n1);
  if (NewHighestBit == n1)
    while ((TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  // create particle at m2
  if ((TmpState & (((unsigned long) (0x1)) << m2))!= 0)
    {
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
  TmpState |= (((unsigned long) (0x1)) << m2);
  // create particle at m1
  if ((TmpState & (((unsigned long) (0x1)) << m1))!= 0)
    {
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
  TmpState |= (((unsigned long) (0x1)) << m1);
  return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
}


// apply a_n1 operator to a given state and search in target space
//
// index = index of the state on which the operator has to be applied
// q = index for annihilation operator
// coefficient = prefactor
// return value =  index in target space

int HardCoreBosonOnLatticeGeneric::A (int index, int q, double &coefficient)
{
  this->ProdATemporaryState = this->StateDescription[index];

  if (((ProdATemporaryState & (((unsigned long) (0x1)) << q)) == 0) )
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }

  this->ProdAHighestBit = this->StateHighestBit[index];

  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << q);

  while ((this->ProdAHighestBit>0)&&(this->ProdATemporaryState >> this->ProdAHighestBit) == 0)
    --this->ProdAHighestBit;
  coefficient = 1.0;
  return this->TargetSpace->FindStateIndex(ProdATemporaryState, ProdAHighestBit);
}


// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on a multiplicative, trivially one here (unreferenced)
// return value = index of the destination state 

int HardCoreBosonOnLatticeGeneric::Ad (int index, int q, double& coefficient)
{
  this->ProdATemporaryState = this->StateDescription[index];
  
  if (((ProdATemporaryState & (((unsigned long) (0x1)) << q)) != 0) )
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  
  this->ProdAHighestBit = this->StateHighestBit[index];
  
  this->ProdATemporaryState |= (((unsigned long) (0x1)) << q);

  while ((this->ProdAHighestBit>0)&&(this->ProdATemporaryState >> this->ProdAHighestBit) == 0)
    --this->ProdAHighestBit;
   coefficient = 1.0;
  return this->TargetSpace->FindStateIndex(ProdATemporaryState, ProdAHighestBit);
}


// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double HardCoreBosonOnLatticeGeneric::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  if (((ProdATemporaryState & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((ProdATemporaryState & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2))
    return 0.0;

  this->ProdAHighestBit = this->StateHighestBit[index];

  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n2);
  
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n1);

  while ((this->ProdAHighestBit>0)&&(this->ProdATemporaryState >> this->ProdAHighestBit) == 0)
    --this->ProdAHighestBit;
  return 1.0;
}


// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on a multiplicative, trivially one here (unreferenced)
// return value = index of the destination state 

int HardCoreBosonOnLatticeGeneric::AdAd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if ((TmpState & (((unsigned long) (0x1)) << m2))!= 0)
    {
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = this->ProdAHighestBit;
  if (m2 > NewHighestBit)
    NewHighestBit = m2;
  TmpState |= (((unsigned long) (0x1)) << m2);
  if ((TmpState & (((unsigned long) (0x1)) << m1))!= 0)
    {
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    NewHighestBit = m1;
  TmpState |= (((unsigned long) (0x1)) << m1);
  coefficient=1.0;
  return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
}


// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored (always 1.0)
// return value = index of the destination state 

int HardCoreBosonOnLatticeGeneric::AdA (int index, int m, int n, double &coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];  
  if (m!=n)
    {
      if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0))
	{
	  return this->HilbertSpaceDimension;
	}
      int NewHighestBit = StateHighestBit;
      unsigned long TmpState = State;
      // perform annihilation operators
      TmpState &= ~(((unsigned long) (0x1)) << n);
      if (NewHighestBit == n)
	while ((NewHighestBit > 0)&&((TmpState >> NewHighestBit) == 0))
	  --NewHighestBit;
      // create particle at m
      if ((TmpState & (((unsigned long) (0x1)) << m))!= 0)
	{
	  return this->HilbertSpaceDimension;
	}
      if (m > NewHighestBit)
	{
	  NewHighestBit = m;
	}
      TmpState |= (((unsigned long) (0x1)) << m);
      coefficient = 1.0;
      return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
    }
  else
    {
      if ((m > StateHighestBit) || ((State & (((unsigned long) (0x1)) << m)) == 0))
	return this->HilbertSpaceDimension;
      else
	{
	  coefficient = 1.0;
	  return index;
	}
    }
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double HardCoreBosonOnLatticeGeneric::AdA (int index, int m)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  if ((m > StateHighestBit) || ((State & (((unsigned long) (0x1)) << m)) == 0))
    return 0.0;
  else
    return 1.0;
}

// apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
// index = index of the state on which the operator has to be applied
// nbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
// qValues = array of quantum numbers where an interaction is present
// interactionPerQ = coefficient U_q of the interaction
//
// no diagonal interaction present, derelict function from inheritance
//
double HardCoreBosonOnLatticeGeneric::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  return 0.0;
}

// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// sublattice = sublattice index
// 
int HardCoreBosonOnLatticeGeneric::EncodeQuantumNumber(int posx, int posy, int sublattice, Complex &translationPhase)
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
void HardCoreBosonOnLatticeGeneric::DecodeQuantumNumber(int q, int &posx, int &posy, int &sublattice)
{
  LatticeGeometry->GetSiteCoordinates(q, TmpCoordinates, sublattice);
  posx=TmpCoordinates[0];
  posy=TmpCoordinates[1];
}

// obtain a list of quantum numbers in state
// index = index of many-body state to be considered
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
// normalization = indicating the multiplicity of the state for bosonic spaces
void HardCoreBosonOnLatticeGeneric::ListQuantumNumbers(int index, int *quantumNumbers, double &normalization)
{
  normalization=1.0;
  this->ListQuantumNumbers(index, quantumNumbers);
}

// obtain a list of quantum numbers in state
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
void HardCoreBosonOnLatticeGeneric::ListQuantumNumbers(int index, int *quantumNumbers)
{
  unsigned long State = this->StateDescription[index];
  int HighestBit = this->StateHighestBit[index];
  int NbrQ=0;
  for (int q=0; q<=HighestBit; ++q)
    if (State&(0x1ul<<q))
      quantumNumbers[NbrQ++]=q;
}

// translate a state by a multiple of the lattice vectors
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// translationPhase = returns phase inccurred by translation
// return value = index of translated state
int HardCoreBosonOnLatticeGeneric::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase)
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

  int ShiftedStateHighestBit=0;
  unsigned long ShiftedState=0x0ul;
  int TemporaryStateHighestBit = this->StateHighestBit[index];
  unsigned long TemporaryState = this->StateDescription[index];
  int BosonsLeft=this->NbrBosons;
  int Q=TemporaryStateHighestBit;
  int NewQ;
  Complex PeriodicPhase;
  Complex CumulatedPhase=1.0;
  while ((Q>-1) && (BosonsLeft>0))
    {
      if (TemporaryState&(0x1ul<<Q))
	{
	  NewQ=CurrentMappings[Q];
	  PeriodicPhase=CurrentTranslationPhases[Q];
	  CumulatedPhase*=PeriodicPhase;
	  if (NewQ>ShiftedStateHighestBit)
	    ShiftedStateHighestBit=NewQ;
	  ShiftedState|=0x1ul<<NewQ;
	  --BosonsLeft;
	}
      --Q;
    }
  translationPhase = Conj(CumulatedPhase);
  return this->FindStateIndex(ShiftedState, ShiftedStateHighestBit);
}


// find whether there is a translation vector from state i to state f
// i = index of initial state
// f = index of final state
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// return value = final state can be reached by translation
bool HardCoreBosonOnLatticeGeneric::IsTranslation(int i, int f, int &shiftX, int &shiftY)
{  
//   int TemporaryStateHighestBit = this->StateHighestBit[i];
//   unsigned long TemporaryState = this->StateDescription[i];
//   int ShiftedStateHighestBit = this->StateHighestBit[f];
//   unsigned long ShiftedState = this->StateDescription[f];
  
  // implementation needed...
  cout << "Implementation of HardCoreBosonOnLatticeGeneric::IsTranslation required"<<endl;
  
  return true;
}


  // evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex HardCoreBosonOnLatticeGeneric::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex HardCoreBosonOnLatticeGeneric::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, 
						     this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex HardCoreBosonOnLatticeGeneric::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
						int firstComponent, int nbrComponent)
{
  return Complex(0.0, 0.0);
}


// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex HardCoreBosonOnLatticeGeneric::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, 
								 int nextCoordinates, int firstComponent, 
								 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void HardCoreBosonOnLatticeGeneric::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& HardCoreBosonOnLatticeGeneric::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrStates; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
  Str << " position = " << this->FindStateIndex(TmpState, this->StateHighestBit[state]) ;// << ", hb="<<this->StateHighestBit[state];
  if (state !=  this->FindStateIndex(TmpState, this->StateHighestBit[state]))
    Str << " error! ";
  return Str;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int HardCoreBosonOnLatticeGeneric::FindStateIndex(unsigned long stateDescription, int highestBit)
{
  // bitset <32> b = stateDescription;
//   cout << "in FindStateIndex: desc=" << b<<", highest bit=" << highestBit << endl;
  long PosMax = stateDescription >> this->LookUpTableShift[highestBit];
  long PosMin = this->LookUpTable[highestBit][PosMax];
  PosMax = this->LookUpTable[highestBit][PosMax + 1];
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
    return PosMin;
}

// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found
int HardCoreBosonOnLatticeGeneric::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
{
//   bitset<20> b=stateDescription;
//   cout << "Searching " << b <<", bitcount="<<bitcount(stateDescription);
  if (bitcount(stateDescription)!=this->NbrBosons)
    {
      cout << "Wrong particle number - state: "<<bitcount(stateDescription)<<" space: "<<this->NbrBosons<<endl;
      return this->HilbertSpaceDimension;
    }
  if (highestBit<0)
    {
      cout << "Negative highest bit"<<endl;
      highestBit = getHighestBit(stateDescription)-1;
    }
//   cout << ", highestBit="<<highestBit;
  if (highestBit >= this->NbrStates)
    {
      cout << "Highest bit out of range"<<endl;
      return this->HilbertSpaceDimension;
    }
//   cout <<", Found index="<<this->FindStateIndex(stateDescription, highestBit)<<endl;
//   this->PrintState(cout, this->FindStateIndex(stateDescription, highestBit));
  return this->FindStateIndex(stateDescription, highestBit);
}



// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// nbrStates = number of elementary states available
//
int HardCoreBosonOnLatticeGeneric::EvaluateHilbertSpaceDimension(int nbrBosons,int nbrStates)
{  
  FactorialCoefficient F;
  F.SetToOne();
  F.FactorialMultiply(nbrStates);
  F.FactorialDivide(nbrStates-nbrBosons);
  F.FactorialDivide(nbrBosons);
  return (int)F.GetIntegerValue();
}

// generate many-body states without any additional symmetries
// nbrBosons = number of bosons
// nbrStates = number of elementary states available
//
// (using routines from UnsignedIntegerTools)
//
int HardCoreBosonOnLatticeGeneric::GenerateStates(int nbrBosons,int nbrStates)
{
  int countdown=this->HilbertSpaceDimension-1;
  int presentHighestBit=nbrBosons-1;    
  this->StateHighestBit[countdown] = presentHighestBit;
  this->StateDescription[countdown--]=smallestOne(nbrBosons);
  while (countdown>-1)
    {      
      this->StateDescription[countdown]=nextone(this->StateDescription[countdown+1]);
      if (this->StateDescription[countdown] & (0x1ul << (presentHighestBit+1)))
	++presentHighestBit;
      this->StateHighestBit[countdown]=presentHighestBit;
      --countdown;
    }
  if (this->StateDescription[0]!=biggestOne(nbrBosons, nbrStates))
    {
      bitset <64> b=this->StateDescription[0];
      cout << "Error in GenerateStates: Last state is: " << b << endl;
      exit(1);
    }
  return 0;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void HardCoreBosonOnLatticeGeneric::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrStates);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrStates)
    this->MaximumLookUpShift = this->NbrStates;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrStates];
  this->LookUpTableShift = new int [this->NbrStates];
  for (int i = 0; i < this->NbrStates; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
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

