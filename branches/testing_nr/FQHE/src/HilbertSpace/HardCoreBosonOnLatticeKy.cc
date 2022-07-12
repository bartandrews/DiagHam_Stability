
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
#include "HilbertSpace/HardCoreBosonOnLatticeKy.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <math.h>

#include <bitset>
using std::bitset;


using std::cout;
using std::endl;


// default constructor
//

HardCoreBosonOnLatticeKy::HardCoreBosonOnLatticeKy ()
{
  this->HilbertSpaceDimension=0;
}

// basic constructor -> yields a square lattice in Landau gauge
// 
// nbrBosons = number of bosons
// lx = length of simulation cell in x-direction
// ly = length of simulation cell in y-direction
// ky = many-body momentum in y-direction
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// memory = memory that can be allocated for precalculations
HardCoreBosonOnLatticeKy::HardCoreBosonOnLatticeKy (int nbrBosons, int lx, int ly, int ky, int nbrFluxQuanta, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->Lx = lx;
  this->Ly = ly;
  this->NbrSublattices = 1;  
  this->NbrStates = Lx*Ly;


#ifdef __64_BITS__  
  if (NbrStates>64)    
#else
  if (NbrStates>32)    
#endif
    {
      cout<<"HardCoreBosonOnLatticeKy: Cannot represent the "<<NbrStates<<" coding bits requested in a single word"<<endl;
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

  this->StateDescription=new unsigned long[this->HilbertSpaceDimension];
  this->StateHighestBit=new int[this->HilbertSpaceDimension];
  
  this->Flag.Initialize();
  
  this->GenerateStates(nbrBosons, NbrStates-1, NbrStates-1, 0, 0, /* debug */ 0);
  
  this->TargetSpace=this;
  this->GenerateLookUpTable(memory);

  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    {
      PrintState(cout, i);
      cout << endl;
    }

  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;

}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

HardCoreBosonOnLatticeKy::HardCoreBosonOnLatticeKy(const HardCoreBosonOnLatticeKy& bosons)
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
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;

  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->TargetSpace = bosons.TargetSpace;
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
}

// destructor
//

HardCoreBosonOnLatticeKy::~HardCoreBosonOnLatticeKy ()
{
  if (this->HilbertSpaceDimension != 0) 
    {
      if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete[] this->StateDescription;
	  if (this->StateHighestBit != 0)
	    delete[] this->StateHighestBit;
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < this->NbrStates; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	}
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

HardCoreBosonOnLatticeKy& HardCoreBosonOnLatticeKy::operator = (const HardCoreBosonOnLatticeKy& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrStates; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
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
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->TargetSpace = bosons.TargetSpace;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* HardCoreBosonOnLatticeKy::Clone()
{
  return new HardCoreBosonOnLatticeKy(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> HardCoreBosonOnLatticeKy::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(this->NbrBosons);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* HardCoreBosonOnLatticeKy::GetQuantumNumber (int index)
{
  return new NumberParticleQuantumNumber(this->NbrBosons);
}

// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id
int HardCoreBosonOnLatticeKy::GetHilbertSpaceAdditionalSymmetry()
{
  return ParticleOnLattice::YTranslations;
}


// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* HardCoreBosonOnLatticeKy::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// it is possible to change the flux through the simulation cell
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
void HardCoreBosonOnLatticeKy::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  cout << "Attention, do not use SetNbrFluxQuanta for Bosons on Lattice with translation symmetry!"<<endl;
  exit(-1);
}

// get maximum possible momentum for this geometry
// return = maximum value of Ky
int HardCoreBosonOnLatticeKy::GetMaximumKy()
{
  return this->Kmax;
}

// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
unsigned long HardCoreBosonOnLatticeKy::Ad (unsigned long state, int q, double &coefficient)
{
  if ((state & (((unsigned long) (0x1)) << q))!= 0)
    {
      return 0x0l;
    }
  coefficient=1.0;
  state |= (((unsigned long) (0x1)) << q);
  return state;
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

int HardCoreBosonOnLatticeKy::AdAdAA (int index, int m1, int m2, int n1, int n2, double &coefficient)
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


// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double HardCoreBosonOnLatticeKy::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  if (((ProdATemporaryState & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((ProdATemporaryState & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2))
    return 0.0;

  this->ProdAHighestBit = this->StateHighestBit[index];

  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n2);
  
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n1);

  while ((this->ProdATemporaryState >> this->ProdAHighestBit) == 0)
    --this->ProdAHighestBit;
  return 1.0;
}


// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on a multiplicative, trivially one here (unreferenced)
// return value = index of the destination state 

int HardCoreBosonOnLatticeKy::AdAd (int m1, int m2, double& coefficient)
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
  return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
}


// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored (always 1.0)
// return value = index of the destination state 

int HardCoreBosonOnLatticeKy::AdA (int index, int m, int n, double &coefficient)
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
double HardCoreBosonOnLatticeKy::AdA (int index, int m)
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
double HardCoreBosonOnLatticeKy::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  return 0.0;
}



// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxQ = maximum value for the quantum number of a boson in the state
// currentMaxQ = current max value for the quantum number of bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

int HardCoreBosonOnLatticeKy::GenerateStates(int nbrBosons, int maxQ, int currentMaxQ, int pos, int currentMomentum, int debugLevel)
{
  // for (int I=0; I<debugLevel; ++I) cout << " "; cout << "Called: GenerateStates("<< nbrBosons <<" ," << maxQ <<", "<< currentMaxQ<<", "<< pos<<", "<< currentMomentum<<")"<<endl;
  bitset<32> b;
  if (nbrBosons == 0)
    {
      // create an entry for a state and set it to zero
      if ((currentMomentum % this->Kmax) == this->Ky)
	{	  
	  this->StateHighestBit[pos] = 0;
	  this->StateDescription[pos] = 0x0l;	  
// 	  for (int I=0; I<debugLevel; ++I) cout << " "; cout << "creating empty record "<<pos<<endl;
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
      if ((nbrBosons == 1)&&((currentMomentum % this->Kmax) == this->Ky))
	{	  
	  this->StateDescription[pos] = 0x1l;
	  this->StateHighestBit[pos] = 0;
// 	  b=this->StateDescription[pos];
// 	  for (int I=0; I<debugLevel; ++I) cout << " "; cout << "creating record "<<pos<<" = " << b << endl;
	  return pos + 1;
	}
      else
	{
	  return pos;
	}
    }

  int TmpNbrBosons = nbrBosons - 1;
  int ReducedCurrentMaxQ = currentMaxQ - 1;
  int TmpPos = pos;
  int currentQMomentum = this->DecodeKy(currentMaxQ);
  // assign a boson to current state and generate corresponding states after return from recursion

  TmpPos = this->GenerateStates(TmpNbrBosons, maxQ, ReducedCurrentMaxQ, pos, currentMomentum + currentQMomentum, debugLevel+1);      
  for (int i = pos; i < TmpPos; i++)
    {      
      this->StateDescription[i] |= 0x1l << currentMaxQ;
      this->StateHighestBit[i] = currentMaxQ;
//       b=this->StateDescription[i];
//       for (int I=0; I<debugLevel; ++I) cout << " "; cout << "updating record "<<i<<" = " << b << " ("<<StateHighestBit[i]<<")"<<endl;
    }
  pos = TmpPos;
  if (maxQ == currentMaxQ)
    return this->GenerateStates(nbrBosons, ReducedCurrentMaxQ, ReducedCurrentMaxQ, pos, currentMomentum, debugLevel+1);
  else
    {
      return this->GenerateStates(nbrBosons, maxQ, ReducedCurrentMaxQ, pos, currentMomentum, debugLevel+1);
    }
}



// recursively evaluate Hilbert space dimension 
//
// nbrBosons = number of bosons
// maxQ = maximum value for the quantum number of a boson in the state
// currentMaxQ = current max value for the quantum number of bosons that are still to be placed
// currentMomentum = current value of the momentum
// return value = Hilbert space dimension

int HardCoreBosonOnLatticeKy::EvaluateHilbertSpaceDimension(int nbrBosons, int maxQ, int currentMaxQ, int currentMomentum)
{
  //cout << "Calling: EvalDim ("<<nbrBosons<<", "<<maxQ<<", "<<currentMaxQ<<",...)"<<endl;
  if ((nbrBosons == 0) || (currentMaxQ == 0) )
    {
      if ((nbrBosons < 2) && ((currentMomentum % this->Kmax) == this->Ky))
	return 1;
      else
	return 0;
    }
    

  int ReducedCurrentMaxQ = currentMaxQ - 1;
  int currentQMomentum = this->DecodeKy(currentMaxQ);
  int rst = this->EvaluateHilbertSpaceDimension(nbrBosons - 1, maxQ, ReducedCurrentMaxQ, currentMomentum + currentQMomentum);

  if (maxQ == currentMaxQ)
    {
      //cout << "reduced MaxQ and CurrentMaxQ="<< ReducedCurrentMaxQ<<endl;
      return rst + this->EvaluateHilbertSpaceDimension(nbrBosons, ReducedCurrentMaxQ, ReducedCurrentMaxQ, currentMomentum);
    }
  else
    {
      //cout << "reduced CurrentMaxQ="<< ReducedCurrentMaxQ<<endl;
      return  rst + this->EvaluateHilbertSpaceDimension(nbrBosons, maxQ, ReducedCurrentMaxQ, currentMomentum);
    }
}




// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// ky = codes momentum in y-direction and 'subcell' index
// sublattice = sublattice index
// 
int HardCoreBosonOnLatticeKy::EncodeQuantumNumber(int posx, int ky, int sublattice, Complex &translationPhase)
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
  //cout<<" tmpPhaseX="<<tmpPhase;
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
void HardCoreBosonOnLatticeKy::DecodeQuantumNumber(int q, int &posx, int &ky, int &sublattice)
{
  int m=this->NbrSublattices;
  sublattice=q%m;
  posx=(q%(m*Lx))/m;
  m*=Lx;
  ky=(q%(m*Ly))/m;
}

// extract the momentum ky from a quantum number q
// return: momentum ky (in range 0...Kmax-1)
int HardCoreBosonOnLatticeKy::DecodeKy(int q)
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
int HardCoreBosonOnLatticeKy::EncodeCompositeMomentum(int ky, int fluxSubLattice)
{
  return ky*TranslationCell+fluxSubLattice;
}

// ky = true momentum in y-direction
// fluxSubLattice = 'sublattice' index remaining after translation symmetry
void HardCoreBosonOnLatticeKy::DecodeCompositeMomentum(int q, int &ky, int &fluxSubLattice)
{
  ky = q/TranslationCell;
  fluxSubLattice = q%TranslationCell;
}


// translate a state by a multiple of the lattice vectors
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// translationPhase = returns phase inccurred by translation
// return value = index of translated state
int HardCoreBosonOnLatticeKy::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase)
{
  int ShiftedStateHighestBit;
  unsigned long ShiftedState;
  int TemporaryStateHighestBit = this->StateHighestBit[index];
  unsigned long TemporaryState = this->StateDescription[index];
  int BosonsLeft=this->NbrBosons;
  int Q=TemporaryStateHighestBit;
  int OldX, OldY, OldSl;
  int NewQ;
  int CountYCoordinates=0; // total phase is shiftX * sum_i y_i in Landau gauge
  Complex PeriodicPhase;
  Complex CumulatedPhase=1.0;
  //cout << "TS:";
  //for (int i=0; i<=TemporaryStateHighestBit; ++i)
  //  cout << " "<<TemporaryState[i];
  //cout << endl;
  ShiftedStateHighestBit=0;
  ShiftedState=0x0l;
  while ((Q>-1) && (BosonsLeft>0))
    {
      if (TemporaryState&(0x1l<<Q))
	{
	  this->DecodeQuantumNumber(Q,OldX, OldY, OldSl);
	  CountYCoordinates+=OldY;
	  NewQ=this->EncodeQuantumNumber(OldX+shiftX, OldY+shiftY, OldSl, PeriodicPhase);
	  CumulatedPhase*=PeriodicPhase;
	  if (NewQ>ShiftedStateHighestBit)
	    {	      
	      ShiftedStateHighestBit=NewQ;
	    }
	  ShiftedState|=0x1ul<<NewQ;
	  --BosonsLeft;
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
  return this->FindStateIndex(ShiftedState, ShiftedStateHighestBit);
}

// find whether there is a translation vector from state i to state f
// i = index of initial state
// f = index of final state
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// return value = final state can be reached by translation
bool HardCoreBosonOnLatticeKy::IsTranslation(int i, int f, int &shiftX, int &shiftY)
{  
//   int TemporaryStateHighestBit = this->StateHighestBit[i];
//   unsigned long TemporaryState = this->StateDescription[i];
//   int ShiftedStateHighestBit = this->StateHighestBit[f];
//   unsigned long ShiftedState = this->StateDescription[f];
  
  // implementation needed...
  cout << "Implementation of HardCoreBosonOnLatticeKy::IsTranslation required"<<endl;
  
  return true;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& HardCoreBosonOnLatticeKy::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrStates; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
  Str << " position = " << this->FindStateIndex(TmpState, this->StateHighestBit[state]);
  if (state !=  this->FindStateIndex(TmpState, this->StateHighestBit[state]))
    Str << " error! ";
  return Str;
}


// find state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int HardCoreBosonOnLatticeKy::FindStateIndex(unsigned long stateDescription, int highestBit)
{
//   bitset <32> b = stateDescription;  
//   cout << "in FindStateIndex: desc=" << b<<", highest bit=" << highestBit << endl; cout.flush();
//   cout << "this="<<this<<endl;
//   cout << "LookUpTable = "<<this->LookUpTableShift<<"   "; cout.flush();
//   cout << "highestBit =" <<highestBit<< endl;
//   cout << "LookUpTable = "<<this->LookUpTableShift[highestBit]<<"   "; cout.flush();
  long PosMax = stateDescription >> this->LookUpTableShift[highestBit];  
//   cout << "PosMax = "<<PosMax<<"   "; cout.flush();
//   cout << "highestBit = "<<highestBit<<endl; cout.flush();
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


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void HardCoreBosonOnLatticeKy::GenerateLookUpTable(unsigned long memory)
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



// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found
int HardCoreBosonOnLatticeKy::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
{
  if (bitcount(stateDescription)!=this->NbrBosons)
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
