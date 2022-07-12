////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Nicolas Regnault                  //
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
#include "HilbertSpace/BosonOnLattice.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/BinomialCoefficients.h"
#include "Architecture/ArchitectureOperation/FQHELatticeParticleEntanglementSpectrumOperation.h"
#include "MathTools/FactorialCoefficient.h" 
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

BosonOnLattice::BosonOnLattice ()
{
  this->HardCoreBasis=0;
}

// basic constructor -> yields a square lattice in Landau gauge
// 
// nbrBosons = number of bosons
// lx = length of simulation cell in x-direction
// ly = length of simulation cell in y-direction
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// memory = memory that can be allocated for precalculations
// solenoidX = solenoid flux through lattice in x-direction (in units of pi)
// solenoidY = solenoid flux through lattice in y-direction (in units of pi)
// landauGaugeAxis = direction of Landau-gauge
// nbrSublattices = number of sublattices to create
BosonOnLattice::BosonOnLattice (int nbrBosons, int lx, int ly, int nbrFluxQuanta, unsigned long memory, double solenoidX, double solenoidY, char landauGaugeAxis, int nbrSublattices)
{
  this->NbrBosons = nbrBosons;
  this->Lx = lx;
  this->Ly = ly;
  this->NbrSublattices = nbrSublattices;  
  this->NbrStates = Lx*Ly*nbrSublattices;
  this->LandauGaugeAxis=landauGaugeAxis;
  
  this->SetNbrFluxQuanta(nbrFluxQuanta, solenoidX, solenoidY);

  this->HardCoreBasis = new HardCoreBoson(nbrBosons, NbrStates + nbrBosons - 1, memory);
  this->HilbertSpaceDimension = this->HardCoreBasis->GetHilbertSpaceDimension();

  this->TemporaryState = new unsigned long [this->NbrStates];
  this->ShiftedState = new unsigned long [this->NbrStates];
  this->ProdATemporaryState = new unsigned long [this->NbrStates];
  this->Flag.Initialize();
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

BosonOnLattice::BosonOnLattice(const BosonOnLattice& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->Lx = bosons.Lx;
  this->Ly = bosons.Ly;
  this->NbrSublattices = bosons.NbrSublattices;
  this->NbrFluxQuanta = bosons.NbrFluxQuanta;
  this->FluxDensity = bosons.FluxDensity;
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

BosonOnLattice::~BosonOnLattice ()
{
  if (this->HardCoreBasis != 0)
    {
      delete this->HardCoreBasis;
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

BosonOnLattice& BosonOnLattice::operator = (const BosonOnLattice& bosons)
{
  if (this->HardCoreBasis != 0)
    {
      delete this->HardCoreBasis;
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLattice::Clone()
{
  return new BosonOnLattice(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnLattice::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(this->NbrBosons);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnLattice::GetQuantumNumber (int index)
{
  return new NumberParticleQuantumNumber(this->NbrBosons);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnLattice::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// get the number of sites
//
// return value = number of sites
int BosonOnLattice::GetNbrSites()
{
  return this->Lx * this->Ly;
}


// it is possible to change the flux through the simulation cell
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
void BosonOnLattice::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->FluxDensity = ((double)NbrFluxQuanta)/this->GetNbrSites();
  switch (this->LandauGaugeAxis)
    {
    case 'x':
      cout << "FluxDensity="<<this->FluxDensity<<endl;
      this->LxTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge...
      this->LyTranslationPhase = Polar(1.0, 2.0*M_PI*FluxDensity*this->Ly);
      cout << "LyTranslationPhase= exp(I*"<<2.0*M_PI*FluxDensity*this->Ly<<")="<<LyTranslationPhase<<endl;
      break;
    case 'y':
      cout << "FluxDensity="<<this->FluxDensity<<endl;
      this->LxTranslationPhase = Polar(1.0, -2.0*M_PI*FluxDensity*this->Lx);
      cout << "LxTranslationPhase= exp(I*"<<-2.0*M_PI*FluxDensity*this->Lx<<")="<<LxTranslationPhase<<endl;
      this->LyTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge...
      break;
    default:
      cout << "Unknown Quantization axis! Exiting BosonOnLattice..."<<endl;
      exit(1);
      break;
    }
}

// change flux through cell and periodic boundary conditions
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
void BosonOnLattice::SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY)
{
  this->SolenoidX = M_PI*solenoidX;
  this->SolenoidY = M_PI*solenoidY;
  this->SetNbrFluxQuanta(nbrFluxQuanta);
}

// request solenoid fluxes
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
//
void BosonOnLattice::GetSolenoidFluxes(double &solenoidX, double &solenoidY)
{
  solenoidX=this->SolenoidX;
  solenoidY=this->SolenoidY;
}

// obtain the current setting of the flux piercing this lattice
int BosonOnLattice::GetNbrFluxQuanta()
{
  return this->NbrFluxQuanta;
}


void print_array2(int length, long unsigned int*array)
{
  if (length>0)
    {
      cout << array[0];
      for (int i=1; i<length; ++i) cout << " " << array[i];
      cout << " (length "<<length<<")"<<endl;
    }
}


// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
unsigned long BosonOnLattice::Ad (unsigned long state, int q, double& coefficient)
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
int BosonOnLattice::AdAdAA (int index, int q1, int q2, int r1, int r2, double& coefficient)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  if ((r1 > this->TemporaryStateHighestBit) || (r2 > this->TemporaryStateHighestBit) || (this->TemporaryState[r1] == 0) || (this->TemporaryState[r2] == 0) || ((r1 == r2) && (this->TemporaryState[r1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = this->TemporaryState[r2];
  --this->TemporaryState[r2];
  if (this->TemporaryState[r2]==0 && TemporaryStateHighestBit==r2 && TemporaryStateHighestBit>0)
    {
      --TemporaryStateHighestBit;
      while (TemporaryState[TemporaryStateHighestBit] == 0 && TemporaryStateHighestBit>0)
        --TemporaryStateHighestBit;
    }
  coefficient *= this->TemporaryState[r1];
  --this->TemporaryState[r1];
  if (this->TemporaryState[r1]==0 && TemporaryStateHighestBit==r1 && TemporaryStateHighestBit>0)
    {
      --TemporaryStateHighestBit;
      while (TemporaryState[TemporaryStateHighestBit] == 0 && TemporaryStateHighestBit>0)
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

double BosonOnLattice::AA (int index, int r1, int r2)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->ProdATemporaryState, this->ProdATemporaryStateHighestBit);
  if ((r1 > this->ProdATemporaryStateHighestBit) || (r2 > this->ProdATemporaryStateHighestBit) || 
      (this->ProdATemporaryState[r1] == 0) || (this->ProdATemporaryState[r2] == 0) || ((r1 == r2) && (this->ProdATemporaryState[r1] == 1)))
    {
      return 0.0;
    }
  for (int i = this->ProdATemporaryStateHighestBit + 1; i < this->NbrStates; ++i)
    this->ProdATemporaryState[i] = 0ul;
  double Coefficient = this->ProdATemporaryState[r2];
  --this->ProdATemporaryState[r2];
  if (ProdATemporaryStateHighestBit == r2 && this->ProdATemporaryState[r2]==0 && ProdATemporaryStateHighestBit > 0)
    {
      --ProdATemporaryStateHighestBit;
      while (ProdATemporaryState[ProdATemporaryStateHighestBit] == 0 && ProdATemporaryStateHighestBit > 0)
	{
	  --ProdATemporaryStateHighestBit;
	  if (ProdATemporaryStateHighestBit<0)
	    {
	      cout << "Error in highest bit:"<<ProdATemporaryStateHighestBit<<endl;
	      exit(1);
	    }
	}
    }
  Coefficient *= this->ProdATemporaryState[r1];
  --this->ProdATemporaryState[r1];
  if (ProdATemporaryStateHighestBit == r1 && this->ProdATemporaryState[r1]==0 && ProdATemporaryStateHighestBit > 0)
    {
      --ProdATemporaryStateHighestBit;
      while (ProdATemporaryState[ProdATemporaryStateHighestBit] == 0 && ProdATemporaryStateHighestBit > 0)
	{
	  --ProdATemporaryStateHighestBit;
	  if (ProdATemporaryStateHighestBit<0)
	    {
	      cout << "Error in highest bit:"<<ProdATemporaryStateHighestBit<<endl;
	      exit(1);
	    }
	}
    }
  return sqrt(Coefficient);
}


// apply a^+_q1 a^+_q2 operator to the state produced using AA method (without destroying it)
//
// q1 = first index for creation operator
// q2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLattice::AdAd (int q1, int q2, double& coefficient)
{
  for (int i = 0; i < this->NbrStates; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[q2];
  TemporaryStateHighestBit=ProdATemporaryStateHighestBit;
  coefficient = this->TemporaryState[q2];
  if (q2 > TemporaryStateHighestBit)
    TemporaryStateHighestBit = q2;
  ++this->TemporaryState[q1];
  coefficient *= this->TemporaryState[q1];
  if (q1 > TemporaryStateHighestBit)
    TemporaryStateHighestBit = q1;
  coefficient = sqrt(coefficient);  
  return this->HardCoreBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
}


// apply a^+_q a_q operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_q a_q

double BosonOnLattice::AdA (int index, int q)
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
int BosonOnLattice::AdA (int index, int q, int r, double& coefficient)
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
double BosonOnLattice::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
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

// calculate non-local density-density interactions \sum q V_{q1,q2} : n_q1 n_q2 :
// index = index of the state on which the operator has to be applied
// nbrInteraction = number of q-values in sum
// interactionPerQ12 = coefficient V_(q1, q2) of the interaction, in order of q12Values
// q12Values = array of quantum numbers of the orbitals in tuples (q1, q2), assuming q1!=q2, 2*nbrInteraction entries in total
//
double BosonOnLattice::RhoRhoDiagonal(int index, int nbrInteraction, double *interactionPerQ12, int *q12Values)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  double result=0.0;
  for (int i=0; i<nbrInteraction; ++i)
    {
      int Q1=q12Values[i<<1];
      int Q2=q12Values[(i<<1)+1];
      result+=interactionPerQ12[i]*TemporaryState[Q1]*TemporaryState[Q2];
    }
//   cout << "state:";
//   for (int i=0; i<NbrStates; ++i)
//     cout << " " << TemporaryState[i];
//   cout << " (hb "<<TemporaryStateHighestBit<<") -> result = " << result<<endl;
    
  return result;
}


// apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
// index = index of the state on which the operator has to be applied
// NbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
// QValues = array of quantum numbers where an interaction is present
// InteractionPerQ = coefficient U_q of the interaction
//
/*
double BosonOnLattice::ProdAdProdADiagonal(int index, int nbrOperators, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  int BosonsLeft=this->NbrBosons;
  double result=0.0;
	double TmpResult = 0.0;
  if (nbrInteraction==NbrStates)
    {
      int Q=TemporaryStateHighestBit;
      while ((Q>-1) && (BosonsLeft>=nbrOperators))
	{
	  if (this->TemporaryState[Q]!= 0)
	    {
				TmpResult = interactionPerQ[Q];
				for (int j = 0 ; j< nbrOperators; j++)
					TmpResult *= (TemporaryState[Q]-j);
	      result+=TmpResult;
	      BosonsLeft-=TemporaryState[Q];
	    }
	  --Q;
	}
    }
  else // cannot assume ordered set of qValues
    {
      for (int i=0; (i<nbrInteraction)&&(BosonsLeft>=nbrOperators); ++i)
	{
	  int Q=qValues[i];
		TmpResult = interactionPerQ[Q];
		for (int j = 0 ; j< nbrOperators; j++)
			TmpResult *= (TemporaryState[Q]-j);
		result+=TmpResult;
		BosonsLeft-=TemporaryState[Q];	  
	}
		}
//   cout << "state:";
//   for (int i=0; i<NbrStates; ++i)
//     cout << " " << TemporaryState[i];
//   cout << " (hb "<<TemporaryStateHighestBit<<") -> result = " << result<<endl;
    
  return result;
}
*/
// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// sublattice = sublattice index
// 
int BosonOnLattice::EncodeQuantumNumber(int posx, int posy, int sublattice, Complex &translationPhase)
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
  while (posy<0)
    {
      posy+=this->Ly;
      ++numYTranslations;      
    }
  while (posy>=this->Ly)
    {
      posy-=this->Ly;
      --numYTranslations;
    }
  int rst = posx + this->Lx*posy;
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
  for (int y=1;y<=posy; ++y)
    translationPhase*=tmpPhase;
  translationPhase*=Polar(1.0, SolenoidX*numXTranslations);
  tmpPhase=1.0;
  if (numYTranslations>0)
    tmpPhase2=LyTranslationPhase;
  else
    tmpPhase2=Conj(LyTranslationPhase);
  for (int i=0; i<abs(numYTranslations); ++i)
    tmpPhase*=tmpPhase2;
  //cout<<" tmpPhaseY="<<tmpPhase;
  for (int x=1;x<=posx; ++x)
    translationPhase*=tmpPhase;
  translationPhase*=Polar(1.0, SolenoidY*numYTranslations);
  //cout << "tX="<<numXTranslations<< ", tY="<<numYTranslations<<", translationPhase= " <<translationPhase<<endl;
  return rst;
}

// decode a single encoded quantum number q to the set of quantum numbers posx, posy
// posx = position along x-direction
// posy = position along y-direction
void BosonOnLattice::DecodeQuantumNumber(int q, int &posx, int &posy, int &sublattice)
{
  int m=this->NbrSublattices;
  sublattice=q%m;
  posx=(q%(m*Lx))/m;
  m*=Lx;
  posy=(q%(m*Ly))/m;  
}

// check whether HilbertSpace implements ordering of operators
//
bool BosonOnLattice::HaveOrder ()
{
  return true;
}


// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
int BosonOnLattice::CheckOrder (int* m, int* n, int nbrIndices)
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
void BosonOnLattice::ListQuantumNumbers(int index, int *quantumNumbers, double &normalization)
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
void BosonOnLattice::ListQuantumNumbers(int index, int *quantumNumbers)
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
int BosonOnLattice::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase)
{
  this->FermionToBoson(this->HardCoreBasis->StateDescription[index], this->HardCoreBasis->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
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
  this->ShiftedStateHighestBit=0;
  this->ShiftedState[0]=0;
  while ((Q>-1) && (BosonsLeft>0))
    {
      if (this->TemporaryState[Q]!=0)
	{
	  this->DecodeQuantumNumber(Q,OldX, OldY, OldSl);
	  CountYCoordinates+=this->TemporaryState[Q]*OldY;
	  NewQ=this->EncodeQuantumNumber(OldX+shiftX, OldY+shiftY, OldSl, PeriodicPhase);
	  //cout << "PeriodicPhase for shift ("<<OldX<<","<<OldY<<")->("<<OldX+shiftX<<","<<OldY+shiftY<<")="<< PeriodicPhase<<endl;
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
  // verify sign of phase!
  //cout << "TranslationPhase for shift by ("<<shiftX<<","<<shiftY<<")="<<Polar(1.0, 2.0*M_PI*FluxDensity*shiftX*CountYCoordinates)<<endl;
  translationPhase = Polar(1.0, 2.0*M_PI*FluxDensity*shiftX*CountYCoordinates)* Conj(CumulatedPhase);
  //cout<<"Cumulated Phase="<<Arg(CumulatedPhase)/M_PI<<"pi"<<endl;
  //cout << "NS:";
  //for (int i=0; i<=ShiftedStateHighestBit; ++i)
  //  cout << " "<<ShiftedState[i];
  //cout << endl;
  return this->HardCoreBasis->FindStateIndex(this->BosonToFermion(this->ShiftedState, this->ShiftedStateHighestBit), this->ShiftedStateHighestBit + this->NbrBosons - 1);
}

// find whether there is a translation vector from state i to state f
// i = index of initial state
// f = index of final state
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// return value = final state can be reached by translation
bool BosonOnLattice::IsTranslation(int i, int f, int &shiftX, int &shiftY)
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


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnLattice::PrintState (ostream& Str, int state)
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
int BosonOnLattice::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
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

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnLattice::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  BosonOnLattice SubsytemSpace (nbrParticleSector, this->Lx, this->Ly, this->NbrFluxQuanta, 10000000, this->SolenoidX, this->SolenoidY);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnLattice ComplementarySpace (ComplementaryNbrParticles, this->Lx, this->Ly, this->NbrFluxQuanta, 10000000, this->SolenoidX, this->SolenoidY);
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

long BosonOnLattice::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnLattice* complementaryHilbertSpace,  ParticleOnLattice* destinationHilbertSpace,
									  ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
   BosonOnLattice* TmpHilbertSpace =  (BosonOnLattice*) complementaryHilbertSpace;
   BosonOnLattice* TmpDestinationHilbertSpace =  (BosonOnLattice*) destinationHilbertSpace;
   int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
   int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
   unsigned long* TmpMonomial2 = new unsigned long [NbrBosonSector];
   unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
   unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
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
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->HardCoreBasis->StateDescription[i], TmpDestinationHilbertSpace->HardCoreBasis->StateHighestBit[i], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateHighestBit);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateHighestBit; ++k)
	TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
   for (; minIndex < MaxIndex; ++minIndex)    
     {
      int Pos = 0;
      TmpHilbertSpace->ConvertToMonomial(TmpHilbertSpace->HardCoreBasis->StateDescription[minIndex], TmpHilbertSpace->HardCoreBasis->StateHighestBit[minIndex], TmpMonomial1);
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->HardCoreBasis->StateDescription[minIndex], TmpHilbertSpace->HardCoreBasis->StateHighestBit[minIndex], TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateHighestBit);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateHighestBit; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->HardCoreBasis->StateDescription[j], TmpDestinationHilbertSpace->HardCoreBasis->StateHighestBit[j], TmpMonomial2);
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < NbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < NbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < NbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->HardCoreBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateHighestBit);
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
   delete[] TmpMonomial2;
   delete[] TmpMonomial1;
   delete[] TmpMonomial3;
   delete[] TmpStatePosition;
   delete[] TmpStatePosition2;
   delete[] TmpStateCoefficient;
   delete[] TmpDestinationLogFactorials;
   return TmpNbrNonZeroElements;
}


// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrozed state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

void BosonOnLattice::SymmetrizeU1U1State (ComplexVector& symmetrizedVector, ComplexVector& leftVector, ComplexVector& rightVector, BosonOnLattice* leftSpace, BosonOnLattice* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
  unsigned long LastComponent = firstComponent + nbrComponents;
  
  FactorialCoefficient Factorial1;
  FactorialCoefficient Factorial2;
  if (unnormalizedBasisFlag == true)
    {
      for (long i = firstComponent; i < LastComponent; ++i)
	{
	  this->FermionToBoson(leftSpace->HardCoreBasis->StateDescription[i], leftSpace->HardCoreBasis->StateHighestBit[i], 
			       leftSpace->TemporaryState, leftSpace->TemporaryStateHighestBit);
	  for (int k = leftSpace->TemporaryStateHighestBit + 1;  k < leftSpace->NbrStates; ++k)
	    leftSpace->TemporaryState[k] = 0;
	  Complex TmpCoefficient = leftVector[i];
	  Factorial1.SetToOne();
	  Factorial1.Power2Divide(leftSpace->NbrBosons);
	  for (int k = 0; k <= leftSpace->TemporaryStateHighestBit; ++k)
	    if (leftSpace->TemporaryState[k] > 1)
	      Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
	  
	  for (long j = 0l; j < rightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      this->FermionToBoson(rightSpace->HardCoreBasis->StateDescription[j], rightSpace->HardCoreBasis->StateHighestBit[j], 
				   rightSpace->TemporaryState, rightSpace->TemporaryStateHighestBit);
	      int k = 0;
	      for (; k <= rightSpace->TemporaryStateHighestBit; ++k)
		this->TemporaryState[k] = leftSpace->TemporaryState[k] + rightSpace->TemporaryState[k];
	      this->TemporaryStateHighestBit = rightSpace->TemporaryStateHighestBit;
	      if (leftSpace->TemporaryStateHighestBit > rightSpace->TemporaryStateHighestBit)
		{
		  for (; k <= leftSpace->TemporaryStateHighestBit; ++k)
		    this->TemporaryState[k] = leftSpace->TemporaryState[k];
		  this->TemporaryStateHighestBit = leftSpace->TemporaryStateHighestBit;
		}
	      int TmpPos = this->HardCoreBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
	      if (TmpPos < this->HilbertSpaceDimension)
		{
		  Factorial2 = Factorial1;
		  for (k = 0; k <= rightSpace->TemporaryStateHighestBit; ++k)
		    if (rightSpace->TemporaryState[k] > 1)
		      Factorial2.FactorialDivide(rightSpace->TemporaryState[k]);
		  for (k = 0; k <= this->TemporaryStateHighestBit; ++k)
		    if (this->TemporaryState[k] > 1)
		      Factorial2.FactorialMultiply(this->TemporaryState[k]);	      
		  symmetrizedVector[TmpPos] += Factorial2.GetNumericalValue() * TmpCoefficient * rightVector[j];
		}
	    }
	}
    }
  else
    {
      for (long i = firstComponent; i < LastComponent; ++i)
	{
	  this->FermionToBoson(leftSpace->HardCoreBasis->StateDescription[i], leftSpace->HardCoreBasis->StateHighestBit[i], 
			       leftSpace->TemporaryState, leftSpace->TemporaryStateHighestBit);
	  for (int k = leftSpace->TemporaryStateHighestBit + 1;  k < leftSpace->NbrStates; ++k)
	    leftSpace->TemporaryState[k] = 0;
	  Complex TmpCoefficient = leftVector[i];
	  Factorial1.SetToOne();
	  Factorial1.Power2Divide(leftSpace->NbrBosons);
	  for (int k = 0; k <= leftSpace->TemporaryStateHighestBit; ++k)
	    if (leftSpace->TemporaryState[k] > 1)
	      Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
	  
	  for (long j = 0l; j < rightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      this->FermionToBoson(rightSpace->HardCoreBasis->StateDescription[j], rightSpace->HardCoreBasis->StateHighestBit[j], 
				   rightSpace->TemporaryState, rightSpace->TemporaryStateHighestBit);
	     
	      int k = 0;
	      for (; k <= rightSpace->TemporaryStateHighestBit; ++k)
		this->TemporaryState[k] = leftSpace->TemporaryState[k] + rightSpace->TemporaryState[k];
	      this->TemporaryStateHighestBit = rightSpace->TemporaryStateHighestBit;
	      if (leftSpace->TemporaryStateHighestBit > rightSpace->TemporaryStateHighestBit)
		{
		  for (; k <= leftSpace->TemporaryStateHighestBit; ++k)
		    this->TemporaryState[k] = leftSpace->TemporaryState[k];
		  this->TemporaryStateHighestBit = leftSpace->TemporaryStateHighestBit;
		}
	      int TmpPos = this->HardCoreBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
	      if (TmpPos < this->HilbertSpaceDimension)
		{
		  Factorial2 = Factorial1;
		  for (k = 0; k <= rightSpace->TemporaryStateHighestBit; ++k)
		    if (rightSpace->TemporaryState[k] > 1)
		      Factorial2.FactorialDivide(rightSpace->TemporaryState[k]);
		  for (k = 0; k <= this->TemporaryStateHighestBit; ++k)
		    if (this->TemporaryState[k] > 1)
		      Factorial2.FactorialMultiply(this->TemporaryState[k]);	      
		  symmetrizedVector[TmpPos] += sqrt(Factorial2.GetNumericalValue()) * TmpCoefficient * rightVector[j];
		}
	    }
	}     
    }
      if ( unnormalizedBasisFlag == false )
   symmetrizedVector /= symmetrizedVector.Norm();
}
