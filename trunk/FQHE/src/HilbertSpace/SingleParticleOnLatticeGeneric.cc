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
#include "HilbertSpace/SingleParticleOnLatticeGeneric.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"

#include <bitset>
#include <iostream>
#include <cstdlib>

using std::bitset;
using std::cout;
using std::endl;

// switch for debugging output
//#define DEBUG_OUTPUT

// default constructor
SingleParticleOnLatticeGeneric::SingleParticleOnLatticeGeneric()
{
  this->NbrBosons=1;
  this->HilbertSpaceDimension=0;
}


// basic constructor -> yields a square lattice in Landau gauge
// 
// latticeGeometry = geometry of lattice system is living on
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// solenoidX = solenoid flux through lattice in x-direction (in units of pi)
// solenoidY = solenoid flux through lattice in y-direction (in units of pi)
SingleParticleOnLatticeGeneric::SingleParticleOnLatticeGeneric (LatticePhases *latticeGeometry, int nbrFluxQuanta, double solenoidX, double solenoidY)
{
  this->NbrBosons = 1;
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
  this->Flag.Initialize();

  this->HilbertSpaceDimension = NbrStates;
  
  this->StateDescription=new int[this->HilbertSpaceDimension];
  for (int i=0; i<NbrStates; ++i)
    this->StateDescription[i]=NbrStates-1-i;
  this->TargetSpace=this;

  this->Flag.Initialize();
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef DEBUG_OUTPUT
  for (int i=0; i<HilbertSpaceDimension; ++i)
    {
      this->PrintState(cout,i);
      cout << endl;
    }
#endif
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension<<endl;
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
SingleParticleOnLatticeGeneric::SingleParticleOnLatticeGeneric(const SingleParticleOnLatticeGeneric& bosons)
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
  this->Flag = bosons.Flag;
  this->TmpTranslations = new int[TrueDimension];
  this->TmpCoordinates = new int[TrueDimension];
  for (int i=0; i<TrueDimension; ++i)
    {
      TmpTranslations[i]=0;
      TmpCoordinates[i]=0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}


// virtual destructor
//

SingleParticleOnLatticeGeneric::~SingleParticleOnLatticeGeneric ()
{
  delete[] this->TmpTranslations;
  delete[] this->TmpCoordinates;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] Length;
    }
}

// assignment (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space
SingleParticleOnLatticeGeneric& SingleParticleOnLatticeGeneric::operator = (const SingleParticleOnLatticeGeneric& bosons)
{
  delete[] this->TmpTranslations;
  delete[] this->TmpCoordinates;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
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
  this->Flag = bosons.Flag;
  this->TmpTranslations = new int[TrueDimension];
  this->TmpCoordinates = new int[TrueDimension];
  for (int i=0; i<TrueDimension; ++i)
    {
      TmpTranslations[i]=0;
      TmpCoordinates[i]=0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}


// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space
AbstractHilbertSpace* SingleParticleOnLatticeGeneric::Clone()
{
  return new SingleParticleOnLatticeGeneric(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void SingleParticleOnLatticeGeneric::SetTargetSpace(ParticleOnLattice* targetSpace)
{
  this->TargetSpace=(SingleParticleOnLatticeGeneric*)targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int SingleParticleOnLatticeGeneric::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->GetHilbertSpaceDimension();
}

// get the particle statistic 
//
// return value = particle statistic
int SingleParticleOnLatticeGeneric::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}

// get the quantization axis 
//
// return value = particle statistic
char SingleParticleOnLatticeGeneric::GetLandauGaugeAxis()
{
  return 'y';
}


// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int SingleParticleOnLatticeGeneric::GetHilbertSpaceAdditionalSymmetry()
{
  return SingleParticleOnLatticeGeneric::NoSymmetry;
}



// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number
List<AbstractQuantumNumber*> SingleParticleOnLatticeGeneric::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(this->NbrBosons);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number
AbstractQuantumNumber* SingleParticleOnLatticeGeneric::GetQuantumNumber (int index)
{
  return new NumberParticleQuantumNumber(this->NbrBosons);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace
AbstractHilbertSpace* SingleParticleOnLatticeGeneric::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// get the number of sites
//
// return value = number of sites
int SingleParticleOnLatticeGeneric::GetNbrSites()
{
  return this->NbrStates;
}

// it is possible to change the flux through the simulation cell
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
void SingleParticleOnLatticeGeneric::SetNbrFluxQuanta(int nbrFluxQuanta)
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
void SingleParticleOnLatticeGeneric::SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY)
{
  this->SolenoidX = M_PI*solenoidX;
  this->SolenoidY = M_PI*solenoidY;
  this->SetNbrFluxQuanta(nbrFluxQuanta);
}

// request solenoid fluxes
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
//
void SingleParticleOnLatticeGeneric::GetSolenoidFluxes(double &solenoidX, double &solenoidY)
{
  solenoidX=this->SolenoidX;
  solenoidY=this->SolenoidY;
}



// obtain the current setting of the flux piercing this lattice
int SingleParticleOnLatticeGeneric::GetNbrFluxQuanta()
{
  return this->NbrFluxQuanta;
}



// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
unsigned long SingleParticleOnLatticeGeneric::Ad (unsigned long state, int q, double& coefficient)
{
  coefficient = 0.0;
  return 0;
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

int SingleParticleOnLatticeGeneric::AdAdAA (int index, int m1, int m2, int n1, int n2, double &coefficient)
{
  return this->HilbertSpaceDimension;
}


// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double SingleParticleOnLatticeGeneric::AA (int index, int n1, int n2)
{
  return 0.0;
}


// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on a multiplicative, trivially one here (unreferenced)
// return value = index of the destination state 

int SingleParticleOnLatticeGeneric::AdAd (int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}


// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored (always 1.0)
// return value = index of the destination state 

int SingleParticleOnLatticeGeneric::AdA (int index, int m, int n, double &coefficient)
{
  int State = this->StateDescription[index];
  if (n!=State)
    return this->HilbertSpaceDimension;
  // new state = m
  coefficient=1.0;
  return this->CarefulFindStateIndex(m);
  return (this->HilbertSpaceDimension-1-m);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double SingleParticleOnLatticeGeneric::AdA (int index, int m)
{
  int State = this->StateDescription[index];
  if (m==State)
    return 1.0;
  else
    return 0.0;
}

// apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
// index = index of the state on which the operator has to be applied
// nbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
// qValues = array of quantum numbers where an interaction is present
// interactionPerQ = coefficient U_q of the interaction
//
// no diagonal interaction present, derelict function from inheritance
//
double SingleParticleOnLatticeGeneric::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  return 0.0;
}

// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// sublattice = sublattice index
// 
int SingleParticleOnLatticeGeneric::EncodeQuantumNumber(int posx, int posy, int sublattice, Complex &translationPhase)
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
void SingleParticleOnLatticeGeneric::DecodeQuantumNumber(int q, int &posx, int &posy, int &sublattice)
{
  LatticeGeometry->GetSiteCoordinates(q, TmpCoordinates, sublattice);
  posx=TmpCoordinates[0];
  posy=TmpCoordinates[1];
}

// obtain a list of quantum numbers in state
// index = index of many-body state to be considered
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
// normalization = indicating the multiplicity of the state for bosonic spaces
void SingleParticleOnLatticeGeneric::ListQuantumNumbers(int index, int *quantumNumbers, double &normalization)
{
  normalization=1.0;
  this->ListQuantumNumbers(index, quantumNumbers);
}

// obtain a list of quantum numbers in state
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
void SingleParticleOnLatticeGeneric::ListQuantumNumbers(int index, int *quantumNumbers)
{
  int State = this->StateDescription[index];
  int NbrQ=0;
  quantumNumbers[NbrQ++]=State;
}

// translate a state by a multiple of the lattice vectors
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// translationPhase = returns phase inccurred by translation
// return value = index of translated state
int SingleParticleOnLatticeGeneric::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase)
{
  cout << "Need to implement BosonOnLatticeGeneric::TranslateState"<<endl;
  return 0;
  /*
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
  */
}

// find whether there is a translation vector from state i to state f
// i = index of initial state
// f = index of final state
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// return value = final state can be reached by translation
bool SingleParticleOnLatticeGeneric::IsTranslation(int i, int f, int &shiftX, int &shiftY)
{  
//   int TemporaryStateHighestBit = this->StateHighestBit[i];
//   unsigned long TemporaryState = this->StateDescription[i];
//   int ShiftedStateHighestBit = this->StateHighestBit[f];
//   unsigned long ShiftedState = this->StateDescription[f];
  
  // implementation needed...
  cout << "Implementation of SingleParticleOnLatticeGeneric::IsTranslation required"<<endl;
  
  return true;
}


  // evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex SingleParticleOnLatticeGeneric::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
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

Complex SingleParticleOnLatticeGeneric::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
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
Complex SingleParticleOnLatticeGeneric::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
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

Complex SingleParticleOnLatticeGeneric::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, 
								 int nextCoordinates, int firstComponent, 
								 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void SingleParticleOnLatticeGeneric::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& SingleParticleOnLatticeGeneric::PrintState (ostream& Str, int state)
{
  int TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrStates; ++i)
    if (TmpState==i) Str << "1 ";
    else Str << "0 ";
  return Str;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int SingleParticleOnLatticeGeneric::FindStateIndex(int stateDescription)
{
  return this->HilbertSpaceDimension-1-stateDescription;
}

// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found
int SingleParticleOnLatticeGeneric::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
{
  return this->CarefulFindStateIndex((int)stateDescription);
}

// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found
int SingleParticleOnLatticeGeneric::CarefulFindStateIndex(int stateDescription)
{
  if ((stateDescription >=0)&&(stateDescription<this->HilbertSpaceDimension))
    {
      int Result=this->HilbertSpaceDimension-1-stateDescription;
      if (this->StateDescription[Result]!=stateDescription)
	cout << "Problem finding states stateDescription="<<stateDescription<<" Result="<<Result<<" this->StateDescription[Result]="<<
	     this->StateDescription[Result]<<endl;
      return Result;
    }
  else return this->HilbertSpaceDimension;
}
