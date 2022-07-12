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
#include "HilbertSpace/SingleParticleOnLattice.h"
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

SingleParticleOnLattice::SingleParticleOnLattice ()
{
  this->HilbertSpaceDimension=0;
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
SingleParticleOnLattice::SingleParticleOnLattice (int nbrBosons, int lx, int ly, int nbrFluxQuanta, unsigned long memory, double solenoidX, double solenoidY, char landauGaugeAxis, int nbrSublattices)
{
  //this->NbrBosons = 1;
  this->Lx = lx;
  this->Ly = ly;
  this->NbrSublattices = nbrSublattices;  
  this->NbrStates = Lx*Ly*nbrSublattices;
  this->LandauGaugeAxis=landauGaugeAxis;
  
  this->SetNbrFluxQuanta(nbrFluxQuanta, solenoidX, solenoidY);

  this->HilbertSpaceDimension = NbrStates;
  
  this->Flag.Initialize();
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

#ifdef DEBUG_OUTPUT
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    {
      PrintState(cout, i);
      cout << endl;
    }
#endif

}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

SingleParticleOnLattice::SingleParticleOnLattice(const SingleParticleOnLattice& bosons)
{
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

SingleParticleOnLattice::~SingleParticleOnLattice ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

SingleParticleOnLattice& SingleParticleOnLattice::operator = (const SingleParticleOnLattice& bosons)
{
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* SingleParticleOnLattice::Clone()
{
  return new SingleParticleOnLattice(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> SingleParticleOnLattice::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(1);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* SingleParticleOnLattice::GetQuantumNumber (int index)
{
  return new NumberParticleQuantumNumber(1);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* SingleParticleOnLattice::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// get the number of sites
//
// return value = number of sites
int SingleParticleOnLattice::GetNbrSites()
{
  return this->Lx * this->Ly;
}


// it is possible to change the flux through the simulation cell
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
void SingleParticleOnLattice::SetNbrFluxQuanta(int nbrFluxQuanta)
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
      cout << "Unknown Quantization axis! Exiting SingleParticleOnLattice..."<<endl;
      exit(1);
      break;
    }
}

// change flux through cell and periodic boundary conditions
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
void SingleParticleOnLattice::SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY)
{
  this->SolenoidX = M_PI*solenoidX;
  this->SolenoidY = M_PI*solenoidY;
  this->SetNbrFluxQuanta(nbrFluxQuanta);
}

// request solenoid fluxes
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
//
void SingleParticleOnLattice::GetSolenoidFluxes(double &solenoidX, double &solenoidY)
{
  solenoidX=this->SolenoidX;
  solenoidY=this->SolenoidY;
}

// obtain the current setting of the flux piercing this lattice
int SingleParticleOnLattice::GetNbrFluxQuanta()
{
  return this->NbrFluxQuanta;
}

// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
unsigned long SingleParticleOnLattice::Ad (unsigned long state, int q, double& coefficient)
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

int SingleParticleOnLattice::AdAdAA (int index, int m1, int m2, int n1, int n2, double &coefficient)
{
  return this->HilbertSpaceDimension;
}


// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double SingleParticleOnLattice::AA (int index, int n1, int n2)
{
  return 0.0;
}


// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on a multiplicative, trivially one here (unreferenced)
// return value = index of the destination state 

int SingleParticleOnLattice::AdAd (int m1, int m2, double& coefficient)
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

int SingleParticleOnLattice::AdA (int index, int m, int n, double &coefficient)
{
  int State = this->NbrStates-1-index; // label states in descending order, for compatibility of BosonOnLattice
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
double SingleParticleOnLattice::AdA (int index, int m)
{
  int State = this->NbrStates-1-index;
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
double SingleParticleOnLattice::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  return 0.0;
}

// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// sublattice = sublattice index
// 
int SingleParticleOnLattice::EncodeQuantumNumber(int posx, int posy, int sublattice, Complex &translationPhase)
{
  //cout << "Encoding " << posx<<", "<<posy<<": ";
  int numXTranslations=0, numYTranslations=0;  
  while (posx<0)
    {
      if (posx>0) cout << "Error in EQN"<<endl;
      posx+=this->Lx;
      ++numXTranslations;      
    }
  while (posx>=this->Lx)
    {
      posx-=this->Lx;
      --numXTranslations;
      if (posx<0) cout << "Error in EQN"<<endl;
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
void SingleParticleOnLattice::DecodeQuantumNumber(int q, int &posx, int &posy, int &sublattice)
{
  int m=this->NbrSublattices;
  sublattice=q%m;
  posx=(q%(m*Lx))/m;
  m*=Lx;
  posy=(q%(m*Ly))/m;  
}

// check whether HilbertSpace implements ordering of operators
//
bool SingleParticleOnLattice::HaveOrder ()
{
  return true;
}


// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
int SingleParticleOnLattice::CheckOrder (int* m, int* n, int nbrIndices)
{
  if (nbrIndices>1 /*this->NbrBosons*/) return 0;
  if (m[0] < n[0])
    return 1;
  else if (m[0] > n[0])
    return -1;
  else return 0;
}

// obtain a list of quantum numbers in state
// index = index of many-body state to be considered
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
// normalization = indicating the multiplicity of the state for bosonic spaces
void SingleParticleOnLattice::ListQuantumNumbers(int index, int *quantumNumbers, double &normalization)
{
  quantumNumbers[0]=this->NbrStates-1-index;
  normalization=1.0;
}

// obtain a list of quantum numbers in state
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
void SingleParticleOnLattice::ListQuantumNumbers(int index, int *quantumNumbers)
{
  quantumNumbers[0]=this->NbrStates-1-index;
}


// translate a state by a multiple of the lattice vectors
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// translationPhase = returns phase inccurred by translation
// return value = index of translated state
int SingleParticleOnLattice::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase)
{
  int q0 = this->NbrStates - 1 - index;
  int posx, posy, sublattice;
  DecodeQuantumNumber(q0, posx, posy, sublattice);
  int qF = EncodeQuantumNumber(posx+shiftX, posy+shiftY, sublattice, translationPhase);
  return this->CarefulFindStateIndex(qF);
}

// find whether there is a translation vector from state i to state f
// i = index of initial state
// f = index of final state
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// return value = final state can be reached by translation
bool SingleParticleOnLattice::IsTranslation(int i, int f, int &shiftX, int &shiftY)
{
  if ((shiftX % this->Lx ==0) && (shiftY % this->Ly ==0))
    return true;
  else
    return false;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& SingleParticleOnLattice::PrintState (ostream& Str, int state)
{
  int q0 = this->NbrStates - 1 - state;
  for (int i=0; i<q0; ++i)
    Str << "0 ";
  Str << "1";
  for (int i=q0+1; i<this->NbrStates; ++i)
    Str << " 0";
  return Str;
}




// carefully test whether state is in Hilbert-space and find corresponding state index
//
// quantumNumber = quantum number of the particle
// return value = corresponding index, or dimension of space, if not found
int SingleParticleOnLattice::CarefulFindStateIndex(int quantumNumber)
{
  if (quantumNumber >=0  && quantumNumber < this->HilbertSpaceDimension)
    return this->NbrStates-1-quantumNumber;
  else
    return this->HilbertSpaceDimension;
}



