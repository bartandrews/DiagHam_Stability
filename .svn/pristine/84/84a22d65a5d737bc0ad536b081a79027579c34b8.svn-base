////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//  class of bosons on linearly indexed system for a system size such that    //
//      NbrStates + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)      //
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


#ifndef BOSONONLATTICEGENERIC_H
#define BOSONONLATTICEGENERIC_H


#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "HilbertSpace/HardCoreBoson.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Tools/FQHESpectrum/LatticePhases.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnLatticeGeneric : public ParticleOnLattice
{

 protected:

  // the pseudo-fermionic Hilbert space associated to the bosonic one
  HardCoreBoson* HardCoreBasis;

  // underlying lattice geometry
  LatticePhases* LatticeGeometry;

  // number of bosons
  int NbrBosons;
  // total number of states
  int NbrStates;
  // lengths in units of lattice-vectors
  int* Length;
  // number of sublattices
  int NbrSublattices;

  // true lattice dimension
  int TrueDimension;
  
  // number of flux quanta piercing the simulation cell
  int NbrFluxQuanta;
  // flux density (flux quanta per unit cell)
  double FluxDensity;

  // solenoid flux through torus around periodic boundary conditions (units of pi)
  double SolenoidX;
  double SolenoidY;

  // phases occurred when translating the system by a full system length (at y=1, or x=1)
  // in the x-direction
  Complex LxTranslationPhase;
  // and in the y-direction
  Complex LyTranslationPhase;

  // internal storage for coordinates/translations for access to LatticePhases
  int* TmpCoordinates;
  int* TmpTranslations;

  // pointer to an integer which indicate which coordinates are kept for the next time step iteration
  int* KeptCoordinates;
  // minors of permanents used for the time coherent wave function evaluation
  Complex** Minors;
  
  // temporary state used when applying operators
  unsigned long* TemporaryState;
  int TemporaryStateHighestBit;
  // state used for translating a state
  unsigned long* ShiftedState;
  int ShiftedStateHighestBit;
  // temporary state used when applying ProdA operator
  unsigned long* ProdATemporaryState;
  int ProdATemporaryStateHighestBit;

  // fields for translations routine
  int *CurrentTranslation;

  int *CurrentMappings;
  Complex* CurrentTranslationPhases;

 public:

  // default constructor
  //
  BosonOnLatticeGeneric ();

  // basic constructor -> yields a square lattice in Landau gauge
  // 
  // nbrBosons = number of bosons
  // latticeGeometry = geometry of lattice system is living on
  // nbrFluxQuanta = number of flux quanta piercing the simulation cell
  // memory = memory that can be allocated for precalculations
  // solenoidX = solenoid flux through lattice in x-direction (in units of pi)
  // solenoidY = solenoid flux through lattice in y-direction (in units of pi)
  BosonOnLatticeGeneric (int nbrBosons, LatticePhases *latticeGeometry, int nbrFluxQuanta, unsigned long memory = 10000000, double solenoidX=0.0, double solenoidY=0.0);


  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnLatticeGeneric(const BosonOnLatticeGeneric& bosons);

  // destructor
  //
  virtual ~BosonOnLatticeGeneric ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLatticeGeneric& operator = (const BosonOnLatticeGeneric& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // get the quantization axis 
  //
  // return value = particle statistic
  virtual char GetLandauGaugeAxis(){return '\0';}

    // get the number of sites
  //
  // return value = number of sites
  virtual int GetNbrSites();

  // get the minimum number of particles
  //
  // return value = smallest number of particles within Hilbert Space
  virtual int GetMinNbrParticles(){return this->NbrBosons;}

  // get the number of sublattices
  //
  // return value = number of sublattices
  virtual int GetNbrSublattices(){return this->NbrSublattices;}

  // check whether HilbertSpace implements ordering of operators
  //
  virtual bool HaveOrder ();
  
  // check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
  virtual int CheckOrder (int* m, int* n, int nbrIndices);

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // it is possible to change the flux through the simulation cell
  // Attention: this does require the Hamiltonian to be recalculated!!
  // nbrFluxQuanta = number of quanta of flux piercing the simulation cell
  virtual void SetNbrFluxQuanta(int nbrFluxQuanta);

  // change flux through cell and periodic boundary conditions
  // Attention: this does require the Hamiltonian to be recalculated!!
  // nbrFluxQuanta = number of quanta of flux piercing the simulation cell
  // solenoidX = new solenoid flux through torus in x-direction
  // solenoidY = new solenoid flux through torus in y-direction
  virtual void SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY);

  // request solenoid fluxes
  // solenoidX = new solenoid flux through torus in x-direction
  // solenoidY = new solenoid flux through torus in y-direction
  //
  virtual void GetSolenoidFluxes(double &solenoidX, double &solenoidY);

  // obtain the current setting of the flux piercing this lattice
  virtual int GetNbrFluxQuanta();


  // apply creation operator to a word, using the conventions
  // for state-coding and quantum numbers of this space
  // state = word to be acted upon
  // q = quantum number of boson to be added
  // coefficient = reference on the double where the multiplicative factor has to be stored
  virtual unsigned long Ad (unsigned long state, int q, double& coefficient);

  // apply a^+_q1 a^+_q2 a_r1 a_r2 operator to a given state (with q1+q2=r1+r2)
  //
  // index = index of the state on which the operator has to be applied
  // q1 = first index for creation operator
  // q2 = second index for creation operator
  // r1 = first index for annihilation operator
  // r2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int q1, int q2, int r1, int r2, double& coefficient);

  // apply a_r1 a_r2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // r1 = first index for annihilation operator
  // r2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int r1, int r2);

  // apply a^+_q1 a^+_q2 operator to the state produced using AA method (without destroying it)
  //
  // q1 = first index for creation operator
  // q2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int q1, int q2, double& coefficient);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (int index, int m);

  // apply a^+_q a_r operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // q = index of the creation operator
  // r = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdA (int index, int q, int r, double& coefficient);

  // apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
  // index = index of the state on which the operator has to be applied
  // nbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
  // qValues = array of quantum numbers where an interaction is present
  // interactionPerQ = coefficient U_q of the interaction
  //
  virtual double AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues);

  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // posy = position along y-direction
  // sublattice = sublattice index
  // translationPhase = returns phase occurred from translating the
  //                    site to the fundamental region [0,Lx-1] x [0,Ly-1]
  virtual int EncodeQuantumNumber(int posx, int posy, int sublattice, Complex &translationPhase);

  // decode a single encoded quantum number q to the set of quantum numbers posx, posy, sublattice
  // posx = position along x-direction
  // posy = position along y-direction
  // sublattice = sublattice index
  virtual void DecodeQuantumNumber(int q, int &posx, int &posy, int &sublattice);

  // obtain a list of quantum numbers in state
  // index = index of many-body state to be considered
  // quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
  // normalization = indicating the multiplicity of the state for bosonic spaces
  virtual void ListQuantumNumbers(int index, int *quantumNumbers, double &normalization);

  // obtain a list of quantum numbers in state
  // index = index of many-body state to be considered
  // quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
  virtual void ListQuantumNumbers(int index, int *quantumNumbers);


  // translate a state by a multiple of the lattice vectors
  // shiftX = length of translation in x-direction
  // shiftY = length of translation in y-direction
  // translationPhase = returns phase inccurred by translation
  // gaugePhases = phases of a gauge transform that was applied previously to the state
  // return value = index of translated state
  virtual int TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase);

  // version with non-standard boundary conditions
  // translationX = array: Exp (i*[k]) gives the phase applied for positive crossing of cell boundary for state with k=y*NbrSubLattices + s
  // translationY = array: Exp (i*[k]) gives the phase applied for positive crossing of cell boundary for state with k=x*NbrSubLattices + s
  virtual int TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase, double *translationX, double *translationY);

  // find whether there is a translation vector from state i to state f
  // i = index of initial state
  // f = index of final state
  // shiftX = length of translation in x-direction
  // shiftY = length of translation in y-direction
  // return value = final state can be reached by translation
  virtual bool IsTranslation(int i, int f, int &shiftX, int &shiftY);

  // apply a gauge transformation
  // phases = phases in array ordered according to the quantum number q
  // input = vector that has to be transformed according to that gauge
  virtual ComplexVector& GaugeTransformVector(double *phases, ComplexVector& input);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // carefully test whether state is in Hilbert-space and find corresponding state index
  //
  // stateDescription = unsigned integer describing the state
  // highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
  // return value = corresponding index, or dimension of space, if not found
  virtual int CarefulFindStateIndex(unsigned long stateDescription, int highestBit);

 protected:

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initialbosonic  state is stored
  // initialStateLzMax = reference on the initial bosonic state maximum Lz value
  // return value = corresponding fermionic state
  unsigned long BosonToFermion(unsigned long*& initialState, int& initialStateLzMax);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax = initial fermionic state maximum Lz value
  // finalState = reference on the array where the bosonic state has to be stored
  // finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored
  void FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax);

};


// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnLatticeGeneric::GetParticleStatistic()
{
  return ParticleOnLattice::BosonicStatistic;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initialbosonic  state is stored
// initialStateLzMax = reference on the initial bosonic state maximum Lz value
// return value = corresponding fermionic state

inline unsigned long BosonOnLatticeGeneric::BosonToFermion(unsigned long*& initialState, int& initialStateLzMax)
{
  unsigned long TmpState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= initialStateLzMax; ++i)
    {
      TmpState |= ((1ul << initialState[i]) - 1ul) << Shift;
      Shift += initialState[i];
      ++Shift;
    }
  return TmpState;
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax = initial fermionic state maximum Lz value
// finalState = reference on the array where the bosonic state has to be stored
// finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored

inline void BosonOnLatticeGeneric::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax)
{
  finalStateLzMax = 0;
  while (initialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState - 1ul) ^ (~initialState);
      TmpState &= ~(TmpState >> 1);
//      cout << hex << initialState << "  " << TmpState << dec << endl;
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
//      cout << TmpPower << endl;
      finalState[finalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState >>= TmpPower;
      ++finalStateLzMax;
      initialStateLzMax -= TmpPower;
    }
  --finalStateLzMax;
}

#endif


