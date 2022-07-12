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


#ifndef BOSONONLATTICE_H
#define BOSONONLATTICE_H


#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "HilbertSpace/HardCoreBoson.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnLattice : public ParticleOnLattice
{

 protected:

  // the pseudo-fermionic Hilbert space associated to the bosonic one
  HardCoreBoson* HardCoreBasis;

  // number of bosons
  int NbrBosons;
  // total number of states
  int NbrStates;
  // length in x-direction
  int Lx;
  // length in y-direction
  int Ly;
  // number of sublattices
  int NbrSublattices;
  // number of flux quanta piercing the simulation cell
  int NbrFluxQuanta;
  // flux density (flux quanta per unit cell)
  double FluxDensity;
  // direction of Landau-gauge (x/y)
  char LandauGaugeAxis;

  // solenoid flux through torus around periodic boundary conditions (units of pi)
  double SolenoidX;
  double SolenoidY;

  // phases occurred when translating the system by a full system length (at y=1, or x=1)
  // in the x-direction
  Complex LxTranslationPhase;
  // and in the y-direction
  Complex LyTranslationPhase;  

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

 public:

  // default constructor
  //
  BosonOnLattice ();

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
  BosonOnLattice (int nbrBosons, int lx, int ly, int nbrFluxQuanta, unsigned long memory = 10000000, double solenoidX=0.0, double solenoidY=0.0, char landauGaugeAxis='y', int nbrSublattices=1);


  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnLattice(const BosonOnLattice& bosons);

  // destructor
  //
  virtual ~BosonOnLattice ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLattice& operator = (const BosonOnLattice& bosons);

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
  virtual char GetLandauGaugeAxis() {return this->LandauGaugeAxis;}

  // get the number of sites
  //
  // return value = number of sites
  virtual int GetNbrSites();

  // get the minimum number of particles
  //
  // return value = smallest number of particles within Hilbert Space
  virtual int GetMinNbrParticles() {return this->NbrBosons;}

  // get the number of sublattices
  //
  // return value = number of sublattices
  virtual int GetNbrSublattices(){return this->NbrSublattices;}
  
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

  // calculate (possibly non-local) density-density interactions \sum q V_{q1,q2} : n_q1 n_q2 :
  // index = index of the state on which the operator has to be applied
  // nbrInteraction = number of q-values in sum
  // interactionPerQ12 = coefficient V_(q1, q2) of the interaction
  // q12Values = array of quantum numbers of the orbitals in tuples (q1, q2), 2*nbrInteraction entries in total
  //
  virtual double RhoRhoDiagonal(int index, int nbrInteraction, double *interactionPerQ12, int *q12Values);

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
  // return value = index of translated state
  virtual int TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase);

  // find whether there is a translation vector from state i to state f
  // i = index of initial state
  // f = index of final state
  // shiftX = length of translation in x-direction
  // shiftY = length of translation in y-direction
  // return value = final state can be reached by translation
  virtual bool IsTranslation(int i, int f, int &shiftX, int &shiftY);


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

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

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

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial  bosonic state
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored
  void ConvertToMonomial(unsigned long* initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial bosonic state in its fermionic representation
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored
  void ConvertToMonomial(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a bosonic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // return value = bosonic state in its fermionic representation
  unsigned long ConvertFromMonomial(unsigned long* initialState);

  // check that the product firstState*secondState is in the lexicographical order
  //
  // firstState = array where the monomial representation  of a state is stored
  // secondState = array where the monomial representation  of another state is stored
  // return value = true if the product firstState*secondState is in the lexicographical order
  virtual bool CheckLexiOrder(int * firstState,unsigned long* secondState,int TailleEgal);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnLattice* complementaryHilbertSpace,  ParticleOnLattice* destinationHilbertSpace,
								  ComplexVector& groundState,  HermitianMatrix* densityMatrix);

public:
  virtual void SymmetrizeU1U1State (ComplexVector& symmetrizedVector, ComplexVector& leftVector, ComplexVector& rightVector, BosonOnLattice* leftSpace, BosonOnLattice* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents);
  
};


// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnLattice::GetParticleStatistic()
{
  return ParticleOnLattice::BosonicStatistic;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initialbosonic  state is stored
// initialStateLzMax = reference on the initial bosonic state maximum Lz value
// return value = corresponding fermionic state

inline unsigned long BosonOnLattice::BosonToFermion(unsigned long*& initialState, int& initialStateLzMax)
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

inline void BosonOnLattice::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax)
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

// convert a bosonic state to its monomial representation
//
// initialState = initial  bosonic state
// initialStateLzMax = initial bosonic state maximum Lz value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnLattice::ConvertToMonomial(unsigned long* initialState, int initialStateLzMax, unsigned long*& finalState)
{
  int Index = 0;
  for (int i = initialStateLzMax; i >= 0; --i)
    for (unsigned long j = 0l; j < initialState[i]; ++j)
      finalState[Index++] = i;
}

// convert a bosonic state to its monomial representation
//
// initialState = initial bosonic state in its fermionic representation
// initialStateLzMax = initial bosonic state maximum Lz value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnLattice::ConvertToMonomial(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState)
{
  int Index = 0;
  int TmpLz = initialStateLzMax - this->NbrBosons + 1;
  while (initialStateLzMax >= 0)
    {
      while ((initialStateLzMax >= 0) && (((initialState >> initialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalState[Index++] = TmpLz;
	  --initialStateLzMax;
	}
      while ((initialStateLzMax >= 0) && (((initialState >> initialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --initialStateLzMax;
	}
    }
}

// convert a bosonic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// return value = bosonic state in its fermionic representation

inline unsigned long BosonOnLattice::ConvertFromMonomial(unsigned long* initialState)
{
  unsigned long Tmp = 0x0ul;
  for (int i = 0; i < this->NbrBosons; ++i)
    Tmp |= 0x1ul << (initialState[i] + ((unsigned long) (this->NbrBosons - i)) - 1ul);
  return Tmp;
}


// check that the product firstState*secondState is in the lexicographical order
//
// firstState = array where the monomial representation  of a state is stored
// secondState = array where the monomial representation  of another state is stored
// return value = true if the product firstState*secondState is in the lexicographical order

inline bool BosonOnLattice::CheckLexiOrder(int * egal,unsigned long* secondState,int TailleEgal)
{
  for (int index = 0; index < TailleEgal; ++index)
    {
      if (secondState[egal[index]]<secondState[egal[index]+1])
	{
	  return false;
	}
    }
  return true;
}


#endif


