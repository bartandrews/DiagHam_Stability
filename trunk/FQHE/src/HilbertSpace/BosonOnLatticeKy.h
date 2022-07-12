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
//                  states indexed by y-momentum and x-position               //
//                    (Landau gauge, A(r) = \alpha x \vec e_y)                //
//                                                                            //
//                        last modification : 29/09/2008                      //
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


#ifndef BOSONONLATTICEKY_H
#define BOSONONLATTICEKY_H

#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;

class BosonOnLattice;


class BosonOnLatticeKy : public ParticleOnLattice
{

 protected:

  // array describing each state
  unsigned long* StateDescription;
  // array giving the highest bit reached for a bososn in a given state
  int* StateHighestBit;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given HighestBit sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // pointer to the target space when an index is require after applying basic operation
  ParticleOnLattice* TargetSpace;

  // number of bosons
  int NbrBosons;
  // total number of states
  int NbrStates;
  // length in x-direction
  int Lx;
  // length in y-direction
  int Ly;
  // many-body momentum k_y constraint
  int Ky;
  // number of sites per translation-cell
  int TranslationCell;
  // maximum many-body momentum (Ly / Kmod)
  int Kmax;
  // number of sublattices
  int NbrSublattices;
  // number of flux quanta piercing the simulation cell
  int NbrFluxQuanta;  
  // flux density (flux quanta per unit cell)
  double FluxDensity;
  // number of bits in internal representation
  int NbrCodingBits;

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

  // temperory items used in conversion to full n-body basis:
  ComplexVector TargetVector;
  BosonOnLattice *FullSpace;

 public:

  // default constructor
  //
  BosonOnLatticeKy ();

  // basic constructor -> yields a square lattice in Landau gauge
  // 
  // nbrBosons = number of bosons
  // lx = length of simulation cell in x-direction
  // ly = length of simulation cell in y-direction
  // ky = many-body momentum in y-direction
  // nbrFluxQuanta = number of flux quanta piercing the simulation cell
  // memory = memory that can be allocated for precalculations
  BosonOnLatticeKy (int nbrBosons, int lx, int ly, int ky, int nbrFluxQuanta, unsigned long memory = 10000000,bool normalTranslation= false);


  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnLatticeKy(const BosonOnLatticeKy& bosons);

  // destructor
  //
  virtual ~BosonOnLatticeKy ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLatticeKy& operator = (const BosonOnLatticeKy& bosons);

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
  virtual char GetLandauGaugeAxis();

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

  // get information about any additional symmetry of the Hilbert space
  //
  // return value = symmetry id  
  int GetHilbertSpaceAdditionalSymmetry();

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
  
  // in presence of translation symmetries, it is NOT possible to change the flux through the simulation cell
  // this function is provided for compatibility with the interface, only
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

  // get maximum possible momentum for this geometry
  // return = maximum value of Ky
  virtual int GetMaximumKy();
  
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
	
  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor   
  virtual double ProdA (int index, int* n, int nbrIndices);
  
  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient);

  // apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
  // index = index of the state on which the operator has to be applied
  // nbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
  // qValues = array of quantum numbers where an interaction is present
  // interactionPerQ = coefficient U_q of the interaction
  //
  virtual double AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues);

  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // ky = momentum in y-direction
  // sublattice = sublattice index
  // translationPhase = returns phase occurred from translating the
  //                    site to the fundamental region [0,Lx-1] x [0,Ly-1]
  virtual int EncodeQuantumNumber(int posx, int ky, int sublattice, Complex &translationPhase); 

  // decode a single encoded quantum number q to the set of quantum numbers posx, posy, sublattice
  // posx = position along x-direction
  // ky = momentum in y-direction
  // sublattice = sublattice index
  virtual void DecodeQuantumNumber(int q, int &posx, int &ky, int &sublattice);

  // ky = true momentum in y-direction
  // fluxSubLattice = 'sublattice' index remaining after translation symmetry
  virtual int EncodeCompositeMomentum(int ky, int fluxSubLattice);

  // decode composite ky-momentum
  // cK = composite momentum (momentum plus flux sublattice)
  // ky = true momentum in y-direction
  // fluxSubLattice = 'sublattice' index remaining after translation symmetry
  virtual void DecodeCompositeMomentum(int cK, int &ky, int &fluxSubLattice);

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

  // extract the momentum ky from a quantum number q
  // return: momentum ky (in range 0...Kmax-1)
  int DecodeKy(int q);
  
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
  
  // conversion to generic (full) many-body representation in real-space basis
  // state: many-body state in Ky-momentum basis
  // nbodyBasis: full Hilbert-space in real-space representation
  // returns: vector in many-body basis of targetSpace
  virtual ComplexVector& ConvertToNbodyBasis(ComplexVector& state, ParticleOnLattice &nbodyBasis);

  // conversion to generic (full) many-body representation in real-space basis
  // state: many-body state in Ky-momentum basis
  // nbodyBasis: full Hilbert-space in real-space representation (should be object of type BosonOnLattice)
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // returns: vector in many-body basis of targetSpace
  virtual ComplexVector& ConvertToNbodyBasis(ComplexVector& state, ParticleOnLattice &nbodyBasis, int firstComponent, int nbrComponent);

  // conversion to generic (full) many-body representation in real-space basis
  // state: many-body state in Ky-momentum basis
  // nbodyBasis: full Hilbert-space in real-space representation
  // returns: vector in many-body basis of targetSpace
  virtual ComplexVector* ConvertFromNbodyBasis(ComplexVector* state, BosonOnLattice &nbodyBasis,int nbrVectors, int NbrComponent, AbstractArchitecture * architecture);

  // conversion to generic (full) many-body representation in real-space basis
  // state: many-body state in Ky-momentum basis
  // nbodyBasis: full Hilbert-space in real-space representation (should be object of type BosonOnLattice)
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // returns: vector in many-body basis of targetSpace
  virtual void ConvertFromNbodyBasis(ComplexVector * initialState, ComplexVector * finalState, ParticleOnLattice &nbodyBasis, int nbrVectors, long firstComponent, long nbrComponent);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)	
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnLattice* complementaryHilbertSpace,  ParticleOnLattice* destinationHilbertSpace, ComplexVector& groundState,  HermitianMatrix* densityMatrix);
  
  inline int GetIndexFromQuantumNumber(int q);
 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // highestBit = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int highestBit);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // maxQ = maximum value for the quantum number of a boson in the state
  // currentMaxQ = current max value for the quantum number of bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrBosons, int maxQ, int currentMaxQ, int pos, int currentMomentum);


  // recursively evaluate Hilbert space dimension 
  //
  // nbrBosons = number of bosons
  // maxQ = maximum value for the quantum number of a boson in the state
  // currentMaxQ = current max value for the quantum number of bosons that are still to be placed
  // currentMomentum = current value of the momentum
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrBosons, int maxQ, int currentMaxQ, int currentMomentum);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table  
  void GenerateLookUpTable(unsigned long memory);


  // recursive expansion of an operator-product of creation operators in k
  // in the full n-body basis
  // nbrOperators = number of operators N remaining to be applied
  // quantumNumbers = array of quantum numbers q1,...qN of creation operators
  // state = state to be acted upon
  // prefactor = previous coefficients applied to state
  // 
  // in last stage of recursion, writes to this->TargetVector using the Hilbert-Space this->FullSpace
  void ExpandBasisState (int nbrOperators, int *quantumNumbers, unsigned long state, Complex prefactor);	
  void ProjectBasisState (int nbrOperators, int *quantumNumbers, unsigned long state, ComplexVector *initialVector,ComplexVector * finalVector, long index,int nbrVectors, Complex prefactor);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initialbosonic  state is stored
  // initialStateLzMax = reference on the initial bosonic state maximum Lz value
  // return value = corresponding fermionic state
  unsigned long BosonToFermion(unsigned long*& initialState, int& initialStateLzMax);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initialbosonic  state is stored
  // initialStateLzMax = reference on the initial bosonic state maximum Lz value
  // finalStateHighestBit = reference on the value where highest bit of fermionic represenation should be stored
  // return value = corresponding fermionic state  
  unsigned long BosonToFermion(unsigned long*& initialState, int& initialStateLzMax, int &finalStateHighestBit);

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

inline int BosonOnLatticeKy::GetParticleStatistic()
{
  return ParticleOnLattice::BosonicStatistic;
}


// get the quantization axis 
//
// return value = particle statistic
inline char BosonOnLatticeKy::GetLandauGaugeAxis()
{
  return 'y';
}


// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initialbosonic  state is stored
// initialStateLzMax = reference on the initial bosonic state maximum Lz value
// finalStateHighestBit = reference on the value where highest bit of fermionic represenation should be stored
// return value = corresponding fermionic state

inline unsigned long BosonOnLatticeKy::BosonToFermion(unsigned long*& initialState, int& initialStateLzMax, int &finalStateHighestBit)
{
  // cout << "BosonToFermion: LzMax="<<initialStateLzMax;
  unsigned long TmpState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= initialStateLzMax; ++i)
    {
      TmpState |= ((1ul << initialState[i]) - 1ul) << Shift;
      Shift += initialState[i];
      ++Shift;
    }
  // cout << ", LzMax2="<<initialStateLzMax<<", Shift="<<Shift<<endl;
  finalStateHighestBit = (Shift>1?Shift-2:0);
  return TmpState;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initialbosonic  state is stored
// initialStateLzMax = reference on the initial bosonic state maximum Lz value
// return value = corresponding fermionic state

inline unsigned long BosonOnLatticeKy::BosonToFermion(unsigned long*& initialState, int& initialStateLzMax)
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


// convert a fermionic state into its bosonic counterpart
//
// initialState = initial fermionic state
// initialStateLzMax = initial fermionic state maximum Lz value
// finalState = reference on the array where the bosonic state has to be stored
// finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored

inline void BosonOnLatticeKy::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax)
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
      //cout << "finalState["<<finalStateLzMax<<"]="<<finalState[finalStateLzMax]<<endl;
      ++TmpPower;
      initialState >>= TmpPower;
      ++finalStateLzMax;
      initialStateLzMax -= TmpPower;
    }
  --finalStateLzMax;
}

inline int BosonOnLatticeKy::GetIndexFromQuantumNumber(int q)
{
  this->TemporaryStateHighestBit = q;
  for(int i = 0; i <q ; i++)
    this->TemporaryState[i] = 0;
  this->TemporaryState[q] = 1;
  int TmpIndex = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1); 
  cout <<q<<" "<<TmpIndex<<endl;
  return TmpIndex;
}

#endif


