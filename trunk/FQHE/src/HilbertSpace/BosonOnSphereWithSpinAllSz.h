////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of bosons on sphere                        //
//                                                                            //
//                        last modification : 05/07/2002                      //
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


#ifndef BOSONONSPHEREWITHSPINALLSZ_H
#define BOSONONSPHEREWITHSPINALLSZ_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

using std::cout;
using std::endl;
using std::hex;
using std::dec;

class BosonOnSphere;
class BosonOnSphereWithSpin;


class BosonOnSphereWithSpinAllSz :  public ParticleOnSphereWithSpin
{

 protected:

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // twice the momentum total value
  int TotalLz;
  // momentum total value shifted by LzMax / 2 * NbrBosons
  int ShiftedTotalLz;
  // twice the maximum Lz value reached by a boson
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;
  // parity flag for Sz (-1=no parity, 0=even, 1=odd)
  int SzParity;
  

  // array describing each state (full storage, temporary use)
  unsigned** StateDescription;
  // array describing each state up / down
  unsigned long* StateDescriptionUp;
  unsigned long* StateDescriptionDown;
  // array giving maximum Lz value reached for a boson in a given state for up and down spin
  unsigned* StateLzMaxUp;
  unsigned* StateLzMaxDown;
  // number of up bosons in state
  unsigned* StateNbrUp;
  // array to store condensed info on states after state generation
  unsigned *StateInfo;

  // Lookup-table for base indices for given sequence of up spins
  unsigned *UpSpinLookUpTable;
  // Lookup-table for offsets for given sequence of down spins
  unsigned *DownSpinLookUpTable;
  // elements of spaces that are being tensored together for both up and down spins (provides index)
  unsigned long* TensoredElements;
  // maximum lz value of TensoredElements
  int *TensoredLzMax;
  // number of tensored elements
  unsigned NbrTensoredElements;

  // maximum shift used for searching a position in the look-up table of array TensoredElements
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // table for coherence factors
  double *CoherenceFactors;
    
  // pointer to an integer which indicate which coordinates are kept for the next time step iteration
  int* KeptCoordinates;
  // minors of permanents used for the time coherent wave function evaluation
  Complex** Minors;

  // temporary state used when applying operators
  unsigned* TemporaryState;
  // temporary state used when applying ProdA operator
  unsigned* ProdATemporaryState;
  // number of up spins in temporary state
  unsigned ProdATemporaryStateNbrUp;

  // flag indicating partial wavefunction calculation
  bool SubspaceRestriction;
  // subspace for this restriction
  int SubspaceSz;

  // target space for operations leaving the Hilbert-space
  BosonOnSphereWithSpinAllSz *TargetSpace;
  
 public:

  // default constructor
  //
  BosonOnSphereWithSpinAllSz ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSpinAllSz (int nbrBosons, int totalLz, int lzMax, unsigned long memory = 10000000);

  // basic constructor with pair parity
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // szParity = parity of the total Sz : (-1)^Sz == 1 (even parity 0 generalized for all N by: N_\up mod 2 == \lgauss N/2 \rgauss mod 2)
  // memory = amount of memory granted for precalculations
  //
  BosonOnSphereWithSpinAllSz (int nbrBosons, int totalLz, int lzMax, int szParity, unsigned long memory = 10000000);


  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSpinAllSz(const BosonOnSphereWithSpinAllSz& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSpinAllSz ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSpinAllSz& operator = (const BosonOnSphereWithSpinAllSz& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

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

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);
  
  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // apply a^+_m1_d a^+_m2_d a_n1_d a_n2_d operator to a given state (with m1+m2=n1+n2, only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m1_u a^+_m2_u a_n1_u a_n2_u operator to a given state (with m1+m2=n1+n2, only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_d_m1 a^+_u_m2 a_d_n1 a_u_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m_d a_m_d operator to a given state (only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AduAu (int index, int m);

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAu (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_u a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient);

  // apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // sigma1 = SU(2) index for the first annihilation operator
  // sigma2 = SU(2) index for the second annihilation operator
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double AuAu (int index, int n1, int n2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AdAd (int index, int n1, int n2);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AuAd (int index, int n1, int n2);

  // apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // sigma1 = SU(2) index for the first creation operator
  // sigma2 = SU(2) index for the second creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient);

  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAdd (int m1, int m2, double& coefficient);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdd (int m1, int m2, double& coefficient);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int* spinIndices, int nbrIndices);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient);

  
  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // print a given State
  //
  // Str = reference on current output stream 
  // myState = explicit form of the state to print
  // return value = reference on current output stream 
  
  ostream& PrintState (ostream& Str, unsigned *myState);


  // set a restriction to evaluate the wavefunction on a given subspace of fixed sz
  // twoSz = subspace to restrict to
  // restriction = flag whether restriction shall be set (true) or deleted (false)
  virtual bool WaveFunctionSubSpace(int twoSz, bool restriction=true);

  // get weight of wavefunction in current subspace
  // state = vector to be considered
  //
  // return = weight
  virtual double GetWeightInSubSpace(RealVector& state);
  
  // evaluate wave function in real space using a given basis
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis);

  // evaluate wave function in real space using a given basis, using time coherence
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // nextCoordinates = index of the coordinate that will be changed during the next time iteration
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							 int nextCoordinates);

  // evaluate wave function in real space using a given basis and only for a given range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
					int firstComponent, int nbrComponent);

  // evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // nextCoordinates = index of the coordinate that will be changed during the next time iteration
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							 int nextCoordinates, int firstComponent, int nbrComponent);

  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  virtual  void InitializeWaveFunctionEvaluation (bool timeCoherence = false);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState);

  // Project the state from the tunneling space (all Sz's)
  // to the space with the fixed projection of Sz (given by SzValue)
  //
  // state = state that needs to be projected
  // su2Space = the subspace onto which the projection is carried out
  // SzValue = the desired value of Sz

  virtual RealVector ForgeSU2FromTunneling(RealVector& state, BosonOnSphereWithSpin& su2Space, int SzValue);

  // Project the state from the tunneling space (all Sz's)
  // to the U(1) space (u1Space)
  //
  // state = state that needs to be projected
  // u1Space = the subspace onto which the projection is carried out
  virtual RealVector ForgeU1FromTunneling(RealVector& state, BosonOnSphere& u1Space);

  // Calculate mean value <Sx> in a given state
  //
  // state = given state
  virtual double MeanSxValue(RealVector& state);

  // Calculate mean value <Sz> in a given state
  //
  // state = given state
  virtual double MeanSzValue(RealVector& state);

  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(char* stateDescription);

  // find state index
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  int FindStateIndex(unsigned* stateDescription, unsigned nbrUp);


  // get Lz component of a component
  //
  // j = index of the component in Hilbert space
  // return value = twice the Lz component
  virtual int GetLzValue(int j=0);

 protected:

  // find state index
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown, int lzMaxUp, int lzMaxDown);

  
  // find index of a tensored configuration
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  unsigned FindTensoredIndex(unsigned long stateDescription, int lzmax);


  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrBosons = number of bosons
  // currentLzMax = two times momentum maximum value for a boson plus one at current stage in recursion
  // totalLz = momentum total value plus nbrBosons * (momentum maximum value for a boson + 1)
  // return value = Hilbert space dimension
  int ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int currentLzMax, int totalLz, int l);

  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrBosons = number of bosons
  // currentLzMax = two times momentum maximum value for a boson plus one at current stage in recursion
  // totalLz = momentum total value plus nbrBosons * (momentum maximum value for a boson + 1)
  // totalNbrUp = number of up-particles placed so far
  // szParity = target parity of totalNbrUp
  // return value = Hilbert space dimension
  int ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int currentLzMax, int totalLz, int totalNbrUp, int szParity, int l);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(int memory);
    
  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual int EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // szParity = parity of the total Sz : (-1)^Sz == 1?
  // return value = Hilbert space dimension
  virtual int EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int szParity);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // currentLzMax = momentum maximum value for bosons that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // currentLzMax = momentum maximum value for bosons that are still to be placed
  // totalLz = momentum total value
  // totalNbrUp = number of up-particles placed so far
  // szParity = target parity of totalNbrUp
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int totalNbrUp, int szParity, int pos);

  // convert a bosonic state into its fermionic counterpart
  //
  // finalStateUp = return value of bit-coded state of up-bosons
  // finalStateDown = return value of bit-coded state of down-bosons
  // finalLzMaxUp = highest bit in fermionic coding finalStateUp
  // finalLzMaxDown = highest bit in fermionic coding finalStateDown
  // initialState = reference on the array where initialbosonic  state is stored
  // initialStateNbrUp = reference on the number of up-bosons in initial state maximum Lz value
  
  void BosonToFermion(unsigned long &finalStateUp, unsigned long &finalStateDown, int &finalLzMaxUp, int &finalLzMaxDown, unsigned*& initialState, unsigned& initialStateNbrUp);


  // convert a fermionic state into its bosonic counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax = initial fermionic state maximum Lz value
  // finalState = reference on the array where the bosonic state has to be stored
  // finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored
  
  void FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown, unsigned initialInfo, unsigned*& finalState, int &finalStateLzMaxUp, int &finalStateLzMaxDown, unsigned& finalStateNbrUp);

  // get LzMax value for a given state
  // index = index of state to analyse
  // return = lzMax value (max of up and down)
  int GetStateLzMax(int index);

  // sort an array and reflect permutations in auxiliary array
  //
  // length = length of arrays
  // sortArray = array to be sorted
  // auxArray = auxiliary array
  //
  void ShellSortAux(unsigned length, unsigned long* sortArray, unsigned *auxArray, int *auxArray2);


};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereWithSpinAllSz::GetParticleStatistic()
{
  return ParticleOnSphere::BosonicStatistic;
}


// convert a bosonic state into its fermionic counterpart
//
// finalStateUp = return value of bit-coded state of up-bosons
// finalStateDown = return value of bit-coded state of down-bosons
// finalLzMaxUp = highest bit in fermionic coding finalStateUp
// finalLzMaxDown = highest bit in fermionic coding finalStateDown
// initialState = reference on the array where initialbosonic  state is stored
// initialStateNbrUp = reference on the number of up-bosons in initial state maximum Lz value

inline void BosonOnSphereWithSpinAllSz::BosonToFermion(unsigned long &finalStateUp, unsigned long &finalStateDown, int &finalLzMaxUp, int &finalLzMaxDown, unsigned*& initialState, unsigned& initialStateNbrUp)
{
/*   cout << "NbrUp="<<initialStateNbrUp <<", InitialState = |"; */
/*   for (int i=0; i<=LzMax; ++i) */
/*     cout << " " << (initialState[i] >> 16)<< "u "<< (initialState[i] & 0xffff)<<"d |"; */
/*   cout << endl; */

  finalStateUp = 0x0ul;
  unsigned ShiftUp = 0;
  finalStateDown = 0x0ul;
  unsigned ShiftDown = 0;
  finalLzMaxUp = 0;
  finalLzMaxDown = 0;
  int RemainingNbrUp=initialStateNbrUp;
  int RemainingNbrDown=NbrBosons-initialStateNbrUp;
  int TmpUp, TmpDown;
  for (int i = 0; (i <= this->LzMax)&&((RemainingNbrUp!=0)||(RemainingNbrDown!=0)); ++i)
    {
      TmpUp = initialState[i]>>16;
      TmpDown = initialState[i]&0xffff;
      finalStateUp |= ((1ul << TmpUp) - 1ul) << ShiftUp;
      ShiftUp += TmpUp;
      RemainingNbrUp -= TmpUp;
      ++ShiftUp;
      finalLzMaxUp = (TmpUp>0)*i+(TmpUp==0)*finalLzMaxUp;
      finalStateDown |= ((1ul << TmpDown) - 1ul) << ShiftDown;
      ShiftDown += TmpDown;
      RemainingNbrDown -= TmpDown;
      ++ShiftDown;
      finalLzMaxDown = (TmpDown>0)*i+(TmpDown==0)*finalLzMaxDown;
    }
  finalLzMaxUp += initialStateNbrUp-(initialStateNbrUp!=0);
  finalLzMaxDown += this->NbrBosons - initialStateNbrUp - 1 + (initialStateNbrUp==(unsigned)this->NbrBosons);
  //  this->PrintState(cout, initialState) << " " << std::hex << finalStateUp << " " << finalStateDown << std::dec << " " << finalLzMaxUp <<" " << finalLzMaxDown<<endl;
  return;
}


// convert a fermionic state into its bosonic counterpart
//
// initialState = initial fermionic state
// initialStateLzMax = initial fermionic state maximum Lz value
// finalState = reference on the array where the bosonic state has to be stored
// finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored

inline void BosonOnSphereWithSpinAllSz::FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown, unsigned initialInfo, unsigned*& finalState, int &finalStateLzMaxUp, int &finalStateLzMaxDown, unsigned& finalStateNbrUp)
{
  //cout << "FermionToBoson :: initialStateUp =" << hex << initialStateUp << ", initialStateDown="<<initialStateDown<<dec<<endl;
  finalStateNbrUp=0;
  int InitialStateLzMax = initialInfo >> 20;
  finalStateLzMaxUp = 0;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateUp - 1ul) ^ (~initialStateUp);
      TmpState &= ~(TmpState >> 1);
      // cout << hex << initialStateUp << "  " << TmpState << dec << endl;
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
      //cout << "initialStateLzMaxUp="<<InitialStateLzMax<<" - setting finalState["<<finalStateLzMaxUp<<"]="<<TmpPower<<"u"<<endl;
      finalState[finalStateLzMaxUp] = TmpPower << 16;
      finalStateNbrUp+=TmpPower;
      ++TmpPower;
      initialStateUp >>= TmpPower;
      ++finalStateLzMaxUp;
      InitialStateLzMax -= TmpPower;
    }
  --finalStateLzMaxUp;

/*   cout << "FinalStateUp ="; */
/*   for (int i=0; i<=LzMax; ++i) */
/*     cout << " " << (finalState[i] >> 16); */
/*   cout << endl; */
  
  for (int i=finalStateLzMaxUp+1; i<NbrLzValue; ++i)
    finalState[i]=0;
  finalStateLzMaxDown = 0;
  InitialStateLzMax = (initialInfo >> 10)&0x03ff;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateDown - 1ul) ^ (~initialStateDown);
      TmpState &= ~(TmpState >> 1);
      // cout << hex << initialStateDown << "  " << TmpState << dec << endl;
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
//      cout << "initialStateLzMaxDown="<<InitialStateLzMax<<" - setting finalState["<<finalStateLzMaxDown<<"]="<<TmpPower<<"d"<<endl;
      finalState[finalStateLzMaxDown] |= TmpPower;
      ++TmpPower;
      initialStateDown >>= TmpPower;
      ++finalStateLzMaxDown;
      InitialStateLzMax -= TmpPower;
    }
  --finalStateLzMaxDown;
/*   cout << "FinalStateDown ="; */
/*   for (int i=0; i<=LzMax; ++i) */
/*     cout << " " << (finalState[i]& 0xffff); */
/*   cout << endl; */
  // this->PrintState(cout, finalState)<<endl;
}


// get LzMax value for a given state
// index = index of state to analyse
// return = lzMax value (max of up and down)

inline int BosonOnSphereWithSpinAllSz::GetStateLzMax(int index)
{
  unsigned Info = StateInfo[index];
  int LzMaxUp = Info >> 20;
  int LzMaxDown = (Info >> 10)&0x03ff;
  unsigned StateNbrUp = Info & 0x03ff;
  LzMaxUp -= StateNbrUp-(StateNbrUp!=0);
  LzMaxDown -= this->NbrBosons - StateNbrUp - 1 + (StateNbrUp==(unsigned)this->NbrBosons);
  return std::max(LzMaxUp, LzMaxDown);
}


// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

inline int BosonOnSphereWithSpinAllSz::FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown, int lzMaxUp, int lzMaxDown)
{
  // cout << "FindStateIndex: "<<hex<<stateDescriptionUp<<" "<<stateDescriptionDown<<" " <<dec << lzMaxUp <<" "<<lzMaxDown;
  int Base = this->UpSpinLookUpTable[this->FindTensoredIndex(stateDescriptionUp, lzMaxUp)];
  int Offset = this->DownSpinLookUpTable[this->FindTensoredIndex(stateDescriptionDown, lzMaxDown)];
  return Base+Offset;
}

// find index of a tensored configuration
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index
inline unsigned BosonOnSphereWithSpinAllSz::FindTensoredIndex(unsigned long elementDescription, int lzmax)
{
  // long PosMax = this->NbrTensoredElements-1;
  long PosMin = elementDescription >> this->LookUpTableShift[lzmax];
  //long PosMin = 0;
  long PosMax = this->LookUpTable[lzmax][PosMin];  
  PosMin = this->LookUpTable[lzmax][PosMin + 1];

  // cout << "Initial propositions searching "<<std::hex<<elementDescription<<std::dec<<" with lzMax="<<lzmax<< " PosMax="<<PosMax<<" PosMin="<<PosMin<<endl;
  long PosMid = (PosMin + PosMax + 1) >> 1;
  unsigned long CurrentElement = this->TensoredElements[PosMid];
  while ((PosMax != PosMid) && (CurrentElement != elementDescription))
    {
      if (CurrentElement > elementDescription)
	{
	  PosMin = PosMid;
	}
      else
	{
	  PosMax = PosMid;
	} 
      PosMid = (PosMin + PosMax + 1) >> 1;
      CurrentElement = this->TensoredElements[PosMid];
    }
  if (CurrentElement == elementDescription)
    {
      return PosMid;
    }
  else
    {
      return PosMin;
    }
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(2) index for the first annihilation operator
// sigma2 = SU(2) index for the second annihilation operator
// return value =  multiplicative factor 

inline double BosonOnSphereWithSpinAllSz::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  cout << "warning : AsigmaAsigma not defined in BosonOnSphereWithSpinAllSz" << endl;
  return 0.0;
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(3) index for the first creation operator
// sigma2 = SU(3) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int BosonOnSphereWithSpinAllSz::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  cout << "warning : AdsigmaAdsigma not defined in BosonOnSphereWithSpinAllSz" << endl;
  return this->HilbertSpaceDimension;
}

#endif
