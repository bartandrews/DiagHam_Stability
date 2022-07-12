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


#ifndef BOSONONSPHEREWITHSPIN_H
#define BOSONONSPHEREWITHSPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

using std::cout;
using std::endl;
using std::hex;
using std::dec;

class BosonOnSphere;
class BosonOnSphereShort;

class BosonOnSphereWithSpin :  public ParticleOnSphereWithSpin
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
    // number of bosons with spin up / down
  int NbrBosonsUp;
  int NbrBosonsDown;
  // twice the total spin value
  int TotalSpin;

  

  // array describing each state (full storage, temporary use)
  unsigned** StateDescription;
  // array describing each state up / down
  unsigned long* StateDescriptionUp;
  unsigned long* StateDescriptionDown;
  // array giving maximum Lz value reached for a boson in a given state for up and down spin
  unsigned* StateLzMaxUp;
  unsigned* StateLzMaxDown;
  // array to store condensed info on states after state generation
  unsigned *StateInfo;

  // Lookup-table for base indices for given sequence of up /down spins
  unsigned long *LookUpTableUp;
  unsigned long *LookUpTableDown;


  // table for coherence factors
  double *CoherenceFactors;
    
  // pointer to an integer which indicate which coordinates are kept for the next time step iteration
  int* KeptCoordinates;
  // minors of permanents used for the time coherent wave function evaluation
  Complex** Minors;

  // temporary state used when applying operators
  unsigned* TemporaryState;  
  // temporary state used when applying operators for type up particles
  unsigned long* TemporaryStateUp;
  // temporary state used when applying operators for type down particles
  unsigned long* TemporaryStateDown;
  // temporary state used when applying ProdA operator
  unsigned* ProdATemporaryState;
  // number of up spins in temporary state
  unsigned ProdATemporaryStateNbrUp;
  // temporary storage for monomial decompositions - should convert to unsigned at some point
  unsigned long* TemporaryMonomials;

    // target space for operations leaving the Hilbert-space
  BosonOnSphereWithSpin *TargetSpace;

  
 public:

  // default constructor
  //
  BosonOnSphereWithSpin ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // totalSpin = twice the total spin
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSpin (int nbrBosons, int totalLz, int lzMax, int totalSpin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSpin(const BosonOnSphereWithSpin& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSpin ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSpin& operator = (const BosonOnSphereWithSpin& bosons);

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

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector,  int szSector, RealVector& groundState);

  // core part of the evaluation density matrix calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  long EvaluatePartialDensityMatrixCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
					 RealVector& groundState, RealSymmetricMatrix* densityMatrix);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrNUpSector = number of spin up  that belong to the subsytem 
  // nbrNDownSector = number of spin down  that belong to the subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
								     int nbrNUpSector, int nbrNDownSector, RealVector& groundState, AbstractArchitecture* architecture);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								  RealVector& groundState, RealSymmetricMatrix* densityMatrix);
  
  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  RealVector ForgeSU2FromU1(RealVector& upState, BosonOnSphere& upStateSpace, RealVector& downState, BosonOnSphere& downStateSpace);

  // Project the state from the su2 space
  // to the U(1) space (u1Space)
  //
  // state = state that needs to be projected
  // u1Space = the subspace onto which the projection is carried out
  virtual RealVector ForgeU1FromSU2(RealVector& state, BosonOnSphere& u1Space);

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
  int FindStateIndex(unsigned* stateDescription);

  // get Lz component of a component
  //
  // j = index of the component in Hilbert space
  // return value = twice the Lz component
  virtual inline int GetLzValue(int j=0);

  //get total Lz of up spins in a state from lookup
  //index = index of state in Hilbert space
  // return value = twice the total lz of up spins in state
  inline int GetTotalLzUp(long index);

  //get total Lz of down spins in a state from lookup
  //index = index of state in Hilbert space
  //return value = twice the total lz of down spins in a state
  inline int GetTotalLzDown(long index);


  //get total Lz of up spins in a state
  //index = index of state in Hilbert space
  // return value = twice the total lz of up spins in state
  inline int CalculateTotalLzUp(long index);

  //get total Lz of down spins in a state
  //index = index of state in Hilbert space
  //return value = twice the total lz of down spins in a state
  inline int CalculateTotalLzDown(long index);
  
  // Compute the product of two states that belong to different Hilbert Spaces
  //
  // firstState = reference on one of the states whose product will be computed
  // polarizedState = reference on the other state whose product will be computed
  // OutputVector = reference on the vector where the result will be stored
  // PolarizedSpace = pointer on the Hilbert Space whose the polarized state belong
  // minIndex = first computed component
  // nbrComponents = Nomber of computed components
  // FinalSpace = pointer on the Hilbert Space whose the final state belong
  void BosonicStateWithSpinTimesBosonicState(RealVector& spinfulState, RealVector& polarizedState, RealVector& outputVector,BosonOnSphereShort * polarizedSpace,int minIndex,int nbrComponents, BosonOnSphereWithSpin * finalSpace);


  //Calculate state and coefficient when symmetrising over two groups of particles
  //PlusStateUp = reference on array where plus state up spin occupations stored
  //MinusStateUp = reference on array where minus state up spin occupations stored
  //PlusStateDown = reference on array where plus state down spin occupations stored
  //MinusStateDown = reference on array where minus state down spin occupations stored
  //coefficient = coefficient of term before symmetrisation
  //OutputVector = vector where resulting term will be added
  void SymmetriseOverGroupsAndAddToVector(unsigned long * & PlusStateUp, unsigned long * & MinusStateUp, unsigned long * & PlusStateDown, unsigned long * &MinusStateDown, double coefficient, RealVector & OutputVector);

  //Compute the geometric correction factor for a given product state when multiplying two monomials and working with second quantised forms on the sphere for fully polarized states
  //
  //firstState = reference on array where monomial representation of first state stored
  //secondState = reference on array where monomial representation of second state stored
  //productState = reference on array where monomial representation of a given final state in the product is stored
  //lzMaxOne = twice maximum lz value for a boson in first state
  //lzMaxTwo = twice maximum lz value for a boson in second state
  double GeometricCorrectionFactor(unsigned long * firstMonomial, unsigned long * secondMonomial, unsigned long * productMonomial, int lzMaxOne, int lzMaxTwo);

  //Compute the geometric correction factor for a given product state when multiplying two monomials and working with second quantised forms on the sphere for for spinful states
  //
  //firstState = reference on array where monomial representation of first state stored
  //secondState = reference on array where monomial representation of second state stored
  //productState = reference on array where monomial representation of a given final state in the product is stored
  //lzMaxOne = twice maximum lz value for a boson in first state
  //lzMaxTwo = twice maximum lz value for a boson in second state
  double GeometricCorrectionFactor(unsigned long * firstMonomial, unsigned long * secondMonomial, unsigned long * productMonomialUp, unsigned long * productMonomialDown, int lzMaxOne, int lzMaxTwo);

  //Compute the occupation correction factor for a given product state when multiplying two monomials and working with second quantised forms on the sphere for spinful states
  //
  //firstState = reference on array where bosonic representation of first state stored
  //secondState = reference on array where bosonic representation of second state stored
  //productState = reference on array where bosonic representation of a given final state in the product is stored
  //lzMaxOne = twice maximum lz value for a boson in first state
  //lzMaxTwo = twice maximum lz value for a boson in second state
  double OccupationCorrectionFactor(unsigned long * firstState, unsigned long * secondState, unsigned long * productState, int lzMaxOne, int lzMaxTwo);

  //Compute the occupation correction factor for a given product state when multiplying two monomials and working with second quantised forms on the sphere for spinful states
  //
  //firstState = reference on array where bosonic representation of first state stored
  //secondState = reference on array where bosonic representation of second state stored
  //productState = reference on array where bosonic representation of a given final state in the product is stored
  //lzMaxOne = twice maximum lz value for a boson in first state
  //lzMaxTwo = twice maximum lz value for a boson in second state
  double OccupationCorrectionFactor(unsigned long * firstPolarizedState, unsigned long * secondUpState, unsigned long * secondDownState, unsigned long * productUpState, unsigned long * productDownState, int lzMaxOne, int lzMaxTwo);


  // convert a fermionic state to its monomial representation
  //
  // index = index of the fermionic state
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void GetMonomial(long index, unsigned long*& finalState);

  //get the bosonic state description of a state
  //index = index of the fermionic state
  //stateDescriptionUp = reference on array where up description will be stored
  //stateDescriptionDown = reference on array where down description will be stored
  void GetBosonicDescription(int index, unsigned long* & stateDescriptionUp, unsigned long* & stateDescriptionDown);

  //convert a vector in the monomial basis to the Fock basis
  //
  //StateInMonomialBasis = state vector components in the basis of symmetric monomials
  //StateInFockBasis = reference on the state vector components in the Fock basis where result is stored
  void MonomialToFockBasis( RealVector & StateInMonomialBasis, RealVector & StateInFockBasis);

  //convert a vector in the Fock basis to the monomial basis
  //
  //StateInFockBasis = state vector componenets in the Fock basis
  //StateInMonomialBasis = reference on the state vector components in the Monomial basis where result is stored
  void FockToMonomialBasis( RealVector & StateInFockBasis, RealVector & StateInMonomialBasis );

  //get the conversion factor to go from a symmetric monomial to the Fock basis
  //
  //index = index of state
  double MonomialToFockConversionFactor( long index );

  //convert a polarized monomial to occupation basis.
  //PolarizedMonomial input monomial
  //PolarizedMonomialOccupationBasis reference on array to store result
  //maxLz twice the maximum angular momentum of the monomial
  void ConvertPolarizedMonomialToOccupationBasis( unsigned long *& PolarizedMonomial, unsigned long *& PolarizedMonomialOccupationBasis, int maxLz );

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component as to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);    

  // convert a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // normalize Jack with respect to cylinder basis
  //
  // state = reference to the Jack state to normalize
  // aspect = aspect ratio of cylinder
  // return value = normalized state
  virtual RealVector& NormalizeJackToCylinder(RealVector& state, double aspect);

  //find state index
  //
  //stateDescriptionUp = fermionic description of the up spins
  //stateDescriptionDown = fermionic description of the down spins
  int FindStateIndex(unsigned long stateDescriptionUp, unsigned long StateDescriptionDown);

  //find state index
  //
  //stateDescriptionUp = reference on array describing up spin occupations
  //stateDescriptionDown = reference on array describing down spin occupations
  int FindStateIndex(unsigned long * stateDescriptionUp, unsigned long * stateDescriptionDown);

 protected:

  //convert a polarized monomial to occupation basis.
  //PolarizedMonomial input monomial
  //PolarizedMonomialOccupationBasis reference on array to store result
  //maxLz twice the maximum angular momentum of the monomial
  void PolarizedMonomialToOccupationBasis( unsigned long *& PolarizedMonomial, unsigned long *& PolarizedMonomialOccupationBasis, int maxLz );

  // find state index
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  //int FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);

  
  // find index of a tensored configuration
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  unsigned FindTensoredIndex(unsigned long stateDescription, int lzmax);

  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrBosons = number of bosons
  // lzMax = two times momentum maximum value for a boson plus one 
  // totalLz = momentum total value plus nbrBosons * (momentum maximum value for a boson + 1)
  // return value = Hilbert space dimension
  long int ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalSpin);


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
  virtual long int EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalSpin);


  // generate all states corresponding to the constraints
  // 
  // nbrBosonsUp = number of bosons with spin up
  // nbrBosonsDown = number of bosons with spin down
  // lzMax = momentum maximum value for a boson in the state
  // totalLz = momentum total value
  // return value = position from which new states have to be stored
  
  long GenerateStates(int nbrBosonsUp, int nbrBosonsDown, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrBosonsUp = number of bosons with spin up
  // nbrBosonsDown = number of bosons with spin down
  // lzMax = momentum maximum value for a boson in the state
  // lzMaxUp = momentum maximum value for a spin-up boson in the state
  // currentLzMax = momentum maximum value for bosons that are still to be placed
  // totalLz = momentum total value
  // totalLzUp = momentum total value for spin-up bosons
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  
  long GenerateStatesWithConstraint(int nbrBosonsUp, int nbrBosonsDown, int lzMax, int lzMaxUp, int currentLzMax, int totalLzUp, int totalLzDown, long pos, int level);

  // convert a bosonic state into its fermionic counterpart
  //
  // finalStateUp = return value of bit-coded state of up-bosons
  // finalStateDown = return value of bit-coded state of down-bosons
  // finalLzMaxUp = highest bit in fermionic coding finalStateUp
  // finalLzMaxDown = highest bit in fermionic coding finalStateDown
  // initialState = reference on the array where initialbosonic  state is stored
  
  void BosonToFermion(unsigned long &finalStateUp, unsigned long &finalStateDown, int &finalLzMaxUp, int &finalLzMaxDown, unsigned*& initialState);


  // convert a fermionic state into its bosonic counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax = initial fermionic state maximum Lz value
  // finalState = reference on the array where the bosonic state has to be stored
  // finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored
  void FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown, unsigned initialInfo, unsigned*& finalState, int &finalStateLzMaxUp, int &finalStateLzMaxDown);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialStateUp = reference on the array where initial bosonic state for the type up particles is stored
  // initialStateDown = reference on the array where initial bosonic state for the type down particles is stored
  // finalStateUp = reference on the corresponding fermionic state for the type up particles
  // finalStateDown = reference on the corresponding fermionic state for the type down particles
  virtual void BosonToFermion(unsigned long*& initialStateUp, unsigned long*& initialStateDown, 
			      unsigned long& finalStateUp, unsigned long& finalStateDown);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialStateUp = initial fermionic state for the type up particles
  // initialStateDown = initial fermionic state for the type down particles
  // finalStateUp = reference on the array where the bosonic state for the type up particles has to be stored
  // finalStateDown = reference on the array where the bosonic state for the type down particles has to be stored
  virtual void FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown, 
			      unsigned long*& finalStateUp, unsigned long*& finalStateDown);

  // convert a bosonic state to its monomial representation
  //
  // initialStateUp = initial spin up bosonic state in its fermionic representation 
  // initialStateDown = initial spin down bosonic state in its fermionic representation
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored  
  void ConvertToMonomial(unsigned long initialStateUp, unsigned long initialStateDown, int initialStateLzMax, unsigned long*& finalState);

  // convert a bosonic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // return value = bosonic state in its fermionic representation  
  unsigned long ConvertFromMonomial(unsigned long* initialState);

  // get LzMax value for a given state
  // index = index of state to analyse
  // return = lzMax value (max of up and down)
  int GetStateLzMax(int index);

  //get LzMax value for a given state up spin
  //index = index of state to analyse
  //return = lzMax value (up)
  int GetStateLzMaxUp(long index);

  //getLzMaxValue for a given state down spin
  //index = index of state to analyse
  //return = lzMax value (down)
  int GetStateLzMaxDown(long index);

  // sort an array and reflect permutations in auxiliary array
  //
  // length = length of arrays
  // sortArray = array to be sorted
  // auxArray = auxiliary array
  //
  void ShellSortAux(unsigned length, unsigned long* sortArray, unsigned *auxArray, int *auxArray2);

  
};


// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the Lz component
int BosonOnSphereWithSpin::GetLzValue(int j)
{
  return this->TotalLz;
}


//get total Lz of up spins in a state using lookup
//index = index of state in Hilbert space
//return value = twice the total lz of up spins in state
inline int BosonOnSphereWithSpin::GetTotalLzUp(long index)
{
  return 2*(this->StateInfo[index]>>20)-this->LzMax*NbrBosonsUp;
}

//get total Lz of up spins in a state using lookup
//index = index of state in Hilbert space
//return value = twice the total lz of up spins in state
inline int BosonOnSphereWithSpin::GetTotalLzDown(long index)
{
  return this->TotalLz - 2*(this->StateInfo[index]>>20)+this->LzMax*NbrBosonsUp;
}


//calculate total Lz of up spins in a state
//index = index of state in Hilbert space
//return value = twice the total lz of up spins in state
inline int BosonOnSphereWithSpin::CalculateTotalLzUp(long index)
{
  this->GetMonomial(index, TemporaryMonomials);
  int totalLzUp=0;
  for(int i=0; i<this->NbrBosonsUp; i++) 
    {
      totalLzUp += (2*TemporaryMonomials[i] - this->LzMax);
    }
  return totalLzUp;
}

//calculate total Lz of down spins in a state
//index = index of state in Hilbert space
//return value = twice the total lz of down spins in a state
inline int BosonOnSphereWithSpin::CalculateTotalLzDown(long index)
{
  this->GetMonomial(index, TemporaryMonomials);
  int totalLzDown=0;
  for(int i=this->NbrBosonsUp; i<this->NbrBosons; i++)
    {
      totalLzDown += (2*TemporaryMonomials[i] - this->LzMax);
    }
  return totalLzDown;
}

//convert a polarized monomial to occupation basis.
//PolarizedMonomial input monomial
//PolarizedMonomialOccupationBasis reference on array to store result
//nbrbosons = number of bosons in the monomial
inline void BosonOnSphereWithSpin::ConvertPolarizedMonomialToOccupationBasis( unsigned long *& PolarizedMonomial, unsigned long *& PolarizedMonomialOccupationBasis, int nbrBosons )
{
  for(int i=0; i<nbrBosons; i++)
    {
      PolarizedMonomialOccupationBasis[ PolarizedMonomial[i] ] ++;
    }
}



// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereWithSpin::GetParticleStatistic()
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

inline void BosonOnSphereWithSpin::BosonToFermion(unsigned long &finalStateUp, unsigned long &finalStateDown, int &finalLzMaxUp, int &finalLzMaxDown, unsigned*& initialState)
{
/*   cout << "NbrUp="<<StateNbrUp <<", InitialState = |"; */
//   for (int i=0; i<=LzMax; ++i) 
//     cout << " " << (initialState[i] >> 16)<< "u "<< (initialState[i] & 0xffff)<<"d |"; 
//   cout << endl; 

  finalStateUp = 0x0ul;
  unsigned ShiftUp = 0;
  finalStateDown = 0x0ul;
  unsigned ShiftDown = 0;
  finalLzMaxUp = 0;
  finalLzMaxDown = 0;
  int RemainingNbrUp=NbrBosonsUp;
  int RemainingNbrDown=NbrBosonsDown;
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
  finalLzMaxUp += this->NbrBosonsUp-(this->NbrBosonsUp!=0);
  finalLzMaxDown += this->NbrBosonsDown -(this->NbrBosonsDown!=0);
  //  this->PrintState(cout, initialState) << " " << std::hex << finalStateUp << " " << finalStateDown << std::dec << " " << finalLzMaxUp <<" " << finalLzMaxDown<<endl;
  return;
}


// convert a fermionic state into its bosonic counterpart
//
// initialState = initial fermionic state
// initialStateLzMax = initial fermionic state maximum Lz value
// finalState = reference on the array where the bosonic state has to be stored
// finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored

inline void BosonOnSphereWithSpin::FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown, unsigned initialInfo, unsigned*& finalState, int &finalStateLzMaxUp, int &finalStateLzMaxDown)
{
  //cout << "FermionToBoson :: initialStateUp =" << hex << initialStateUp << ", initialStateDown="<<initialStateDown<<dec<<endl;
  int StateNbrUp=0;
  int InitialStateLzMax = (initialInfo >> 10)&0x3ffu;
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
      StateNbrUp+=TmpPower;
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
  InitialStateLzMax = initialInfo&0x3ff;
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


// convert a bosonic state into its fermionic counterpart
//
// initialStateUp = reference on the array where initial bosonic state for the type up particles is stored
// initialStateDown = reference on the array where initial bosonic state for the type down particles is stored
// finalStateUp = reference on the corresponding fermionic state for the type up particles
// finalStateDown = reference on the corresponding fermionic state for the type down particles

inline void BosonOnSphereWithSpin::BosonToFermion(unsigned long*& initialStateUp, unsigned long*& initialStateDown,
						  unsigned long& finalStateUp, unsigned long& finalStateDown)
{
  finalStateUp = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateUp |= ((1ul << initialStateUp[i]) - 1ul) << Shift;
      Shift += initialStateUp[i];
      ++Shift;
    }
  finalStateDown = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateDown |= ((1ul << initialStateDown[i]) - 1ul) << Shift;
      Shift += initialStateDown[i];
      ++Shift;
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialStateUp = initial fermionic state for the type up particles
// initialStateDown = initial fermionic state for the type down particles
// finalStateUp = reference on the array where the bosonic state for the type up particles has to be stored
// finalStateDown = reference on the array where the bosonic state for the type down particles has to be stored

inline void BosonOnSphereWithSpin::FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown,
						  unsigned long*& finalStateUp, unsigned long*& finalStateDown)
{
  int FinalStateLzMax = 0;
  int InitialStateLzMax = this->LzMax + NbrBosonsUp - 1;
  while ((InitialStateLzMax >= 0) && ((initialStateUp >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateUp - 1ul) ^ (~initialStateUp);
      TmpState &= ~(TmpState >> 1);
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
      finalStateUp[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateUp >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateUp[FinalStateLzMax] = 0x0ul;

  FinalStateLzMax = 0;
  InitialStateLzMax = this->LzMax + NbrBosonsDown - 1;
  while ((InitialStateLzMax >= 0) && ((initialStateDown >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateDown - 1ul) ^ (~initialStateDown);
      TmpState &= ~(TmpState >> 1);
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
      finalStateDown[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateDown >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateDown[FinalStateLzMax] = 0x0ul;
}

// convert a fermionic state to its monomial representation
//
// index = index of the fermionic state
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnSphereWithSpin::GetMonomial(long index, unsigned long*& finalState)
{
  int Index = 0;
  unsigned long InitialStateUp = this->StateDescriptionUp[index];
  unsigned long InitialStateDown = this->StateDescriptionDown[index];
  int InitialStateLzMax  = this->LzMax + this->NbrBosonsUp -1;
  int TmpLz = this->LzMax;
  while (InitialStateLzMax >= 0)
    {
      while ((InitialStateLzMax >= 0) && (((InitialStateUp >> InitialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalState[Index++] = (unsigned long) TmpLz;
	  --InitialStateLzMax;
	}
      while ((InitialStateLzMax >= 0) && (((InitialStateUp >> InitialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --InitialStateLzMax;
	}
    }
  InitialStateLzMax = this->LzMax + this->NbrBosonsDown - 1;
  TmpLz = this->LzMax;
  while (InitialStateLzMax >= 0)
    {
      while ((InitialStateLzMax >= 0) && (((InitialStateDown >> InitialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalState[Index++] = (unsigned long) TmpLz;
	  --InitialStateLzMax;
	}
      while ((InitialStateLzMax >= 0) && (((InitialStateDown >> InitialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --InitialStateLzMax;
	}
    }
  
}

// convert a bosonic state to its monomial representation
//
// initialStateUp = initial spin up bosonic state in its fermionic representation 
// initialStateDown = initial spin down bosonic state in its fermionic representation
// initialStateLzMax = initial bosonic state maximum Lz value
// finalState = reference on the array where the monomial representation has to be stored

//let's try with just the one initialStateLzMax

inline void BosonOnSphereWithSpin::ConvertToMonomial(unsigned long initialStateUp, unsigned long initialStateDown, int initialStateBosonLzMax, unsigned long*& finalState) 
{
  int initialStateLzMax = initialStateBosonLzMax + this->NbrBosonsUp -1;
  int Index = 0;
  int TmpLz = initialStateLzMax - this->NbrBosonsUp + 1;
  while (initialStateLzMax >= 0)
    {
      while ((initialStateLzMax >= 0) && (((initialStateUp >> initialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalState[Index++] = TmpLz;
	  --initialStateLzMax;
	}
      while ((initialStateLzMax >= 0) && (((initialStateUp >> initialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --initialStateLzMax;
	}
    }
  initialStateLzMax = initialStateBosonLzMax + this->NbrBosonsDown - 1;
  TmpLz = initialStateLzMax - this->NbrBosonsDown + 1;
  while (initialStateLzMax >= 0)
    {
      while ((initialStateLzMax >= 0) && (((initialStateDown >> initialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalState[Index++] = TmpLz;
	  --initialStateLzMax;
	}
      while ((initialStateLzMax >= 0) && (((initialStateDown >> initialStateLzMax) & 0x1ul) == 0x0ul))
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

inline unsigned long BosonOnSphereWithSpin::ConvertFromMonomial(unsigned long* initialState)
{
  unsigned long Tmp = 0x0ul;
  for (int i = 0; i < this->NbrBosonsUp; ++i)
    Tmp |= 0x1ul << (initialState[i] + ((unsigned long) (this->NbrBosons - i)) - 1ul);
  for (int i = this->NbrBosonsUp; i < this->NbrBosons; ++i)
    Tmp |= 0x1ul << (initialState[i] + ((unsigned long) (this->NbrBosons - i)) - 1ul);
  return Tmp;
}

// get LzMax value for a given state
// index = index of state to analyse
// return = lzMax value (max of up and down)
inline int BosonOnSphereWithSpin::GetStateLzMax(int index)
{
  unsigned Info = StateInfo[index];
  int LzMaxUp = (Info >> 10)&0x3ffu;
  int LzMaxDown = Info&0x3ff;
  LzMaxUp -= NbrBosonsUp-(NbrBosonsUp!=0);
  LzMaxDown -= this->NbrBosonsDown - (NbrBosonsDown!=0);
  return std::max(LzMaxUp, LzMaxDown);
}

//get LzMax value for a given state up spin
//index = index of state to analyse
//return = lzMax value of up spins
inline int BosonOnSphereWithSpin::GetStateLzMaxUp(long index)
{
  /*
  unsigned Info = StateInfo[index];
  int LzMaxUp = (Info >> 10)&0x3ffu;
  LzMaxUp -= NbrBosonsUp-(NbrBosonsUp!=0);
  return LzMaxUp;  
  */
  this->GetMonomial(index, TemporaryMonomials);
  unsigned long lzmaxup = 0;
  if(this->NbrBosonsUp != 0)
    {
      for(int i = 0; i<this->NbrBosonsUp; i++)
	{
	  if(TemporaryMonomials[i] > lzmaxup)
	    lzmaxup = TemporaryMonomials[i];
	}
    }
  else
    {
      cout << "Error NbrBosons == 0\n";
      lzmaxup = -1;
    }
  /*
  int myLzMax = (((StateInfo[index])>>10)&0x3ffu) - this->NbrBosonsUp;
  if (NbrBosonsUp!=0)
    ++myLzMax;
  cout << "myLzMax = "<<myLzMax<<" lzmaxup="<<lzmaxup<<" state: ";
  PrintState(cout,index)<<endl;
  */
      
  return (int) lzmaxup;
}

//get LzMax value for a given state up spin
//index = index of state to analyse
//return = lzMax value of down spins
inline int BosonOnSphereWithSpin::GetStateLzMaxDown(long index)
{
  /*
  unsigned Info = StateInfo[index];
  int LzMaxDown = (Info & 0x3ffu);
  LzMaxDown -= NbrBosonsDown-(NbrBosonsDown!=0);
  return LzMaxDown;
  */
  this->GetMonomial(index, TemporaryMonomials);
  unsigned long lzmaxdown = 0;
  if(this->NbrBosonsDown != 0)
    {
      for(int i = this->NbrBosonsUp; i<this->NbrBosons; i++)
	{
	  if(TemporaryMonomials[i] > lzmaxdown)
	    lzmaxdown =  TemporaryMonomials[i];
	}
    }
  else
    {
      cout << "Error lzmaxdown\n";
      lzmaxdown = -1;

    }

  return (int) lzmaxdown;

}


// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

inline int BosonOnSphereWithSpin::FindStateIndex(unsigned* stateDescription)
{
  unsigned long StateDescriptionUp;
  unsigned long StateDescriptionDown;
  int LzMaxUp, LzMaxDown;
  this->BosonToFermion(StateDescriptionUp, StateDescriptionDown, LzMaxUp, LzMaxDown, stateDescription);
  //  cout << "up: "<< hex << FinalStateUp << dec << " deduced lzmax="<<LzMaxUp<<endl;
  //  cout << "down: "<< hex << FinalStateDown << dec << " deduced lzmax="<<LzMaxDown<<endl;
  return this->LookUpTableUp[StateDescriptionUp]+this->LookUpTableDown[StateDescriptionDown];
}


// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index
inline int BosonOnSphereWithSpin::FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  // cout << "FindStateIndex: "<<hex<<stateDescriptionUp<<" "<<stateDescriptionDown<<" " <<dec << lzMaxUp <<" "<<lzMaxDown;
  int Base = this->LookUpTableUp[stateDescriptionUp];
  int Offset = this->LookUpTableDown[stateDescriptionDown];
  return Base+Offset;
}

inline long factorial(int n) {
  if(n==0 || n==1)
    return 1l;
  else
    return n*factorial(n-1l);
}

inline int BosonOnSphereWithSpin::FindStateIndex(unsigned long * stateDescriptionUp, unsigned long * stateDescriptionDown)
{
  for(int i=0; i<this->NbrLzValue; i++)
    {
      TemporaryState[i] = ((unsigned) stateDescriptionDown[i]) | (((unsigned) stateDescriptionUp[i])<<16);
    }
  unsigned long fermionUp, fermionDown;
  int fermionLzMaxUp, fermionLzMaxDown;
  this->BosonToFermion(fermionUp, fermionDown, fermionLzMaxUp, fermionLzMaxDown, TemporaryState);
  //cout << "Fermion up description " << fermionUp << " fermion down description  " << fermionDown << "\n";
  int index = this->FindStateIndex( fermionUp, fermionDown );
  
  //int index = this->FindStateIndex( TemporaryState );

  //  cout << "Index found to be " << index << " with state description ";
  //  this->PrintState(cout, index);
  //  cout << "\n";
  return index;
}


//convert a vector in the monomial basis to the Fock basis
//
//StateInMonomialBasis = state vector components in the basis of symmetric monomials
//StateInFockBasis = reference on the state vector components in the Fock basis where result is stored

inline void BosonOnSphereWithSpin::MonomialToFockBasis( RealVector & StateInMonomialBasis, RealVector & StateInFockBasis )
{
  for(long i=0x0l; i<this->GetHilbertSpaceDimension(); i++)
    {
      (StateInFockBasis)[i] = StateInMonomialBasis[i] * this->MonomialToFockConversionFactor(i);
    }
}

//convert a vector in the Fock basis to the monomial basis
//
//StateInFockBasis = state vector components in the Fock basis
//StateInMonomialBasis = reference on the state vector components in the basis of symmetric monomials where result is stored

inline void BosonOnSphereWithSpin::FockToMonomialBasis( RealVector & StateInFockBasis, RealVector & StateInMonomialBasis )
{
  for(long i=0x0l; i<this->GetHilbertSpaceDimension(); i++)
    {
      (StateInMonomialBasis)[i] = StateInFockBasis[i] / this->MonomialToFockConversionFactor(i);
    }
}

//get the conversion factor to go from a symmetric monomial to the Fock basis
//
//index = index of state

inline double BosonOnSphereWithSpin::MonomialToFockConversionFactor( long index )
{
  long SphericalGeometricFactorSquared = 1l;
  long OccupationFactorSquared = 1l;
  this->GetMonomial(index, TemporaryMonomials);
  unsigned long * OccupationBasisUp = new unsigned long[this->LzMax+1];
  unsigned long * OccupationBasisDown = new unsigned long[this->LzMax+1];
  for(int i=0; i < this->LzMax+1; i++)
    {
      OccupationBasisUp[i] = 0x0ul;
      OccupationBasisDown[i] = 0x0ul;
    }

  int N_phi = this->LzMax;
  for( int i=0; i<this->NbrBosonsUp; i++ )
    {
      ++OccupationBasisUp[ TemporaryMonomials[i] ];
      SphericalGeometricFactorSquared *= (factorial( (int)TemporaryMonomials[i] ))*factorial(N_phi - (int)TemporaryMonomials[i]);
      //cout << "UpBoson " << i << " 2lz " << TemporaryMonomials[i] << " N_phi - 2lz " << N_phi - TemporaryMonomials[i] << "\n";
    }
  for( int i=this->NbrBosonsUp; i<this->NbrBosons; i++ )
    {
      OccupationBasisDown[ TemporaryMonomials[i] ]++;
      SphericalGeometricFactorSquared *= (factorial( (int)TemporaryMonomials[i] ))*factorial(N_phi - (int)TemporaryMonomials[i]);
      //cout << "DownBoson " << i << " 2lz " << TemporaryMonomials[i] << " N_phi - 2lz " << N_phi - TemporaryMonomials[i] << "\n";
    }
  for( int i=0; i <= this->LzMax; i++)
    {
      OccupationFactorSquared *= factorial( OccupationBasisUp[i] )*factorial( OccupationBasisDown[i] );
      //cout << "OccupationBasisUp " << OccupationBasisUp[i] << " OccupationBasisDown " << OccupationBasisDown[i] << " for lz " << i << "\n";
    }

  double conversionFactor = sqrt( (double)factorial(this->NbrBosons) * (double)SphericalGeometricFactorSquared / (double) OccupationFactorSquared );
  cout << "conversionFactor " << conversionFactor << "\n";

  delete [] OccupationBasisUp;
  delete [] OccupationBasisDown;

  return conversionFactor;
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(2) index for the first annihilation operator
// sigma2 = SU(2) index for the second annihilation operator
// return value =  multiplicative factor 

inline double BosonOnSphereWithSpin::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  cout << "warning : AsigmaAsigma not defined in BosonOnSphereWithSpin" << endl;
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

inline int BosonOnSphereWithSpin::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  cout << "warning : AdsigmaAdsigma not defined in BosonOnSphereWithSpin" << endl;
  return this->HilbertSpaceDimension;
}


#endif
