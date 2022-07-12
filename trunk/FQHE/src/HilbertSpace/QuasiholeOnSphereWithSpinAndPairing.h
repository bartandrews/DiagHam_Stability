////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                        Class author Cecile Repellin                        //
//                                                                            //
//                                                                            //
//              class of quasiholes on sphere with spin and pairing           //
//      (i.e. Sz conservation but no conservation of the particle number)     //
//                                                                            //
//                        last modification : 05/05/2016                      //
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


#ifndef QUASIHOLEONSPHEREWITHSPINANDPAIRING_H
#define QUASIHOLEONSPHEREWITHSPINANDPAIRING_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Matrix/SparseRealMatrix.h"

#include <iostream>


class QuasiholeOnSphereWithSpinAndPairing :  public ParticleOnSphereWithSpin
{


 protected:

  // number of fermions
  int NbrFermions;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a fermion
  int LzMax;
  // number of fermions with spin up / down
  int NbrFermionsUp;
  int NbrFermionsDown;
  // number of Lz values in a stat
  int NbrLzValue;
  // twice the total spin value
  int TotalSpin;
  

   // first index of the (k, r) exclusion principle in one layer
   int KValue;
   // second index of the (k, r) exclusion principle in one layer
   int RValue;
   // fermion factor (0 if boson, 1 if fermion)
   int FermionFactor;
   // number of admissible (NbrParticles, Lz) in one layer
   int NbrQuasiholeEntriesSingleLayer;
   // array containing the dimension of each quasihole subspace for a single layer
   int* NbrQuasiholesPerNPerLzSingleLayer;
   // array containing the number of fermions in the up layer for each Hilbert space index
   int* NbrFermionUpFullSpace;
   // array containing the value of Lz in the up layer for each Hilbert space index
   int* LzValueUpFullSpace;
   // array containing the value of the lowest Hilbert space index with given quantum numbers N and Lz
   int* FirstIndexWithNbrParticlesUpLzValueUp;
   // maximal number of fermions in the upper layer
   int NbrFermionsUpMax;
   // minimal number of fermions in the upper layer
   int NbrFermionsUpMin;
   // maximal number of fermions in the lower layer
   int NbrFermionsDownMax;
   // minimal number of fermions in the lower layer
   int NbrFermionsDownMin;
   // maximal number of states that any state can be coupled to
   int MaximalNumberCouplingElements;

   // maximal total Lz value that can be reached in a single layer for a given number of particles
   int* MaximalLzSingleLayer;
   // minimal amplitude to consider for the Jack coupling elements
   double Error;

   
   // array where the all the a^+a matrix elements are stored
   RealMatrix*** SingleLayerAdAMatrices;
   // array where all the auad matrix elements are stored
   RealMatrix*** SingleLayerAnnihilationMatrices;

   
    // linearized indices for a single layer
    int** SingleLayerLinearIndices;

 public:

  // default constructor
  //
  QuasiholeOnSphereWithSpinAndPairing();

  // basic constructor
  // 
  // kExclusionPrinciple = k value of the exclusion principle
  // rExclusionPrinciple = r value of the exclusion principle
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twice the total spin value
  // directory = optional path to data files
  // filePrefix = prefix for all input file (should include everything related to the statistics and the geometry)
  // discardPairing = if true, do not load the matrix element required to compute the pairing term
  QuasiholeOnSphereWithSpinAndPairing (int kExclusionPrinciple, int rExclusionPrinciple, int totalLz, int lzMax, int totalSpin, 
				       const char* directory, const char* filePrefix, bool discardPairing = false);

  // copy constructor (without duplicating datas)
  //
  // quasiholes = reference on the hilbert space to copy to copy
  QuasiholeOnSphereWithSpinAndPairing(const QuasiholeOnSphereWithSpinAndPairing& quasiholes);

  // destructor
  //
  ~QuasiholeOnSphereWithSpinAndPairing ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  QuasiholeOnSphereWithSpinAndPairing& operator = (const QuasiholeOnSphereWithSpinAndPairing& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();
  
  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);
  
  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // apply a_u_m a_d_m to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for destruction operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AuAd (int index, int m, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_u_m a^\dagger_d_m to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for creation operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AduAdd (int index, int m, int*& leftIndices, double*& interactionElements);
  
  // apply a_u_m a_d_n to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for the left annihilation operator
  // n = index for the right annihilation operator
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AuAd (int index, int m, int n, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_u_m a^\dagger_d_n to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for the left annihilation operator
  // n = index for the right annihilation operator
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AduAdd (int index, int m, int n, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_u_m a_u_m to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for destruction operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AduAu (int index, int m, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_d_m a_d_m to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for destruction operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AddAd (int index, int m, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_u_m a_u_n to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for creation  operators
  // n = index for annihilation operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AduAu (int index, int m, int n, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_d_m a_d_n to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for creation  operators
  // n = index for annihilation operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AddAd (int index, int m, int n, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_u_m to a given state defined only in the up layer
  //
  // m = index for destruction operators
  // inputState = reference of the state to act on
  // outputState  = reference of the state where the result will be stored
  // nbrParticlesUp = number of particles for the up layer input state
  // lzUp = momentum of the up layer input state
  virtual void Adu (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp);
  
  // apply a^\dagger_d_m to a given state defined only in the up layer
  //
  // m = index for destruction operators
  // inputState = reference of the state to act on
  // outputState  = reference of the state where the result will be stored
  // nbrParticlesDown = number of particles for the down layer input state
  // lzDown = momentum of the down layer input state  
  virtual void Add (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp);
  
  // apply a_u_m to a given state defined only in the up layer
  //
  // m = index for destruction operators
  // inputState = reference of the state to act on
  // outputState  = reference of the state where the result will be stored
  // nbrParticlesUp = number of particles for the up layer input state
  // lzUp = momentum of the up layer input state  
  virtual void Au (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp);
  
  // apply a_d_m to a given state defined only in the up layer
  //
  // m = index for destruction operators
  // inputState = reference of the state to act on
  // outputState  = reference of the state where the result will be stored
  // nbrParticlesDown = number of particles for the down layer input state
  // lzDown = momentum of the down layer input state
  virtual void Ad (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp);

  // apply a^\dagger_u_m a_u_m to a given state defined only in the up layer
  //
  // m = index for destruction operators
  // inputState = reference of the state to act on
  // outputState  = reference of the state where the result will be stored
  // nbrParticlesUp = number of particles in the up layer
  // lzUp = momentum of the up layer eigenstate
  virtual void AduAu (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp);
  
  // apply a^\dagger_u_m a_u_m to a given state defined only in the down layer
  //
  // m = index for destruction operators
  // inputState = reference of the state to act on
  // outputState  = reference of the state where the result will be stored
  // nbrParticlesDown = number of particles in the down layer
  // lzDown = momentum of the down layer eigenstate
  virtual void AddAd (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesDown, int lzDown);
  
  // get the number of particles in a given state
  //
  // index =index of the state whose number of particles has to be returned
  // return value = number of particles
  virtual int GetTotalNumberOfParticles (int index);
  
  //get the maximal number of states that any state can be coupled to
  //
  //return value = number of coupling elements
  virtual int GetMaximalNumberCouplingElements();

  // convert a given state from a given  n-body basis basis to another one
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis where state is defined
  // return value = converted vector
  virtual RealVector ConvertToNbodyBasis(RealVector& state, QuasiholeOnSphereWithSpinAndPairing* nbodyBasis);

  // create a state from two single eigenstates
  //
  // eigenstateUp = reference of the eigenstate for the up layer
  // nbrParticlesUp = number of particles in the up layer
  // lzUp = momentum of the up layer eigenstate
  // eigenstateDown = reference of the eigenstate for the down layer
  // nbrParticlesDown = number of particles in the down layer
  // lzDown = momentum of the down layer eigenstate
  // return value = state built from the tensor product of the two single layer eigenstates
  RealVector BuildFromTwoSingleLayerEigenstates(RealVector& eigenstateUp, int nbrParticlesUp, int lzUp,
						RealVector& eigenstateDown, int nbrParticlesDown, int lzDown);

 protected:

  // evaluate Hilbert space dimension
  //
  // return value = Hilbert space dimension      
  virtual long EvaluateHilbertSpaceDimension();


  // generate all states corresponding to the constraints
  // 
  // return value = Hilbert space dimension   
  virtual long GenerateStates();
  
  
  // get the linear index corresponding to a set of number of fermions and momentum
  //
  // NbrParticles = number of fermions
  // totalLz = value of the angular momentum
  // return value = linearized index
  virtual int GetLinearIndexSingleLayer(int nbrParticles, int totalLz);
  
  // get the maximal value of Lz in one layer for a given number of particles in this layer
  //
  // nbrParticles = number of particles
  // return value = maximal Lz
  virtual int GetMaximalLzSingleLayer(int nbrParticles);
  
  // find state index from the value of each number in one layer
  //
  // nbrParticlesUp = number of particles with up spin
  // lzValueUp = value of the angular momentum for up spins
  // alpha = index of state with up spin
  // beta = index of state with down spin
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(int nbrParticlesUp, int lzValueUp, int alpha, int beta);
  
  // find the values of beta in each layer for a given state
  //
  // index = state index
  // alpha = reference on the value of alpha (up spin)
  // beta = reference on the value of beta (down spsin)
  virtual void FindBetaIndices(int index, int& alpha, int& beta);
  
};


// get the particle statistic 
//
// return value = particle statistic

inline int QuasiholeOnSphereWithSpinAndPairing::GetParticleStatistic()
{
  return ParticleOnSphereWithSpin::FermionicStatistic;
}

// get the linear index corresponding to a set of number of fermions and momentum
//
// NbrParticles = number of fermions
// totalLz = value of the angular momentum
// return value = linearized index

inline int QuasiholeOnSphereWithSpinAndPairing::GetLinearIndexSingleLayer(int nbrParticles, int totalLz)
{  
  return this->SingleLayerLinearIndices[nbrParticles][(totalLz + this->GetMaximalLzSingleLayer(nbrParticles)) / 2];
}
  

// get the maximal value of Lz in one layer for a given number of particles in this layer
//
// nbrParticles = number of particles
// return value = maximal Lz

inline int QuasiholeOnSphereWithSpinAndPairing::GetMaximalLzSingleLayer(int nbrParticles)
{
  return this->MaximalLzSingleLayer[nbrParticles];
}

  
// get the number of particles in a given state
//
// index =index of the state whose number of particles has to be returned
// return value = number of particles

inline int QuasiholeOnSphereWithSpinAndPairing::GetTotalNumberOfParticles (int index)
{
  return (2 * this->NbrFermionUpFullSpace[index] - this->TotalSpin);
}
  
// get the maximal number of states that any state can be coupled to
//
// return value = number of coupling elements

inline int QuasiholeOnSphereWithSpinAndPairing::GetMaximalNumberCouplingElements()
{
  return this->MaximalNumberCouplingElements;
}

// apply a_u_m a_d_n to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for the left annihilation operator
// n = index for the right annihilation operator
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

inline int QuasiholeOnSphereWithSpinAndPairing::AuAd (int index, int m, int n, int*& leftIndices, double*& interactionElements)
{
  if (m == n)
    return this->AuAd(index, m, leftIndices, interactionElements);
  else
    return 0;
}
  
// apply a^\dagger_u_m a^\dagger_d_n to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for the left annihilation operator
// n = index for the right annihilation operator
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

inline int QuasiholeOnSphereWithSpinAndPairing::AduAdd (int index, int m, int n, int*& leftIndices, double*& interactionElements)
{
  if (m == n)
    return this->AduAdd(index, m, leftIndices, interactionElements);
  else
    return 0;
}
  
// apply a^\dagger_u_m a_u_n to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation  operators
// n = index for annihilation operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

inline int QuasiholeOnSphereWithSpinAndPairing::AduAu (int index, int m, int n, int*& leftIndices, double*& interactionElements)
{
  if (m == n)
    return this->AduAu(index, m, leftIndices, interactionElements);
  else
    return 0;
}
  
// apply a^\dagger_d_m a_d_n to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation  operators
// n = index for annihilation operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

inline int QuasiholeOnSphereWithSpinAndPairing::AddAd (int index, int m, int n, int*& leftIndices, double*& interactionElements)
{
  if (m == n)
    return this->AddAd(index, m, leftIndices, interactionElements);
  else
    return 0;
}
  
#endif
