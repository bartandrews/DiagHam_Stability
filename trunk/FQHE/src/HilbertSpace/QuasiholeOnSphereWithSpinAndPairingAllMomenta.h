////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of quasiholes on sphere with spin,                //
//                       pairing and all momentum sectors                     //
//                                                                            //
//                        last modification : 26/08/2016                      //
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


#ifndef QUASIHOLEONSPHEREWITHSPINANDPAIRINGALLMOMENTA_H
#define QUASIHOLEONSPHEREWITHSPINANDPAIRINGALLMOMENTA_H


#include "config.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "Matrix/SparseRealMatrix.h"

#include <iostream>


using std::cout;
using std::endl;


class QuasiholeOnSphereWithSpinAndPairingAllMomenta :  public QuasiholeOnSphereWithSpinAndPairing
{

 protected:

   // array containing the value of Lz in the down layer for each Hilbert space index
   int* LzValueDownFullSpace;
   // array containing the value of the lowest Hilbert space index with given quantum numbers N, Lzup and Lzdown
   int** FirstIndexWithNbrParticlesUpLzValueUpLzValueDown;

   // array where the all the a^+a matrix elements are stored
   RealMatrix**** SingleLayerFullAdAMatrices;

 public:

  // default constructor
  //
  QuasiholeOnSphereWithSpinAndPairingAllMomenta();

  // basic constructor
  // 
  // kExclusionPrinciple = k value of the exclusion principle
  // rExclusionPrinciple = r value of the exclusion principle
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twice the total spin value
  // maxMomentumTransfer = maximum momentum transfer allowed in the density operators
  // directory = optional path to data files
  // filePrefix = prefix for all input file (should include everything related to the statistics and the geometry)
  // memory = amount of memory granted for precalculations
  QuasiholeOnSphereWithSpinAndPairingAllMomenta (int kExclusionPrinciple, int rExclusionPrinciple, int lzMax, int totalSpin, 
						 int maxMomentumTransfer, const char* directory, const char* filePrefix, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // quasiholes = reference on the hilbert space to copy to copy
  QuasiholeOnSphereWithSpinAndPairingAllMomenta(const QuasiholeOnSphereWithSpinAndPairingAllMomenta& quasiholes);

  // destructor
  //
  ~QuasiholeOnSphereWithSpinAndPairingAllMomenta ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  QuasiholeOnSphereWithSpinAndPairingAllMomenta& operator = (const QuasiholeOnSphereWithSpinAndPairingAllMomenta& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();
  
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
  //  virtual int GetTotalNumberOfParticles (int index);
  
  // get the maximal number of states that any state can be coupled to
  //
  // return value = number of coupling elements
  //  virtual int GetMaximalNumberCouplingElements();

 protected:

  // evaluate Hilbert space dimension
  //
  // return value = Hilbert space dimension      
  virtual long EvaluateHilbertSpaceDimension();


  // generate all states corresponding to the constraints
  // 
  // return value = Hilbert space dimension   
  virtual long GenerateStates();
    
  // find state index from the value of each number in one layer
  //
  // nbrParticlesUp = number of particles with up spin
  // nbrParticlesDown = number of particles with down spin
  // lzValueUp = value of the angular momentum for up spins
  // lzValueDown = value of the angular momentum for down spins
  // alpha = index of state with up spin
  // beta = index of state with down spin
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(int nbrParticlesUp, int nbrParticlesDown, int lzValueUp, int lzValueDown, int alpha, int beta);
  
  // find the values of beta in each layer for a given state
  //
  // index = state index
  // alpha = reference on the value of alpha (up spin)
  // beta = reference on the value of beta (down spsin)
  virtual void FindBetaIndices(int index, int& alpha, int& beta);
  
};

// apply a^\dagger_u_m a^\dagger_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

inline int QuasiholeOnSphereWithSpinAndPairingAllMomenta::AduAdd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  return this->AduAdd(index, m, m, leftIndices, interactionElements);
}

  
#endif
