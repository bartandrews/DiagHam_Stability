////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                        Class author Cecile Repellin                        //
//                                                                            //
//                                                                            //
//                    class of quasiholes on sphere with spin                 //
//               (i.e. both  Sz and the particle number conservation)         //
//                                                                            //
//                        last modification : 28/05/2016                      //
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


#ifndef QUASIHOLEONSPHEREWITHSPIN_H
#define QUASIHOLEONSPHEREWITHSPIN_H


#include "config.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "Matrix/SparseRealMatrix.h"

#include <iostream>


using std::cout;
using std::endl;


class QuasiholeOnSphereWithSpin :  public QuasiholeOnSphereWithSpinAndPairing
{

 public:

  // default constructor
  //
  QuasiholeOnSphereWithSpin();

  // basic constructor
  // 
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // nbrParticles = number of particles
  // totalSpin = twice the total spin value
  // directory = optional path to data files
  // filePrefix = prefix for all input file (should include everything related to the statistics and the geometry)
  // memory = amount of memory granted for precalculations
  QuasiholeOnSphereWithSpin (int kExclusionPrinciple, int rExclusionPrinciple, int totalLz, int lzMax, int nbrParticles, int totalSpin, 
			     const char* directory, const char* filePrefix, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // quasiholes = reference on the hilbert space to copy to copy
  QuasiholeOnSphereWithSpin(const QuasiholeOnSphereWithSpin& quasiholes);

  // destructor
  //
  ~QuasiholeOnSphereWithSpin ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  QuasiholeOnSphereWithSpin& operator = (const QuasiholeOnSphereWithSpin& fermions);

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
  
  // get the number of particles in a given state
  //
  // index =index of the state whose number of particles has to be returned
  // return value = number of particles
  virtual int GetTotalNumberOfParticles (int index);
  

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


// get the linear index corresponding to a set of number of fermions and momentum
//
// NbrParticles = number of fermions
// totalLz = value of the angular momentum
// return value = linearized index

inline int QuasiholeOnSphereWithSpin::GetLinearIndexSingleLayer(int nbrParticles, int totalLz)
{
  return this->SingleLayerLinearIndices[nbrParticles][(totalLz + this->GetMaximalLzSingleLayer(nbrParticles)) / 2];
}
  

// get the maximal value of Lz in one layer for a given number of particles in this layer
//
// nbrParticles = number of particles
// return value = maximal Lz

inline int QuasiholeOnSphereWithSpin::GetMaximalLzSingleLayer(int nbrParticles)
{
  return this->MaximalLzSingleLayer[nbrParticles];
}

  
// get the number of particles in a given state
//
// index =index of the state whose number of particles has to be returned
// return value = number of particles

inline int QuasiholeOnSphereWithSpin::GetTotalNumberOfParticles (int index)
{
  return this->NbrFermions;
}
  
// apply a_u_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

inline int QuasiholeOnSphereWithSpin::AuAd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  return 0;
}


// apply a_u_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

inline int QuasiholeOnSphereWithSpin::AduAdd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  return 0;    
}

#endif
