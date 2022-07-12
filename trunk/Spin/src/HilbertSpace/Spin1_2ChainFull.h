////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of spin 1/2 chain without any Sz contraint           //
//                                                                            //
//                        last modification : 11/12/2013                      //
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

 
#ifndef SPIN1_2CHAINFULL_H
#define SPIN1_2CHAINFULL_H


#include "config.h"
#include "HilbertSpace/Spin1_2Chain.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainFull : public Spin1_2Chain
{

  friend class Spin1_2ChainFullAnd2DTranslation;
  friend class Spin1_2ChainFullInversionAnd2DTranslation;
  friend class ComplexPEPSPBC;

 protected:


 public:


  // default constructor
  //
  Spin1_2ChainFull ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  Spin1_2ChainFull (int chainLength);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainFull (const Spin1_2ChainFull& chain);

  // destructor
  //
  ~Spin1_2ChainFull ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainFull& operator = (const Spin1_2ChainFull& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSmj (int i, int j, int state, double& coefficient);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);
	
  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
	
  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // 
  // sites = list of sites that define the A subsystem 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int* sites, int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // 
  // sites = list of sites that define the A subsystem 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int* sites, int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseComplexMatrix* bMatrices, ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
					      long initialIndex = 0l, long nbrComponents = 0l);

  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (ComplexMatrix* bMatrices, ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
					      long initialIndex = 0l, long nbrComponents = 0l);

 protected:

  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long state);

};

// find state index
//
// state = state description
// return value = corresponding index

inline int Spin1_2ChainFull::FindStateIndex(unsigned long state)
{
  return ((int) state);    
}

#endif


