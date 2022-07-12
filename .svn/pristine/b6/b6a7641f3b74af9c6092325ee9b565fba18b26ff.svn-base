////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of abstract doubled spin chain                  //
//                                                                            //
//                        last modification : 21/01/2016                      //
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


#ifndef ABSTRACTDOUBLEDSPINCHAIN_H
#define ABSTRACTDOUBLEDSPINCHAIN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include <iostream>


using std::ostream;

class AbstractDoubledSpinChain : public AbstractSpinChain
{
  
  friend class ComplexPEPSTransfertMatrixPBC;
  
 protected:
  
  // flag to indicate if the total sz component is fixed
  bool FixedSpinProjectionFlag;

  //  total difference sz component (if fixed)
  int DiffSz;

  // shift to apply to a state to obtain an index to the look-up table 
  int LookUpTableShift;
  // look-up table (LookUpTable[i] gives the index of the smallest state that greater than i <<  LookUpTableShift)
  long* LookUpTable;

  // array describing each state
  unsigned long* ChainDescriptionBra;
  unsigned long* ChainDescriptionKet;


  // sorted array that contains each unique configuration for the type up particles
  unsigned long* UniqueStateDescriptionBra;
  // number of time each unique configuration for the type up particles appears in StateDescriptionUp
  int* UniqueStateDescriptionSubArraySizeBra;
  // number of unique configurations for the type up-plus particles
  long NbrUniqueStateDescriptionBra;
  // first time a type up appears in the Hilbert space
  int* FirstIndexUniqueStateDescriptionBra;


 public:

  // destructor
  //
  ~AbstractDoubledSpinChain ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  AbstractDoubledSpinChain & operator = (const AbstractDoubledSpinChain & chain);

  // return Hilbert space dimension
  //
  // return value = Hilbert space dimension
  int GetHilbertSpaceDimension();

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  int TotalSz (int index);
  
  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateBra, unsigned long stateKet);

  virtual int FindStateIndexFromLinearizedIndex(unsigned long linearizedState);
  
  virtual void GetChainDescriptionInCondensedForm(unsigned long * HilbertSpaceDescription);

  virtual void AddConvertFromGeneralSpace(ComplexVector vSource,ComplexVector & vDestination){cout<<"using undefined   virtual void AddConvertFromGeneralSpace(ComplexVector vSource,ComplexVector & vDestination) in AbstractDoubledSpinChain"<<endl;}
  virtual void ConvertToGeneralSpace(ComplexVector vSource,ComplexVector & vDestination){cout<<"using undefined   virtual void ConvertToGeneralSpace(ComplexVector vSource,ComplexVector & vDestination) in AbstractDoubledSpinChain"<<endl;};


  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int szSector, RealVector& groundState){ return RealSymmetricMatrix(); };
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int szSector, int momentumSector, RealVector& groundState){ return RealSymmetricMatrix(); };
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int szSector, ComplexVector& groundState){return HermitianMatrix(); } ;
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int szSector, int momentumSector, ComplexVector& groundState){return HermitianMatrix(); } ;
  

 protected:

  // generate all states corresponding to the constraints
  // 
  // lengthBra = length of the chain to be decided for bra spins
  // lengthBra = length of the chain to be decided for ket spins
  // diffSz = difference of spin projection between bra and ket chain
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos);
  
  // create look-up table used to speed up index search
  //
  void GenerateLookUpTable(unsigned long memory);


  virtual int GetTotalSz (unsigned long stateDescriptionBra,unsigned long stateDescriptionKet){ return 0;};
  
/*  double TotalSzSz (int index);
  double SziSzj (int i, int j, int state);
  int SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int SpiSpi (int i, int state, double& coefficient, int& nbrTranslation);
  int SmiSmi (int i, int state, double& coefficient, int& nbrTranslation);
  int SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int Spi (int i, int state, double& coefficient, int& nbrTranslation);
  int Smi (int i, int state, double& coefficient, int& nbrTranslation);*/

  // return index of resulting state from application of S-_i operator on a given state
  //
  // i = position of S- operator
  // state = index of the state to be applied on S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Smi (int i, int state, double& coefficient);
  
  // return index of resulting state from application of Sz_i operator on a given state
  //
  // i = position of Sz operator
  // state = index of the state to be applied on Sz_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state 
  virtual int Szi (int i, int state, double& coefficient);
  
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
  
  // return index of resulting state from application of S+_i Sz_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of Sz operator
  // state = index of the state to be applied on S+_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSzj (int i, int j, int state, double& coefficient);
  
  // return index of resulting state from application of S-_i Sz_j operator on a given state
  //
  // i = position of S- operator
  // j = position of Sz operator
  // state = index of the state to be applied on S-_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSzj (int i, int j, int state, double& coefficient);
  
  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state 
  virtual int SmiSpj (int i, int j, int state, double& coefficient);
  
  // return index of resulting state from application of S+_i operator on a given state
  //
  // i = position of S+ operator
  // state = index of the state to be applied on S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Spi (int i, int state, double& coefficient);

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  virtual double SziSzj (int i, int j, int state) ;
  

};



#endif


