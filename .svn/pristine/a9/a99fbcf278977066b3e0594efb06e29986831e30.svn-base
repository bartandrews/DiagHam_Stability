////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of abstract spin chain with translations             //
//                                                                            //
//                        last modification : 04/03/2002                      //
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


#ifndef ABSTRACTSPINCHAINWITHTRANSLATIONS_H
#define ABSTRACTSPINCHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"


class Matrix;


class AbstractSpinChainWithTranslations : public AbstractSpinChain
{

 protected:
  
  // momentum 
  int Momentum;

 public:

  // virtual destructor
  //
  virtual ~AbstractSpinChainWithTranslations ();

  // return value of the value of the sum of the square of spin projection on (Oz) 
  //
  // index = index of the state to test
  // return value = twice spin projection on (Oz)
  virtual double TotalSzSz (int index) = 0;

  // get the momentum of each state in the current Hilbert space
  //
  // return value = momentum value
  virtual int GetMomentum() {return this->Momentum;};

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  virtual double SziSzj (int i, int j, int state) = 0;

  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation) = 0;

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation) = 0;

  // return index of resulting state from application of S+_i S+_i operator on a given state
  //
  // i = position of first S+ operator
  // state = index of the state to be applied on S+_i S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SpiSpi (int i, int state, double& coefficient, int& nbrTranslation) = 0;

  // return index of resulting state from application of S-_i S-_i operator on a given state
  //
  // i = position of the S- operator
  // state = index of the state to be applied on S-_i S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SmiSmi (int i, int state, double& coefficient, int& nbrTranslation) = 0;

  // return index of resulting state from application of S+_i S-_i operator on a given state
  //
  // i = position of first S+ operator
  // state = index of the state to be applied on S+_i S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SpiSmi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i S+_i operator on a given state
  //
  // i = position of first S+ operator
  // state = index of the state to be applied on S-_i S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SmiSpi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i Sz_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of Sz operator
  // state = index of the state to be applied on S+_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation) = 0;

  // return index of resulting state from application of S-_i Sz_j operator on a given state
  //
  // i = position of S- operator
  // j = position of Sz operator
  // state = index of the state to be applied on S-_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation) = 0;

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state (orbit index)
  virtual int SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation) = 0;

  // return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of leftmost S- operator
  // j1 = position of leftmost S+ operator
  // i2 = position of rightmost S- operator
  // j2 = position of rightmost S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state (orbit index)
  virtual int SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i1 S-_j1 S+_i2 S+_j2 operator on a given state
  //
  // i1 = position of leftmost S- operator
  // j1 = position of rightmost S- operator
  // i2 = position of leftmost S+ operator
  // j2 = position of rightmost S+ operator
  // state = index of the state to be applied on  S-_i1 S-_j1 S+_i2 S+_j2 operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state (orbit index)
  virtual int SmiSmjSpiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of Sz_i1 Sz_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of first Sz operator
  // j1 = position of second Sz operator
  // i2 = position of S- operator
  // j2 = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state (orbit index)
  virtual int SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i operator on a given state (only valid if there is no constraint on total Sz)
  //
  // i = operator position
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int Spi (int i, int state, double& coefficient, int& nbrTranslation) = 0;

  // return index of resulting state from application of S-_i operator on a given state (only valid if there is no constraint on total Sz)
  //
  // i = operator position
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int Smi (int i, int state, double& coefficient, int& nbrTranslation) = 0;

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

  // compute all the states connected to a single one by a two site spin 0 projector
  //
  // i = index of the first site
  // j = index of the first site
  // state = index of state on whcih the projector should be applied
  // indices = pointer to the array where the connected state indeices will be stored
  // coefficients = pointer to the array where coefficients related to each connected states will be stored
  // nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
  // return value = number of connected states
  virtual int Spin0Projector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations);

  // compute all the states connected to a single one by a two site spin 1 projector
  //
  // i = index of the first site
  // j = index of the first site
  // state = index of state on whcih the projector should be applied
  // indices = pointer to the array where the connected state indeices will be stored
  // coefficients = pointer to the array where coefficients related to each connected states will be stored
  // nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
  // return value = number of connected states
  virtual int Spin1Projector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations);

  // compute all the states connected to a single one by a two site spin 2 projector
  //
  // i = index of the first site
  // j = index of the first site
  // state = index of state on whcih the projector should be applied
  // indices = pointer to the array where the connected state indeices will be stored
  // coefficients = pointer to the array where coefficients related to each connected states will be stored
  // nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
  // return value = number of connected states
  virtual int Spin2Projector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations);

  // compute all the states connected to a single one by a two site spin 3 projector
  //
  // i = index of the first site
  // j = index of the first site
  // state = index of state on whcih the projector should be applied
  // indices = pointer to the array where the connected state indeices will be stored
  // coefficients = pointer to the array where coefficients related to each connected states will be stored
  // nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
  // return value = number of connected states
  virtual int Spin3Projector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations);

  // compute all the states connected to a single one by a two site spin 4 projector
  //
  // i = index of the first site
  // j = index of the first site
  // state = index of state on whcih the projector should be applied
  // indices = pointer to the array where the connected state indeices will be stored
  // coefficients = pointer to the array where coefficients related to each connected states will be stored
  // nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
  // return value = number of connected states
  virtual int Spin4Projector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations);

  // find the canonical form of a state
  //
  // state = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // return value = canonical form of the state
  virtual unsigned long FindCanonicalForm(unsigned long state, int& nbrTranslation);

  // return the scaling factor when going from state i to state j
  //
  virtual double GetRescalingFactor(int i,int j) const {return 0;};
};

#endif


