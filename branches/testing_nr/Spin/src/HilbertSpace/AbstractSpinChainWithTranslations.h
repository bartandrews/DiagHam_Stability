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
#include "HilbertSpace/AbstractHilbertSpace.h"


class Matrix;


class AbstractSpinChainWithTranslations : public AbstractHilbertSpace
{

 public:

  // virtual destructor
  //
  virtual ~AbstractSpinChainWithTranslations ();

  // return value of spin projection on (Oz) for a given state
  //
  // Str = reference on current output stream 
  // return value = spin projection on (Oz)
  virtual int TotalSz (int state) = 0;

  // return value of the value of the sum of the square of spin projection on (Oz) 
  //
  // index = index of the state to test
  // return value = twice spin projection on (Oz)
  virtual double TotalSzSz (int index) = 0;

  // get the momentum of each state in the current Hilbert space
  //
  // return value = momentum value
  virtual int GetMomentum() = 0;

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
    
  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long state) = 0;

};

#endif

