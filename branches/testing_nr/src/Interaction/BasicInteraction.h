////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract interaction                      //
//                                                                            //
//                        last modification : 05/04/2001                      //
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


#ifndef BASICINTERACTION_H
#define BASICINTERACTION_H


#include "config.h"
#include "Interaction/AbstractInteraction.h"


class BasicInteraction : public AbstractInteraction
{

 protected:

  double* CouplingConstants;
  int NbrCouplingConstant;

 public:
  
  // default constructor
  //
  BasicInteraction();

  // constructor from partial datas
  //
  // couplingConstants = pointer to array containing coupling constants
  // nbrCouplingConstant = number of coupling constants
  BasicInteraction(double* couplingConstants, int nbrCouplingConstant);

  // constructor from complete datas
  //
  // couplingConstants = pointer to array containing coupling constants
  // nbrCouplingConstant = number of coupling constants
  // leftSpaceIndex = index of space where left interaction acts
  // rightSpaceIndex = index of space where right interaction acts
  // struture = reference on tensor product structure
  BasicInteraction(double* couplingConstants, int nbrCouplingConstant, int leftSpaceIndex, 
		   int rightSpaceIndex, AbstractTensorProductStructure* structure);

  // copy constructor
  //
  // interaction = reference to interaction to copy
  BasicInteraction (const BasicInteraction& interaction);

  // destructor
  //
  ~BasicInteraction();

  // assignment
  //
  // interaction = reference to interaction to assign
  // return value = reference to current interaction
  BasicInteraction& operator = (const BasicInteraction& interaction);

  // evaluate interaction between two systems from operators of each system
  //
  // firstSystemOperators = list of operators associated to the first system
  // secondSystemOperators = list of operators associated to the second system
  // return value = tensor corresponding to the interaction
  TwoSpaceTensor Interaction (List<Matrix*>& firstSystemOperators, 
			      List<Matrix*>& secondSystemOperators);

  // evaluate interaction between two systems from operators of each system with a given
  // space decomposition for each space
  //
  // firstSystemOperators = list of operators associated to the first system
  // secondSystemOperators = list of operators associated to the second system
  // firstSpaceDecomposition = space decomposition of the space associated to the first system
  // secondSpaceDecomposition = space decomposition of the space associated to the second system
  // return value = tensor corresponding to the interaction
  TwoSpaceTensor Interaction (List<Matrix*>& firstSystemOperators, 
			      List<Matrix*>& secondSystemOperators, 
			      SpaceDecomposition& firstSpaceDecomposition,
			      SpaceDecomposition& secondSpaceDecomposition);

};

#endif
