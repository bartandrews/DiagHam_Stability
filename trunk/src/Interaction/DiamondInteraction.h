////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of diamond interaction                       //
//                                                                            //
//                        last modification : 14/12/2001                      //
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


#ifndef DIAMONDINTERACTION_H
#define DIAMONDINTERACTION_H


#include "config.h"
#include "Interaction/AbstractInteraction.h"


class DiamondInteraction : public AbstractInteraction
{

 protected:

  double UpperCouplingConstant;
  double LowerCouplingConstant;

 public:
  
  // default constructor
  //
  DiamondInteraction();

  // constructor from partial datas
  //
  // upperCouplingConstant = couling constant for upper link
  // lowerCoupligConstant = couling constant for lower link
  DiamondInteraction(double upperCouplingConstant, double lowerCouplingConstant);

  // constructor from complete datas
  //
  // upperCouplingConstant = couling constant for upper link
  // lowerCoupligConstant = couling constant for lower link
  // leftSpaceIndex = index of space where left interaction acts
  // rightSpaceIndex = index of space where right interaction acts
  // struture = reference on tensor product structure
  DiamondInteraction(double upperCouplingConstant, double lowerCouplingConstant, int leftSpaceIndex, 
		     int rightSpaceIndex, AbstractTensorProductStructure* structure);

  // copy constructor
  //
  // interaction = reference to interaction to copy
  DiamondInteraction (const DiamondInteraction& interaction);

  // destructor
  //
  ~DiamondInteraction();

  // assignment
  //
  // interaction = reference to interaction to assign
  // return value = reference to current interaction
  DiamondInteraction& operator = (const DiamondInteraction& interaction);

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
