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


#ifndef ABSTRACTINTERACTION_H
#define ABSTRACTINTERACTION_H


#include "config.h"
#include "GeneralTools/List.h"
#include "TensorProduct/AbstractTensorProductStructure.h"


class TwoSpaceTensor;
class Matrix;
class SpaceDecomposition;

class AbstractInteraction
{

 protected:

  int LeftSpaceIndex;
  int RightSpaceIndex;
  AbstractTensorProductStructure* Structure;
  
 public:


  // virtual destructor
  //
  virtual ~AbstractInteraction();

  // evaluate interaction between two systems from operators of each system
  //
  // firstSystemOperators = list of operators associated to the first system
  // secondSystemOperators = list of operators associated to the second system
  // return value = tensor corresponding to the interaction
  virtual TwoSpaceTensor Interaction (List<Matrix*>& firstSystemOperators, 
				      List<Matrix*>& secondSystemOperators) = 0;

  // evaluate interaction between two systems from operators of each system with a given
  // space decomposition for each space
  //
  // firstSystemOperators = list of operators associated to the first system
  // secondSystemOperators = list of operators associated to the second system
  // firstSpaceDecomposition = space decomposition of the space associated to the first system
  // secondSpaceDecomposition = space decomposition of the space associated to the second system
  // return value = tensor corresponding to the interaction
  virtual TwoSpaceTensor Interaction (List<Matrix*>& firstSystemOperators, 
				      List<Matrix*>& secondSystemOperators, 
				      SpaceDecomposition& firstSpaceDecomposition,
				      SpaceDecomposition& secondSpaceDecomposition) = 0;

  // set index of space where left interaction acts
  //
  // index = space index
  virtual void SetLeftSpaceIndex (int index);

  // set index of space where right interaction acts
  //
  // index = space index
  virtual void SetRightSpaceIndex (int index);

  // set structure of tensor space where interactions act
  //
  // structure = tensor space structure
  virtual void SetTensorProductStructure (AbstractTensorProductStructure* structure);

};

#endif
