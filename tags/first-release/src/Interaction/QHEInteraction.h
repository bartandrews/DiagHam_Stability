////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quantum hall effect interaction                 //
//                                                                            //
//                        last modification : 25/11/2002                      //
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


#ifndef QHEINTERACTION_H
#define QHEINTERACTION_H


#include "config.h"
#include "Interaction/AbstractInteraction.h"


class QHEInteraction : public AbstractInteraction
{

 protected:

  // maximum momentum in the system
  int MaxMomentum;

  // flag to indicate if particles are fermions
  bool FermionFlag;

  // precision on interaction coefficient
  double Precision;

  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;

  double* CouplingConstants;
  int NbrCouplingConstant;

 public:
  
  // constructor
  //
  // fermion = flag to indicate if particles are fermions
  // maxMomentum = maximum momentum that can be reached in the system
  // ratio = ratio between the width in the x direction and the width in the y direction
  // precision = precision on interaction coefficient 
  QHEInteraction(bool fermion, int maxMomentum, double ratio, double precision = MACHINE_PRECISION);

  // copy constructor
  //
  // interaction = reference to interaction to copy
  QHEInteraction (const QHEInteraction& interaction);

  // destructor
  //
  ~QHEInteraction();

  // assignment
  //
  // interaction = reference to interaction to assign
  // return value = reference to current interaction
  QHEInteraction& operator = (const QHEInteraction& interaction);

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
