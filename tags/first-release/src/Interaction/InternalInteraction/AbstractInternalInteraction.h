////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of abstract internal interaction                  //
//                                                                            //
//                        last modification : 08/11/2002                      //
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


#ifndef ABSTRACTINTERNALINTERACTION_H
#define ABSTRACTINTERNALINTERACTION_H


#include "config.h"
#include "GeneralTools/List.h"


class Matrix;


class AbstractInternalInteraction
{

 protected:

  // real system size (needed when interactions depends on system real size)
  int SystemSize;

 public:


  // virtual destructor
  //
  virtual ~AbstractInternalInteraction();

  // evaluate amd add to hamitonian internal interaction terms
  //
  // hamiltonian = reference on hamiltonian matrix representation
  // operators = list of operators used to construct interaction terms
  // return value = reference on hamiltonian matrix representation
  virtual Matrix& AddInteraction (Matrix& hamiltonian, List<Matrix*>& operators) = 0;

  // set  real system size (needed when interactions depends on system real size)
  //
  // size = system size
  virtual void SetSystemSize (int size);

};

#endif
