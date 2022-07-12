////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                                                                            //
//                       class of abstract spin collection                    //
//                                                                            //
//                        last modification : 18/04/2001                      //
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


#ifndef ABSTRACTSUNSPINCOLLECTION_H
#define ABSTRACTSUNSPINCOLLECTION_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"


class AbstractSUNSpinCollection : public AbstractHilbertSpace
{

 public:

  // virtual destructor
  //
  virtual ~AbstractSUNSpinCollection ();
  
  // get diagonal terms for an S*S interaction (counting instances for connections with same prefactor)
  // index = number of state to be considered
  virtual int S2DiagonalElements(int index) = 0;

  // get off-diagonal terms for an S*S interaction (counting instances for connections with same prefactor)
  // index = number of state to be considered
  // spin = number of spin associated with first spin operator
  // targetIndices = states connected to (size: N-1)
  //
  virtual void S2OffDiagonalElements(int index, int spin, int *targetIndices) = 0;
  
  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long state) = 0;

};

#endif


