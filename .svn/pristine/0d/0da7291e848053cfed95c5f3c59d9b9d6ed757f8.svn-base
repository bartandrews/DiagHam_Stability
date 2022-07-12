////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of numerical interval representing R                //
//                                                                            //
//                        last modification : 07/07/2004                      //
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


#ifndef RNUMERICALINTERVAL_H
#define RNUMERICALINTERVAL_H


#include "config.h"
#include "MathTools/NumericalAnalysis/AbstractNumericalInterval.h"


class RNumericalInterval : public AbstractNumericalInterval
{

 public:

  // default contructor
  //
  RNumericalInterval();

  // destructor
  //
  ~RNumericalInterval();

  // clone interval
  //
  // pointer value = pointer to a clone interval
  AbstractNumericalInterval* Clone();

  // get intersection of two intervals
  //
  // interval1 = first interval
  // interval2 = second interval
  // return value = pointer to the interval correponding to the inetrsection
  friend AbstractNumericalInterval* operator & (RNumericalInterval& interval1, AbstractNumericalInterval& interval2);

  // check if two intervals are identical (same mathematical interval and same subdivision)
  //
  // interval1 = first interval
  // interval2 = second interval
  // return value =  true if the two intervals are identical
  friend bool operator == (RNumericalInterval& interval1, AbstractNumericalInterval& interval2);

  // check if two intervals are different (different mathematical interval or different subdivision)
  //
  // interval1 = first interval
  // interval2 = second interval
  // return value = true if the two intervals are different
  friend bool operator != (RNumericalInterval& interval1, AbstractNumericalInterval& interval2);

};

#endif
