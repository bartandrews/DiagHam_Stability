////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of abstract numerical interval                   //
//                    (also representing the empty interval)                  //
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


#ifndef ABSTRACTNUMERICALINTERVAL_H
#define ABSTRACTNUMERICALINTERVAL_H


#include "config.h"


class AbstractNumericalInterval
{

 protected:

  // interval type including its subdivision scheme
  int Type;


 public:

  enum IntervalType
    {
      Empty = 0x00,
      Finite = 0x01,
      FullRInterval = 0x02,
      SemiInfiniteLeftBounded = 0x03,
      SemiInfiniteRightBounded = 0x04
    };

  enum SubdivisionScheme
    {
      NoSubdivision = 0x0000,
      RegularSubdivision = 0x0100,
    };

  // default contructor
  //
  AbstractNumericalInterval();

  // virtual destructor
  //
  virtual ~AbstractNumericalInterval();

  // clone interval
  //
  // pointer value = pointer to a clone interval
  virtual AbstractNumericalInterval* Clone();

  // get subdivision step around a given point
  // 
  // x = point where to look at
  // return value = local subdivision step
  virtual double GetSubdivisionStep (double x = 0.0);

  // get the total number of subdivisions
  // 
  // return value = total number of subdivisions
  virtual long GetNbrSubdivision ();

  // get interval type and its subdivision scheme
  //
  // return value = interval type
  int GetType();

  // get intersection of two intervals
  //
  // interval1 = first interval
  // interval2 = second interval
  // return value = pointer to the interval correponding to the inetrsection
  friend AbstractNumericalInterval* operator & (AbstractNumericalInterval& interval1, AbstractNumericalInterval& interval2);

  // check if two intervals are identical (same mathematical interval and same subdivision)
  //
  // interval1 = first interval
  // interval2 = second interval
  // return value =  true if the two intervals are identical
  friend bool operator == (AbstractNumericalInterval& interval1, AbstractNumericalInterval& interval2);

  // check if two intervals are different (different mathematical interval or different subdivision)
  //
  // interval1 = first interval
  // interval2 = second interval
  // return value = true if the two intervals are different
  friend bool operator != (AbstractNumericalInterval& interval1, AbstractNumericalInterval& interval2);

};

// get subdivision step around a given point
// 
// x = point where to look at
// return value = local subdivision step

inline double AbstractNumericalInterval::GetSubdivisionStep (double x)
{
  return 0.0;
}

// get the total number of subdivisions
// 
// return value = total number of subdivisions

inline long AbstractNumericalInterval::GetNbrSubdivision ()
{
  return (long) 0;
}

// get interval type and its subdivision scheme
//
// return value = interval type

inline int AbstractNumericalInterval::GetType()
{
  return this->Type;
}

#endif
