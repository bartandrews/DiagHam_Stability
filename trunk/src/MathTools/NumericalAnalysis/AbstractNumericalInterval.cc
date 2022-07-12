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


#include "config.h"
#include "MathTools/NumericalAnalysis/AbstractNumericalInterval.h"


// default contructor
//

AbstractNumericalInterval::AbstractNumericalInterval()
{
  this->Type = AbstractNumericalInterval::Empty | AbstractNumericalInterval::NoSubdivision;
}

// virtual destructor
//

AbstractNumericalInterval::~AbstractNumericalInterval()
{
}

// clone interval
//
// pointer value = pointer to a clone interval

AbstractNumericalInterval* AbstractNumericalInterval::Clone()
{
  return new AbstractNumericalInterval;
}

// get intersection of two intervals
//
// interval1 = first interval
// interval2 = second interval
// return value = pointer to the interval correponding to the inetrsection

AbstractNumericalInterval* operator & (AbstractNumericalInterval& interval1, AbstractNumericalInterval& interval2)
{
  if ((interval1.Type == (AbstractNumericalInterval::Empty | AbstractNumericalInterval::NoSubdivision)) || 
      (interval2.Type == (AbstractNumericalInterval::Empty | AbstractNumericalInterval::NoSubdivision)))
    {
      return new AbstractNumericalInterval;      
    }
  return new AbstractNumericalInterval;
}

// check if two intervals are identical (same mathematical interval and same subdivision)
//
// interval1 = first interval
// interval2 = second interval
// return value =  true if the two intervals are identical

bool operator == (AbstractNumericalInterval& interval1, AbstractNumericalInterval& interval2)
{
  if (interval1.Type != interval2.Type)
    return false;
  return true;
}

// check if two intervals are different (different mathematical interval or different subdivision)
//
// interval1 = first interval
// interval2 = second interval
// return value = true if the two intervals are different

bool operator != (AbstractNumericalInterval& interval1, AbstractNumericalInterval& interval2)
{
  return (!(interval1 == interval2));
}
