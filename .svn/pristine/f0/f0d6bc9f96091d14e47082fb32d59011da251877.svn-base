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


#include "config.h"
#include "MathTools/NumericalAnalysis/RNumericalInterval.h"


// default contructor
//

RNumericalInterval::RNumericalInterval()
{
  this->Type = AbstractNumericalInterval::FullRInterval | AbstractNumericalInterval::NoSubdivision;
}

// destructor
//

RNumericalInterval::~RNumericalInterval()
{
}

// clone interval
//
// pointer value = pointer to a clone interval

AbstractNumericalInterval* RNumericalInterval::Clone()
{
  return new RNumericalInterval;
}

// get intersection of two intervals
//
// interval1 = first interval
// interval2 = second interval
// return value = pointer to the interval correponding to the inetrsection

AbstractNumericalInterval* operator & (RNumericalInterval& interval1, AbstractNumericalInterval& interval2)
{
  if (interval2.GetType() == (AbstractNumericalInterval::Empty | AbstractNumericalInterval::NoSubdivision))
    {
      return new AbstractNumericalInterval;
    }
  else
    {
      return interval2.Clone();
    }
}

// check if two intervals are identical (same mathematical interval and same subdivision)
//
// interval1 = first interval
// interval2 = second interval
// return value =  true if the two intervals are identical

bool operator == (RNumericalInterval& interval1, AbstractNumericalInterval& interval2)
{
  if (interval2.GetType() == (AbstractNumericalInterval::FullRInterval | AbstractNumericalInterval::NoSubdivision))
    {
      return true;
    }
  else
    {
      return false;
    }
}

// check if two intervals are different (different mathematical interval or different subdivision)
//
// interval1 = first interval
// interval2 = second interval
// return value = true if the two intervals are different

bool operator != (RNumericalInterval& interval1, AbstractNumericalInterval& interval2)
{
  if (interval2.GetType() != (AbstractNumericalInterval::FullRInterval | AbstractNumericalInterval::NoSubdivision))
    {
      return true;
    }
  else
    {
      return false;
    }
}
