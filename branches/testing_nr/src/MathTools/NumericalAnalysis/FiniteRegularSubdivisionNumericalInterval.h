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


#ifndef FINITEREGULARSUBDIVISIONNUMERICALINTERVAL_H
#define FINITEREGULARSUBDIVISIONNUMERICALINTERVAL_H


#include "config.h"
#include "MathTools/NumericalAnalysis/AbstractNumericalInterval.h"


class FiniteRegularSubdivisionNumericalInterval : public AbstractNumericalInterval
{

 protected:

  // interval lower bound value
  double Minimum;
  // interval upper bound value
  double Maximum;

  // number of subdivisions
  unsigned long NbrStep;
  // setp value
  double StepValue;

 public:

  // contructor knowing number of subdivision
  //
  // minimum = interval lower bound value
  // maximum = interval upper bound value
  // nbrStep = number of subdivisions
  FiniteRegularSubdivisionNumericalInterval(double minimum, double maximum, unsigned long nbrStep);

  // contructor knowing step value
  //
  // minimum = interval lower bound value
  // maximum = interval upper bound value
  // stepValue = step value
  FiniteRegularSubdivisionNumericalInterval(double minimum, double maximum, double stepValue);

  // copy contructor
  //
  // interval = reference on the interval to copy
  FiniteRegularSubdivisionNumericalInterval(const FiniteRegularSubdivisionNumericalInterval& interval);

  // destructor
  //
  ~FiniteRegularSubdivisionNumericalInterval();

  // assignment
  //
  // interval = reference on the interval to assign
  //return value = reference on the current interval
  FiniteRegularSubdivisionNumericalInterval& operator = (const FiniteRegularSubdivisionNumericalInterval& interval);

  // clone interval
  //
  // pointer value = pointer to a clone interval
  AbstractNumericalInterval* Clone();

  // get intersection of two intervals
  //
  // interval1 = first interval
  // interval2 = second interval
  // return value = pointer to the interval correponding to the inetrsection
  friend AbstractNumericalInterval* operator & (FiniteRegularSubdivisionNumericalInterval& interval1, AbstractNumericalInterval& interval2);

  // check if two intervals are identical (same mathematical interval and same subdivision)
  //
  // interval1 = first interval
  // interval2 = second interval
  // return value =  true if the two intervals are identical
  friend bool operator == (FiniteRegularSubdivisionNumericalInterval& interval1, AbstractNumericalInterval& interval2);

  // check if two intervals are different (different mathematical interval or different subdivision)
  //
  // interval1 = first interval
  // interval2 = second interval
  // return value = true if the two intervals are different
  friend bool operator != (FiniteRegularSubdivisionNumericalInterval& interval1, AbstractNumericalInterval& interval2);

};

#endif
