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
#include "MathTools/NumericalAnalysis/FiniteRegularSubdivisionNumericalInterval.h"


#include <math.h>


// contructor knowing number of subdivision
//
// minimum = interval lower bound value
// maximum = interval upper bound value
// nbrStep = number of subdivisions

FiniteRegularSubdivisionNumericalInterval::FiniteRegularSubdivisionNumericalInterval(double minimum, double maximum, unsigned long nbrStep)
{
  this->Type = AbstractNumericalInterval::Finite | AbstractNumericalInterval::RegularSubdivision;
  this->Minimum = minimum;
  this->Maximum = maximum;
  this->NbrStep = nbrStep;
  this->StepValue = (this->Maximum - this->Minimum) / ((double) this->NbrStep);
}

// contructor knowing step value
//
// minimum = interval lower bound value
// maximum = interval upper bound value
// stepValue = step value

FiniteRegularSubdivisionNumericalInterval::FiniteRegularSubdivisionNumericalInterval(double minimum, double maximum, double stepValue)
{
  this->Type = AbstractNumericalInterval::Finite | AbstractNumericalInterval::RegularSubdivision;
  this->Minimum = minimum;
  this->Maximum = maximum;
  double Tmp = (this->Maximum - this->Minimum) / stepValue;
  if (Tmp != floor(Tmp))
    {
      this->NbrStep = (unsigned long) Tmp ;  
    }
  else
    {
      this->NbrStep = (unsigned long) (Tmp + 1.0);  
    }
  this->StepValue = (this->Maximum - this->Minimum) / ((double) this->NbrStep);
}

// copy contructor
//
// interval = reference on the interval to copy

FiniteRegularSubdivisionNumericalInterval::FiniteRegularSubdivisionNumericalInterval(const FiniteRegularSubdivisionNumericalInterval& interval)
{
  this->Type = AbstractNumericalInterval::Finite | AbstractNumericalInterval::RegularSubdivision;
  this->Minimum = interval.Minimum;
  this->Maximum = interval.Maximum;
  this->NbrStep = interval.NbrStep;
  this->StepValue = (this->Maximum - this->Minimum) / ((double) this->NbrStep);
}

// destructor
//

FiniteRegularSubdivisionNumericalInterval::~FiniteRegularSubdivisionNumericalInterval()
{
}

// assignment
//
// interval = reference on the interval to assign
//return value = reference on the current interval

FiniteRegularSubdivisionNumericalInterval& FiniteRegularSubdivisionNumericalInterval::operator = (const FiniteRegularSubdivisionNumericalInterval& interval)
{
  this->Type = AbstractNumericalInterval::Finite | AbstractNumericalInterval::RegularSubdivision;
  this->Minimum = interval.Minimum;
  this->Maximum = interval.Maximum;
  this->NbrStep = interval.NbrStep;
  this->StepValue = (this->Maximum - this->Minimum) / ((double) this->NbrStep);
  return *this;
}

// clone interval
//
// pointer value = pointer to a clone interval

AbstractNumericalInterval* FiniteRegularSubdivisionNumericalInterval::Clone()
{
  return new FiniteRegularSubdivisionNumericalInterval(*this);
}

// get intersection of two intervals
//
// interval1 = first interval
// interval2 = second interval
// return value = pointer to the interval correponding to the inetrsection

AbstractNumericalInterval* operator & (FiniteRegularSubdivisionNumericalInterval& interval1, 
				       AbstractNumericalInterval& interval2)
{
  if (interval2.GetType() == (AbstractNumericalInterval::Empty | AbstractNumericalInterval::NoSubdivision))
    return new AbstractNumericalInterval;
  if (interval2.GetType() == (AbstractNumericalInterval::Finite | AbstractNumericalInterval::RegularSubdivision))
    {
      FiniteRegularSubdivisionNumericalInterval& TmpInterval = (FiniteRegularSubdivisionNumericalInterval&) interval2;
      if (interval1.Minimum > TmpInterval.Minimum)
	{
	  FiniteRegularSubdivisionNumericalInterval& TmpInterval2 = interval1;
	  interval1 = TmpInterval;
	  TmpInterval = TmpInterval2;
	}
      if (interval1.Maximum < TmpInterval.Minimum)
	return new AbstractNumericalInterval;
      if ((interval1.Minimum == TmpInterval.Minimum) && (interval1.Maximum == TmpInterval.Maximum))
	{
	  if (interval1.NbrStep >= TmpInterval.NbrStep)
	    {
	      return new FiniteRegularSubdivisionNumericalInterval (interval1);
	    }
	  else
	    {
	      return new FiniteRegularSubdivisionNumericalInterval (TmpInterval);	
	    }
	}
      if (interval1.Maximum < TmpInterval.Maximum)
	{
	  if (interval1.StepValue >= TmpInterval.StepValue)
	    {
	      return new FiniteRegularSubdivisionNumericalInterval (TmpInterval.Minimum, interval1.Maximum, TmpInterval.StepValue);      
	    }
	  else
	    {
	      return new FiniteRegularSubdivisionNumericalInterval (TmpInterval.Minimum, interval1.Maximum, interval1.StepValue);      
	    }
	}
      else
	{
	  if (interval1.StepValue >= TmpInterval.StepValue)
	    {
	      return new FiniteRegularSubdivisionNumericalInterval (TmpInterval.Minimum, TmpInterval.Maximum, TmpInterval.StepValue);      
	    }
	  else
	    {
	      return new FiniteRegularSubdivisionNumericalInterval (TmpInterval.Minimum, TmpInterval.Maximum, interval1.StepValue);      
	    }
	}
    }
  return new AbstractNumericalInterval;
}

// check if two intervals are identical (same mathematical interval and same subdivision)
//
// interval1 = first interval
// interval2 = second interval
// return value =  true if the two intervals are identical

bool operator == (FiniteRegularSubdivisionNumericalInterval& interval1, AbstractNumericalInterval& interval2)
{
  if (interval2.GetType() != (AbstractNumericalInterval::Finite | AbstractNumericalInterval::RegularSubdivision))
    {
      return false;
    }
  if ((interval1.Minimum != ((FiniteRegularSubdivisionNumericalInterval&) interval2).Minimum) ||
      (interval1.Maximum != ((FiniteRegularSubdivisionNumericalInterval&) interval2).Maximum) ||
      (interval1.NbrStep != ((FiniteRegularSubdivisionNumericalInterval&) interval2).NbrStep))
    return false;
  return true;
}

// check if two intervals are different (different mathematical interval or different subdivision)
//
// interval1 = first interval
// interval2 = second interval
// return value = true if the two intervals are different

bool operator != (FiniteRegularSubdivisionNumericalInterval& interval1, AbstractNumericalInterval& interval2)
{
  if (interval2.GetType() != (AbstractNumericalInterval::Finite | AbstractNumericalInterval::RegularSubdivision))
    {
      return true;
    }
  if ((interval1.Minimum != ((FiniteRegularSubdivisionNumericalInterval&) interval2).Minimum) ||
      (interval1.Maximum != ((FiniteRegularSubdivisionNumericalInterval&) interval2).Maximum) ||
      (interval1.NbrStep != ((FiniteRegularSubdivisionNumericalInterval&) interval2).NbrStep))
    return true;
  return false;
}
