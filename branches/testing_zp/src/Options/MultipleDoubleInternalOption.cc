////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of multiple double option                      //
//                                                                            //
//                        last modification : 16/09/2001                      //
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


#include "Options/MultipleDoubleInternalOption.h"

#include <iostream>
#include <stdio.h>
#include <string.h>


using std::cout;
using std::endl;


// constructor from default datas
//
// optionCode = character associated to the option
// optionName = string corresponding to option name
// optionDescription = string describing option (used for -h option)
// separator = character used to separate single entries on the command line
// external = flag to indicate whether visible to the external world (in help output)
// altSeparator = character used to separate entries in output
// minValueFlag = flag to indicates an double minimum value
// minValue = double minimum value 
// maxValueFlag = flag to indicates an double maximum value
// maxValue = double maximum value (no maximum value if lower or equal to minValue) 
MultipleDoubleInternalOption::MultipleDoubleInternalOption(char optionCode, const char* optionName,
							   const char* optionDescription, char separator,
		       bool external, char altSeparator ,
		       char* defaultValues, bool minValueFlag, double minValue, 
		       bool maxValueFlag, double maxValue):
  MultipleDoubleOption( optionCode, optionName, optionDescription, separator, altSeparator,
			defaultValues, minValueFlag, minValue, maxValueFlag, maxValue )
{
  this->extFlag=external;
}

// destructor
//
MultipleDoubleInternalOption::~MultipleDoubleInternalOption()
{
}


// print the current option and its values
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream

ostream& MultipleDoubleInternalOption::DisplayOption (ostream& output, bool shortVersion)
{
   if (extFlag) return this->MultipleDoubleOption::DisplayOption(output, shortVersion);
  else return output;
}

// print help concerning current option
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& MultipleDoubleInternalOption::DisplayHelp (ostream& output)
{
   if (extFlag) return this->MultipleDoubleOption::DisplayHelp(output);
  else return output;
}


// Set doubles[pos] 
// to newVal

void MultipleDoubleInternalOption::SetSingleDouble(int pos, double newVal)
{
  if (pos < this->Length)
    this->Doubles[pos] = newVal;
}


// Set new array of doubles with given length
// 
void MultipleDoubleInternalOption::SetAllDoubles(int length, double *newVals)
{
  if (this->Doubles!=NULL) delete [] this->Doubles;
  this->Length=length;
  this->Doubles=new double [Length];
  for (int i=0; i<Length; ++i) this->Doubles[i]=newVals[i];
}
  

// multiply with factor
// returns: new value of internal variable
//
double* MultipleDoubleInternalOption::Multiply(double factor)
{
  for (int i=0; i<this->Length; ++i) this->Doubles[i] *= factor;
  return this->Doubles;
}

