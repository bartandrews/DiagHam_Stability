////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of single double option                      //
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


#include "Options/SingleDoubleInternalOption.h"

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
// defaultValue = double default value 
// minValueFlag = flag to indicates an double minimum value
// minValue = double minimum value (no minimum value if greater or equal to maxValue) 
// maxValueFlag = flag to indicates an double maximum value
// maxValue = double maximum value (no maximum value if lower or equal to minValue) 

SingleDoubleInternalOption::SingleDoubleInternalOption(char optionCode, char* optionName, char* optionDescription, double defaultValue, bool external, bool minValueFlag, double minValue, 
						       bool maxValueFlag, double maxValue):
  SingleDoubleOption( optionCode, optionName, optionDescription, defaultValue, 
		      minValueFlag,  minValue, maxValueFlag, maxValue)
{
  this->extFlag=external;
}

// destructor
//
SingleDoubleInternalOption::~SingleDoubleInternalOption()
{
}


// print the current option and its values
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream

ostream& SingleDoubleInternalOption::DisplayOption (ostream& output, bool shortVersion)
{
   if (extFlag) return this->SingleDoubleOption::DisplayOption(output, shortVersion);
  else return output;
}

// print help concerning current option
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& SingleDoubleInternalOption::DisplayHelp (ostream& output)
{
   if (extFlag) return this->SingleDoubleOption::DisplayHelp(output);
  else return output;
}


// Set double 
// set to newVal

void SingleDoubleInternalOption::SetDouble(double newVal)
{
  this->Double = newVal;
}

// multiply with factor
// returns: new value of internal variable
//
double SingleDoubleInternalOption::Multiply(double factor)
{
  this->Double *= factor;
  return this->Double;
}

