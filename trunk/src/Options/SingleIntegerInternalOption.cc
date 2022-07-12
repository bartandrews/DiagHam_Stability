////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of single integer option                      //
//                                                                            //
//                        last modification : 19/08/2001                      //
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


#include "Options/SingleIntegerInternalOption.h"

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
// defaultValue = integer default value 
// minValueFlag = flag to indicates an integer minimum value
// minValue = integer minimum value (no minimum value if greater or equal to maxValue) 
// maxValueFlag = flag to indicates an integer maximum value
// maxValue = integer maximum value (no maximum value if lower or equal to minValue) 

SingleIntegerInternalOption::SingleIntegerInternalOption(char optionCode, const char* optionName,
							 const char* optionDescription,
							 int defaultValue, bool external, 
							 bool minValueFlag, int minValue, 
							 bool maxValueFlag, int maxValue):
  SingleIntegerOption( optionCode,  optionName,  optionDescription,
		       defaultValue,  minValueFlag,  minValue, 
		       maxValueFlag,  maxValue)
    
{
  this->extFlag=external;
}

// destructor
//

SingleIntegerInternalOption::~SingleIntegerInternalOption()
{
}
 

// Get read integer
//
// return value = integer value

int SingleIntegerInternalOption::GetInteger()
{
  return this->Integer;
}

// Set new integer
//

void SingleIntegerInternalOption::SetInteger(int newValue)
{
  this->Integer=newValue;
}

void SingleIntegerInternalOption::Increase()
{
  this->Integer++;
}


// print the current option and its values
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream

ostream& SingleIntegerInternalOption::DisplayOption (ostream& output, bool shortVersion)
{
  if (extFlag) return this->SingleIntegerOption::DisplayOption(output, shortVersion);
  else return output;    
}

// print help concerning current option
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& SingleIntegerInternalOption::DisplayHelp (ostream& output)
{
  if (extFlag) return this->SingleIntegerOption::DisplayHelp(output);
  else return output;    
}
