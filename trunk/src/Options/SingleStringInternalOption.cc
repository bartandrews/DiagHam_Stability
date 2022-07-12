////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of single string option                       //
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


#include "Options/SingleStringInternalOption.h"

#include <iostream>
#include <string.h>


using std::cout;
using std::endl;


// constructor from default datas
//
// optionCode = character associated to the option
// optionName = string corresponding to option name
// optionDescription = string describing option (used for -h option)
// string = default value of the string (null pointer if none)
// same as SingleStringOption:
SingleStringInternalOption::SingleStringInternalOption(char optionCode, const char* optionName, const char* optionDescription, const char* string, bool external):SingleStringOption(optionCode, optionName, optionDescription, string)
{
  extFlag=external;
}

// destructor
//
// same as SingleStringOption:
SingleStringInternalOption::~SingleStringInternalOption()
{
}
 
// Test if an argument corresponds to the current option and read its content
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// argumentPosition = position of the first argument to read
// return value = number of arguments that have been read (-1 if an error occured)

int SingleStringInternalOption::ReadOption(char** argumentValues, int nbrArgument, int argumentPosition)
{
  return this->SingleStringOption::ReadOption(argumentValues, nbrArgument, argumentPosition);
}

// print error message on output stream
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& SingleStringInternalOption::PrintError (ostream& output)
{
  return this->SingleStringOption::PrintError(output);
}

// Get read string
//
// return value = pointer to string

char* SingleStringInternalOption::GetString()
{
  return this->String;
}

// Set new string
//

void SingleStringInternalOption::SetString(char* NewValue)
{
  if (this->String != NULL)
    {
      delete [] this->String;
      this->String=0;
    }
  this->String = NewValue;
}


// print the current option and its values
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream
// Internal Option is made invisible -> no help output!
ostream& SingleStringInternalOption::DisplayOption (ostream& output, bool shortVersion)
{
  if (extFlag) return this->SingleStringOption::DisplayOption(output, shortVersion);
  else return output;    
}

// print help concerning current option
//
// output = reference on output stream;
// return value = reference on current output stream
// Internal Option is made invisible -> no help output!
ostream& SingleStringInternalOption::DisplayHelp (ostream& output)
{
   if (extFlag) return this->SingleStringOption::DisplayHelp(output);
  else return output;
}
