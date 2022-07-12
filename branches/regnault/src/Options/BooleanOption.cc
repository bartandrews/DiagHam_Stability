////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of boolean option                         //
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


#include "Options/BooleanOption.h"

#include <iostream.h>
#include <stdio.h>
#include <string.h>


// constructor from default datas
//
// optionCode = character associated to the option
// optionName = string corresponding to option name
// optionDescription = string describing option (used for -h option)
// defaultValue = boolean default value 

BooleanOption::BooleanOption(char optionCode, char* optionName, char* optionDescription, bool defaultValue)
{
  this->OptionCode = optionCode;
  this->OptionName = new char [strlen(optionName) + 2];
  this->OptionName[0] = '-';
  strcpy (&(this->OptionName[1]), optionName);
  this->OptionDescription = new char [strlen(optionDescription) + 1];
  strcpy (this->OptionDescription, optionDescription);
  this->Boolean = defaultValue;
}

// destructor
//

BooleanOption::~BooleanOption()
{
  delete[] this->OptionName;
  delete[] this->OptionDescription;
}
 
// Test if an argument corresponds to the current option and read its content
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// argumentPosition = position of the first argument to read
// return value = number of arguments that have been read (-1 if an error occured)

int BooleanOption::ReadOption(char** argumentValues, int nbrArgument, int argumentPosition)
{
  char* Argument = argumentValues[argumentPosition];
  if (this->OptionCode == 0)
    {
      return 0;
    }
  if (Argument[0] != '-')
    return 0;
  if ((Argument[1] != this->OptionCode) && (strcmp(&(Argument[1]), this->OptionName) != 0))
    return 0;
  this->Boolean = true;
  return 1;
}

// print error message on output stream
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& BooleanOption::PrintError (ostream& output)
{
  if ((this->ErrorCode == BooleanOption::NoError) || (this->OptionCode == 0))
    return output;
  return output;
}

// Get read boolean
//
// return value = boolean value

bool BooleanOption::GetBoolean()
{
  return this->Boolean;
}

// print error message on output stream
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& BooleanOption::DisplayHelp (ostream& output)
{
  if ((this->OptionCode != '\0') && (this->OptionCode != '\n'))
    {
      output << "-" << this->OptionCode << ", ";
    }
  output << "-" << this->OptionName << " : " <<  this->OptionDescription;
  return output;
}

