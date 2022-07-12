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


#include "Options/SingleStringOption.h"

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

SingleStringOption::SingleStringOption(char optionCode, const char* optionName, const char* optionDescription, const char* string)
{
  this->OptionCode = optionCode;
  this->OptionType = AbstractOption::OTString;
  this->OptionName = new char [strlen(optionName) + 2];
  this->OptionName[0] = '-';
  strcpy (&(this->OptionName[1]), optionName);
  this->OptionDescription = new char [strlen(optionDescription) + 1];
  strcpy (this->OptionDescription, optionDescription);
  if (string != 0)
    {
      this->String = new char [strlen(string) + 1];
      strcpy (this->String, string);
    }
  else
    this->String = 0;
  this->ErrorCode = 0;
}

// destructor
//

SingleStringOption::~SingleStringOption()
{
  delete[] this->OptionName;
  delete[] this->OptionDescription;
  if (this->String != 0)
    delete[] this->String;
}
 
// Test if an argument corresponds to the current option and read its content
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// argumentPosition = position of the first argument to read
// return value = number of arguments that have been read (-1 if an error occured)

int SingleStringOption::ReadOption(char** argumentValues, int nbrArgument, int argumentPosition)
{
  char* Argument = argumentValues[argumentPosition];
  if (this->OptionCode == 0)
    if (Argument[0] != '-')
      {
	this->String = new char [strlen (Argument) + 1];
	strcpy (this->String, Argument);
	return 1;
      }
    else
      return 0;
  if (Argument[0] != '-')
    return 0;
  if ((Argument[1] != this->OptionCode) && (strncmp(&(Argument[1]), this->OptionName, 
						    strlen (this->OptionName)) != 0))
    return 0;
  int Pos = 2;
  if (strncmp(&(Argument[1]), this->OptionName, strlen (this->OptionName)) == 0)
    Pos = 1 + strlen(this->OptionName);
  int Lim = strlen (Argument);
  while ((Pos < Lim) && (Argument[Pos] == ' '))
    Pos++;
  if (Pos != Lim)
    {
      this->String = new char [Lim - Pos + 1];
      strcpy (this->String, &(Argument[Pos]));
      return 1;
    }
  if ((argumentPosition + 1) == nbrArgument)
    {
      this->ErrorCode = SingleStringOption::NoString;
      return -1;
    }
  Argument = argumentValues[argumentPosition + 1];
  this->String = new char [strlen (Argument) + 1];
  strcpy (this->String, Argument);
  return 2;
}

// print error message on output stream
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& SingleStringOption::PrintError (ostream& output)
{
  if ((this->ErrorCode == SingleStringOption::NoError) || (this->OptionCode == 0))
    return output;
  if (this->ErrorCode == SingleStringOption::NoString)
    {
	output << "option -" << this->OptionName;
	if ((this->OptionCode != '\n') && (this->OptionCode != '\0'))
	  output << " (-" << this->OptionCode <<")";
	output << " needs one string as argument" << endl;
      return output;
    }
  return output;
}

// Get read string
//
// return value = pointer to string

char* SingleStringOption::GetString()
{
  return this->String;
}

// get option value as a string
// 
// return value = corresponding string (deallocation has to be done manually, 0 if an error occured)

char* SingleStringOption::GetAsAString()
{
  char* TmpString = new char [strlen(this->String) + 1];
  strcpy (TmpString, this->String);
  return TmpString;
}

// print the current option and its values
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream

ostream& SingleStringOption::DisplayOption (ostream& output, bool shortVersion)
{
  if (this->String != 0)
    {
      if (shortVersion)
	{
	  output << "-" << this->OptionName << " " <<  this->String;
	  return output;    
	}
      else
	{
	  output << "-" << this->OptionName << " : " << this->OptionDescription << " : " << this->String;
	  return output;   
	}
    }
  return output; 
}

// print help concerning current option
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& SingleStringOption::DisplayHelp (ostream& output)
{
  if ((this->OptionCode != '\0') && (this->OptionCode != '\n'))
    {
      output << "-" << this->OptionCode << ", ";
    }
  output << "-" << this->OptionName << " : " <<  this->OptionDescription;
  return output;
}
