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
// defaultValue = boolean default value 

BooleanOption::BooleanOption(char optionCode, const char* optionName, const char* optionDescription, bool defaultValue)
{
  this->OptionCode = optionCode;
  this->OptionType = AbstractOption::OTBoolean;
  this->OptionName = new char [strlen(optionName) + 2];
  this->OptionName[0] = '-';
  strcpy (&(this->OptionName[1]), optionName);
  this->OptionDescription = new char [strlen(optionDescription) + 1];
  strcpy (this->OptionDescription, optionDescription);
  this->Boolean = defaultValue;
  this->FalseString=NULL;
  this->TrueString=NULL;
}

// constructor from default datas
//
// optionCode = character associated to the option
// optionName = string corresponding to option name
// optionDescription = string describing option (used for -h option)
// trueString = string output by GetAsString() if optionvalue is true
// falseString = string output by GetAsString() if optionvalue is false
// defaultValue = boolean default value 
BooleanOption::BooleanOption(char optionCode, const char* optionName, const char* optionDescription,
			     const char* trueString, const char* falseString, bool defaultValue)
{
  this->OptionCode = optionCode;
  this->OptionType = AbstractOption::OTBoolean;
  this->OptionName = new char [strlen(optionName) + 2];
  this->OptionName[0] = '-';
  strcpy (&(this->OptionName[1]), optionName);
  this->OptionDescription = new char [strlen(optionDescription) + 1];
  strcpy (this->OptionDescription, optionDescription);
  this->Boolean = defaultValue;
  if (falseString!=NULL)
    {
      this->FalseString=new char [strlen(falseString) + 1];
      strcpy (this->FalseString, falseString);
    }
  else this->FalseString=NULL;
  if (trueString!=NULL)
    {
      this->TrueString=new char [strlen(trueString) + 1];
      strcpy (this->TrueString, trueString);
    }
  else this->TrueString=NULL;
}

// destructor
//

BooleanOption::~BooleanOption()
{
  delete[] this->OptionName;
  delete[] this->OptionDescription;
  if (FalseString!=NULL) delete [] FalseString;
  if (TrueString!=NULL) delete [] TrueString;
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

// get option value as a string
// 
// return value = corresponding string (deallocation has to be done manually, 0 if an error occured)

char* BooleanOption::GetAsAString()
{
  char* TmpString;
  if ((FalseString==NULL)&&(TrueString==NULL))
    {
      TmpString = new char [2];
      if (this->Boolean == false)
	TmpString[0] = '0';
      else
	TmpString[0] = '1';
      TmpString[1] = '\0';
    }
  else
    {
      if (this->Boolean == false)
	{
	  if (FalseString!=NULL)
	    {
	      TmpString = new char[strlen(FalseString)+1];
	      strcpy(TmpString,FalseString);
	    }
	  else
	    {
	      TmpString = new char[1];
	      TmpString[0]='\0';
	    }
	}
      else
	{	  
	  if (TrueString!=NULL)
	    {
	      TmpString = new char[strlen(TrueString)+1];
	      strcpy(TmpString,TrueString);
	    }
	  else
	    {
	      TmpString = new char[1];
	      TmpString[0]='\0';
	    }
	}
    }
  return TmpString;
}

// print the current option and its values
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream

ostream& BooleanOption::DisplayOption (ostream& output, bool shortVersion)
{
  if (shortVersion)
    {
      if (this->Boolean)	
	output << "-" << this->OptionName;      
      return output;     
    }
  else
    {
      if (this->Boolean) 
	output << "-" << this->OptionName << " : " << this->OptionDescription << " : true";
      else
	output << "-" << this->OptionName << " : " << this->OptionDescription << " : false";
      return output;   
    }
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

