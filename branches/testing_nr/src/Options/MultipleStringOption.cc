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


#include "Options/MultipleStringOption.h"

#include <iostream>
#include <string.h>


using std::cout;
using std::endl;


// constructor from default datas
//
// optionCode = character associated to the option
// optionName = string corresponding to option name
// optionDescription = string describing option (used for -h option)

MultipleStringOption::MultipleStringOption(char optionCode, const char* optionName, const char* optionDescription)
{
  this->OptionCode = optionCode;
  this->OptionType = AbstractOption::OTStrings;
  this->OptionName = new char [strlen(optionName) + 2];
  this->OptionName[0] = '-';
  strcpy (&(this->OptionName[1]), optionName);
  this->OptionDescription = new char [strlen(optionDescription) + 1];
  strcpy (this->OptionDescription, optionDescription);
  this->NbrStrings = 0;
  this->Strings = 0;
  this->ErrorCode = 0;
}

// destructor
//

MultipleStringOption::~MultipleStringOption()
{
  delete[] this->OptionName;
  delete[] this->OptionDescription;
  for (int i=0; i<this->NbrStrings; ++i)
    delete[] this->Strings[i];
  if (this->NbrStrings>0)
    delete[] this->Strings;
}
 
// Test if an argument corresponds to the current option and read its content
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// argumentPosition = position of the first argument to read
// return value = number of arguments that have been read (-1 if an error occured)

int MultipleStringOption::ReadOption(char** argumentValues, int nbrArgument, int argumentPosition)
{
  char* Argument = argumentValues[argumentPosition];
  char* TmpStrings[100];
  int StringsRead=0;
  int ArgumentPos=argumentPosition;
  if (this->OptionCode == 0)
    if (Argument[0] != '-')
      {
	TmpStrings[StringsRead] = new char [strlen (Argument) + 1];
	strcpy (TmpStrings[StringsRead], Argument);
	++ArgumentPos;
	++StringsRead;
	while ((ArgumentPos<nbrArgument)&&(argumentValues[ArgumentPos][0] != '-'))
	  {
	    TmpStrings[StringsRead] = new char [strlen (argumentValues[ArgumentPos]) + 1];
	    strcpy (TmpStrings[StringsRead], argumentValues[ArgumentPos]);
	    ++ArgumentPos;
	    ++StringsRead;
	  }
	char **TmpStrings2=new char*[this->NbrStrings+StringsRead];
	for (int i=0; i<NbrStrings; ++i)
	  TmpStrings2[i]=Strings[i];
	if (this->NbrStrings>0) delete [] this->Strings;
	this->Strings=TmpStrings2;
	for (int i=0; i<StringsRead; ++i,this->NbrStrings++)
	  this->Strings[NbrStrings]=TmpStrings[i];
	return StringsRead;
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
      TmpStrings[StringsRead] = new char [Lim - Pos + 1];
      strcpy (TmpStrings[StringsRead], &(Argument[Pos]));
      ++StringsRead;      
    }
  if ((argumentPosition + 1) == nbrArgument)
    {
      this->ErrorCode = MultipleStringOption::NoString;
      return -1;
    }
  ++ArgumentPos;

  while ((ArgumentPos<nbrArgument)&&(argumentValues[ArgumentPos][0] != '-'))
    {
      TmpStrings[StringsRead] = new char [strlen (argumentValues[ArgumentPos]) + 1];
      strcpy (TmpStrings[StringsRead], argumentValues[ArgumentPos]);
      ++ArgumentPos;
      ++StringsRead;
    }
  char **TmpStrings2=new char*[this->NbrStrings+StringsRead];
  for (int i=0; i<NbrStrings; ++i)
    TmpStrings2[i]=Strings[i];
  if (this->NbrStrings>0) delete [] this->Strings;
  this->Strings=TmpStrings2;
  for (int i=0; i<StringsRead; ++i,this->NbrStrings++)
    this->Strings[NbrStrings]=TmpStrings[i];
  return (ArgumentPos-argumentPosition);
}

// print error message on output stream
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& MultipleStringOption::PrintError (ostream& output)
{
  if ((this->ErrorCode == MultipleStringOption::NoError) || (this->OptionCode == 0))
    return output;
  if (this->ErrorCode == MultipleStringOption::NoString)
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
// index = position of string read
// return value = pointer to string

char* MultipleStringOption::GetString(int index)
{
  if ((index>0)&&(index<NbrStrings))
    return this->Strings[index];
  else
    return 0;
}

// Get array of read string
//
// return value = pointer to string
char** MultipleStringOption::GetStrings()
{
  return this->Strings;
}

// Get array of read string
//
// nbrStrings = returns number of available strings
// return value = pointer to string
char** MultipleStringOption::GetStrings(int &nbrStrings)
{
  nbrStrings = NbrStrings;
  return this->Strings;
}

// get option value as a string
//
// index = position of string read
// return value = corresponding string (deallocation has to be done manually, 0 if an error occured)

char* MultipleStringOption::GetAsAString(int index)
{
  if ((index>0)&&(index<NbrStrings))
    {
      char* TmpString = new char [strlen(this->Strings[index]) + 1];
      strcpy (TmpString, this->Strings[index]);
      return TmpString;
    }
  else
    return 0;
}

// print the current option and its values
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream

ostream& MultipleStringOption::DisplayOption (ostream& output, bool shortVersion)
{
  if (shortVersion)
    {
      output << "-" << this->OptionName;
      for (int i=0; i<NbrStrings; ++i)
	cout << " " <<  this->Strings[i];      
      return output;    
    }
  else
    {
      output << "-" << this->OptionName << " : " << this->OptionDescription << " :";
      for (int i=0; i<NbrStrings; ++i)
	cout << " " <<  this->Strings[i];      
      return output;   
    }
}

// print help concerning current option
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& MultipleStringOption::DisplayHelp (ostream& output)
{
  if ((this->OptionCode != '\0') && (this->OptionCode != '\n'))
    {
      output << "-" << this->OptionCode << ", ";
    }
  output << "-" << this->OptionName << " : " <<  this->OptionDescription << " (multiple arguments)";
  return output;
}
