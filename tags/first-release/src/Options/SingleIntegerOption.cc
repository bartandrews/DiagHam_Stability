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


#include "Options/SingleIntegerOption.h"

#include <iostream.h>
#include <stdio.h>
#include <string.h>


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

SingleIntegerOption::SingleIntegerOption(char optionCode, char* optionName, char* optionDescription, int defaultValue, 
					 bool minValueFlag, int minValue, 
					 bool maxValueFlag, int maxValue)
{
  this->OptionCode = optionCode;
  this->OptionName = new char [strlen(optionName) + 2];
  this->OptionName[0] = '-';
  strcpy (&(this->OptionName[1]), optionName);
  this->OptionDescription = new char [strlen(optionDescription) + 1];
  strcpy (this->OptionDescription, optionDescription);
  this->MinValueFlag = minValueFlag;
  this->MaxValueFlag = maxValueFlag;
  this->MaxValue = maxValue;
  this->MinValue = minValue;
  this->Integer = defaultValue;
  this->DefaultValue = defaultValue;
}

// destructor
//

SingleIntegerOption::~SingleIntegerOption()
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

int SingleIntegerOption::ReadOption(char** argumentValues, int nbrArgument, int argumentPosition)
{
  char* Argument = argumentValues[argumentPosition];
  if (this->OptionCode == 0)
    if (Argument[0] != '-')
      {
	if (sscanf (Argument, "%d", &(this->Integer)) != 1)
	  {	    
	    this->ErrorCode = SingleIntegerOption::NotAnInteger;
	    return -1;
	  }
	if ((this->MinValueFlag == true) && (this->MinValue > this->Integer))
	  {
	    this->ErrorCode = SingleIntegerOption::Lower;
	    return -1;
	  }
	if ((this->MaxValueFlag == true) && (this->MaxValue < this->Integer))
	  {
	    this->ErrorCode = SingleIntegerOption::Greater;
	    return -1;
	  }	  
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
      if (sscanf (&(Argument[Pos]), "%d", &(this->Integer)) != 1)
	{	    
	  this->ErrorCode = SingleIntegerOption::NotAnInteger;
	  return -1;
	}
      if ((this->MinValueFlag == true) && (this->MinValue > this->Integer))
	{
	  this->ErrorCode = SingleIntegerOption::Lower;
	  return -1;
	}
      if ((this->MaxValueFlag == true) && (this->MaxValue < this->Integer))
	{
	  this->ErrorCode = SingleIntegerOption::Greater;
	  return -1;
	}	  
      return 1;
    }
  if ((argumentPosition + 1) == nbrArgument)
    {
      this->ErrorCode = SingleIntegerOption::NoInteger;
      return -1;
    }
  Argument = argumentValues[argumentPosition + 1];
  if (sscanf (Argument, "%d", &(this->Integer)) != 1)
    {	    
      this->ErrorCode = SingleIntegerOption::NotAnInteger;
      return -1;
    }
  if ((this->MinValueFlag == true) && (this->MinValue > this->Integer))
    {
      this->ErrorCode = SingleIntegerOption::Lower;
      return -1;
    }
  if ((this->MaxValueFlag == true) && (this->MaxValue < this->Integer))
    {
      this->ErrorCode = SingleIntegerOption::Greater;
      return -1;
    }	  
  return 2;
}

// print error message on output stream
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& SingleIntegerOption::PrintError (ostream& output)
{
  if ((this->ErrorCode == SingleIntegerOption::NoError) || (this->OptionCode == 0))
    return output;
  switch (this->ErrorCode)
    {
    case SingleIntegerOption::NotAnInteger:
      {
	output << "option -" << this->OptionName << " (-" << this->OptionCode
	       <<") needs an integer as argument" << endl;
      }
      break;
    case SingleIntegerOption::NoInteger:
      {
	output << "option -" << this->OptionName << " (-" << this->OptionCode
	       <<") needs one integer as argument" << endl;
      }
      break;
    case SingleIntegerOption::Greater:
      {
	output << "option -" << this->OptionName << " (-" << this->OptionCode
	       <<") needs an integer lower than " << this->MaxValue << " as argument" << endl;
      }
      break;
    case SingleIntegerOption::Lower:
      {
	output << "option -" << this->OptionName << " (-" << this->OptionCode
	       <<") needs an integer greater than " << this->MinValue << " as argument" << endl;
      }
      break;
    }
  return output;
}

// Get read integer
//
// return value = integer value

int SingleIntegerOption::GetInteger()
{
  return this->Integer;
}

// print help concerning current option
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& SingleIntegerOption::DisplayHelp (ostream& output)
{
  if ((this->OptionCode != '\0') && (this->OptionCode != '\n'))
    {
      output << "-" << this->OptionCode << ", ";
    }
  output << "-" << this->OptionName << " : " <<  this->OptionDescription << " (";
  if ((this->MinValueFlag == false) && (this->MaxValueFlag == true))
    {
      output << "value must be lower than " << this->MaxValue << ", ";
    }
  else
    if ((this->MinValueFlag == true) && (this->MaxValueFlag == false))
      {
	output << "value must be greater than " << this->MinValue << ", ";
      }
    else
      if ((this->MinValueFlag == true) && (this->MaxValueFlag == true))
	{
	  output << "value must be greater than " << this->MinValue << " and lower than " << this->MaxValue << ", ";
	}
  output << "default value = " << this->DefaultValue << ")";
  return output;
}
