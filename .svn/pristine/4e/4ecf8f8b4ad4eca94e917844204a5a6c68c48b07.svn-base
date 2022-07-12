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


#include "Options/SingleDoubleOption.h"

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
// minValueFlag = flag to indicates a double minimum value
// minValue = double minimum value (no minimum value if greater or equal to maxValue) 
// maxValueFlag = flag to indicates a double maximum value
// maxValue = double maximum value (no maximum value if lower or equal to minValue) 

SingleDoubleOption::SingleDoubleOption(char optionCode, const char* optionName, const char* optionDescription, double defaultValue, 
				       bool minValueFlag, double minValue, 
				       bool maxValueFlag, double maxValue)
{
  this->OptionCode = optionCode;
  this->OptionType = AbstractOption::OTDouble;
  this->OptionName = new char [strlen(optionName) + 2];
  this->OptionName[0] = '-';
  strcpy (&(this->OptionName[1]), optionName);
  this->OptionDescription = new char [strlen(optionDescription) + 1];
  strcpy (this->OptionDescription, optionDescription);
  this->MinValueFlag = minValueFlag;
  this->MaxValueFlag = maxValueFlag;
  this->MaxValue = maxValue;
  this->MinValue = minValue;
  this->Double = defaultValue;
  this->DefaultValue = defaultValue;
  this->OutputFormat = new char[7];
  strcpy(this->OutputFormat,"%.14f");
}

// destructor
//

SingleDoubleOption::~SingleDoubleOption()
{
  delete[] this->OptionName;
  delete[] this->OptionDescription;
  delete[] this->OutputFormat;
}
 
// Test if an argument corresponds to the current option and read its content
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// argumentPosition = position of the first argument to read
// return value = number of arguments that have been read (-1 if an error occured)

int SingleDoubleOption::ReadOption(char** argumentValues, int nbrArgument, int argumentPosition)
{
  char* Argument = argumentValues[argumentPosition];
  if (this->OptionCode == 0)
    {
      if (Argument[0] != '-')
	{
	  if (sscanf (Argument, "%lf", &(this->Double)) != 1)
	    {	    
	      this->ErrorCode = SingleDoubleOption::NotAnDouble;
	      return -1;
	    }
	  if ((this->MinValueFlag == true) && (this->MinValue > this->Double))
	    {
	      this->ErrorCode = SingleDoubleOption::Lower;
	      return -1;
	    }
	  if ((this->MaxValueFlag == true) && (this->MaxValue < this->Double))
	    {
	      this->ErrorCode = SingleDoubleOption::Greater;
	      return -1;
	    }	  
	  return 1;
	}
      else
	{
	  return 0;
	}
    }

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
      if (sscanf (&(Argument[Pos]), "%lf", &(this->Double)) != 1)
	{	    
	  this->ErrorCode = SingleDoubleOption::NotAnDouble;
	  return -1;
	}
      if ((this->MinValueFlag == true) && (this->MinValue > this->Double))
	{
	  this->ErrorCode = SingleDoubleOption::Lower;
	  return -1;
	}
      if ((this->MaxValueFlag == true) && (this->MaxValue < this->Double))
	{
	  this->ErrorCode = SingleDoubleOption::Greater;
	  return -1;
	}	  
      return 1;
    }
  if ((argumentPosition + 1) == nbrArgument)
    {
      this->ErrorCode = SingleDoubleOption::NoDouble;
      return -1;
    }
  Argument = argumentValues[argumentPosition + 1];
  if (sscanf (Argument, "%lf", &(this->Double)) != 1)
    {	    
      this->ErrorCode = SingleDoubleOption::NotAnDouble;
      return -1;
    }

  if ((this->MinValueFlag == true) && (this->MinValue > this->Double))
    {
      this->ErrorCode = SingleDoubleOption::Lower;
      return -1;
    }
  if ((this->MaxValueFlag == true) && (this->MaxValue < this->Double))
    {
      this->ErrorCode = SingleDoubleOption::Greater;
      return -1;
    }	  
  return 2;
}

// print error message on output stream
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& SingleDoubleOption::PrintError (ostream& output)
{
  if ((this->ErrorCode == SingleDoubleOption::NoError) || (this->OptionCode == 0))
    return output;
  switch (this->ErrorCode)
    {
    case SingleDoubleOption::NotAnDouble:
      {
	output << "option -" << this->OptionName;
	if ((this->OptionCode != '\n') && (this->OptionCode != '\0'))
	  output << " (-" << this->OptionCode <<")";
	output << " needs a double as argument" << endl;
      }
      break;
    case SingleDoubleOption::NoDouble:
      {
	output << "option -" << this->OptionName;
	if ((this->OptionCode != '\n') && (this->OptionCode != '\0'))
	  output << " (-" << this->OptionCode <<")";
	output <<" needs one double as argument" << endl;
      }
      break;
    case SingleDoubleOption::Greater:
      {
	output << "option -" << this->OptionName;
	if ((this->OptionCode != '\n') && (this->OptionCode != '\0'))
	  output << " (-" << this->OptionCode <<")";
	output << " needs a double lower than " << this->MaxValue << " as argument" << endl;
      }
      break;
    case SingleDoubleOption::Lower:
      {
	if ((this->OptionCode != '\n') && (this->OptionCode != '\0'))
	  output << " (-" << this->OptionCode <<")";
	output << " needs a double greater than " << this->MinValue << " as argument" << endl;
      }
      break;
    }
  return output;
}

// Get read double
//
// return value = double value

double SingleDoubleOption::GetDouble()
{
  return this->Double;
}

// get option value as a string
// 
// return value = corresponding string (deallocation has to be done manually, 0 if an error occured)

char* SingleDoubleOption::GetAsAString()
{
  char* TmpString = new char [32];
  sprintf (TmpString, this->OutputFormat, this->Double);
  return TmpString;
}

// set output format used by GetAsAString()
// 
// format = format for a double, using conventions of printf
void SingleDoubleOption::SetStringFormat(const char *format)
{
  delete [] this->OutputFormat;
  this->OutputFormat = new char[strlen(format)+1];
  strcpy(this->OutputFormat,format);
}

// print the current option and its values
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream

ostream& SingleDoubleOption::DisplayOption (ostream& output, bool shortVersion)
{
  if (shortVersion)
    {
      output << "-" << this->OptionName << " " <<  this->Double;
      return output;    
    }
  else
    {
      output << "-" << this->OptionName << " : " << this->OptionDescription << " : " << this->Double;
      return output;   
    }
}

// print help concerning current option
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& SingleDoubleOption::DisplayHelp (ostream& output)
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
