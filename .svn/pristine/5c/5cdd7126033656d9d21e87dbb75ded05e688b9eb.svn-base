////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DiagHam version  0.05                           //
//                                                                            //
//                    Copyright (C) 1998-2007 Gunnar Moller                   //
//                                                                            //
//                                                                            //
//                         class of multiple double option                    //
//                                                                            //
//                        last modification : 06/12/2006                      //
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


#include "Options/MultipleDoubleOption.h"

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
// altSeparator = character used to separate entries in output
// minValueFlag = flag to indicates a double minimum value
// minValue = double minimum value (no minimum value if greater or equal to maxValue) 
// maxValueFlag = flag to indicates a double maximum value
// maxValue = double maximum value (no maximum value if lower or equal to minValue) 

MultipleDoubleOption::MultipleDoubleOption(char optionCode, const char* optionName, const char* optionDescription,
					   char separator, char altSeparator, const char* defaultValues,
					   bool minValueFlag, double minValue, 
					   bool maxValueFlag, double maxValue)
{
  this->OptionCode = optionCode;
  this->OptionType = AbstractOption::OTDoubles;
  this->OptionName = new char [strlen(optionName) + 2];
  this->OptionName[0] = '-';
  strcpy (&(this->OptionName[1]), optionName);
  this->OptionDescription = new char [strlen(optionDescription) + 1];
  strcpy (this->OptionDescription, optionDescription);
  this->Separator=separator;
  if (altSeparator != 0) 
    this->AltSeparator = altSeparator;
  else 
    this->AltSeparator = separator;  
  this->Length = 0;
  this->MinValueFlag = minValueFlag;
  this->MaxValueFlag = maxValueFlag;
  this->MaxValue = maxValue;
  this->MinValue = minValue;
  this->Doubles = NULL;
  if (defaultValues!=NULL) this->AnalyzeString(defaultValues);
}

// destructor
//

MultipleDoubleOption::~MultipleDoubleOption()
{
  delete[] this->OptionName;
  delete[] this->OptionDescription;
  if (Doubles!=NULL) 
    delete[] Doubles;
}
 
// Test if an argument corresponds to the current option and read its content
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// argumentPosition = position of the first argument to read
// return value = number of arguments that have been read (-1 if an error occured)

int MultipleDoubleOption::ReadOption(char** argumentValues, int nbrArgument, int argumentPosition)
{
  char* Argument = argumentValues[argumentPosition];
  int val;
  char *String;
  if (this->OptionCode == 0)
    {
      if (Argument[0] != '-')
	{
	  String = new char [strlen (Argument) + 1];
	  strcpy (String, Argument);
	  val=this->AnalyzeString(String);
	  if (val!=0) return val;
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
      String = new char [Lim - Pos + 1];
      strcpy (String, &(Argument[Pos]));
      val=this->AnalyzeString(String);
      if (val!=0) return val;
      delete [] String;
      return 1;
    }
  if ((argumentPosition + 1) == nbrArgument)
    {
      this->ErrorCode = MultipleDoubleOption::NoDouble;
      return -1;
    }
  Argument = argumentValues[argumentPosition + 1];
  String = new char [strlen (Argument) + 1];
  strcpy (String, Argument);
  val=this->AnalyzeString(String);
  if (val!=0) return val;
  delete [] String;
  return 2;
}

// print error message on output stream
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& MultipleDoubleOption::PrintError (ostream& output)
{
  if ((this->ErrorCode == MultipleDoubleOption::NoError) || (this->OptionCode == 0))
    return output;
  switch (this->ErrorCode)
    {
    case MultipleDoubleOption::NotAnDouble:
      {
	output << "option -" << this->OptionName;
	if ((this->OptionCode != '\n') && (this->OptionCode != '\0'))
	  output << " (-" << this->OptionCode <<")";
	output << " needs a list of doubles separated by '" <<this->Separator<<"' as argument" << endl;
      }
      break;
    case MultipleDoubleOption::NoDouble:
      {
	output << "option -" << this->OptionName;
	if ((this->OptionCode != '\n') && (this->OptionCode != '\0'))
	  output << " (-" << this->OptionCode <<")";
	output <<" needs at least one double as argument" << endl;
      }
      break;
    case MultipleDoubleOption::Greater:
      {
	output << "option -" << this->OptionName;
	if ((this->OptionCode != '\n') && (this->OptionCode != '\0'))
	  output << " (-" << this->OptionCode <<")";
	output << " needs doubles lower than " << this->MaxValue << " as argument" << endl;
      }
      break;
    case MultipleDoubleOption::Lower:
      {
	if ((this->OptionCode != '\n') && (this->OptionCode != '\0'))
	  output << " (-" << this->OptionCode <<")";
	output << " needs doubles greater than " << this->MinValue << " as argument" << endl;
      }
      break;
    }
  return output;
}

// Get read double
//
// return value = double value

double* MultipleDoubleOption::GetDoubles()
{
  double *tmp=0;
  if (Length>0)
    {
      tmp= new double[Length];
      for (int i=0;i<Length;++i) tmp[i]=Doubles[i];
    }
  return tmp;
}


// get option value as a string
// 
// return value = corresponding string (deallocation has to be done manually, 0 if an error occured)

char* MultipleDoubleOption::GetAsAString()
{
  char* TmpString;
  if (Length>0)
    {
      TmpString = new char [32*Length];
      sprintf (TmpString, "%.14g", this->Doubles[0]);
      for (int i=1;i<Length;++i) 
	sprintf (TmpString, "%s%c%.14g", TmpString, this->AltSeparator, this->Doubles[i]);
    }
  else
    {
      TmpString= new char[7];
      sprintf(TmpString,"%s","empty");
    }
  return TmpString;
}

// print the current option and its values
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream

ostream& MultipleDoubleOption::DisplayOption (ostream& output, bool shortVersion)
{
  if (shortVersion)
    {
      output << "-" << this->OptionName;
      return output;    
    }
  else
    {
      output << "-" << this->OptionName << " : " << this->OptionDescription;
      return output;   
    }
}

// print help concerning current option
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& MultipleDoubleOption::DisplayHelp (ostream& output)
{
  if ((this->OptionCode != '\0') && (this->OptionCode != '\n'))
    {
      output << "-" << this->OptionCode << ", ";
    }
  output << "-" << this->OptionName << " : " <<  this->OptionDescription;
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
  return output;
}



int MultipleDoubleOption::AnalyzeString(const char *String)
{
  char *tmpC, *tmpC2, *token;
  char sep[2];
  sprintf(sep,"%c",this->Separator);
  double tmp;
  tmpC = new char [strlen(String)+2];
  tmpC2 = new char [strlen(String)+2];
  strcpy (tmpC, String);
  
  //got string value of parameter: now count the number of tokens:
  if (strchr(tmpC,this->Separator)==NULL)
    Length=1;
  else
    {
      token=strtok(tmpC,sep);
      Length=0;
      while (token!=NULL)
	{
	  Length++;
	  token=strtok(NULL,sep);
	}
    }
  // know the count now, but tmpC garbled -> get new copy of string
  if (this->Doubles!=NULL) delete [] this->Doubles;
  this->Doubles=new double[Length];
  strcpy (tmpC, String);
  int n=0;
  while (n<Length)
    {
      if (tmpC[0]==this->Separator)
	{
	  for (unsigned i=0;i<strlen(tmpC);++i) tmpC[i]=tmpC[i+1];
	}
      else //this should now be a number
	{
	  if( sscanf(tmpC,"%lf%s",&tmp,tmpC2) >0)
	    {	      
	      this->Doubles[n++]=tmp;
	      if (n<Length)
		strcpy(tmpC,tmpC2);
	    }
	  else  {	    
	    this->ErrorCode = MultipleDoubleOption::NotAnDouble;
	    return -1;
	  }
	  if ((this->MinValueFlag == true) && (this->MinValue > tmp))
	    {
	      this->ErrorCode = MultipleDoubleOption::Lower;
	      return -1;
	    }
	  if ((this->MaxValueFlag == true) && (this->MaxValue < tmp))
	    {
	      this->ErrorCode = MultipleDoubleOption::Greater;
	      return -1;
	    }	  
	}
    }
  delete[] tmpC;
  delete[] tmpC2;
  return 0;
}
  
