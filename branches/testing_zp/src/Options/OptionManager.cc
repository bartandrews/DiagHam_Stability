////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of option manager                         //
//                                                                            //
//                        last modification : 19/05/2004                      //
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


#include "config.h"
#include "GeneralTools/ListIterator.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/OptionManager.h"
#include "Options/Options.h"

#include <iostream>
#include <string.h>


using std::ostream;
using std::endl;
using std::cout;

// constructor
//
// groupName = group full name

OptionManager::OptionManager(const char* programName, const char* programVersion, const char* programAdditionalInformations)
{
  this->ProgramName = new char [strlen(programName) + 1];
  strcpy (this->ProgramName, programName);      
  if (programVersion == 0)
    {
      this->ProgramVersion = 0;
    }
  else
    {
      this->ProgramVersion = new char [strlen(programVersion) + 1];
      strcpy (this->ProgramVersion, programVersion);      
    }
  if (programAdditionalInformations == 0)
    {
      this->ProgramAdditionalInformations = 0;
    }
  else
    {
      this->ProgramAdditionalInformations = new char [strlen(programAdditionalInformations) + 1];
      strcpy (this->ProgramAdditionalInformations, programAdditionalInformations);      
    }
}

// destructor
//

OptionManager::~OptionManager()
{
  delete[] this->ProgramName;
  if (this->ProgramVersion != 0)
    delete[] this->ProgramVersion;
  if (this->ProgramAdditionalInformations != 0)
    delete[] this->ProgramAdditionalInformations;
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while ((TmpGroup = IterGroup()))
    {
      delete *TmpGroup;
    }
}

// add an option group to the manager
// 
// group = pointer to the option group to add 
// return value = reference on the current option manager

OptionManager& OptionManager::operator += (OptionGroup* group)
{
  this->Groups += group;
  return *this;
}

// get option from its name
//
// optionName = string containing option name
// return value = poitner to the option if it has been found, 0 either

AbstractOption* OptionManager::operator[] (const char* optionName)
{
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  AbstractOption* TmpOption = 0;
  while ((TmpOption == 0) && ((TmpGroup = IterGroup())))
    {
      TmpOption = (**TmpGroup)[optionName];
    }
  return TmpOption; 
}

// get an option group from its name
//
// optionGroupName = string containing option group name
// return value = pointer to the option group if it has been found, 0 either

OptionGroup* OptionManager::GetOptionGroup(const char* optionGroupName)
{
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while ((TmpGroup = IterGroup()))
    {
      if ((*TmpGroup)->IsGroupName(optionGroupName) == true)
	return (*TmpGroup);
    }
  return 0;  
}

// Proceed running options from command line arguments
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// output = reference on output stream used to display errors
// return value = true if proceeding succeded, false if an error occurs

bool OptionManager::ProceedOptions (char** argumentValues, int nbrArgument, ostream& output)
{
  int Pos = 1;
  int Inc;
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while (Pos < nbrArgument)
    {
      Inc = 0;
      IterGroup.DefineList(this->Groups);
      while ((Inc == 0) && ((TmpGroup = IterGroup())) )
	{
	  Inc = (*TmpGroup)->ReadOption(argumentValues, nbrArgument, Pos);
	  if (Inc == -1)
	    {
	      (*TmpGroup)->PrintError(output);
	      return false; 
	    }
	}
      if (Inc == 0)
	{
	  output << "unknown option " <<  argumentValues[Pos] << endl;
	  return false;
	}
      Pos += Inc;
    }
  return true;
}

// ProceedOptions and test if help should be displayed
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// output = reference on output stream used to display errors  
void OptionManager::StandardProceedings(char** argumentValues, int nbrArgument, ostream& output)
{
  if (this->ProceedOptions(argumentValues, nbrArgument, output) == false)
    {
      cout << "see man page for option syntax or type " << this->ProgramName << " -h" << endl;
      exit(-1);
    }
  
  if (this->GetBoolean("help") == true)
    {
      this->DisplayHelp (output);
      exit(0);
    }
}

// print the options and their values in the current group
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// comment = if different from the null character, add it in front of each line
// return value = reference on current output stream

ostream& OptionManager::DisplayOption (ostream& output, bool shortVersion, char comment)
{
  if (shortVersion)
    {
      if (comment != '\0')
	output << comment << " ";
      output << this->ProgramName << "  ";
      ListIterator<OptionGroup*> IterGroup(this->Groups);
      OptionGroup** TmpGroup;
      while ((TmpGroup = IterGroup()))
	{
	  (*TmpGroup)->DisplayOption(output, shortVersion, comment);
	}
      return output;
    }
  else
    {
      if ((this->ProgramVersion != 0) || (this->ProgramAdditionalInformations != 0))
	{
	  if (comment != '\0')
	    output << comment << " ";
	  output << this->ProgramName;
	  if (this->ProgramVersion != 0)
	    {
	      output << ", version " << this->ProgramVersion << endl;
	    }
	  else
	    {
	      output << endl;
	    }
	  if (this->ProgramAdditionalInformations != 0)
	    {
	      if (comment != '\0')
		output << comment << " ";
	      output << this->ProgramAdditionalInformations << endl;
	    }
	}  
      if (comment != '\0')
	output << comment << " ";
      output << endl;
      if (comment != '\0')
	output << comment << " ";
      output << "Options:" << endl;
      if (comment != '\0')
	output << comment << " ";
      cout << endl;
      ListIterator<OptionGroup*> IterGroup(this->Groups);
      OptionGroup** TmpGroup;
      while ((TmpGroup = IterGroup()))
	{
	  (*TmpGroup)->DisplayOption(output, shortVersion, comment);
	  if (comment != '\0')
	    output << comment << " ";	  
	  cout << endl;
	}
      return output;
    }
}

// print help concerning current option group
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& OptionManager::DisplayHelp (ostream& output)
{
  if ((this->ProgramVersion != 0) || (this->ProgramAdditionalInformations != 0))
    {
      output << this->ProgramName;
      if (this->ProgramVersion != 0)
	{
	   output << ", version " << this->ProgramVersion << endl;
	}
      else
	{
	  output << endl;
	}
      if (this->ProgramAdditionalInformations != 0)
	{
	  output << this->ProgramAdditionalInformations << endl;
	}
      output << endl;      
    }  
  output << "Usage: " << this->ProgramName << " [options] " << endl;
  output << endl << "Options:" << endl << endl;
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while ((TmpGroup = IterGroup()))
    {
      (*TmpGroup)->DisplayHelp(output) << endl;
    }
  return output;
}

// dump some of the options into a formatted string
//
// format = string describing the format to use (each symbol %optionname% is replaced by the value associated to the option referred as optionname)
// return value = formatted string (deallocation has to be done manually, 0 if an error occured)

char* OptionManager::GetFormattedString (const char* format)
{
  char* TmpScratch = new char [strlen(format)];
  int StringLength = 0;
  List<char*> Values;
  char* TmpFormat = (char*)format;
  char** TmpValue;
  AbstractOption* TmpOption;
  while ((*TmpFormat) != '\0')
    {
      while (((*TmpFormat) != '%') && ((*TmpFormat) != '\0'))
	{
	  ++TmpFormat;
	  ++StringLength;
	}
      if ((*TmpFormat) != '\0')
	{
	  ++TmpFormat;
	  int TmpSize = 0;
	  while (((*TmpFormat) != '%') && ((*TmpFormat) != '\0'))
	    {
	      TmpScratch[TmpSize] = (*TmpFormat);
	      ++TmpSize;
	      ++TmpFormat;
	    }
	  if ((*TmpFormat) == '\0')
	    {
	      ListIterator<char*> IterValues(Values);
	      while ((TmpValue = IterValues()))
		{
		  delete[] (*TmpValue);
		} 
	      delete[] TmpScratch;
	      return 0l;
	    }
	  TmpScratch[TmpSize] = '\0';
	  TmpOption = (*this)[TmpScratch];
	  if (TmpOption == 0)
	    {
	      ListIterator<char*> IterValues(Values);
	      while ((TmpValue = IterValues()))
		{
		  delete[] (*TmpValue);
		} 
	      delete[] TmpScratch;
	      return 0l;
	    }
	  else
	    {
	      char* TmpString = TmpOption->GetAsAString();
	      StringLength += strlen (TmpString);
	      Values += TmpString;
	    }
	  ++TmpFormat;
	}
    }
  TmpFormat = (char*)format;
  ListIterator<char*> IterValues(Values);
  char* TmpFormattedString = new char [StringLength + 1];
  char* TmpFormattedString2 = TmpFormattedString;
  while ((*TmpFormat) != '\0')
    {
      while (((*TmpFormat) != '%') && ((*TmpFormat) != '\0'))
	{
	  (*TmpFormattedString2) = (*TmpFormat);
	  ++TmpFormattedString2;
	  ++TmpFormat;
	}
      if ((*TmpFormat) != '\0')
	{
	  ++TmpFormat;
	  while ((*TmpFormat) != '%')
	    {
	      ++TmpFormat;	  
	    }
	  ++TmpFormat;
	  TmpValue = IterValues();
	  strcpy (TmpFormattedString2, *TmpValue);
	  TmpFormattedString2 += strlen(*TmpValue);
	  delete[] (*TmpValue);
	}
    }
  delete[] TmpScratch;
  (*TmpFormattedString2) = '\0';
  return TmpFormattedString;
}


//accessor routine for Boolean value
bool OptionManager::GetBoolean(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTBoolean)
	return ((BooleanOption*)OptionPointer)->GetBoolean();
      else
	{
	  cout << "Boolean value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//accessor routine for Double value
double OptionManager::GetDouble(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTDouble)
	return ((SingleDoubleOption*)OptionPointer)->GetDouble();
      else
	{
	  cout << "Double value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//accessor routine for Multiple Double value
double* OptionManager::GetDoubles(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTDoubles)
	return ((MultipleDoubleOption*)OptionPointer)->GetDoubles();
      else
	{
	  cout << "MultipleDouble value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

// alternative accessor routine for Multiple Double value
double* OptionManager::GetDoubles(const char *optionName, int & length)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTDoubles)
	{
	  length = ((MultipleDoubleOption*)OptionPointer)->GetLength();
	  return ((MultipleDoubleOption*)OptionPointer)->GetDoubles();
	}
      else
	{
	  cout << "MultipleDouble value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//accessor routine for Integer value
long OptionManager::GetInteger(const char *optionName)
  {
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTInteger)
	return ((SingleIntegerOption*)OptionPointer)->GetInteger();
      else
	{
	  cout << "Integer value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//accessor routine for Multiple Integer value
int* OptionManager::GetIntegers(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTIntegers)
	return ((MultipleIntegerOption*)OptionPointer)->GetIntegers();
      else
	{
	  cout << "MultipleInteger value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

// alternative accessor routine for Multiple Integer value
int* OptionManager::GetIntegers(const char *optionName, int & length)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTIntegers)
	{
	  length = ((MultipleIntegerOption*)OptionPointer)->GetLength();
	  return ((MultipleIntegerOption*)OptionPointer)->GetIntegers();
	}
      else
	{
	  cout << "MultipleInteger value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}


//accessor routine for String value
char* OptionManager::GetString(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTString)
	return ((SingleStringOption*)OptionPointer)->GetString();
      else
	{
	  cout << "Boolean value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//accessor routine for Multiple String value
char** OptionManager::GetStrings(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTStrings)
	return ((MultipleStringOption*)OptionPointer)->GetStrings();
      else
	{
	  cout << "MultipleString value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

// alternative accessor routine for Multiple String value
char** OptionManager::GetStrings(const char *optionName, int & length)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTStrings)
	{
	  length = ((MultipleStringOption*)OptionPointer)->GetNbrStrings();
	  return ((MultipleStringOption*)OptionPointer)->GetStrings();
	}
      else
	{
	  cout << "MultipleString value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}
