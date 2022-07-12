////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of running option                          //
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


#include "Options/AbstractOption.h"
#include "GeneralTools/ListIterator.h"

#include <iostream>
#include <string.h>


using std::cout;
using std::endl;


// virtual destructor
//

AbstractOption::~AbstractOption()
{
}

// test if a string matches the option name
//
// optionName = string to test
// return value = true if the string matches the option name

bool AbstractOption::IsOptionName (const char* optionName)
{
  if (strlen(&(this->OptionName[1])) != strlen(optionName))
    return false;
  if (strncmp(optionName, &(this->OptionName[1]), strlen(this->OptionName) - 1) != 0)
    {
      return false;
    }
  else
    {
      return true;
    }
}

// Proceed running options 
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// options = option list
// return value = true if proceeding succeded

bool ProceedOptions (char** argumentValues, int nbrArgument, List<AbstractOption*>& option)
{
  int Pos = 1;
  int Inc;
  AbstractOption** TmpOption;
  ListIterator<AbstractOption*> OptionIter(option);
  while (Pos < nbrArgument)
    {
      Inc = 0;
      OptionIter.DefineList(option);
      while ((Inc == 0) && ((TmpOption = OptionIter())) )
	{
	  Inc = (*TmpOption)->ReadOption(argumentValues, nbrArgument, Pos);
	  if (Inc == -1)
	    {
	      (*TmpOption)->PrintError(cout);
	      return false; 
	    }
	}
      if (Inc == 0)
	{
	  cout << "unknown option " <<  argumentValues[Pos] << endl;
	  return false;
	}
      Pos += Inc;
    }
  return true;
}

// Proceed running options from a file 
//
// inputFile = name of the file from which options have to be retrieved
// options = option list
// str = reference on output stream to use if any error occurs while parsing
// return value = true if proceeding succeded

bool ProceedOptions (char* inputFile, List<AbstractOption*>& option, ostream& str)
{
  return true;
}

// parse a line 
//
// line = pointer to the line to parse
// options = option list
// str = reference on output stream to use if any error occurs while parsing
// lineNumber = number of the line which has to be parse
// return value = true if parsing succeded

bool LineParse (char* line, char*& optionName, char*& optionValue)
{
  return true;
}

// display help
//
// options = option list
// str = reference on output stream to use
// return value = true if proceeding succeded

void DisplayHelp (List<AbstractOption*>& option, ostream& str)
{
  AbstractOption** TmpOption;
  ListIterator<AbstractOption*> OptionIter(option);
  str << "usage : " << endl << endl;
  while ((TmpOption = OptionIter()))
    {
      str << "  ";
      (*TmpOption)->DisplayHelp(str) << endl;
    }
}

// display help
//
// options = option list
// str = reference on output stream to use
// programName = string containing the program name
// return value = true if proceeding succeded

void DisplayHelp (List<AbstractOption*>& option, ostream& str, char* programName)
{
  AbstractOption** TmpOption;
  ListIterator<AbstractOption*> OptionIter(option);
  str << "Usage: " << programName << " [options] ";
  while ((TmpOption = OptionIter()))
    {
      if ((*TmpOption)->OptionCode == '\0')
	{
	  str << &((*TmpOption)->OptionName[1]) << endl << endl;
	  str << programName << " has to be called with at least the parameter " << &((*TmpOption)->OptionName[1]) 
	      << ". This parameter can be passed as an argument or using the --" << &((*TmpOption)->OptionName[1]) << " option.";
	}
    }  
  str << endl << endl << "Options:" << endl;
  OptionIter.DefineList(option);
  while ((TmpOption = OptionIter()))
    {
      str << "  ";
      (*TmpOption)->DisplayHelp(str) << endl;
    }
}
