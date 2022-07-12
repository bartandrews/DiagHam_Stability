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


#ifndef OPTIONMANAGER_H
#define OPTIONMANAGER_H


#include "config.h"
#include "GeneralTools/List.h"

#include <iostream>


using std::ostream;


class OptionGroup;
class AbstractOption;

class OptionManager
{

 protected:
  
  // ordered list containing all option groups which have to be handled by the option manager
  List<OptionGroup*> Groups;

  // program name 
  char* ProgramName;
 
  // string containing program version
  char* ProgramVersion;

  // an optional string that can be displayed with help
  char* ProgramAdditionalInformations;

 public:

  // constructor
  //
  // groupName = group full name
  OptionManager(const char* programName, const char* programVersion = 0, const char* programAdditionalInformations = 0);

  // destructor
  //
  ~OptionManager();

  // add an option group to the manager
  // 
  // group = pointer to the option group to add 
  // return value = reference on the current option manager
  OptionManager& operator += (OptionGroup* group);

  // get option from its name
  //
  // optionName = string containing option name
  // return value = pointer to the option if it has been found, 0 either
  AbstractOption* operator[] (const char* optionName);

  // get an option group from its name
  //
  // optionGroupName = string containing option group name
  // return value = pointer to the option group if it has been found, 0 either
  OptionGroup* GetOptionGroup(const char* optionGroupName);

  // Proceed running options from command line arguments
  //
  // argumentValues = string array of arguments
  // nbrArgument = number of arguments in argumentValues array
  // output = reference on output stream used to display errors
  // return value = true if proceeding succeded, false if an error occurs
  bool ProceedOptions (char** argumentValues, int nbrArgument, ostream& output);

  // ProceedOptions and test if help should be displayed
  //
  // argumentValues = string array of arguments
  // nbrArgument = number of arguments in argumentValues array
  // output = reference on output stream used to display errors  
  void StandardProceedings(char** argumentValues, int nbrArgument, ostream& output);

  // print the options and their values in the current group
  //  
  // output = reference on output stream;
  // shortVersion = true if return only option code and the option value, false if return option description in addition
  // comment = if different from the null character, add it in front of each line
  // return value = reference on current output stream
  ostream& DisplayOption (ostream& output, bool shortVersion = false, char comment = '\0');

  // print help concerning current option group
  //
  // output = reference on output stream;
  // return value = reference on current output stream
  ostream& DisplayHelp (ostream& output);

  // dump some of the options into a formatted string
  //
  // format = string describing the format to use (each symbol %optionname% is replaced by the value associated to the option referred as optionname)
  // return value = formatted string (deallocation has to be done manually)
  char* GetFormattedString (const char* format);

  // accessor routines for managed options:
  // these are wrappers of the operator[] routine
  // returning the value of the option, if found
  // or exiting with a warning message, otherwise  
  bool GetBoolean(const char *optionName);
  double GetDouble(const char *optionName);
  // from MultipleDoubleOption: possibility to request the length of the vector, if wanted
  double* GetDoubles(const char *optionName);
  double* GetDoubles(const char *optionName, int &length);
  long GetInteger(const char *optionName);
  // from MultipleIntegerOption: possibility to request the length of the vector, if wanted
  int* GetIntegers(const char *optionName);
  int* GetIntegers(const char *optionName, int &length);
  char* GetString(const char *optionName);
  //accessor routine for Multiple String value
  char** GetStrings(const char *optionName);
  // alternative accessor routine for Multiple String value
  char** GetStrings(const char *optionName, int & length);


  
};

#endif
