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


#ifndef SINGLEINTEGEROPTION_H
#define SINGLEINTEGEROPTION_H


#include "config.h"
#include "Options/AbstractOption.h"


class SingleIntegerOption : public AbstractOption
{

 protected:
  
  // value assigned to the option
  int Integer;
  
  // default value
  double DefaultValue;

  // maximum value allowed
  int MaxValue;
  // flag indicating there is a maximum value allowed
  bool MaxValueFlag;
  // minimum value allowed
  int MinValue;
  // flag indicating there is a minimum value allowed
  bool MinValueFlag;

 public:

  enum
  {
    NoError = 0,
    NotAnInteger = 1,
    Greater = 2,
    Lower = 3,
    NoInteger = 4
  };

  // constructor from default datas
  //
  // optionCode = character associated to the option
  // optionName = string corresponding to option name
  // optionDescription = string describing option (used for -h option)
  // defaultValue = integer default value 
  // minValueFlag = flag to indicates an integer minimum value
  // minValue = integer minimum value 
  // maxValueFlag = flag to indicates an integer maximum value
  // maxValue = integer maximum value (no maximum value if lower or equal to minValue) 
  SingleIntegerOption(char optionCode, const char* optionName, const char* optionDescription, int defaultValue = 0, 
		      bool minValueFlag = false, int minValue = 0, 
		      bool maxValueFlag = false, int maxValue = 0);

  // destructor
  //
  ~SingleIntegerOption();

  // Test if an argument corresponds to the current option and read its content
  //
  // argumentValues = string array of arguments
  // nbrArgument = number of arguments in argumentValues array
  // argumentPosition = position of the first argument to read
  // return value = number of arguments that have been read (-1 if an error occured)
  int ReadOption(char** argumentValues, int nbrArgument, int argumentPosition);

  // print error message on output stream
  //
  // output = reference on output stream;
  // return value = reference on current output stream
  ostream& PrintError (ostream& output);

  // print the current option and its values
  //  
  // output = reference on output stream;
  // shortVersion = true if return only option code and the option value, false if return option description in addition
  // return value = reference on current output stream
  virtual ostream& DisplayOption (ostream& output, bool shortVersion = false);

  // print help concerning current option
  //
  // output = reference on output stream;
  // return value = reference on current output stream
  ostream& DisplayHelp (ostream& output);

  // Get read integer
  //
  // return value = integer value
  int GetInteger();

  // get option value as a string
  // 
  // return value = corresponding string (deallocation has to be done manually, 0 if an error occured)
  char* GetAsAString();

};

#endif
