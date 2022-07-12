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


#ifndef SINGLEDOUBLEOPTION_H
#define SINGLEDOUBLEOPTION_H


#include "config.h"
#include "Options/AbstractOption.h"


class SingleDoubleOption : public AbstractOption
{

 protected:
  
  // value assigned to the option
  double Double;
  
  // default value
  double DefaultValue;

  // maximum value allowed
  double MaxValue;
  // flag indicating there is a maximum value allowed
  bool MaxValueFlag;
  // minimum value allowed
  double MinValue;
  // flag indicating there is a minimum value allowed
  bool MinValueFlag;

  // format option for output
  char *OutputFormat;

 public:

  enum
  {
    NoError = 0,
    NotAnDouble = 1,
    Greater = 2,
    Lower = 3,
    NoDouble = 4
  };

  // constructor from default datas
  //
  // optionCode = character associated to the option
  // optionName = string corresponding to option name
  // optionDescription = string describing option (used for -h option)
  // defaultValue = double default value 
  // minValueFlag = flag to indicates an double minimum value
  // minValue = double minimum value 
  // maxValueFlag = flag to indicates an double maximum value
  // maxValue = double maximum value (no maximum value if lower or equal to minValue) 
  SingleDoubleOption(char optionCode, const char* optionName, const char* optionDescription, double defaultValue = 0, 
		     bool minValueFlag = false, double minValue = 0.0, 
		     bool maxValueFlag = false, double maxValue = 0.0);

  // destructor
  //
  ~SingleDoubleOption();

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

  // Get read double
  //
  // return value = double value
  double GetDouble();

  // get option value as a string
  // 
  // return value = corresponding string (deallocation has to be done manually, 0 if an error occured)
  char* GetAsAString();

  // set output format used by GetAsAString()
  // 
  // format = format for a double, using conventions of printf
  void SetStringFormat(const char *format);

};

#endif
