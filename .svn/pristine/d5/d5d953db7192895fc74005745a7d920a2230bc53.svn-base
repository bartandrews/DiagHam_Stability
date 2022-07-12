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


#ifndef BOOLEANINTERNALOPTION_H
#define BOOLEANINTERNALOPTION_H


#include "config.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include <string>

class BooleanInternalOption : public BooleanOption
{

 protected:
  
  bool ExtFlag;

  
 public:

  // constructor from default datas
  //
  // optionCode = character associated to the option
  // optionName = string corresponding to option name
  // optionDescription = string describing option (used for -h option)
  // defaultValue = boolean default value
  // external = flag to indicate whether visible to the external world (in help output);
  BooleanInternalOption(char optionCode, const char* optionName, const char* optionDescription, bool defaultValue = false, bool external=false);

  // constructor from default datas
  //
  // optionCode = character associated to the option
  // optionName = string corresponding to option name
  // optionDescription = string describing option (used for -h option)
  // trueString = string output by Manager.GetFormattedString if optionvalue is true
  // falseString = string output by Manager.GetFormattedString if optionvalue is false
  // defaultValue = boolean default value
  // external = flag to indicate whether visible to the external world (in help output);
  BooleanInternalOption(char optionCode, const char* optionName, const char* optionDescription,
			const char* trueString, const char* falseString, bool defaultValue = false, bool external=false);

  // destructor
  //
  ~BooleanInternalOption();

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

  // Set Boolean value
  // set to newVal
  // 
  void SetBoolean(bool newVal);


};

#endif
