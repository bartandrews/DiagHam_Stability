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


#ifndef ABSTRACTOPTION_H
#define ABSTRACTOPTION_H


#include "config.h"
#include "GeneralTools/List.h"

#include <ostream.h>


class ostream;

class AbstractOption
{

 protected:
  
  // standalone character used for otion
  char OptionCode;

  // option full name
  char* OptionName;

  // string describing option (used for -h option)
  char* OptionDescription;

  // error code used dor print error method
  int ErrorCode;

 public:

  // virtual destructor
  //
  virtual ~AbstractOption();

  // Test if an argument corresponds to the current option and read its content
  //
  // argumentValues = string array of arguments
  // nbrArgument = number of arguments in argumentValues array
  // argumentPosition = position of the first argument to read
  // return value = number of arguments that have been read (-1 if an error occured)
  virtual int ReadOption(char** argumentValues, int nbrArgument, int argumentPosition) = 0;

  // print error message on output stream
  //
  // output = reference on output stream;
  // return value = reference on current output stream
  virtual ostream& PrintError (ostream& output) = 0;

  // print help concerning current option
  //
  // output = reference on output stream;
  // return value = reference on current output stream
  virtual ostream& DisplayHelp (ostream& output) = 0;

};

// Proceed running options 
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// options = option list
// return value = true if proceeding succeded
bool ProceedOptions (char** argumentValues, int nbrArgument, List<AbstractOption*>& option);

// display help
//
// options = option list
// str = reference on output stream to use
// return value = true if proceeding succeded
void DisplayHelp (List<AbstractOption*>& option, ostream& str);

#endif
