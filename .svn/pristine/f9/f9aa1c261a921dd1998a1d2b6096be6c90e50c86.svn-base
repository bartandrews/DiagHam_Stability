////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of abstract output                         //
//                                                                            //
//                        last modification : 25/01/2001                      //
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


#ifndef ABSTRACTOUTPUT_H
#define ABSTRACTOUTPUT_H


#include <iostream>
#include <strstream>


using std::ostream;
using std::strstream;


class AbstractOutput
{

 protected:

  strstream File;

 public:

  // default constructor
  //
  AbstractOutput();

  // destructor
  //
  virtual ~AbstractOutput();

  // Save output in a file
  //
  // File = string corresponding to file name
  // return value = true if operation succeded
  virtual bool Save(char* File);

  // write a string of char
  //
  // c = pointer to string first character 
  // n = number of character to write
  // return value = reference on current output stream
  virtual AbstractOutput& write(char* c, int n);

  // push a null-ended string on output
  //
  // Out = reference on current output stream
  // c = string pointer
  // return value = reference on current output stream
  friend AbstractOutput& operator << (AbstractOutput& Out, const char* c);

  // push a character on output
  //
  // Out = reference on current output stream
  // c = character to push
  // return value = reference on current output stream
  friend AbstractOutput& operator << (AbstractOutput& Out, const char c);

  // push an integer on output
  //
  // Out = reference on current output stream
  // x = integer to push
  // return value = reference on current output stream
  friend AbstractOutput& operator << (AbstractOutput& Out, const int x);

  // push a double on output
  //
  // Out = reference on current output stream
  // d = double to push
  // return value = reference on current output stream
  friend AbstractOutput& operator << (AbstractOutput& Out, const double x);

  // push an output stream on output
  //
  // Out = reference on current output stream
  // Str = reference on output stream to push
  // return value = reference on current output stream
  friend AbstractOutput& operator << (AbstractOutput& Out, ostream& Str);

  // push an abstract output on output stream
  //
  // Str = reference on current output stream
  // Out = reference on abstract output to push
  // return value = reference on current output stream
  friend ostream& operator << (ostream& Str, AbstractOutput& Out);

  // push an abstract ouput on output
  //
  // Out1 = reference on current output stream
  // Out2 = abstract ouput to push
  // return value = reference on current output stream
  friend AbstractOutput& operator << (AbstractOutput& Out1, AbstractOutput& Out2);

};

#endif
