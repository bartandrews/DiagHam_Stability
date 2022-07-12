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


#include "Output/AbstractOutput.h"

#include <fstream>


using std::ofstream;
using std::ios;


// default constructor
//

AbstractOutput::AbstractOutput() 
{
}

// destructor
//

AbstractOutput::~AbstractOutput() 
{
}

// Save output in a file
//
// FileName = string corresponding to file name
// return value = true if operation succeded

bool AbstractOutput::Save(char* FileName) 
{
  if (FileName == 0)
    return false;
  ofstream RealFile;
  RealFile.open(FileName, ios::binary | ios::out);
//  RealFile << this->File;
  RealFile.write (this->File.str(), this->File.pcount()); 
  RealFile.close();
  return true;
}

// write a string of char
//
// c = pointer to string first character 
// n = number of character to write
// return value = reference on current output stream

AbstractOutput& AbstractOutput::write(char* c, int n) 
{
  this->File.write(c, n);
  return *this;
}

// push a null-ended string on output
//
// Out = reference on current output stream
// c = string pointer
// return value = reference on current output stream

AbstractOutput& operator << (AbstractOutput& Out, const char* c) 
{
  Out.File << c;
  return Out;
}

// push a character on output
//
// Out = reference on current output stream
// c = character to push
// return value = reference on current output stream

AbstractOutput& operator << (AbstractOutput& Out, const char c) 
{
  Out.File << c;
  return Out;
}

// push an integer on output
//
// Out = reference on current output stream
// x = integer to push
// return value = reference on current output stream

AbstractOutput& operator << (AbstractOutput& Out, const int x) 
{
  Out.File << x;
  return Out;
}

// push a double on output
//
// Out = reference on current output stream
// d = double to push
// return value = reference on current output stream

AbstractOutput& operator << (AbstractOutput& Out, const double x) 
{
  Out.File << x;
  return Out;
}

// push an output stream on output
//
// Out = reference on current output stream
// Str = reference on output stream to push
// return value = reference on current output stream

AbstractOutput& operator << (AbstractOutput& Out, ostream& Str) 
{
  Out.File << Str;
  return Out;
}

// push an abstract output on output stream
//
// Str = reference on current output stream
// Out = reference on abstract output to push
// return value = reference on current output stream

ostream& operator << (ostream& Str, AbstractOutput& Out)
{
  Str.write (Out.File.str(), Out.File.pcount());
  return Str;
}

// push an abstract ouput on output
//
// Out1 = reference on current output stream
// Out2 = abstract ouput to push
// return value = reference on current output stream

AbstractOutput& operator << (AbstractOutput& Out1, AbstractOutput& Out2) 
{
  Out1.File.write(Out2.File.str(), Out2.File.pcount());
  return Out1;
}
