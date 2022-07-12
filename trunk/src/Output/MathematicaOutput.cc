////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of Mathematica output                        //
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


#include "Output/MathematicaOutput.h"


// default constructor
//

MathematicaOutput::MathematicaOutput() 
{
}

// destructor
//

MathematicaOutput::~MathematicaOutput() 
{
}

// push a null-ended string on output
//
// Out = reference on current output stream
// c = string pointer
// return value = reference on current output stream

MathematicaOutput& operator << (MathematicaOutput& Out, const char* c) 
{
  Out.File << c;
  return Out;
}

// push a character on output
//
// Out = reference on current output stream
// c = character to push
// return value = reference on current output stream

MathematicaOutput& operator << (MathematicaOutput& Out, const char& c) 
{
  Out.File << c;
  return Out;
}

// push an integer on output
//
// Out = reference on current output stream
// x = integer to push
// return value = reference on current output stream

MathematicaOutput& operator << (MathematicaOutput& Out, const int& x) 
{
  Out.File << x;
  return Out;
}

// push a long on output
//
// Out = reference on current output stream
// x = long to push
// return value = reference on current output stream

MathematicaOutput& operator << (MathematicaOutput& Out, const long& x) 
{
  Out.File << x;
  return Out;
}

// push a double on output
//
// Out = reference on current output stream
// d = double to push
// return value = reference on current output stream

MathematicaOutput& operator << (MathematicaOutput& Out, const double& x) 
{
  Out.File << x;
  return Out;
}

// push a long double on output
//
// Out = reference on current output stream
// d = long double to push
// return value = reference on current output stream

MathematicaOutput& operator << (MathematicaOutput& Out, const long double& x) 
{
  Out.File << x;
  return Out;
}

// push an output stream on output
//
// Out = reference on current output stream
// Str = reference on output stream to push
// return value = reference on current output stream

MathematicaOutput& operator << (MathematicaOutput& Out, ostream& Str) 
{
//  Out.File << Str;
  return Out;
}

// push an Mathematica output on output stream
//
// Str = reference on current output stream
// Out = reference on Mathematica output to push
// return value = reference on current output stream

ostream& operator << (ostream& Str, MathematicaOutput& Out)
{
//  Str << Out.File;
  return Str;
}

// push an Mathematica ouput on output
//
// Out1 = reference on current output stream
// Out2 = Mathematica ouput to push
// return value = reference on current output stream

MathematicaOutput& operator << (MathematicaOutput& Out1, MathematicaOutput& Out2) 
{
#ifdef __SSTREAM_STYLE__
  Out1 << Out2;
#else
  Out1.File.write(Out2.File.str(), Out2.File.pcount());
#endif
  return Out1;
}

