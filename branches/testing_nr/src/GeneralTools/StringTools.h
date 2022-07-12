////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        various tools to handle strings                     //
//                                                                            //
//                        last modification : 08/08/2008                      //
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


#ifndef STRINGTOOLS_H
#define STRINGTOOLS_H

#include "config.h"

#include <ostream>
using std::ostream;

// split a line using a given separator
//
// string = string to split
// array = reference on the array where elements have to be stored (allocation is done by the method itself, de-allocation has to be done by hand)
// separator = character which is used as separator between columns 
//             if \s (i.e \t or space) is used, then any number of consecutive \s or \t are identify as one separator
// return value = number of elements in the line (zero if an error occured)
int SplitLine(char* string, char**& array, char separator);

// split a line using a given separator and requesting a fixed number of elements
//
// string = string to split
// array = pointer on the array to use to store elements
// nbrElements = number of elements to retrieve 
// separator = character which is used as separator between columns 
//             if \s (i.e \t or space) is used, then any number of consecutive \s or \t are identify as one separator
// return value = number of elements in the line (should be equal if no error occured)
int FixedSplitLine(char* string, char** array, int nbrElements, char separator);

// clean a line from useless comments and spaces
// 
// line = pointer to the line to clean
// return value = true if the line still contains usefull information
bool CleanLine (char* line);


// print the given memory size in b, kb, Mb, or Gb
// str = stream to write to
// bytes = size in bytes
// return = reference on stream
ostream& PrintMemorySize(ostream &str, int bytes);




#endif
