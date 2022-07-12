////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class that provide multicolumn ascii file                //
//                                                                            //
//                        last modification : 31/12/2006                      //
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


#ifndef MULTICOLUMNSCIIFILE_H
#define MULTICOLUMNSCIIFILE_H

#include "config.h"
#include "GeneralTools/List.h"
#include "MathTools/Complex.h"

#include <iostream>


using std::ostream;


class Rational;
class LongRational;


class MultiColumnASCIIFile
{

 protected:

  // number of columns
  int NbrColumns;

  // number of lines (must be the same for each column)
  int NbrLines;

  // column separator
  char Separator;

  // data columns
  char*** Data;

  // list of all error messages 
  List<char*> ErrorLog;

 public:

  // default constructor
  //
  // separator = character which is used as separator between columns 
   //             if \s (i.e \t or space) is used, then any number of consecutive \s or \t are identified as one separator
  MultiColumnASCIIFile(char separator = ' ');

  // destructor
  //
  ~MultiColumnASCIIFile();

  // parse a file
  //
  // filename = name of the file to parse
  // return value = true if no error occurs
  bool Parse(char* filename);

  // get the number of columns
  //
  // return value = number of columns
  int GetNbrColumns();

  // number of lines (same for each column)
  //
  // return value = number of lines
  int GetNbrLines();

  // retrieve a given value (warning the string is not duplicated, thus should not be modified nor de-allocated)
  //
  // column = index of the column (zero based)
  // line = index of the line (zero based)
  // return value = string corresponding to the requested element (null pointer if out of range)
  char* operator ()(int column, int line);


  // get a column converting it to strings
  //
  // column = column index
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by hand)
  char** GetAsStringArray (int column);

  // get a column converting it to integer
  //
  // column = column index
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by hand)
  int* GetAsIntegerArray (int column);

  // get a column converting it to long
  //
  // column = column index
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself,  de-allocation has to be done by hand)
  long* GetAsLongArray (int column);

  // get a column converting it to double
  //
  // column = column index
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by hand)
  double* GetAsDoubleArray (int column);

  // get a column converting it to complex
  //
  // column = column index
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by hand)
  Complex* GetAsComplexArray (int column);

  // get a column converting it to rational
  //
  // column = column index
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself,  de-allocation has to be done by hand)
  Rational* GetAsRationalArray (int column);

  // get a column converting it to long rational
  //
  // column = column index
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself,  de-allocation has to be done by hand)
  LongRational* GetAsLongRationalArray (int column);

  
  // get a line converting it to integer
  //
  // line = line index
  // firstColumn = index of the first column to store
  // lastColumn = index of the last column to store (go up to last column if lastColumn <= firstColumn)
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by hand)
  int* GetLineAsIntegerArray (int line, int firstColumn = 0, int lastColumn = 0);

  // get a line converting it to long
  //
  // line = line index
  // firstColumn = index of the first column to store
  // lastColumn = index of the last column to store (go up to last column if lastColumn <= firstColumn)
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself,  de-allocation has to be done by hand)
  long* GetLineAsLongArray (int line, int firstColumn = 0, int lastColumn = 0);

  // get a line converting it to double
  //
  // line = line index
  // firstColumn = index of the first column to store
  // lastColumn = index of the last column to store (go up to last column if lastColumn <= firstColumn)
  // return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by hand)
  double* GetLineAsDoubleArray (int line, int firstColumn = 0, int lastColumn = 0);

  // print last error encountered during parsing operation
  //
  // str = reference on the output stream
  // return value = reference on the output stream
  ostream& PrintLastError (ostream& str);

  // dump all errors encountered during parsing operation
  //
  // str = reference on the output stream
  // return value = reference on the output stream
  ostream& DumpErrors (ostream& str);

};

// get the number of columns
//
// return value = number of columns

inline int MultiColumnASCIIFile::GetNbrColumns()
{
  return this->NbrColumns;
}


// number of lines (same for each column)
//
// return value = number of lines

inline int MultiColumnASCIIFile::GetNbrLines()
{
  return this->NbrLines;
}

#endif
