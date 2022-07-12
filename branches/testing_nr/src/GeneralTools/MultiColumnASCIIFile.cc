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


#include "config.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/StringTools.h"

#include <string.h>
#include <fstream>
#include <iostream>


using std::ifstream;
using std::ios;
using std::endl;
using std::cout;


// default constructor
//
// separator = character which is used as separator between integer values in the string 
//             if \s (i.e \t or space) is used, then any number of consecutive \s or \t are identify as one separator

MultiColumnASCIIFile::MultiColumnASCIIFile(char separator)
{
  this->Separator = separator;
  this->Data = 0;
}

// destructor
//

MultiColumnASCIIFile::~MultiColumnASCIIFile()
{
  if (this->Data != 0)
    {
      for (int i = 0; i < this->NbrColumns; ++i)
	{
	  for (int j = 0; j < this->NbrLines; ++j)
	    delete[] this->Data[i][j];
	  delete[] this->Data[i];
	}
      delete[] this->Data;
    }
}

// parse a file
//
// filename = name of the file to parse
// return value = true if no error occurs

bool MultiColumnASCIIFile::Parse(char* filename)
{
  ifstream File;
  File.open(filename, ios::binary | ios::in);
  if (!File.is_open())
    {
      char* TmpString = new char [strlen (filename) + 32]; 
      sprintf (TmpString, "cannot open file: %s\n", filename);
      this->ErrorLog += TmpString;
      return false;
    }
  File.seekg(0, ios::end);
  unsigned int Size = File.tellg();
  File.seekg(0, ios::beg);
  char* TmpBuffer = new char [Size + 1];
  if (TmpBuffer == 0)
    {
      char* TmpString = new char [strlen (filename) + 32]; 
      sprintf (TmpString, "%s is too big\n", filename);
      this->ErrorLog += TmpString;
      return false;
    }
  File.read(TmpBuffer, Size);
  TmpBuffer[Size] = '\0';
  File.close();
  unsigned int Pos = 0;
  int MaxNbrLines = 0;
  while (Pos < Size)
    {
      if (TmpBuffer[Pos] == '\n')
	++MaxNbrLines;
      ++Pos;
    }
  Pos = 0;
  char* Start = TmpBuffer;
  int LineNumber = 1;
  // bool Flag = true;
  this->NbrColumns = 0;
  this->NbrLines = 0;
  char** TmpArray = 0;
  while (Pos < Size)
    {
      Start = TmpBuffer + Pos;
      while ((Pos < Size) && (TmpBuffer[Pos] != '\n'))
	++Pos;
      TmpBuffer[Pos] = '\0';
      if (CleanLine (Start) == true)
	{
	  if (this->NbrColumns == 0)
	    {
	      this->NbrColumns = SplitLine(Start, TmpArray, this->Separator);
	      if (this->NbrColumns <= 0)
		{
		  char* TmpString = new char [strlen (filename) + 256]; 
		  sprintf (TmpString, "fatal error at line %d in file %s : can't retrieve number of columns\n", LineNumber, filename);
		  this->ErrorLog += TmpString;
		  return false;
		}
	      this->Data = new char** [this->NbrColumns];
	      for (int i = 0; i < this->NbrColumns; ++i)
		{
		  this->Data[i] = new char* [MaxNbrLines];
		  this->Data[i][0] = TmpArray[i];
		}	      
	      ++this->NbrLines;
	    }
	  else
	    {
	      if (FixedSplitLine(Start, TmpArray, this->NbrColumns, this->Separator) != this->NbrColumns)
		{
		  char* TmpString = new char [strlen (filename) + 256]; 
		  sprintf (TmpString, "fatal error at line %d in file %s : the number of columns is different from %d\n", LineNumber, filename, this->NbrColumns);
		  this->ErrorLog += TmpString;
		  return false;		  
		}
	      for (int i = 0; i < this->NbrColumns; ++i)
		this->Data[i][this->NbrLines] = TmpArray[i];
	      ++this->NbrLines;
	    }
	}
      ++LineNumber;
      ++Pos;
    }
  delete[] TmpArray;
  if (this->NbrLines != MaxNbrLines)
    {
      for (int i = 0; i < this->NbrColumns; ++i)
	{
	  char** TmpColumn = new char* [this->NbrLines];
	  for (int j = 0; j < this->NbrLines; ++j)
	    TmpColumn[j] = this->Data[i][j];
	  delete[] this->Data[i];
	  this->Data[i] = TmpColumn;
	}
    }
  return true;
}

// retrieve a given value (warning the string is not duplicated, thus should not be modified nor de-allocated)
//
// column = index of the column (zero based)
// line = index of the line (zero based)
// return value = string corresponding to the requested element (null pointer if out of range)

char* MultiColumnASCIIFile::operator ()(int column, int line)
{
  if ((column < 0) || (column >= this->NbrColumns) || (line < 0) || (line >= this->NbrLines))
    return 0;
  else
    return this->Data[column][line];
}

// get a column converting it to integer
//
// column = column index
// return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by han)

int* MultiColumnASCIIFile::GetAsIntegerArray (int column)
{
  int* TmpColumn = new int [this->NbrLines];
  char** TmpASCIIColumn = this->Data[column];
  char* TmpError;
  for (int i = 0; i < this->NbrLines; ++i)
    {
      long Tmp = strtol(TmpASCIIColumn[i], &TmpError, 0);
      if ((TmpError == TmpASCIIColumn[i]) || ((*TmpError) != '\0'))
	{
	  char* TmpString = new char [256 + strlen(TmpASCIIColumn[i])]; 
	  sprintf (TmpString, "element in column %d and line %d is an invalid integer value (%s)\n", column, i, TmpASCIIColumn[i]);
	  this->ErrorLog += TmpString;
	  delete[] TmpColumn;
	  return 0;
	}
      else 
	TmpColumn[i] = (int) Tmp;
    } 
  return TmpColumn;
}

// get a column converting it to long
//
// column = column index
// return value = reference on the array where the read values have to be stored (allocation is done by the method itself,  de-allocation has to be done by hand)

long* MultiColumnASCIIFile::GetAsLongArray (int column)
{
  long* TmpColumn = new long [this->NbrLines];
  char** TmpASCIIColumn = this->Data[column];
  char* TmpError;
  for (int i = 0; i < this->NbrLines; ++i)
    {
      TmpColumn[i] = strtol(TmpASCIIColumn[i], &TmpError, 0);
      if ((TmpError == TmpASCIIColumn[i]) || ((*TmpError) != '\0'))
	{
	  char* TmpString = new char [256 + strlen(TmpASCIIColumn[i])]; 
	  sprintf (TmpString, "element in column %d and line %d is an invalid long integer value (%s)\n", column, i, TmpASCIIColumn[i]);
	  this->ErrorLog += TmpString;
	  delete[] TmpColumn;
	  return 0;
	}
    } 
  return TmpColumn;
}

// get a column converting it to double
//
// column = column index
// return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by han)

double* MultiColumnASCIIFile::GetAsDoubleArray (int column)
{
  double* TmpColumn = new double [this->NbrLines];
  char** TmpASCIIColumn = this->Data[column];
  char* TmpError;
  for (int i = 0; i < this->NbrLines; ++i)
    {
      TmpColumn[i] = strtod(TmpASCIIColumn[i], &TmpError);
      if ((TmpError == TmpASCIIColumn[i]) || ((*TmpError) != '\0'))
	{
	  char* TmpString = new char [256 + strlen(TmpASCIIColumn[i])]; 
	  sprintf (TmpString, "element in column %d and line %d is an invalid double value (%s)\n", column, i, TmpASCIIColumn[i]);
	  this->ErrorLog += TmpString;
	  delete[] TmpColumn;
	  return 0;
	}
    } 
  return TmpColumn;
}

// print last error encountered during parsing operation
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& MultiColumnASCIIFile::PrintLastError (ostream& str)
{
  if (this->ErrorLog.GetNbrElement() > 0)
    str << this->ErrorLog[this->ErrorLog.GetNbrElement() - 1];
  return str;
}


// dump all errors encountered during parsing operation
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& MultiColumnASCIIFile::DumpErrors (ostream& str)
{
  ListIterator<char*> IterError(this->ErrorLog);
  char** TmpError;
  while ((TmpError = IterError()))
    {
      str << *TmpError;
    }
  return str;
}
