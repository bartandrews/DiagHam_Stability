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


#include "config.h"
#include "GeneralTools/StringTools.h"

#include <string.h>

// split a line using a given separator
//
// string = string to split
// array = reference on the array where elements have to be stored (allocation is done by the method itself, de-allocation has to be done by hand)
// separator = character which is used as separator between columns 
//             if \s (i.e \t or space) is used, then any number of consecutive \s or \t are identify as one separator
// return value = number of elements in the line (zero if an error occured)

int SplitLine(char* string, char**& array, char separator)
{
  char* End = string;
  int NbrElements = 1;
  if ((separator == ' ') || (separator == '\t'))
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if (((*End) == ' ') || ((*End) == '\t'))
	  {
	    ++NbrElements;
	    while ((((*End) != '\0') || ((*End) != '\n')) && (((*End) == ' ') || ((*End) == '\t')))
	      ++End;
	  }
	else
	  ++End;
      }
  else
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if ((*End) == separator)
	  ++NbrElements;
	++End;
      }
  array = new char*[NbrElements];  
  NbrElements = 0;
  long TmpLength;
  End = string;
  if ((separator == ' ') || (separator == '\t'))
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if (((*End) == ' ') || ((*End) == '\t'))
	  {
	    TmpLength = End - string;
	    array[NbrElements] = new char [TmpLength + 1];
	    strncpy (array[NbrElements], string, TmpLength);
	    array[NbrElements][TmpLength] = '\0';
	    ++NbrElements;
	    while ((((*End) != '\0') && ((*End) != '\n')) && (((*End) == ' ') || ((*End) == '\t')))
	      ++End;
	    string = End;
	  }
	else
	  ++End;
      }
  else
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if ((*End) == separator)
	  {
	    TmpLength = End - string;
	    array[NbrElements] = new char [TmpLength + 1];
	    strncpy (array[NbrElements], string, TmpLength);
	    array[NbrElements][TmpLength] = '\0';
	    ++NbrElements;
	    string = End + 1;
	  }
	++End;
      }
  TmpLength = End - string;
  array[NbrElements] = new char [TmpLength + 1];
  strncpy (array[NbrElements], string, TmpLength);
  array[NbrElements][TmpLength] = '\0';
  ++NbrElements;	
  return NbrElements;
}

// split a line using a given separator and requesting a fixed number of elements
//
// string = string to split
// array = pointer on the array to use to store elements
// nbrElements = number of elements to retrieve 
// separator = character which is used as separator between columns 
//             if \s (i.e \t or space) is used, then any number of consecutive \s or \t are identify as one separator
// return value = number of elements in the line (should be equal if no error occured)

int FixedSplitLine(char* string, char** array, int nbrElements, char separator)
{
  char* End = string;
  int NbrElements = 0;
  long TmpLength;
  if ((separator == ' ') || (separator == '\t'))
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if (((*End) == ' ') || ((*End) == '\t'))
	  {
	    if (NbrElements < nbrElements)
	      {
		TmpLength = End - string;
		array[NbrElements] = new char [TmpLength + 1];
		strncpy (array[NbrElements], string, TmpLength);
		array[NbrElements][TmpLength] = '\0';
	      }
	    ++NbrElements;
	    while ((((*End) != '\0') && ((*End) != '\n')) && (((*End) == ' ') || ((*End) == '\t')))
	      ++End;
	    string = End;
	  }
	else
	  ++End;
      }
  else
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if ((*End) == separator)
	  {
	    if (NbrElements < nbrElements)
	      {
		TmpLength = End - string;
		array[NbrElements] = new char [TmpLength + 1];
		strncpy (array[NbrElements], string, TmpLength);
		array[NbrElements][TmpLength] = '\0';
	      }
	    ++NbrElements;
	    string = End + 1;
	  }
	++End;
      }
  if (NbrElements < nbrElements)
    {
      TmpLength = End - string;
      array[NbrElements] = new char [TmpLength + 1];
      strncpy (array[NbrElements], string, TmpLength);
      array[NbrElements][TmpLength] = '\0';
    }
  ++NbrElements;	
  return NbrElements;
}

// clean a line from useless comments and spaces
// 
// line = pointer to the line to clean
// return value = true if the line still contains usefull information

bool CleanLine (char* line)
{
  int NbrCharacters = strlen(line);
  if (NbrCharacters == 0)
    return false;
  int Index = 0;
  while ((Index < NbrCharacters) && ((line[Index] == ' ') || (line[Index] == '\t')))
    Index ++;
  if (Index == NbrCharacters)
    return false;
  char* Start = line + Index;
  NbrCharacters -= Index;
  Index = 0;
  while ((Index < NbrCharacters) && ((Start[Index] != '#') || ((Index > 0) && (Start[Index - 1] == '\\'))))
    Index ++;
  while ((Index > 0) && ((line[Index - 1] == ' ') || (line[Index - 1] == '\t')))
    --Index;
  if (Index == 0) 
    return false;
  for (int i = 0; i < Index; ++i)
    line[i] = Start[i];
  line[Index] = '\0';
  return true;
}


// print the given memory size in b, kb, Mb, or Gb
// str = stream to write to
// bytes = size in bytes
// return = reference on stream
ostream& PrintMemorySize(ostream &str, int bytes)
{
  if (bytes >= 1024)
    if (bytes >= 1048576)
      if (bytes >= 1073741824)
	str << (bytes >> 30) << "Go";
      else
	str << (bytes >> 20) << "Mo";
    else
      str << (bytes >> 10) << "ko";
  else
    str << bytes;
  return str;
}

// print the given memory size in b, kb, Mb, or Gb
// str = stream to write to
// bytes = size in bytes
// return = reference on stream
ostream& PrintMemorySize(ostream &str, long bytes)
{
  if (bytes >= 1024)
    if (bytes >= 1048576)
      if (bytes >= 1073741824)
#ifdef __64_BITS__
	if (bytes >= (0x1l<<40))
	  str << (bytes >> 40) << "To";
	else
#endif
	  str << (bytes >> 30) << "Go";
      else
	str << (bytes >> 20) << "Mo";
    else
      str << (bytes >> 10) << "ko";
  else
    str << bytes;
  return str;
}
