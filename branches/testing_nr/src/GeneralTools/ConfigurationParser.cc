////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class that provide parser for configuration files             //
//                                                                            //
//                        last modification : 06/01/2005                      //
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
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/ListIterator.h"

#include <string.h>
#include <fstream>
#include <iostream>


using std::ifstream;
using std::ios;
using std::endl;
using std::cout;


// default constructor
//

ConfigurationParser::ConfigurationParser()
{
}

// destructor
//

ConfigurationParser::~ConfigurationParser()
{
  if (this->ParameterNames.GetNbrElement() > 0)
    {
      ListIterator<char*> IterNames(this->ParameterNames);
      ListIterator<char*> IterValues(this->ParameterValues);
      char** TmpName;
      char** TmpValue;
      while ((TmpName = IterNames()))
	{
	  TmpValue = IterValues();
	  delete[] *TmpName;
	  delete[] *TmpValue;
	}
    }
  if (this->ErrorLog.GetNbrElement() > 0)
    {
      ListIterator<char*> IterError(this->ErrorLog);
      char** TmpError;
      while ((TmpError = IterError()))
	{
	  delete[] (*TmpError);
	}
    }
}

// parse configuration from a file
//
// filename = name of the file to parse
// return value = true if no error occurs

bool ConfigurationParser::Parse(const char* filename)
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
  char* Start = TmpBuffer;
  int LineNumber = 1;
  bool Flag = true;
  while (Pos < Size)
    {
      Start = TmpBuffer + Pos;
      while ((Pos < Size) && (TmpBuffer[Pos] != '\n'))
	{
	  Pos++;
	}
      TmpBuffer[Pos] = '\0';
      if (this->CleanLine (Start) == true)
	{
	  unsigned int Pos2 = 0;
	  while ((Start[Pos2] != '\0') && (Start[Pos2] != '='))
	    Pos2++;
	  if (Start[Pos2] == '\0')
	    {
	      char* TmpString = new char [64]; 
	      sprintf (TmpString, "syntax error at line %d\n", LineNumber);
	      this->ErrorLog += TmpString;
	      Flag = false;
	    }
	  else
	    {
	      unsigned int Pos3 = Pos2;
	      if (Pos2 > 0)
		{
		  --Pos2;
		  while ((Pos2 > 0) && ((Start[Pos2] == ' ') || (Start[Pos2] == '\t')))
		    --Pos2;
		}
	      if ((Pos2 == 0)  && ((Start[Pos2] == ' ') || (Start[Pos2] == '\t') || (Start[Pos2] == '=')))
		{
		  char* TmpString = new char [64]; 
		  sprintf (TmpString, "no parameter name defined at line %d\n", LineNumber);
		  this->ErrorLog += TmpString;
		  Flag = false;
		}
	      else
		{
		  char* TmpName = new char [Pos2 + 2];
		  strncpy (TmpName, Start, Pos2 + 1);
		  TmpName[Pos2 + 1] = '\0';
		  Pos2 = Pos3 + 1;
		  while ((Start[Pos2] != '\0') && ((Start[Pos2] == ' ') || (Start[Pos2] == '\t')))
		    ++Pos2;
		  if (Start[Pos2] == '\0')
		    {
		      char* TmpString = new char [64]; 
		      sprintf (TmpString, "no parameter value defined at line %d\n", LineNumber);
		      this->ErrorLog += TmpString;
		      delete[] TmpName;
		      Flag = false;
		    }
		  else
		    {
		      char* TmpValue = new char [strlen(Start + Pos2) + 1];
		      strcpy (TmpValue, Start + Pos2);
		      this->ParameterNames += TmpName;
		      this->ParameterValues += TmpValue;
		    }
		}
	    }
	}
      Pos++;
      LineNumber++;
    }
  delete[] TmpBuffer;
  return Flag;
}

// get the string corresponding to a configuration parameter 
//
// parameterName = string corresponding to a parameter name
// return value = string (must not be altered) corresponding to a configuration parameter, null if the parameter is not defined

char* ConfigurationParser::operator [] (const char* parameterName)
{
  ListIterator<char*> IterNames(this->ParameterNames);
  ListIterator<char*> IterValues(this->ParameterValues);
  char** TmpName;
  char** TmpValue;
  while ((TmpName = IterNames()))
    {
      TmpValue = IterValues();
      if (strcmp (*TmpName, parameterName) == 0)
	{
	  return *TmpValue;
	}
    }
  char* TmpString = new char [64 + strlen (parameterName)]; 
  sprintf (TmpString, "parameter %s is not defined\n", parameterName);
  this->ErrorLog += TmpString;
  return 0;
}

// dump configuration file
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& ConfigurationParser::DumpConfiguration (ostream& str)
{
  ListIterator<char*> IterNames(this->ParameterNames);
  ListIterator<char*> IterValues(this->ParameterValues);
  char** TmpName;
  char** TmpValue;
  while ((TmpName = IterNames()))
    {
      TmpValue = IterValues();
      str << *TmpName << "=" << *TmpValue << endl;
    }
  return str;
}

// print last error encountered during parsing operation
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& ConfigurationParser::PrintLastError (ostream& str)
{
  if (this->ErrorLog.GetNbrElement() > 0)
    str << this->ErrorLog[this->ErrorLog.GetNbrElement() - 1];
  return str;
}

// dump all errors encountered during parsing operation
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& ConfigurationParser::DumpErrors (ostream& str)
{
  ListIterator<char*> IterError(this->ErrorLog);
  char** TmpError;
  while ((TmpError = IterError()))
    {
      str << *TmpError;
    }
  return str;
}

// clean a line from useless comments and spaces
// 
// line = pointer to the line to clean
// return value = true if the line still contains usefull information

bool ConfigurationParser::CleanLine (char* line)
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
  if (Index == 0)
    return false;
  for (int i = 0; i < Index; ++i)
    line[i] = Start[i];
  line[Index] = '\0';
  return true;
}


// get boolean value corresponding to a configuration parameter 
//
// parameterName = string corresponding to a parameter name
// value = reference on the boolean where the read value has to be stored
// trueString = string which means true (case insensitive)
// falseString = string which means false (case insensitive)
// return value = true if no error occured

bool ConfigurationParser::GetAsBoolean (const char* parameterName, bool& value, const char* trueString, const char* falseString)
{
  char* TmpValue = (*this)[parameterName];
  if (TmpValue == 0)
    return false;
  if (strcasecmp(trueString, TmpValue) == 0)
    {
      value = true;
      return true;
    }
  if (strcasecmp(falseString, TmpValue) == 0)
    {
      value = false;
      return true;
    }
  char* TmpString = new char [128 + strlen(parameterName) + strlen(TmpValue) + strlen (trueString) + strlen(falseString)]; 
  sprintf (TmpString, "parameter %s is fixed to an invalid boolean value (%s). Must be set to %s (true) or %s (false)\n", parameterName, TmpValue,
	   trueString, falseString);
  this->ErrorLog += TmpString;
  return false;
}


// get integer value corresponding to a configuration parameter 
//
// parameterName = string corresponding to a parameter name
// value = reference on the integer where the read value has to be stored
// return value = true if no errro occured

bool ConfigurationParser::GetAsSingleInteger (const char* parameterName, int& value)
{
  char* TmpValue = (*this)[parameterName];
  if (TmpValue == 0)
    return false;
  char* TmpError;
  value = strtol(TmpValue, &TmpError, 0);
  if ((TmpError == TmpValue) || ((*TmpError) != '\0'))
    {
      char* TmpString = new char [64 + strlen(parameterName) + strlen(TmpValue)]; 
      sprintf (TmpString, "parameter %s is fixed to an invalid integer value (%s)\n", parameterName, TmpValue);
      this->ErrorLog += TmpString;
      return false;
    }
  return true;
}

// get integer value corresponding to a configuration parameter 
//
// parameterName = string corresponding to a parameter name
// separator = character which is used as separator between integer values in the string 
//             (if \s is used, then any number of consecutive \s or \t are identify as one separator)
// array = reference on the array where the read values have to be stored (allocation is done by the method itself)
// nbrValues = reference on the integer where the number of read values has to be stored
// return value = true if no errro occured

bool ConfigurationParser::GetAsIntegerArray (const char* parameterName, char separator, int*& array, int& nbrValues)
{
  char* TmpValue = (*this)[parameterName];
  if (TmpValue == 0)
    return false;

  nbrValues = 0;
  if ((separator != ' ')  && (separator != '\t'))
    {
      char* Start = TmpValue;
      while ((*Start) != '\0')
	{
	  if ((*Start) == separator)
	    ++nbrValues;
	  ++Start;
	}
      cout << nbrValues << endl;
      array = new int [nbrValues + 1];
      nbrValues = 0;
      Start = TmpValue;
      char* End = TmpValue;
      char* Error;
      while ((*Start) != '\0')
	{
	  while (((*End) != '\0') && ((*End) != separator))
	    {
	      ++End;
	    }	  
	  array[nbrValues] = strtol (Start, &Error, 0);
	  ++nbrValues;
	  if (Error != End)
	    {
	      char* TmpString = new char [64 + strlen(parameterName) + strlen(TmpValue)]; 
	      sprintf (TmpString, "parameter %s contains an invalid integer value (%s)\n", parameterName, TmpValue);
	      this->ErrorLog += TmpString;
	      return false;
	    }
	  if ((*End) == separator)
	    ++End;
	  while ((((*End) == ' ') || ((*End) == '\t')) && ((*End) != '\0'))
	    {
	      ++End;
	    }
	  Start = End;
	}
    }
  else
    {
      char* Start = TmpValue;
      while ((*Start) != '\0')
	{
	  if (((*Start) == ' ') || ((*Start) == '\t'))
	    {
	      ++nbrValues;
	      while (((*Start) != '\0') && (((*Start) == ' ') || ((*Start) == '\t')))
		++Start;		
	    }
	  else
	    ++Start;
	}      
      array = new int [nbrValues + 1];
      nbrValues = 0;
      Start = TmpValue;
      char* End = TmpValue;
      char* Error;
      while ((*Start) != '\0')
	{
	  while (((*End) != '\0') && ((*End) != ' ') && ((*End) != '\t'))
	    {
	      ++End;
	    }	  
	  array[nbrValues] = strtol (Start, &Error, 0);
	  ++nbrValues;
	  if (Error != End)
	    {
	      char* TmpString = new char [64 + strlen(parameterName) + strlen(TmpValue)]; 
	      sprintf (TmpString, "parameter %s contains an invalid integer value (%s)\n", parameterName, TmpValue);
	      this->ErrorLog += TmpString;
	      return false;
	    }
	  while ((((*End) == ' ') || ((*End) == '\t')) && ((*End) != '\0'))
	    {
	      ++End;
	    }
	  Start = End;
	}
    }
  return true;
}

// get double value corresponding to a configuration parameter 
//
// parameterName = string corresponding to a parameter name
// value = reference on the double where the read value has to be stored
// return value = true if no errro occured

bool ConfigurationParser::GetAsSingleDouble (const char* parameterName, double& value)
{
  char* TmpValue = (*this)[parameterName];
  if (TmpValue == 0)
    return false;
  char* TmpError;
  value = strtod(TmpValue, &TmpError);
  if ((TmpError == TmpValue) || ((*TmpError) != '\0'))
    {
      char* TmpString = new char [64 + strlen(parameterName) + strlen(TmpValue)]; 
      sprintf (TmpString, "parameter %s is fixed to an invalid double value (%s)\n", parameterName, TmpValue);
      this->ErrorLog += TmpString;
      return false;
    }
  return true;
}

// get double value corresponding to a configuration parameter 
//
// parameterName = string corresponding to a parameter name
// separator = character which is used as separator between double values in the string 
//             (if \s is used, then any number of consecutive \s or \t are identify as one separator)
// array = reference on the array where the read values have to be stored (allocation is done by the method itself)
// nbrValues = reference on the double where the number of read values has to be stored
// return value = true if no errro occured

bool ConfigurationParser::GetAsDoubleArray (const char* parameterName, char separator, double*& array, int& nbrValues)
{
  char* TmpValue = (*this)[parameterName];
  if (TmpValue == 0)
    return false;

  nbrValues = 0;
  if ((separator != ' ')  && (separator != '\t'))
    {
      char* Start = TmpValue;
      while ((*Start) != '\0')
	{
	  if ((*Start) == separator)
	    ++nbrValues;
	  ++Start;
	}
      cout << nbrValues << endl;
      array = new double [nbrValues + 1];
      nbrValues = 0;
      Start = TmpValue;
      char* End = TmpValue;
      char* Error;
      while ((*Start) != '\0')
	{
	  while (((*End) != '\0') && ((*End) != separator))
	    {
	      ++End;
	    }	  
	  array[nbrValues] = strtod (Start, &Error);
	  ++nbrValues;
	  if (Error != End)
	    {
	      char* TmpString = new char [64 + strlen(parameterName) + strlen(TmpValue)]; 
	      sprintf (TmpString, "parameter %s contains an invalid double value (%s)\n", parameterName, TmpValue);
	      this->ErrorLog += TmpString;
	      return false;
	    }
	  if ((*End) == separator)
	    ++End;
	  while ((((*End) == ' ') || ((*End) == '\t')) && ((*End) != '\0'))
	    {
	      ++End;
	    }
	  Start = End;
	}
    }
  else
    {
      char* Start = TmpValue;
      while ((*Start) != '\0')
	{
	  if (((*Start) == ' ') || ((*Start) == '\t'))
	    {
	      ++nbrValues;
	      while (((*Start) != '\0') && (((*Start) == ' ') || ((*Start) == '\t')))
		++Start;		
	    }
	  else
	    ++Start;
	}      
      array = new double [nbrValues + 1];
      nbrValues = 0;
      Start = TmpValue;
      char* End = TmpValue;
      char* Error;
      while ((*Start) != '\0')
	{
	  while (((*End) != '\0') && ((*End) != ' ') && ((*End) != '\t'))
	    {
	      ++End;
	    }	  
	  array[nbrValues] = strtod (Start, &Error);
	  ++nbrValues;
	  if (Error != End)
	    {
	      char* TmpString = new char [64 + strlen(parameterName) + strlen(TmpValue)]; 
	      sprintf (TmpString, "parameter %s contains an invalid double value (%s)\n", parameterName, TmpValue);
	      this->ErrorLog += TmpString;
	      return false;
	    }
	  while ((((*End) == ' ') || ((*End) == '\t')) && ((*End) != '\0'))
	    {
	      ++End;
	    }
	  Start = End;
	}
    }
  return true;
}

// get string value corresponding to a configuration parameter 
//
// parameterName = string corresponding to a parameter name
// separator = character which is used as separator between integer values in the string 
//             (if \s is used, then any number of consecutive \s or \t are identify as one separator)
// array = reference on the array where the read values have to be stored (allocation is done by the method itself)
// nbrValues = reference on the integer where the number of read values has to be stored
// return value = true if no errro occured

bool ConfigurationParser::GetAsStringArray (const char* parameterName, char separator, char**& array, int& nbrValues)
{
  char* TmpValue = (*this)[parameterName];
  if (TmpValue == 0)
    return false;

  nbrValues = 0;
  if ((separator != ' ')  && (separator != '\t'))
    {
      char* Start = TmpValue;
      while ((*Start) != '\0')
	{
	  if ((*Start) == separator)
	    ++nbrValues;
	  ++Start;
	}
      array = new char* [nbrValues + 1];
      nbrValues = 0;
      Start = TmpValue;
      char* End = TmpValue;
      // char* Error;
      while ((*Start) != '\0')
	{
	  while (((*End) != '\0') && ((*End) != separator))
	    {
	      ++End;
	    }	  
	  array[nbrValues] = new char [End - Start + 1];
	  strncpy(array[nbrValues], Start, End - Start);
	  array[nbrValues][End - Start] = '\0';
	  ++nbrValues;
	  if ((*End) == separator)
	    ++End;
	  while ((((*End) == ' ') || ((*End) == '\t')) && ((*End) != '\0'))
	    {
	      ++End;
	    }
	  Start = End;
	}
    }
  else
    {
      char* Start = TmpValue;
      while ((*Start) != '\0')
	{
	  if (((*Start) == ' ') || ((*Start) == '\t'))
	    {
	      ++nbrValues;
	      while (((*Start) != '\0') && (((*Start) == ' ') || ((*Start) == '\t')))
		++Start;		
	    }
	  else
	    ++Start;
	}      
      array = new char* [nbrValues + 1];
      nbrValues = 0;
      Start = TmpValue;
      char* End = TmpValue;
      // char* Error;
      while ((*Start) != '\0')
	{
	  while (((*End) != '\0') && ((*End) != ' ') && ((*End) != '\t'))
	    {
	      ++End;
	    }	  
	  array[nbrValues] = new char [End - Start + 1];
	  strncpy(array[nbrValues], Start, End - Start);
	  array[nbrValues][End - Start] = '\0';
	  ++nbrValues;
	  while ((((*End) == ' ') || ((*End) == '\t')) && ((*End) != '\0'))
	    {
	      ++End;
	    }
	  Start = End;
	}
    }
  return true;
}
