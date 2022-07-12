////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             various generic tools about file name manipulations            //
//                                                                            //
//                        last modification : 23/02/2006                      //
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
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ArrayTools.h"

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>


using std::cout;
using std::endl;
using std::ifstream;


// list all files or directories that obey a given pattern (that can include relative/absolute path) /to/directory/patternxxxsuffix where xxx is an integer
//
// pattern = string that corresponds to the pattern  (i.e. /to/directory/pattern)
// matchedFileArray = reference on the sorted array (with respect to xxx) of files or directories names (with the optional relative/absolute path), 
//                    memory allocation isd one by the function itself
// suffix = optional suffix  to test
// return value = number of matched files

int GetAllFilesDirectories(const char* pattern, char**& matchedFileArray, const char* suffix)
{
  const char* Path = strrchr(pattern, '/');
  char* TmpPattern;
  long PatternLength;
  DIR* TmpDirectory;
  long PathLength = 0;
  if (Path == 0)
    {
      TmpDirectory = opendir(".");
      PatternLength = strlen(pattern);
      TmpPattern = new char [PatternLength + 1];
      strcpy(TmpPattern, pattern);
    }
  else
    {
      PathLength = (Path - pattern) + 1;
      PatternLength = strlen(Path) - 1;
      TmpPattern = new char [PatternLength + 1];
      strcpy(TmpPattern, Path + 1); 
      char* TmpPath = new char [PathLength + 1];
      strncpy (TmpPath, pattern, PathLength);
      TmpPath[PathLength] = '\0';
      TmpDirectory = opendir(TmpPath); 
      delete[] TmpPath;
    }
  dirent* DirectoryContent;
  List<char*> MatchedFiles;
  if (suffix == 0)
    {
      while ((DirectoryContent = readdir(TmpDirectory)))
	{
	  if ((strncmp(DirectoryContent->d_name, TmpPattern, PatternLength) == 0) && ((*(DirectoryContent->d_name + PatternLength)) >= '0') && 
	  ((*(DirectoryContent->d_name + PatternLength)) <= '9'))
	    {
	      char* TmpName  = new char [strlen(DirectoryContent->d_name) + 1];
	      strcpy (TmpName, DirectoryContent->d_name);
	      MatchedFiles += TmpName;
	    }
	}
    }
  else
    {
      while ((DirectoryContent = readdir(TmpDirectory)))
	{
	  if (strncmp(DirectoryContent->d_name, TmpPattern, PatternLength) == 0)
	    {
	      char* EndPos = DirectoryContent->d_name + PatternLength;
	      while (((*EndPos) != '\0') && ((*EndPos) >= '0') && ((*EndPos) <= '9'))
		++EndPos;
	      if (((*EndPos) != '\0') && (EndPos != (DirectoryContent->d_name + PatternLength)) && (strcmp(suffix, EndPos) == 0))
		{
		  char* TmpName  = new char [strlen(DirectoryContent->d_name) + 1];
		  strcpy (TmpName, DirectoryContent->d_name);
		  MatchedFiles += TmpName;		
		}
	    }
	}
    }
  closedir(TmpDirectory);
  if (MatchedFiles.GetNbrElement() == 0)
    {
      matchedFileArray = 0;
      return 0;
    }
  ListIterator<char*> MatchedFileIterator(MatchedFiles);
  char** TmpName2;
  int* FileIndices = new int [MatchedFiles.GetNbrElement()];
  matchedFileArray = new char* [MatchedFiles.GetNbrElement()];
  int Pos = 0;
  while ((TmpName2 = MatchedFileIterator()))
    {
      FileIndices[Pos] = atoi ((*TmpName2) + PatternLength); 
      if (PathLength == 0)
	{
	  matchedFileArray[Pos] = (*TmpName2);
	}
      else
	{
	  char* TmpName3 = new char[PathLength + 2 + strlen(*TmpName2)];
	  strncpy (TmpName3, pattern, PathLength);
	  TmpName3[PathLength] = '/';
	  strcpy (TmpName3 + PathLength + 1,(*TmpName2));	  
	  matchedFileArray[Pos] = TmpName3;
	  delete[] (*TmpName2);
	}
      ++Pos;
    }
  SortArrayUpOrdering<char*>(FileIndices, matchedFileArray, Pos);
  delete[] FileIndices;
  delete[] TmpPattern;
  return Pos;
}

// test if a file exist and can be opened
//
// fileName = name of the file (must include relative/absolute path)
// return value = true if the file exists

bool IsFile (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      return false;
    }
  File.close();
  return true;
}

// concatenate path and file name
//
// path = string corresponding to the path 
// fileName = string corresponding to the file name (with optional relative path)
// return value = corresponding string

char* ConcatenatePathAndFileName (char* path, char* fileName)
{
  long TmpLength = strlen(path);
  long TmpLength2 = TmpLength + strlen(fileName) + 1;
  char* TmpFileName;
  if (path[TmpLength - 1] == '/')
    {
      TmpFileName = new char[TmpLength2];
      strcpy(TmpFileName, path);
    }
  else
    {
      ++TmpLength2;
      TmpFileName = new char[TmpLength2];
      strcpy(TmpFileName, path);
      TmpFileName[TmpLength] = '/';
      ++TmpLength;
   }
  strcpy(TmpFileName + TmpLength, fileName);
  return TmpFileName;
}


// concatenate path and file name
//
// input = input string
// path = string corresponding to the path 
// fileName = string corresponding to the file name
//
void ExtractPathAndFileName (const char* input, char* &path, char* &fileName)
{
  const char* TmpFile = strrchr(input,'/');
  if (TmpFile==NULL)
    {
      path = new char[3];
      sprintf(path,"./");
      fileName = new char[strlen(input)+1];
      strcpy(fileName,input);
    }
  else
    {
      TmpFile+=1;
      long TmpLength = strlen(input) - strlen(TmpFile);
      path = new char[TmpLength + 1];
      strncpy(path, input, TmpLength);
      fileName = new char[strlen(TmpFile) + 1];
      strcpy(fileName,TmpFile);
    }
}

// add a given extension to a file name
//
// fileName = string corresponding to the file name (with optional relative path)
// extension = extension (without initial dot)
// return value = corresponding string

char*  AddExtensionToFileName(char* fileName, const char* extension)
{
  long TmpLength = strlen(fileName);
  long TmpLength2 = TmpLength + strlen(extension) + 1;
  char* TmpFileName;
  if (fileName[TmpLength - 1] == '.')
    {
      TmpFileName = new char[TmpLength2];
      strcpy(TmpFileName, fileName);
    }
  else
    {
      ++TmpLength2;
      TmpFileName = new char[TmpLength2];
      strcpy(TmpFileName, fileName);
      TmpFileName[TmpLength] = '.';
      ++TmpLength;
   }
  strcpy(TmpFileName + TmpLength, extension);
  return TmpFileName;
}

// get the existing extension of a file name
//
// fileName = string corresponding to the file name (with optional relative path)
// maxL = maximum length of the extension to be searched for
// return value = corresponding string
char*  GetExtensionFromFileName(char* fileName, int maxL)
{
  int MaxL;  
  if (maxL>0)
    MaxL=maxL;
  else
    MaxL=strlen(fileName);
  char* Extension = fileName+strlen(fileName)-1;
  int Length=1;
  while ((*Extension != '.')&&(Length<MaxL))
    {
      Extension--;
      Length++;
    }
  if (*Extension == '.')
    {
      char *rst = new char[strlen(Extension)+1];
      strcpy(rst,Extension);
      return rst;
    }
  else
    return 0;
}

  


// replace extension to a file name
//
// fileName = string corresponding to the file name (with optional relative path)
// oldExtension = extension to replace (without initial dot)
// newExtension = new extension (without initial dot)
// return value = corresponding string (0 if the old extension was not found)

char* ReplaceExtensionToFileName(char* fileName, const char* oldExtension, const char* newExtension)
{
  char* ExtensionPosition = strstr(fileName, oldExtension);
  if (ExtensionPosition == 0)
    return 0;
  long TmpLength = strlen(fileName);
  long TmpLength2 = TmpLength - strlen(oldExtension);
  char* TmpFileName = new char[TmpLength2 + strlen(newExtension) + 1];
  strncpy (TmpFileName, fileName, TmpLength2);
  strcpy (TmpFileName + TmpLength2, newExtension);
  return TmpFileName;
}

// remove extension from a file name
//
// fileName = string corresponding to the file name (with optional relative path)
// oldExtension = extension to remove (with initial dot)
// return value = corresponding string (0 if the old extension was not found)

char* RemoveExtensionFromFileName(char* fileName, const char* oldExtension)
{
  char* ExtensionPosition = strstr(fileName, oldExtension);
  if (ExtensionPosition == 0)
    return 0;
  int TmpLength = strlen(fileName);
  int TmpLength2 = TmpLength - strlen(oldExtension);  
  char* TmpFileName = new char[TmpLength2 + 1];
  strncpy (TmpFileName, fileName, TmpLength2);
  TmpFileName[TmpLength2]='\0';
  return TmpFileName;
}




// get unique filename by appending a counter to a requested name, if necessary
// baseName = main part of the file name (initial part, excluding a final dot)
// optExtension = optional extension to add after the counter (including the initial dot)
// minCounter = optional minimum value of the counter
// return value = unique file name formed as baseName.[Count]optExtension

char* GetUniqueFileName(const char* baseName, const char* optExtension, int minCounter)
{
  return GetUniqueFileName(baseName, minCounter, optExtension);
}


// get unique filename by appending a counter to a requested name, if necessary
// baseName = main part of the file name (initial part, excluding a final dot)
// minCounter = provide minimum value of the counter and return the actual counter used in filename
// optExtension = optional extension to add after the counter (including the initial dot)
// return value = unique file name formed as baseName.[Count]optExtension
//
char* GetUniqueFileName(const char* baseName, int & minCounter, const char* optExtension)
{
  char *TheExtension;
  if (optExtension!=NULL)
    {
      TheExtension = new char[strlen(optExtension)+1];
      strcpy(TheExtension, optExtension);
    }
  else
    {
      TheExtension = new char[2];
      strcpy(TheExtension, "");
    }
  char *UniqueFileName = new char[strlen(baseName)+strlen(TheExtension)+10];
  sprintf(UniqueFileName,"%s.%d%s", baseName, minCounter++, TheExtension);
  std::ifstream testExistant(UniqueFileName,std::ios::in);
  while (testExistant.is_open())
    {
      testExistant.close();
      sprintf(UniqueFileName,"%s.%d%s", baseName, minCounter++, TheExtension);
      testExistant.open(UniqueFileName,std::ios::in);
    }
  delete [] TheExtension;
  std::ofstream TouchIt(UniqueFileName,std::ios::out);
  TouchIt.close();
  --minCounter;
  return UniqueFileName;
}


// test if file with given filename already exists, and if so, back it up
// fileName = filename to save
// optExtension = optional extension to add to backups (with leading dot, augmented by counter)
//
// returns = number of existing backup files (-1 if no backup required)
//
int BackUpFile(const char* fileName, const char* optExtension)
{
  int BackupCounter=-1;
  std::ifstream testExistant(fileName,std::ios::in);
  if (testExistant.is_open())
    {
      testExistant.close();
      char *BackupFileName;
      char *Extension;
      BackupCounter = 0;
      if (optExtension!=NULL)
	{
	  BackupFileName = new char[strlen(fileName)+strlen(optExtension)+5];
	  Extension = new char[strlen(optExtension)+5];
	  strcpy(Extension,optExtension);
	}
      else
	{
	  BackupFileName = new char[strlen(fileName)+5];
	  Extension = new char[1];
	  Extension[0]='\0';
	}
      sprintf(BackupFileName, "%s%s", fileName, Extension);
      testExistant.open(BackupFileName,std::ios::in);
      while (testExistant.is_open())
	{
	  testExistant.close();
	  ++BackupCounter;
	  sprintf(BackupFileName, "%s%s%d", fileName, Extension, BackupCounter);
	  testExistant.open(BackupFileName,std::ios::in);
	}
      testExistant.close();
      if ( rename(fileName, BackupFileName) != 0)
	{
	  cout << "Problem backing up existing file "<<fileName<<endl;
	}
      delete [] Extension;
      delete [] BackupFileName;
    }
  return BackupCounter;
}


// compute the number of lines in a text file
//
// fileName = text file name 
// return value = number of lines 

long GetFileNbrLines (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "error : cannot open file " <<  fileName << endl;
      return -1;
    }
  File.seekg(0, ios::end);
  long Size = File.tellg();
  long Size2 = Size & (~255l);
  File.seekg(0, ios::beg);
  char* TmpBuffer = new char [256];
  
  long NbrLines = 0l;
  for (long i = 0; i < Size2; i += 256l)
    {
      File.read(TmpBuffer, 256l);
      for (int j = 0; j < 256; ++ j)
	if (TmpBuffer[j] == '\n')
	  ++NbrLines;
    }
  Size &= 255l;
  if (Size > 0)
    {
      File.read(TmpBuffer, Size);
      for (int j = 0; j < Size; ++ j)
	if (TmpBuffer[j] == '\n')
	  ++NbrLines;
    }
  else
    Size = 256;
  if (TmpBuffer[Size - 1] != '\n')
    ++NbrLines;
  delete[] TmpBuffer;
  return NbrLines;
}

// get the requested line number from file as a string
// fileName = file to read
// nbrLine = line number to return (first line starts at zero)
// return = line, upon success
char* GetLineFromFile (char* fileName, int nbrLine)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "error : cannot open file " <<  fileName << endl;
      return 0;
    }
  File.seekg(0, ios::end);
  long Size = File.tellg();
  long Size2 = Size & (~255l);
  File.seekg(0, ios::beg);
  char* TmpBuffer = new char [256];
  
  long NbrLines = 0l;
  int j;
  long Start=0;
  long End=0;
  for (long i = 0; (i < Size2)&&(NbrLines<nbrLine); i += 256l)
    {
      File.read(TmpBuffer, 256l);
      for (j = 0; (j < 256)&&(NbrLines<nbrLine); ++ j)
	if (TmpBuffer[j] == '\n')
	  {
	    ++NbrLines;
	    if (NbrLines==nbrLine-1)
	      {
		Start = File.tellg();
		Start += j;
	      }
	    if (NbrLines==nbrLine)
	      {
		End = File.tellg();
		End +=j;
	      }
	  }
    }  
  Size &= 255l;
  if ((End==0)&&(Size > 0))
    {
      File.read(TmpBuffer, Size);
      for (int j = 0; j < Size; ++ j)
	if (TmpBuffer[j] == '\n')
	  {
	    ++NbrLines;
	    if (NbrLines==nbrLine-1)
	      {
		Start = File.tellg();
		Start += j;
	      }
	    if (NbrLines==nbrLine)
	      {
		End = File.tellg();
		End += j;
	      }
	  }
    }
  if ((Start!=0)&&(End==0)&&(TmpBuffer[Size - 1] != '\n')&&(NbrLines==nbrLine-1))
    End = Size;
  delete [] TmpBuffer;

  if ((Start!=0)&&(End!=0))
    {
      File.seekg(0, ios::beg);
      TmpBuffer = new char[End-Start+2];
      File.read(TmpBuffer, End-Start+1);
      TmpBuffer[End-Start+1]='\0';
      return TmpBuffer;
    }
  else
    return 0;
}



// get unique filename by appending a counter to a requested name, if necessary
// inputName = previous file name
// insertion = string to insert after element
// element = segment to be searched for
// HaveIntValue = element has optional integer argument
// return value = file name with inserted string
//
char* AddSegmentInFileName(char* inputName, const char* insertion, const char* element, bool HaveIntValue)
{
  char* StrElement;
  char* Result = new char[strlen(inputName)+strlen(insertion)+1];
  StrElement = strstr(inputName, element);
  if (StrElement != 0)
    {      
      StrElement += strlen(element);
      if (HaveIntValue)
	{
	  int SizeString = 0;
	  if (StrElement[SizeString] == '-')
	    ++SizeString;
	  while ((StrElement[SizeString] != '\0') && (StrElement[SizeString] != '.') && (StrElement[SizeString] >= '0') 
		 && (StrElement[SizeString] <= '9'))
	    ++SizeString;	  
	  StrElement += SizeString;
	}
      if (*StrElement!='\0')
	{
	  char* StringEnd=NULL;
	  StringEnd = new char[strlen(StrElement)+1];
	  strcpy(StringEnd, StrElement);
	  StrElement[0]='\0';
	  sprintf(Result,"%s%s%s",inputName,insertion,StringEnd);
	  delete [] StringEnd;
	}
      else
	sprintf(Result,"%s%s",inputName,insertion);
    }
  else
    sprintf(Result,"%s%s",inputName,insertion);
  return Result;
}

// create a directory if it does not exist
//
// directory = name of the directory to create
// return value = true if the directory has been created wothout error

bool CreateDirectory (const char* directory)
{
  struct stat Tmp = {0};
  if (stat(directory, &Tmp) != -1)
    {
      return true;
    }
  if (mkdir(directory, 0700) < 0)
    {
      return false;
    }
  return true;
}


// search for occurances of "fermion" and "boson" in a file name, then report false (i.e. bosons), true (i.e. fermions), or cannot be determined from file name
// OutputBool[in] = if the boolean is not true on input, do not parse and return immediately
// OutputBool[out] = boolean reporting the result of the statistics, as parsed from file name
// MyFilename = file name to search
// return = true, if the statistics was checked in the file name and a boolean is returned
//
bool FilenameStatisticsCheck (bool& OutputBool, char* MyFilename)
{
  if (OutputBool == true)
    {
      if (strstr(MyFilename, "fermion") == 0)
	{
	  if (strstr(MyFilename, "boson") == 0)
	    {
	      cout << "can't guess particle statistics from file name " << MyFilename << endl;
	      return false;	  
	    }
	  else
	    {
	      OutputBool = false; //for bosons
	    }
	}
      else
	{
	  OutputBool = true; //for fermions
	}
    }
  return true;
}

// extract the integer following the search string from a file name
// OutputInt[in] = if the integer is not zero on input, do not parse and return immediately
// OutputInt[out] = integer following search string, as parsed from file name
// MyFilename = file name to examine
// SearchString = string to search in file name
// return = true, if the search string was located and an integer is returned
// 
bool FilenameIntegerSearch (int& OutputInt, char* MyFilename, const char* SearchString)
{
  if (OutputInt==0)
    {
      char* MyPointer = strstr(MyFilename, SearchString); //set pointer to beginning of substring
      if (MyPointer != 0)
	{
	  MyPointer += strlen(SearchString); //move pointer to the end of substring
	  OutputInt = atoi(MyPointer); //convert subsequent number to integer
	}
      else
	{
	  cout << "can't guess " << SearchString << " from file name " << MyFilename << endl;
	  return false;            
	}
    }
  return true;
}

// extract the double following the search string from a file name
// OutputDouble[in] = if the double is not zero on input, do not parse and return immediately
// OutputDouble[out] = double following search string, as parsed from file name
// MyFilename = file name to examine
// SearchString = string to search in file name
// return = true, if the search string was located and a double is returned
// 
bool FilenameDoubleSearch (double& OutputDouble, char* MyFilename, const char* SearchString)
{
  if (OutputDouble==0)
    {
      char* MyPointer = strstr(MyFilename, SearchString); //set pointer to beginning of substring
      if (MyPointer != 0)
	{
	  MyPointer += strlen(SearchString); //move pointer to the end of substring
	  OutputDouble = strtod(MyPointer,NULL); //convert subsequent number to double
	}
      else
	{
	  cout << "can't guess " << SearchString << " from file name " << MyFilename << endl;
	  return false;            
	}
    }
  return true;
}

// search for an occurance of a string in a file name, then report true (i.e. present) or false (i.e. absent)
// OutputBool[in] = the derisred boolean variable, which needs updating
// OutputBool[out] = boolean reporting the result of the string search in the file name 
// MyFilename = file name to search
// SearchString = string to search in file name
// return = true, if the string was searched in the file name and a boolean is returned
//
bool FilenameBooleanSearch (bool& OutputBool, char* MyFilename, const char* SearchString)
{
  if (strstr(MyFilename, SearchString) != 0) //if there is an occurence of the substring, assign MyParameter=true
    {
      OutputBool=true;
    }
  return true;
}

// extract the character following the search string from a file name
// OutputChar[in] = if the character is not zero on input, do not parse and return immediately
// OutputChar[out] = character following search string, as parsed from file name
// MyFilename = file name to examine
// SearchString = string to search in file name
// return = true, if the search string was located and a character is returned
// 
bool FilenameCharacterSearch (char& OutputChar, char* MyFilename,  char const* SearchString)
{
  if (OutputChar==0)
    {
      char* MyPointer = strstr(MyFilename, SearchString); //set pointer to beginning of substring
      if (MyPointer != 0)
	{
	  MyPointer += strlen(SearchString); //move pointer to the end of substring	  
	  OutputChar = *MyPointer; // return the subsequent character
	}
      else
	{
	  cout << "can't find the character search string " << SearchString << " in file name " << MyFilename << endl;
	  return false;
	}
    }
  return true;
}

// extract the integer following the penultimate dot in a file name
// OutputInt[in] = if the integer is not zero on input, do not parse and return immediately
// OutputInt[out] = int following penultimate dot, as parsed from file name 
// MyFilename = file name to search
// return = true, if the penultimate dot was located and the following integer returned
//
bool FilenamePenultimateDotIntegerSearch (int& OutputInt, char* MyFilename)
{
  if (OutputInt==0)
    {
      int len=strlen(MyFilename);
      char* MyFilenameTmp = new char[len+1]; //create a dynamic-sized array on the heap
      strcpy(MyFilenameTmp, MyFilename);

      char* d1=NULL; char* d2=NULL;
      int i=0;
      for (char* c = MyFilenameTmp; *c != '\0'&& i <= len; c++, i++) //go through each character, terminate when you reach the end of the string (i.e. the null character)
	{
	  if(*c == '.')
	    {
	      d1 = d2; //penultimate dot
	      d2 = c; //final dot
	    }
	}
	
      OutputInt = atoi(d1+1); //position after the penultimate dot

      *d2 = '\0'; //terminate the stripped char array (e.g. without the ".0.vec" ending)
	
      delete[] MyFilenameTmp; //delete array copy from heap
    }
  return true;
}
