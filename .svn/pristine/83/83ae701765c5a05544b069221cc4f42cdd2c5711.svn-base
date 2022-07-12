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


// list all files or directories that obey a given pattern (that can include relative/absolute path) /to/directory/patternxxxsuffix where xxx is an integer
//
// pattern = string that corresponds to the pattern  (i.e. /to/directory/pattern)
// matchedFileArray = reference on the sorted array (with respect to xxx) of files or directories names (with the optional relative/absolute path), 
//                    memory allocation isd one by the function itself
// suffix = optional suffix  to test
// return value = number of matched files

int GetAllFilesDirectories(const char* pattern, char**& matchedFileArray, const char* suffix)
{
  char* Path = strrchr(pattern, '/');
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

