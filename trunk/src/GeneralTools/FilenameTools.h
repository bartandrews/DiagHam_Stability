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

#ifndef FILENAMETOOLS_H
#define FILENAMETOOLS_H

#include "config.h"

#include <fstream>
#include <string>
#include <dirent.h>

using std::ifstream;
using std::ofstream;
using std::ios;


// list all files or directories that obey a given pattern (that can include relative/absolute path) /to/directory/patternxxxsuffix where xxx is an integer
//
// pattern = string that corresponds to the pattern  (i.e. /to/directory/pattern)
// matchedFileArray = reference on the sorted array (with respect to xxx) of files or directories names (with the optional relative/absolute path), 
//                    memory allocation is done by the function itself
// suffix = optional suffix  to test
// return value = number of matched files
int GetAllFilesDirectories(const char* pattern, char**& matchedFileArray, const char* suffix = 0);

// test if a file exist and can be opened
//
// fileName = name of the file (must include relative/absolute path)
// return value = true if the file exists
bool IsFile (const char* fileName);

// concatenate path and file name
//
// path = string corresponding to the path 
// fileName = string corresponding to the file name (with optional relative path)
// return value = corresponding string
char* ConcatenatePathAndFileName (char* path, char* fileName);

// concatenate path and file name
//
// input = input string
// path = string corresponding to the path 
// fileName = string corresponding to the file name
//
void ExtractPathAndFileName (const char* input, char* &path, char* &fileName);

// add a given extension to a file name
//
// fileName = string corresponding to the file name (with optional relative path)
// extension = extension (without initial dot)
// return value = corresponding string
char*  AddExtensionToFileName(char* fileName, const char* extension);

// get the existing extension of a file name
//
// fileName = string corresponding to the file name (with optional relative path)
// maxL = maximum length of the extension to be searched for
// return value = corresponding string
char*  GetExtensionFromFileName(char* fileName, int maxL=0);

// replace extension to a file name
//
// fileName = string corresponding to the file name (with optional relative path)
// oldExtension = extension to replace (without initial dot)
// newExtension = new extension (without initial dot)
// return value = corresponding string (0 if the old extension was not found)
char* ReplaceExtensionToFileName(char* fileName, const char* oldExtension, const char* newExtension);


// remove extension from a file name
//
// fileName = string corresponding to the file name (with optional relative path)
// oldExtension = extension to remove (with initial dot)
// return value = corresponding string (0 if the old extension was not found)
char* RemoveExtensionFromFileName(char* fileName, const char* oldExtension);


// get unique filename by appending a counter to a requested name, if necessary
//
// baseName = main part of the file name (initial part, excluding a final dot)
// optExtension = optional extension to add after the counter (including the initial dot)
// minCounter = optional minimum value of the counter
// return value = unique file name formed as baseName.[Count]optExtension
char* GetUniqueFileName(const char* baseName, const char* optExtension = NULL, int minCounter = 0);

// get unique filename by appending a counter to a requested name, if necessary
//
// baseName = main part of the file name (initial part, excluding a final dot)
// minCounter = provide minimum value of the counter and return the actual counter used in filename
// optExtension = optional extension to add after the counter (including the initial dot)
// return value = unique file name formed as baseName.[Count]optExtension
char* GetUniqueFileName(const char* baseName, int & minCounter, const char* optExtension = NULL);

// test if file with given filename already exists, and if so, back it up
//
// fileName = filename to save
// optExtension = optional extension to add to backups (with leading dot, augmented by counter)
// return = number of existing backup files (-1 if no backup required)
int BackUpFile(const char* fileName, const char* optExtension = ".bak");

// compute the number of lines in a text file
//
// fileName = text file name 
// return value = number of lines 
long GetFileNbrLines (char* fileName);

// get the requested line number from file as a string
// fileName = file to read
// nbrLine = line number to return (first line starts at zero)
// return = line, upon success
char* GetLineFromFile (char* fileName, int nbrLine);

// get unique filename by appending a counter to a requested name, if necessary
//
// inputName = previous file name
// insertion = string to insert after element
// element = segment to be searched for
// HaveIntValue = element has optional integer argument
// return value = file name with inserted string
char* AddSegmentInFileName(char* inputName, const char* insertion, const char* element, bool HaveIntValue=false);

// create a directory if it does not exist
//
// directory = name of the directory to create
// return value = true if the directory has been created wothout error
bool CreateDirectory (const char* directory);


// tools for parsing parts of filenames

// extract the integer following the search string from a file name
//
// OutputInt[in] = if the integer is not zero on input, do not parse and return immediately
// OutputInt[out] = integer following search string, as parsed from file name
// MyFilename = file name to examine
// SearchString = string to search in file name
// return value = true if the search string was located and an integer is returned
bool FilenameIntegerSearch (int& OutputInt, char* MyFilename, const char* SearchString);

// extract the double following the search string from a file name
//
// OutputDouble[in] = if the double is not zero on input, do not parse and return immediately
// OutputDouble[out] = double following search string, as parsed from file name
// MyFilename = file name to examine
// SearchString = string to search in file name
// return value = true if the search string was located and a double is returned
bool FilenameDoubleSearch (double& OutputDouble, char* MyFilename, const char* SearchString);

// search for an occurance of a string in a file name, then report true (i.e. present) or false (i.e. absent)
//
// OutputBool[in] = the derisred boolean variable, which needs updating
// OutputBool[out] = boolean reporting the result of the string search in the file name 
// MyFilename = file name to search
// SearchString = string to search in file name
// return value = true if the string was searched in the file name and a boolean is returned
bool FilenameBooleanSearch (bool& OutputBool, char* MyFilename, const char* SearchString);

// extract the character following the search string from a file name
//
// OutputChar[in] = if the character is not zero on input, do not parse and return immediately
// OutputChar[out] = character following search string, as parsed from file name
// MyFilename = file name to examine
// SearchString = string to search in file name
// return value = true if the search string was located and a character is returned
bool FilenameCharacterSearch (char& OutputChar, char* MyFilename, char const* SearchString);

// search for occurances of "fermion" and "boson" in a file name, then report false (i.e. bosons), true (i.e. fermions), or cannot be determined from file name
//
// OutputBool[in] = if the boolean is not true on input, do not parse and return immediately
// OutputBool[out] = boolean reporting the result of the statistics, as parsed from file name
// MyFilename = file name to search
// return value = true if the statistics was checked in the file name and a boolean is returned
bool FilenameStatisticsCheck (bool& OutputBool, char* MyFilename);

// extract the integer following the penultimate dot in a file name
//
// OutputInt[in] = if the integer is not zero on input, do not parse and return immediately
// OutputInt[out] = int following penultimate dot, as parsed from file name 
// MyFilename = file name to search
// return value = true if the penultimate dot was located and the following integer returned
bool FilenamePenultimateDotIntegerSearch (int& OutputInt, char* MyFilename);


#endif
