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

#endif
