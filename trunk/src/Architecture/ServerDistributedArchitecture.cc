////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of server for ditributed architecture                //
//                                                                            //
//                        last modification : 09/04/2003                      //
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
#include "Architecture/ServerDistributedArchitecture.h"
#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"


// read node network configuration
//
// fileName = name of the file containing node network configuration
// return value = true if no error occurs

bool ServerDistributedArchitecture::ReadNodeNetworkConfiguration(char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  File.seekg(0, ios::end);
  int FileSize = File.tellg();
  char* TmpBuffer = new char [FileSize];
  File.seekg(0, ios::beg); 
  File.read (TmpBuffer, FileSize);
  delete[] TmpBuffer;
  File.close();
  int StartPosition = 0;
  int CurrentPosition = 0;
  while (StartPosition < FileSize)
    {
      while (TmpBuffer[CurrentPosition] != '\n')
	++CurrentPosition;
      if (TmpBuffer[StartPosition] == '#')
	{
	  
	}
      ++CurrentPosition;
      StartPosition = CurrentPosition;
    }
#ifdef __DEBUG__
#endif
}
