////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.08                           //
//                                                                            //
//                  Copyright (C) 1998-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         a set of functions used for little<->big endian convertion         //
//                                                                            //
//                          first release : 05/06/2003                        //
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


#include <fstream>


using std::ofstream;
using std::ifstream;


// function to read Little Endian encoded variable from a file
//
// file = reference on the input file stream
// var = reference on the variable to store the result

template<class ClassName>
void ReadLittleEndian (ifstream& file, ClassName& var)
{
  file.read ((char*) &var, sizeof(ClassName));
#ifdef __BIGENDIAN__
  ClassName TmpVar = var;
  unsigned char* TmpBin1 = (unsigned char*) &var;
  unsigned char* TmpBin2 = (unsigned char*) &TmpVar;
  int max = sizeof(ClassName) >> 1;
  for (int i = 0; i < max; i++)
    {
      TmpBin1[i] = TmpBin2[sizeof(ClassName) - i -1];
      TmpBin1[sizeof(ClassName) - i -1] = TmpBin2[i];
    }
#endif
};

// function to write Little Endian encoded variable from a file using 
//
// file = reference on the output file stream
// var = reference on the variable to store the result

template<class ClassName>
void WriteLittleEndian (ofstream& file, ClassName& var)
{
#ifdef __BIGENDIAN__
  ClassName TmpVar = var;
  unsigned char* TmpBin2 = (unsigned char*) &var;
  unsigned char* TmpBin1 = (unsigned char*) &TmpVar;
  int max = sizeof(ClassName) >> 1;
  for (int i = 0; i < max; i++)
    {
      TmpBin1[i] = TmpBin2[sizeof(ClassName) - i -1];
      TmpBin1[sizeof(ClassName) - i -1] = TmpBin2[i];
    }
  file.write ((char*) &TmpVar, sizeof(ClassName));
#else
  file.write ((char*) &var, sizeof(ClassName));
#endif
};

