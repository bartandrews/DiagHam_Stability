////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of Abstract Architecture                      //
//                                                                            //
//                        last modification : 30/04/2002                      //
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

#include "Architecture/AbstractArchitecture.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include <sys/time.h>
#include <string.h>


// destructor
//

AbstractArchitecture::~AbstractArchitecture()
{
}

// get the amount of memory available for the local architecture
//
// return value = amount of memory in byte (negative if the information is not available)

long AbstractArchitecture::GetLocalMemory()
{
  return -1l;
}
  
// get typical range of indices on which the local architecture acts
//
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)

void AbstractArchitecture::GetTypicalRange (long& minIndex, long& maxIndex)
{
  minIndex = 0;
  maxIndex = this->HilbertSpaceDimension - 1;
}
  
// get a new real vector with memory alloaction depending on the architecture
//
// return value = pointer to the requested vector (zero if an error occurs)

RealVector* AbstractArchitecture::GetNewRealVector ()
{
  return new RealVector;
}
  
// get a new real vector with memory alloaction depending on the architecture
//
// dimension = dimension of the requested vector
// zeroFlag = true if all vector entries has to be set to zero
// return value = pointer to the requested vector (zero if an error occurs)

RealVector* AbstractArchitecture::GetNewRealVector (long dimension, bool zeroFlag)
{
  return new RealVector(dimension, zeroFlag);
}
  
// get a new complex vector with memory alloaction depending on the architecture
//
// return value = pointer to the requested vector (zero if an error occurs)

ComplexVector* AbstractArchitecture::GetNewComplexVector ()
{
  return new ComplexVector;
}
  
// get a new complex vector with memory alloaction depending on the architecture
//
// dimension = dimension of the requested vector
// zeroFlag = true if all vector entries has to be set to zero
// return value = pointer to the requested vector (zero if an error occurs)

ComplexVector* AbstractArchitecture::GetNewComplexVector (long dimension, bool zeroFlag)
{
  return new ComplexVector (dimension, zeroFlag);
}
  
// set dimension of the Hilbert space on which the architecture has to work
// 
// dimension = dimension of the Hilbert space

void AbstractArchitecture::SetDimension (long dimension)
{
  this->HilbertSpaceDimension = dimension;
}

// get a temporary file name
//
// return value = string corresponding to a temporary file name

char* AbstractArchitecture::GetTemporaryFileName()
{
  timeval Time;
  gettimeofday (&Time, 0);
  char* TmpString = new char [32];
  sprintf (TmpString, "diagam%d%d.tmp",(int)  Time.tv_sec, (int)  Time.tv_usec);
  return TmpString;
}
  
