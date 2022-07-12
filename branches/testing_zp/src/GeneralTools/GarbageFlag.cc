////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                             class of garbage flag                          //
//                                                                            //
//                        last modification : 06/07/2001                      //
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


#include "GeneralTools/GarbageFlag.h"

#include <iostream>


// default constructor
//

GarbageFlag::GarbageFlag()
{
  this->Counter = 0;
}

// copy constructor
//
// flag = garbage flag to copy

GarbageFlag::GarbageFlag(const GarbageFlag& flag)
{
  this->Counter = flag.Counter;
  if (this->Counter != 0)
    (*(this->Counter))++;
}

// destructor
//

GarbageFlag::~GarbageFlag()
{
  if (this->Counter != 0)
    {
      if ((*(this->Counter)) == 1)
	{
	  (*(this->Counter))--;
	  delete this->Counter;
	}
      else
	(*(this->Counter))--;
    }
}
 
// assignement
// 
// flag = garbage flag to assign
// return value = reference on current garbage flag

GarbageFlag& GarbageFlag::operator = (const GarbageFlag& flag)
{
  if (this->Counter != 0)
    if ((*(this->Counter)) == 1)
      delete this->Counter;
    else
      (*(this->Counter))--;
  this->Counter = flag.Counter;
  if (this->Counter != 0)
    (*(this->Counter))++;
  return *this;
}

// initialize garbage flag
//

void GarbageFlag::Initialize ()
{
  if (this->Counter != 0)
    if ((*(this->Counter)) == 1)
      delete this->Counter;
    else
      (*(this->Counter))--;
    this->Counter = new int;
  (*(this->Counter)) = 1;    
}

// test if datas are shared
//
// return value = true if datas are shared

bool GarbageFlag::Shared ()
{
  if (this->Counter == 0)
    return false;
  else
    if ((*(this->Counter)) == 1)
      return false;
    else
      return true;
}

// test if datas are used
//
// return value = true if datas are used

bool GarbageFlag::Used ()
{
  if (this->Counter == 0)
    return false;
  else
    return true;
}
