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


#ifndef GARBAGEFLAG_H
#define GARBAGEFLAG_H

#include "config.h"

#include <pthread.h>

class GarbageFlag
{

 private:
 
  int* Counter;

#ifdef __SMP__
  pthread_mutex_t *FlagMutex;
#endif


 public:

  // default constructor
  //
  GarbageFlag();

  // copy constructor
  //
  // flag = garbage flag to copy
  GarbageFlag(const GarbageFlag& flag);

  // destructor
  //
  ~GarbageFlag();
  
  // assignement
  // 
  // flag = garbage flag to assign
  // return value = reference on current garbage flag
  GarbageFlag& operator = (const GarbageFlag& flag);

  // initialize garbage flag
  //
  void Initialize ();

  // test if datas are shared
  //
  // return value = true if datas are shared
  bool Shared ();

  // test if datas are used
  //
  // return value = true if datas are used
  bool Used ();

  // test routine to check for SMP compatibility
  //
  bool ThreadSafe ();

};

#endif
