////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//      class for a basic Monte Carlo algorith for particles on a sphere      //
//                                                                            //
//                        last modification : 23/01/2008                      //
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


#ifndef ABSTRACTBLOCKOBSERVABLE_H
#define ABSTRACTBLOCKOBSERVABLE_H

#include "config.h"

#include <iostream>
#include "AbstractObservable.h"
#include "MathTools/Complex.h"

class AbstractBlockObservable : public AbstractObservable
{
  
 public:
  
  // destructor
  virtual ~AbstractBlockObservable();

  // call to make an observation
  // weight = relative weight of this sample
  virtual void RecordRealValue(Complex *weights) = 0;

  // call to make an observation
  // weight = relative weight of this sample
  virtual void RecordValue(Complex *weights) = 0;
  
};

#endif
