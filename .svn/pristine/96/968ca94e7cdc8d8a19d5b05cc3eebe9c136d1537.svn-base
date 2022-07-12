////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of generic 2D tight binding model                 //                                                                          //
//                        last modification : 03/10/2012                      //
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


#ifndef GENERIC2DTIGHTBINDINGMODEL_H
#define GENERIC2DTIGHTBINDINGMODEL_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class Generic2DTightBindingModel : public Abstract2DTightBindingModel
{

 public:

  // default constructor
  //
  // fileName = name of the binary file that contains the band structure information
  Generic2DTightBindingModel(char* fileName);
  
    // destructor
  //
  ~Generic2DTightBindingModel();

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);
  
 };

#endif
