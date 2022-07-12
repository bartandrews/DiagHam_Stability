////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of Architecture Manager                      //
//                                                                            //
//                        last modification : 28/05/2004                      //
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


#ifndef ARCHITECTUREMANAGER_H
#define ARCHITECTUREMANAGER_H


#include "config.h"


class AbstractArchitecture;
class OptionManager;


class ArchitectureManager
{

 protected:

  // pointer to the best avalaible architecture in agreement with the option constraints
  AbstractArchitecture* Architecture;

  // pointer to the option manager
  OptionManager* Options;

 public:

  // default constructor
  //
  ArchitectureManager();

  // destructor
  //
  ~ArchitectureManager();
 
  // add an option group containing all options related to the architecture 
  //
  // manager = pointer to the option manager
  void AddOptionGroup(OptionManager* manager);

  // get the best avalaible architecture in agrement with the option constraints
  //
  // return value = pointer to the architecture
  AbstractArchitecture* GetArchitecture();

};

#endif

