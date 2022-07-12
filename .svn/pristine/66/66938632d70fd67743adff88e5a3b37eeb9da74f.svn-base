////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of abstract architecture operation manager              //
//                                                                            //
//                        last modification : 13/09/2005                      //
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


#ifndef ABSTRACTARCHITECTUREOPERATIONMANAGER_H
#define ABSTRACTARCHITECTUREOPERATIONMANAGER_H


#include "config.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"


class AbstractArchitectureOperationManager
{

 protected:

  // pointer to the architecture (must provide a communicator)
  SimpleMPIArchitecture* Architecture;

 public:
  
  // destructor
  //
  virtual ~AbstractArchitectureOperationManager();

  // test if an operation can be handled by the manager
  // 
  // operationI = ID of the operation to get
  // return value = true if the operation can be handled
  virtual bool IsHandled(int operationID) = 0;
  
  // retrieve an operation from its ID, and initialize using information communicated through the archiecture communicator
  //
  // operationI = ID of the operation to get
  // return value = pointer to the operation (null if not handled by the manager, or if error occured during initialization)
  virtual AbstractArchitectureOperation* GetOperation(int operationID) = 0;
  
};


#endif
