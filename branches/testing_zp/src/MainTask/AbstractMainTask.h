////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract main task                        //
//                                                                            //
//                        last modification : 09/06/2004                      //
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


#ifndef ABSTRACTMAINTASK_H
#define ABSTRACTMAINTASK_H


#include "config.h"
#include "GeneralTools/List.h"


class AbstractArchitecture;
class AbstractArchitectureOperationManager;


class AbstractMainTask
{

 protected:

  // pointer to the architecture 
  AbstractArchitecture* Architecture;

  // list of managers of architecture operations that can be used by the main task
  List<AbstractArchitectureOperationManager*> OperationManagers;

 public:

  // virtual destructor
  //  
  virtual ~AbstractMainTask();
  
  // execute the main task
  // 
  // return value = 0 if no error occurs, else return error code
  virtual int ExecuteMainTask() = 0;

  // set architecture binded to the task
  // 
  // architecture = pointer to the architecture to use
  virtual void SetArchitecture(AbstractArchitecture* architecture);

  // execute a given architecture-dependent operation requested by the main task
  //
  // operationID = architecture operation ID
  // return value = true if no error occured
  virtual bool ExecuteOperation(int operationID);

};

#endif
