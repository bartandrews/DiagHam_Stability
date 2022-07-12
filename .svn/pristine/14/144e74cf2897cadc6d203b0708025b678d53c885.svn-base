////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of matrix main task operation                   //
//                                                                            //
//                        last modification : 10/06/2004                      //
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


#ifndef MAINTASKOPERATION_H
#define MAINTASKOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"


class AbstractMainTask;


class MainTaskOperation: public AbstractArchitectureOperation
{

 protected:

  // pointer to the main task
  AbstractMainTask* Task;

 public:
  
  // constructor 
  //
  // task = pointer to the main task
  MainTaskOperation(AbstractMainTask* task);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MainTaskOperation(const MainTaskOperation& operation);
  
  // destructor
  //
  ~MainTaskOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // get the main task
  //
  // return value = pointer to the main task
  AbstractMainTask* GetMainTask();

  // apply operation for a given architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ApplyOperation(AbstractArchitecture* architecture);

  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();

 protected:

  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
};

// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

inline void MainTaskOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
}

// get the main task
//
// return value = pointer to the main task

inline AbstractMainTask* MainTaskOperation::GetMainTask()
{
  return this->Task;
}

#endif
