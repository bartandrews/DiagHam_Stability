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


#include "config.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "MainTask/AbstractMainTask.h"

#include <iostream>


using std::cout;
using std::endl;


// constructor 
//
// task = pointer to the main task

MainTaskOperation::MainTaskOperation(AbstractMainTask* task)
{
  this->Task = task;
}


// copy constructor 
//
// operation = reference on operation to copy

MainTaskOperation::MainTaskOperation(const MainTaskOperation& operation)
{
  this->Task = operation.Task;
}

// destructor
//
MainTaskOperation::~MainTaskOperation()
{
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* MainTaskOperation::Clone()
{
  return new MainTaskOperation(*this);
}
  
// apply operation for a given architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool MainTaskOperation::ApplyOperation(AbstractArchitecture* architecture)
{
  switch (architecture->GetArchitectureID())
    {
    case AbstractArchitecture::SimpleMPI:
      return this->ArchitectureDependentApplyOperation((SimpleMPIArchitecture*) architecture);
    case AbstractArchitecture::MixedMPISMP:
      return this->ArchitectureDependentApplyOperation((SimpleMPIArchitecture*) architecture);  
    default:
      this->Task->SetArchitecture(architecture);
      return this->RawApplyOperation();
    }
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool MainTaskOperation::RawApplyOperation()
{
  if (this->Task->ExecuteMainTask() != 0)
    return false;
  else
    return true;
}
 
// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool MainTaskOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
  this->Task->SetArchitecture(architecture);
#ifdef __MPI__
  int TmpOperationID = 0x0;
  if (architecture->IsMasterNode())
    {
      bool Flag = true;
      architecture->RequestOperation(SimpleMPIArchitecture::SynchronizeSignal);
      if (this->Task->ExecuteMainTask() != 0)
  	Flag = false;
      architecture->RequestOperation(SimpleMPIArchitecture::FreeSlaveSignal);
      return Flag;
    }
  else
    {
      while (architecture->WaitOperation(TmpOperationID) == true)
	{
	  architecture->SendAcknowledge();
	  this->Task->ExecuteOperation(TmpOperationID);
	}
      architecture->SendAcknowledge();
      return true;
    }
#else
    return false;
#endif
}
