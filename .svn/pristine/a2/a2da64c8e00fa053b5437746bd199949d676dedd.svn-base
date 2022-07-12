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



#include "config.h"
#include "MainTask/AbstractMainTask.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperationManager.h"
#include "GeneralTools/ListIterator.h"


// virtual destructor
//  

AbstractMainTask::~AbstractMainTask()
{
  if (this->OperationManagers.GetNbrElement()>0)
    {
      ListIterator<AbstractArchitectureOperationManager*> Operations(this->OperationManagers);
      AbstractArchitectureOperationManager** TmpOperation;
      while ( (TmpOperation = Operations()) != NULL )
	{
	  delete *TmpOperation;
	}
    }
}

// set architecture binded to the task
// 
// architecture = pointer to the architecture to use

void AbstractMainTask::SetArchitecture(AbstractArchitecture* architecture)
{
  this->Architecture = architecture;
}

// execute a given architecture-dependent operation requested by the main task
//
// operationID = architecture operation ID
// return value = true if no error occured

bool AbstractMainTask::ExecuteOperation(int operationID)
{
  AbstractArchitectureOperationManager** TmpOperationManager;
  ListIterator<AbstractArchitectureOperationManager*> OperationManagerIterator(this->OperationManagers);
  while ((TmpOperationManager = OperationManagerIterator()))
    {
      if ((*TmpOperationManager)->IsHandled(operationID))
	{
	  AbstractArchitectureOperation* TmpOperation = (*TmpOperationManager)->GetOperation(operationID);
	  if (TmpOperation == 0)
	    return false;
	  else
	    {
	      bool Rst = TmpOperation->ApplyOperation(this->Architecture);
	      delete TmpOperation;
	      return Rst;
	    }
	}
    }
  return false;
}

