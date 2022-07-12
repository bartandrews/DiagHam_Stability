////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of  architecture base operation manager                 //
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


#include "config.h"
#include "Architecture/ArchitectureOperation/ArchitectureBaseOperationManager.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/MultipleVectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/HamiltonianFullDiagonalizeOperation.h"


// constructor
//
// architecture = pointer to the architecture
// hamiltonian = pointer to the hamiltonian that will be provide to operation

ArchitectureBaseOperationManager::ArchitectureBaseOperationManager(SimpleMPIArchitecture* architecture, AbstractHamiltonian* hamiltonian)
{
  this->Architecture = architecture;
  this->Hamiltonian = hamiltonian;
}

// destructor
//

ArchitectureBaseOperationManager::~ArchitectureBaseOperationManager()
{
}

// test if an operation can be handled by the manager
// 
// operationI = ID of the operation to get
// return value = true if the operation can be handled

bool ArchitectureBaseOperationManager::IsHandled(int operationID)
{
  switch (operationID)
    {
    case AbstractArchitectureOperation::VectorHamiltonianMultiply:
      if (this->Hamiltonian != 0)
	return true;
      else
	return false;
    case AbstractArchitectureOperation::MultipleVectorHamiltonianMultiply:
      if (this->Hamiltonian != 0)
	return true;
      else
	return false;
    case AbstractArchitectureOperation::HamiltonianFullDiagonalize:
      if (this->Hamiltonian != 0)
	return true;
      else
	return false;
    default:
      return false; 
    }
  return false; 
}
  
// retrieve an operation from its ID, and initialize using information communicated through the archiecture communicator
//
// operationI = ID of the operation to get
// return value = pointer to the operation (null if not handled by the manager, or if error occured during initialization)

AbstractArchitectureOperation* ArchitectureBaseOperationManager::GetOperation(int operationID)
{
  if (this->Architecture == 0)
    return 0;
  switch (operationID)
    {
    case AbstractArchitectureOperation::VectorHamiltonianMultiply:
      if (this->Hamiltonian != 0)
	return new VectorHamiltonianMultiplyOperation(this->Hamiltonian, this->Architecture);
      else
	return 0;
    case AbstractArchitectureOperation::MultipleVectorHamiltonianMultiply:
      if (this->Hamiltonian != 0)
	return new MultipleVectorHamiltonianMultiplyOperation(this->Hamiltonian, this->Architecture);
      else
	return 0;
    case AbstractArchitectureOperation::HamiltonianFullDiagonalize:
      if (this->Hamiltonian != 0)
	return new HamiltonianFullDiagonalizeOperation(this->Hamiltonian, this->Architecture);
      else
	return 0;
    default:
      return 0;
    }
  return 0; 
}
  
