////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of Abstract Architecture                      //
//                                                                            //
//                        last modification : 30/04/2002                      //
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


#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MatrixMatrixMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
//#include "Architecture/ArchitectureOperation/GenericOperation.h"


// destructor
//

AbstractArchitecture::~AbstractArchitecture()
{
}

// execute an architecture-dependent operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (AbstractArchitectureOperation* operation)
{
  switch (operation->GetOperationType())
    {
    case AbstractArchitectureOperation::VectorHamiltonianMultiply:
      return this->ExecuteOperation((VectorHamiltonianMultiplyOperation*) operation);
      break;
    case AbstractArchitectureOperation::AddRealLinearCombination:
      return this->ExecuteOperation((AddRealLinearCombinationOperation*) operation);
      break;
    case AbstractArchitectureOperation::MultipleRealScalarProduct:
      return this->ExecuteOperation((MultipleRealScalarProductOperation*) operation);
      break;
    case AbstractArchitectureOperation::MatrixMatrixMultiply:
      return this->ExecuteOperation((MatrixMatrixMultiplyOperation*) operation);
      break;
//    case AbstractArchitectureOperation::Generic:
//      return this->ExecuteOperation((GenericOperation*) operation);
//      break;
    default:
      return false;
    }
  return false;
}

// execute an architecture-dependent vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (VectorHamiltonianMultiplyOperation* operation)
{
  return false;
}

// execute an architecture-dependent add real linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (AddRealLinearCombinationOperation* operation)
{
  return false;
}

// execute an architecture-dependent multiple real scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (MultipleRealScalarProductOperation* operation)
{
  return false;
}

// execute an architecture-dependent matrix matrix multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (MatrixMatrixMultiplyOperation* operation)
{
  return false;
}

// execute an architecture-dependent QHE particle hamiltonian precalculation operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (QHEParticlePrecalculationOperation* operation)
{
  return false;
}
    
