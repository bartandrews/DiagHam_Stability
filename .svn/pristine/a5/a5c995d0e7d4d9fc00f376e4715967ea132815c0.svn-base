////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of MonoProcessor Architecture                    //
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


#include "config.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MatrixMatrixMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"


// constructor
//

MonoProcessorArchitecture::MonoProcessorArchitecture()
{
}
  
// destructor
//

MonoProcessorArchitecture::~MonoProcessorArchitecture()
{
}
  
// multiply a vector by an hamiltonian and store the result in another vector
//
// hamiltonian = pointer to the hamiltonian to use
// vSource = vector to multiply 
// vDestination = vector where result has to be stored 

void MonoProcessorArchitecture::Multiply (AbstractHamiltonian* hamiltonian, Vector& vSource, Vector& vDestination)
{
  hamiltonian->Multiply(vSource, vDestination);
}

// execute an architecture-dependent vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool MonoProcessorArchitecture::ExecuteOperation (VectorHamiltonianMultiplyOperation* operation)
{
  return operation->ApplyOperation();
}
  

// execute an architecture-dependent add real linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool MonoProcessorArchitecture::ExecuteOperation (AddRealLinearCombinationOperation* operation)
{
  return operation->ApplyOperation();
}  

// execute an architecture-dependent multiple real scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool MonoProcessorArchitecture::ExecuteOperation (MultipleRealScalarProductOperation* operation)
{
  return operation->ApplyOperation();
}  

// execute an architecture-dependent matrix matrix multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool MonoProcessorArchitecture::ExecuteOperation (MatrixMatrixMultiplyOperation* operation)
{
  bool Flag = operation->ApplyOperation();
  return Flag;
}
 
// execute an architecture-dependent QHE particle hamiltonian precalculation operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool MonoProcessorArchitecture::ExecuteOperation (QHEParticlePrecalculationOperation* operation)
{
  return operation->ApplyOperation();
}
    
