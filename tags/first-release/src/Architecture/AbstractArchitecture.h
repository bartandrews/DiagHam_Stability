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
//                        last modification : 10/04/2002                      //
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


#ifndef ABSTRACTARCHITECTURE_H
#define ABSTRACTARCHITECTURE_H


#include "config.h"


class Vector;
class AbstractHamiltonian;

class AbstractArchitectureOperation;
class VectorHamiltonianMultiplyOperation;
class AddRealLinearCombinationOperation;
class MultipleRealScalarProductOperation;
class MatrixMatrixMultiplyOperation;
class QHEParticlePrecalculationOperation;


class AbstractArchitecture
{

public:
  
  // destructor
  //
  virtual ~AbstractArchitecture();
  
  // multiply a vector by an hamiltonian and store the result in another vector
  //
  // hamiltonian = pointer to the hamiltonian to use
  // vSource = vector to multiply 
  // vDestination = vector where result has to be stored 
  virtual void Multiply (AbstractHamiltonian* hamiltonian, Vector& vSource, Vector& vDestination) = 0;

  // execute an architecture-dependent operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (AbstractArchitectureOperation* operation);

  // execute an architecture-dependent vector hamiltonian multiplication operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (VectorHamiltonianMultiplyOperation* operation);
  
  // execute an architecture-dependent add real linear combination operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (AddRealLinearCombinationOperation* operation);
  
  // execute an architecture-dependent multiple real scalar product operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (MultipleRealScalarProductOperation* operation);
  
  // execute an architecture-dependent matrix matrix multiplication operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (MatrixMatrixMultiplyOperation* operation);
    
  // execute an architecture-dependent QHE particle hamiltonian precalculation operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (QHEParticlePrecalculationOperation* operation);
    
};

#endif
