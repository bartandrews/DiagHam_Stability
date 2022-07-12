////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of cluster architecture using MPI and for which          //
//    each vector has its associated full datas only located on one process   //
//                                                                            //
//                        last modification : 22/06/2004                      //
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
#include "Architecture/ClusterArchitecture/DelocalizedVectorMPIArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MatrixMatrixMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"

#ifdef __MPI__
#include <mpi.h>
#endif

// constructor
// 

DelocalizedVectorMPIArchitecture::DelocalizedVectorMPIArchitecture()
{
  this->LocalArchitecture = new MonoProcessorArchitecture;
}

// destructor
//

DelocalizedVectorMPIArchitecture::~DelocalizedVectorMPIArchitecture()
{
#ifdef __MPI__
  MPI::Finalize();
#endif
  delete this->LocalArchitecture;
  if (this->ClusterPerformanceArray != 0)
    delete[] this->ClusterPerformanceArray;
}
  
// get typical range of indices on which the local architecture acts
//
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)

void DelocalizedVectorMPIArchitecture::GetTypicalRange (long& minIndex, long& maxIndex)
{
  minIndex = this->MinimumIndex;
  maxIndex = this->MaximumIndex;
}
  
// execute an architecture-dependent vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool DelocalizedVectorMPIArchitecture::ExecuteOperation (VectorHamiltonianMultiplyOperation* operation)
{
}
  
// execute an architecture-dependent add real linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool DelocalizedVectorMPIArchitecture::ExecuteOperation (AddRealLinearCombinationOperation* operation)
{
}
  
// execute an architecture-dependent add complex linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool DelocalizedVectorMPIArchitecture::ExecuteOperation (AddComplexLinearCombinationOperation* operation)
{
}

// execute an architecture-dependent multiple real scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool DelocalizedVectorMPIArchitecture::ExecuteOperation (MultipleRealScalarProductOperation* operation)
{
}
  
// execute an architecture-dependent multiple complex scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool DelocalizedVectorMPIArchitecture::ExecuteOperation (MultipleComplexScalarProductOperation* operation)
{
}
  
// execute an architecture-dependent matrix matrix multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool DelocalizedVectorMPIArchitecture::ExecuteOperation (MatrixMatrixMultiplyOperation* operation)
{
}

// execute an architecture-dependent abstract hamiltonian precalculation operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool DelocalizedVectorMPIArchitecture::ExecuteOperation (AbstractPrecalculationOperation* operation)
{
}

// get a new real vector with memory alloaction depending on the architecture
//
// return value = pointer to the requested vector (zero if an error occurs)

RealVector* DelocalizedVectorMPIArchitecture::GetNewRealVector ()
{
}
  
// get a new real vector with memory alloaction depending on the architecture
//
// dimension = dimension of the requested vector
// zeroFlag = true if all vector entries has to be set to zero
// return value = pointer to the requested vector (zero if an error occurs)

RealVector* DelocalizedVectorMPIArchitecture::GetNewRealVector (long dimension, bool zeroFlag = false)
{
}
  
// get a new complex vector with memory alloaction depending on the architecture
//
// return value = pointer to the requested vector (zero if an error occurs)

ComplexVector* DelocalizedVectorMPIArchitecture::GetNewComplexVector ()
{
}
  
// get a new complex vector with memory alloaction depending on the architecture
//
// dimension = dimension of the requested vector
// zeroFlag = true if all vector entries has to be set to zero
// return value = pointer to the requested vector (zero if an error occurs)

ComplexVector* DelocalizedVectorMPIArchitecture::GetNewComplexVector (long dimension, bool zeroFlag = false)
{
}
  
// get a vector element from a real vector
//
// vectorId = id of vector where the lement is located
// index = element index in the vector
// return value = corresponding vector element

double DelocalizedVectorMPIArchitecture::RequestRealVectorElement(int vectorId, int index)
{
}
  
// set a vector element into a real vector
//
// component = value of the vector element
// vectorId = id of vector where the element is located
// index = element index in the vector

void DelocalizedVectorMPIArchitecture::SetRealVectorElement(const double& component, int vectorId, int index)
{
}
 
// get a vector element from a complex vector
//
// vectorId = id of vector where the lement is located
// index = element index in the vector
// return value = corresponding vector element

Complex DelocalizedVectorMPIArchitecture::RequestComplexVectorElement(int vectorId, int index)
{
}
  
// set a vector element into a real vector
//
// component = value of the vector element
// vectorId = id of vector where the lement is located
// index = element index in the vector

void DelocalizedVectorMPIArchitecture::SetComplexVectorElement(const Complex& component, int vectorId, int index)
{
}
  
// get a real vector 
// 
// vectorId = id of vector to get
// return value = corresponding vector

RealVector DelocalizedVectorMPIArchitecture::GetRealVector(int vectorId)
{
}
  
// set a real vector 
// 
// vector = reference on the vector which contains thedatas 
// vectorId = id of vector where to store the datas

void DelocalizedVectorMPIArchitecture::SetRealVector(RealVector& vector, int vectorId)
{
}
  
// indicate an allocation of memory to the architecture
//
// pointer = pointer to the memory zone which will be allocated
// memory = amount of requested memory in bytes

void DelocalizedVectorMPIArchitecture::AllocateMemory (void* pointer, unsigned long memory)
{
  this->
}

// indicate an deallocation of memory to the architecture
//
// pointer = pointer to the memory zone which will be free

void DelocalizedVectorMPIArchitecture::DeallocateMemory (void* pointer)
{
}
