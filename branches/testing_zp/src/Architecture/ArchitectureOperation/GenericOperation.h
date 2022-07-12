////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of hamiltonian vector multiplication operation           //
//                                                                            //
//                        last modification : 10/03/2003                      //
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


#ifndef GENERICOPERATION_H
#define GENERICOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"

template<class ClassName>
class GenericOperation: public AbstractArchitectureOperation
{

 protected:

  // method to call for each operation (first argument is job index and the second argument is the total number of jobs)
  void (ClassName::*Method) (int, int);
  // object from whose method belongs to
  ClassName* Object;

  // job index
  int JobIndex;
  // job total number
  int NbrJob;

 public:
  
  // constructor 
  //
  // object = object 
  // method = method to call for each operation (first argument is job index and the second argument is the total number of jobs)
 GenericOperation (ClassName* object, void (ClassName::*method) (int, int));

  // copy constructor 
  //
  // operation = reference on operation to copy
 GenericOperation(const GenericOperation<ClassName>& operation);
  
  // destructor
  //
  ~GenericOperation();
  
  // set job index and job total number
  //
  // jobIndex = job index
  // nbrJob = job total number
  void SetJobIndices(int jobIndex, int nbrJob);

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation
  //
  // return value = true if no error occurs
  bool ApplyOperation();
  
};

// constructor 
//
// object = object 
// method = method to call for each operation (first argument is job index and the second argument is the total number of jobs)

template <class ClassName>
GenericOperation<ClassName>::GenericOperation (ClassName* object, void (ClassName::*method) (int, int))
{
  this->Object = object;
  this->Method = method;
  this->JobIndex = 0;
  this->NbrJob = 0;
  this->OperationType = AbstractArchitectureOperation::Generic;
}

// copy constructor 
//
// operation = reference on operation to copy

template <class ClassName>
GenericOperation<ClassName>::GenericOperation(const GenericOperation<ClassName>& operation)
{
  this->Object = this->Object;
  this->Method = this->Method;
  this->JobIndex = this->JobIndex;
  this->NbrJob = this->NbrJob;
  this->OperationType = AbstractArchitectureOperation::Generic;
}
  
// destructor
//

template <class ClassName>
GenericOperation<ClassName>::~GenericOperation()
{
}
  
// set job index and job total number
//
// jobIndex = job index
// nbrJob = job total number

template <class ClassName>
void GenericOperation<ClassName>::SetJobIndices(int jobIndex, int nbrJob)
{
  this->JobIndex = jobIndex;
  this->NbrJob = nbrJob;
}

// clone operation
//
// return value = pointer to cloned operation

template <class ClassName>
AbstractArchitectureOperation* GenericOperation<ClassName>::Clone()
{
  return new GenericOperation<ClassName>(*this);
}
  
// apply operation
//
// return value = true if no error occurs

template <class ClassName>
bool GenericOperation<ClassName>::ApplyOperation()
{
  ((this->Object)->*(this->Method))(this->JobIndex, this->NbrJob);
  return true;
}

#endif
