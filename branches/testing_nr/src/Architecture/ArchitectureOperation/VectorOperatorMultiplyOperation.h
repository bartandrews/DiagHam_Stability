////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of operator vector multiplication operation           //
//                                                                            //
//                        last modification : 23/10/2002                      //
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


#ifndef VECTOROPERATORMULTIPLYOPERATION_H
#define VECTOROPERATORMULTIPLYOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"


class AbstractOperator;
class Vector;


class VectorOperatorMultiplyOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the operator
  AbstractOperator* Operator;

  // vector to be multiplied by the operator
  Vector* SourceVector;
  // vector where the result has to be stored
  Vector* DestinationVector;  

 public:
  
  // constructor 
  //
  // operator = pointer to the operator to use
  // sourceVector = vector to be multiplied by the operator
  // destinationVector = vector where the result has to be stored
  VectorOperatorMultiplyOperation(AbstractOperator* theOperator, Vector* sourceVector, Vector* destinationVector);

  // copy constructor 
  //
  // operation = reference on operation to copy
  VectorOperatorMultiplyOperation(const VectorOperatorMultiplyOperation& operation);

  // constructor from a master node information
  //
  // operator = pointer to the operator to use
  // architecture = pointer to the distributed architecture to use for communications
  VectorOperatorMultiplyOperation(AbstractOperator* theOperator, SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~VectorOperatorMultiplyOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetDestinationVector (Vector* DestinationVector);

  // get destination vector 
  // 
  // return value = pointer to destination vector
  Vector* GetDestinationVector ();

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();

 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
};

// get destination vector 
// 
// return value = pointer to destination vector

inline Vector* VectorOperatorMultiplyOperation::GetDestinationVector ()
{
  return this->DestinationVector;
}

#endif
