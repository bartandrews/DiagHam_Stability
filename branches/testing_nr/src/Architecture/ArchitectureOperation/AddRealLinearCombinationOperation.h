////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of add linear combination operation                //
//                                                                            //
//                        last modification : 24/10/2002                      //
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


#ifndef ADDREALLINEARCOMBINATIONOPERATION_H
#define ADDREALLINEARCOMBINATIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"


class AddRealLinearCombinationOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // vector array containing the vectors that have to be added
  RealVector* SourceVector;
  // vector array containing pointers to the vectors that have to be added
  RealVector** SourceVectorByPointers;
  // matrix containing the vectors that have to be added (can be used instead of SourceVector)
  RealMatrix SourceVectorMatrix;
  // number of vector that have to be added
  int NbrVector;
  // coefficient of the linear combination
  double* Coefficients;

  // vector where the result has to be stored
  RealVector* DestinationVector;  

 public:
  
  // constructor 
  //
  // destinationVector = vector to which the linear combination has to be added
  // sourceVector = array containing the vectors that have to be added
  // nbrVector = number of vector that have to be added
  // coefficients = coefficient of the linear combination
  AddRealLinearCombinationOperation(RealVector* destinationVector, RealVector* sourceVector, int nbrVector, double* coefficients);

  // constructor 
  //
  // destinationVector = vector to which the linear combination has to be added
  // sourceVector = array containing the pointers to the vectors that have to be added
  // nbrVector = number of vector that have to be added
  // coefficients = coefficient of the linear combination
  AddRealLinearCombinationOperation(RealVector* destinationVector, RealVector** sourceVector, int nbrVector, double* coefficients);

  // constructor 
  //
  // destinationVector = vector to which the linear combination has to be added
  // sourceVector = matrix containing the vectors that have to be added
  // nbrVector = number of vector that have to be added
  // coefficients = coefficient of the linear combination
  AddRealLinearCombinationOperation(RealVector* destinationVector, RealMatrix& sourceVector, int nbrVector, double* coefficients);

  // copy constructor 
  //
  // operation = reference on operation to copy
  AddRealLinearCombinationOperation(const AddRealLinearCombinationOperation& operation);
  
  // constructor from a master node information
  //
  // architecture = pointer to the distributed architecture to use for communications
  AddRealLinearCombinationOperation(SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~AddRealLinearCombinationOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

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

inline Vector* AddRealLinearCombinationOperation::GetDestinationVector ()
{
  return this->DestinationVector;
}

#endif
