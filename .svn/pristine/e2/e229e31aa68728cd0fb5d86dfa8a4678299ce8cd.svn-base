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
//                        last modification : 26/05/2003                      //
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


#ifndef ADDCOMPLEXLINEARCOMBINATIONOPERATION_H
#define ADDCOMPLEXLINEARCOMBINATIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"


class AddComplexLinearCombinationOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // vector array containing the vectors that have to be added
  ComplexVector* SourceVector;
  // vector array containing pointers to the vectors that have to be added
  ComplexVector** SourceVectorByPointers;
  // matrix containing the vectors that have to be added (can be used instead of SourceVector)
  ComplexMatrix SourceVectorMatrix;
  // number of vector that have to be added
  int NbrVector;
  // coefficient of the linear combination (null if coefficients are real)
  Complex* Coefficients;
  // coefficient of the linear combination (null if coefficients are complex)
  double* RealCoefficients;

  // vector where the result has to be stored
  ComplexVector* DestinationVector;  

 public:
  
  // constructor 
  //
  // destinationVector = vector to which the linear combination has to be added
  // sourceVector = array containing the vectors that have to be added
  // nbrVector = number of vector that have to be added
  // coefficients = coefficient of the linear combination
  AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexVector* sourceVector, int nbrVector, Complex* coefficients);

  // constructor 
  //
  // destinationVector = vector to which the linear combination has to be added
  // sourceVector = array containing pointers to the vectors that have to be added
  // nbrVector = number of vector that have to be added
  // coefficients = coefficient of the linear combination
  AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexVector** sourceVector, int nbrVector, Complex* coefficients);

  // constructor 
  //
  // destinationVector = vector to which the linear combination has to be added
  // sourceVector = matrix containing the vectors that have to be added
  // nbrVector = number of vector that have to be added
  // coefficients = coefficient of the linear combination
  AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexMatrix& sourceVector, int nbrVector, Complex* coefficients);

  // constructor 
  //
  // destinationVector = vector to which the linear combination has to be added
  // sourceVector = array containing the vectors that have to be added
  // nbrVector = number of vector that have to be added
  // coefficients = coefficient of the linear combination
  AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexVector* sourceVector, int nbrVector, double* coefficients);

  // constructor 
  //
  // destinationVector = vector to which the linear combination has to be added
  // sourceVector = array containing poiinters to the vectors that have to be added
  // nbrVector = number of vector that have to be added
  // coefficients = coefficient of the linear combination
  AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexVector** sourceVector, int nbrVector, double* coefficients);

  // constructor 
  //
  // destinationVector = vector to which the linear combination has to be added
  // sourceVector = matrix containing the vectors that have to be added
  // nbrVector = number of vector that have to be added
  // coefficients = coefficient of the linear combination
  AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexMatrix& sourceVector, int nbrVector, double* coefficients);

  // copy constructor 
  //
  // operation = reference on operation to copy
  AddComplexLinearCombinationOperation(const AddComplexLinearCombinationOperation& operation);
  
  // destructor
  //
  ~AddComplexLinearCombinationOperation();
  
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
  
};

// get destination vector 
// 
// return value = pointer to destination vector

inline Vector* AddComplexLinearCombinationOperation::GetDestinationVector ()
{
  return this->DestinationVector;
}

#endif
