////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of operation to evaluate the characteristic polynomial       //
//                           of a long integer matrix                         //
//                                                                            //
//                        last modification : 13/01/2022                      //
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


#ifndef LONGINTEGERMATRIXCHARACTERISTICPOLYNOMIALOPERATION_H
#define LONGINTEGERMATRIXCHARACTERISTICPOLYNOMIALOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Matrix/LongIntegerMatrix.h"


class LongIntegerMatrixCharacteristicPolynomialOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the matrix
  LongIntegerMatrix* SourceMatrix;

  // number of non-zero matrix element per row
  int* NbrNonZeroMatrixElements;
  // positions of the non-zero matrix element per row
  int** NonZeroMatrixElementPositions;

  // pointer to characteristic polynomial
#ifdef __GMP__
  mpz_t* CharacteristicPolynomial;
#else
  LONGLONG* CharacteristicPolynomial;
#endif

  // two temporary matrices used during the evaluation
  LongIntegerMatrix TemporaryMatrix1;
  LongIntegerMatrix TemporaryMatrix2;

  
 public:
  
  // constructor 
  //
  // sourceMatrix = pointer to the matrix for which the characteric polynomial should be evaluated
  LongIntegerMatrixCharacteristicPolynomialOperation(LongIntegerMatrix* sourceMatrix);

  // constructor providing the sparse structure of the input matrix
  //
  // sourceMatrix = pointer to the matrix for which the characteric polynomial should be evaluated
  // nbrNonZeroMatrixElements = number of non-zero matrix element per row
  // nonZeroMatrixElementPositions = positions of the non-zero matrix element per row
  LongIntegerMatrixCharacteristicPolynomialOperation (LongIntegerMatrix* sourceMatrix, int* nbrNonZeroMatrixElements, int** nonZeroMatrixElementPositions);

  // copy constructor 
  //
  // operation = reference on operation to copy
  LongIntegerMatrixCharacteristicPolynomialOperation(const LongIntegerMatrixCharacteristicPolynomialOperation& operation);
  
  // destructor
  //
  ~LongIntegerMatrixCharacteristicPolynomialOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
  // get the characteristic polynomial
  // 
  // return value = pointer to characteristic polynomial
#ifdef __GMP__
  mpz_t* GetCharacteristicPolynomial();
#else
  LONGLONG* GetCharacteristicPolynomial();
#endif
  
 protected:

  // apply operation for mono processor architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture);
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // common code for the algorithm initialization
  //
  // trace = pointer to the trace
#ifdef __GMP__
  void AlgorithmInitialization(mpz_t* trace);
#else
  void AlgorithmInitialization(LONGLONG* trace);
#endif
  
};

// get the characteristic polynomial
// 
// return value = pointer to characteristic polynomial

#ifdef __GMP__
inline mpz_t* LongIntegerMatrixCharacteristicPolynomialOperation::GetCharacteristicPolynomial()
#else
inline LONGLONG* LongIntegerMatrixCharacteristicPolynomialOperation::GetCharacteristicPolynomial()
#endif  
{
  return this->CharacteristicPolynomial;
}

#endif
