////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of multiple complex scalar product operation            //
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


#ifndef MULTIPLECOMPLEXSCALARPRODUCTOPERATION_H
#define MULTIPLECOMPLEXSCALARPRODUCTOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"


class MultipleComplexScalarProductOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component of each partial scalar product
  int FirstComponent;
  // number of component to take into account for each partial scalar product
  int NbrComponents;

  // number of scalar products
  int NbrScalarProduct;

  // strategy used to do the scalar products
  int Strategy;

  // array containing the scalar product
  Complex* ScalarProducts;

  // array of vectors to use for the right hand side of the scalar product
  ComplexVector* RightVectors;
  // array of pointers to the vectors to use for the right hand side of the scalar product
  ComplexVector** RightVectorsByPointers;
  // real matrix where vectors to use for the right hand side of the scalar product are stored (can be used instead of RightVectors)
  ComplexMatrix RightVectorMatrix;
  
  // pointer to the vector to use for the left hand side of the scalar product
  ComplexVector* LeftVector;  

  // execution time measured in RawApply
  double ExecutionTime;

public:
  
  enum StrategyType
    {
      // each process achieves part of each scalar product
      VectorSubdivision = 0x01,
      // each process achieves full scalar product for a given group of vectors
      GroupSubdivision = 0x02
    };
  
  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = array of vectors to use for the right hand side of the scalar product
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  MultipleComplexScalarProductOperation(ComplexVector* leftVector, ComplexVector* rightVectors, int nbrScalarProduct, Complex* scalarProducts);

  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = array of  pointers to the vectors to use for the right hand side of the scalar product
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  MultipleComplexScalarProductOperation(ComplexVector* leftVector, ComplexVector** rightVectors, int nbrScalarProduct, Complex* scalarProducts);

  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = real matrix where vectors to use for the right hand side of the scalar product are stored
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  MultipleComplexScalarProductOperation(ComplexVector* leftVector, ComplexMatrix& rightVectors, int nbrScalarProduct, Complex* scalarProducts);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MultipleComplexScalarProductOperation(const MultipleComplexScalarProductOperation& operation);
  
  // destructor
  //
  ~MultipleComplexScalarProductOperation();
  
  // set the array where scalar products have to be stored
  //
  // scalarProducts = array where scalar products have to be stored
  void SetScalarProducts (Complex* scalarProducts);

  // get the array where scalar products have to be stored
  //
  // return value = array where scalar products have to be stored
  Complex* GetScalarProducts ();

  // set the strategy used to do the scalar products (per vector subdivision or per group)
  //
  // strategy = flag corresponding to the strategy
  void SetStrategy (int strategy);

  // get the strategy used to do the scalar products (per vector subdivision or per group)
  //
  // return value = flag corresponding to the strategy
  int GetStrategy ();

  // get the vector used as the left hand side of the scalar product
  //
  // return value = pointer to the vector
  ComplexVector* GetLeftVector();

  // set index range of scalar product that have to be calculated
  // 
  // firstComponent = index of the first component of each partial scalar product for per vector subdivision stategy, 
  //                  or index of the first scalar product to evaluate for per group subdivision stategy
  // nbrComponent = number of component to take into account for each partial scalar product for per vector subdivision stategy,
  //                or number of scalar products that have to be evaluated for per group subdivision stategy
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // get the number of scalar products that have to be evaluated
  // 
  // return value = number of scalar products
  int GetNbrScalarProduct ();

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

// get the number of scalar products that have to be evaluated
// 
// return value = number of scalar products

inline int MultipleComplexScalarProductOperation::GetNbrScalarProduct ()
{
  return this->NbrScalarProduct;
}

// get the array where scalar products have to be stored
//
// return value = array where scalar products have to be stored

inline Complex* MultipleComplexScalarProductOperation::GetScalarProducts ()
{
  return this->ScalarProducts;
}

// get the strategy used to do the scalar products (per vector subdivision or per group)
//
// return value = flag corresponding to the strategy

inline int MultipleComplexScalarProductOperation::GetStrategy ()
{
  return this->Strategy;
}

// get the vector used as the left hand side of the scalar product
//
// return value = pointer to the vector

inline ComplexVector* MultipleComplexScalarProductOperation::GetLeftVector()
{
  return this->LeftVector;
}

#endif
