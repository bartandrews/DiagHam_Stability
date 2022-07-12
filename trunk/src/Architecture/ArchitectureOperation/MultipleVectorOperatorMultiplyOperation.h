////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of multiple operator-vector multiplication operation         //
//                                                                            //
//                        last modification : 19/07/2016                      //
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


#ifndef MULTIPLEVECTOROPERATORMULTIPLYOPERATION_H
#define MULTIPLEVECTOROPERATORMULTIPLYOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"


class AbstractOperator;
class Vector;
class RealVector;
class ComplexVector;
class PartialRealVector;
class PartialComplexVector;


class MultipleVectorOperatorMultiplyOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the operator
  AbstractOperator* Oper;

  // array of real vectors to be multiplied by the operator
  RealVector* RealSourceVectors;
  // array of real vectors where the result has to be stored
  RealVector* RealDestinationVectors;  
  // array of complex vectors to be multiplied by the operator
  ComplexVector* ComplexSourceVectors;
  // array of complex vectors where the result has to be stored
  ComplexVector* ComplexDestinationVectors;  

  // array of real partial vectors to be multiplied by the operator
  PartialRealVector* RealSourcePartialVectors;
  // array of real partial vectors where the result has to be stored
  PartialRealVector* RealDestinationPartialVectors;  
  // array of complex partial vectors to be multiplied by the operator
  PartialComplexVector* ComplexSourcePartialVectors;
  // array of complex partial vectors where the result has to be stored
  PartialComplexVector* ComplexDestinationPartialVectors;  

  // number of vectors that have to be evaluated together
  int NbrVectors;

  // execution time measured in RawApply
  double ExecutionTime;

 public:
  
  // constructor for real vectors
  //
  // oper = pointer to the operator to use
  // sourceVectors = array of vectors to be multiplied by the operator
  // destinationVectors = array of vectors where the result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  MultipleVectorOperatorMultiplyOperation(AbstractOperator* oper, RealVector* sourceVectors, RealVector* destinationVectors, int nbrVectors);

  // constructor for complex vectors
  //
  // oper = pointer to the operator to use
  // sourceVectors = array of vectors to be multiplied by the operator
  // destinationVectors = array of vectors where the result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  MultipleVectorOperatorMultiplyOperation(AbstractOperator* oper, ComplexVector* sourceVectors, ComplexVector* destinationVectors, int nbrVectors);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MultipleVectorOperatorMultiplyOperation(const MultipleVectorOperatorMultiplyOperation& operation);
  
  // constructor from a master node information
  //
  // oper = pointer to the operator to use
  // architecture = pointer to the distributed architecture to use for communications
  MultipleVectorOperatorMultiplyOperation(AbstractOperator* oper, SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~MultipleVectorOperatorMultiplyOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // set destination vectors 
  // 
  // destinationVectors = array of vector where the result has to be stored
  void SetDestinationVectors (RealVector* destinationVectors);

  // set destination vectors 
  // 
  // destinationVectors = array of vector where the result has to be stored
  void SetDestinationVectors (ComplexVector* destinationVectors);

  // get destination real vector array
  // 
  // return value = array of destination vectors
  RealVector* GetDestinationRealVectors ();

  // get destination complex vector array
  // 
  // return value = array of destination vectors
  ComplexVector* GetDestinationComplexVectors ();

  // get number of vectors that have to be evaluated together
  //
  // return value = number of vectors
  int GetNbrVectors ();

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

// get number of vectors that have to be evaluated together
//
// return value = number of vectors

inline int MultipleVectorOperatorMultiplyOperation::GetNbrVectors()
{
  return this->NbrVectors;
}


// get destination real vector array
// 
// return value = array of destination vectors

inline RealVector* MultipleVectorOperatorMultiplyOperation::GetDestinationRealVectors ()
{
  return this->RealDestinationVectors;
}
     
// get destination complex vector array
// 
// return value = array of destination vectors

inline ComplexVector* MultipleVectorOperatorMultiplyOperation::GetDestinationComplexVectors ()
{
  return this->ComplexDestinationVectors;
}

#endif
