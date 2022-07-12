////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of multiple hamiltonian-vector multiplication operation        //
//                                                                            //
//                        last modification : 15/03/2005                      //
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


#ifndef MULTIPLEVECTORHAMILTONIANMULTIPLYOPERATION_H
#define MULTIPLEVECTORHAMILTONIANMULTIPLYOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"


class AbstractHamiltonian;
class Vector;
class RealVector;
class ComplexVector;


class MultipleVectorHamiltonianMultiplyOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the hamiltonian
  AbstractHamiltonian* Hamiltonian;

  // array of real vectors to be multiplied by the hamiltonian
  RealVector* RealSourceVectors;
  // array of real vectors where the result has to be stored
  RealVector* RealDestinationVectors;  
  // array of complex vectors to be multiplied by the hamiltonian
  ComplexVector* ComplexSourceVectors;
  // array of complex vectors where the result has to be stored
  ComplexVector* ComplexDestinationVectors;  
  // number of vectors that have to be evaluated together
  int NbrVectors;

 public:
  
  // constructor for real vectors
  //
  // hamiltonian = pointer to the hamiltonian to use
  // sourceVectors = array of vectors to be multiplied by the hamiltonian
  // destinationVectors = array of vectors where the result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  MultipleVectorHamiltonianMultiplyOperation(AbstractHamiltonian* hamiltonian, RealVector* sourceVectors, RealVector* destinationVectors, int nbrVectors);

  // constructor for complex vectors
  //
  // hamiltonian = pointer to the hamiltonian to use
  // sourceVectors = array of vectors to be multiplied by the hamiltonian
  // destinationVectors = array of vectors where the result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  MultipleVectorHamiltonianMultiplyOperation(AbstractHamiltonian* hamiltonian, ComplexVector* sourceVectors, ComplexVector* destinationVectors, int nbrVectors);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MultipleVectorHamiltonianMultiplyOperation(const MultipleVectorHamiltonianMultiplyOperation& operation);
  
  // constructor from a master node information
  //
  // hamiltonian = pointer to the hamiltonian to use
  // architecture = pointer to the distributed architecture to use for communications
  MultipleVectorHamiltonianMultiplyOperation(AbstractHamiltonian* hamiltonian, SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~MultipleVectorHamiltonianMultiplyOperation();
  
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

inline int MultipleVectorHamiltonianMultiplyOperation::GetNbrVectors()
{
  return this->NbrVectors;
}


// get destination real vector array
// 
// return value = array of destination vectors

inline RealVector* MultipleVectorHamiltonianMultiplyOperation::GetDestinationRealVectors ()
{
  return this->RealDestinationVectors;
}
     
// get destination complex vector array
// 
// return value = array of destination vectors

inline ComplexVector* MultipleVectorHamiltonianMultiplyOperation::GetDestinationComplexVectors ()
{
  return this->ComplexDestinationVectors;
}

#endif
