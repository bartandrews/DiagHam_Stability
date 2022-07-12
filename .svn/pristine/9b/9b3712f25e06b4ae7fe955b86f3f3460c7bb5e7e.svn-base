////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hamiltonian full diagonalization operation           //
//                                                                            //
//                        last modification : 06/01/2012                      //
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


#ifndef HAMILTONIANFULLDIAGONALIZEOPERATION_H
#define HAMILTONIANFULLDIAGONALIZEOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"


class AbstractHamiltonian;
class Vector;


class HamiltonianFullDiagonalizeOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the hamiltonian
  AbstractHamiltonian* Hamiltonian;
  // use full hermitian structure of the hamiltonian
  bool UseHermitianFlag;

  // execution time measured in RawApply
  double ExecutionTime;

  // true if the hamiltonian is complex
  bool ComplexFlag;

  // true if the eigenstates have to be computed
  bool EigenstateFlag ;
  // number of eigenstates that have to be computed
  int NbrEigenstates; 

  // matrix where the eigenvalues will be stored
  RealDiagonalMatrix DiagonalizedMatrix;

  // matrix where the eigenstates will be stored (for a real hamiltonian)
  RealMatrix RealEigenstates;
  // matrix where the eigenstates will be stored (for a complex hamiltonian)
  ComplexMatrix ComplexEigenstates;


 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the hamiltonian to use
  // complexFlag = true if the hamiltonian is complex
  // eigenstateFlag = true if the eigenstates have to be computed
  // nbrEigenstates = number of eigenstates that have to be computed (<=0 if all eigenstates have to be computed)
  HamiltonianFullDiagonalizeOperation(AbstractHamiltonian* hamiltonian, bool complexFlag = false, bool eigenstateFlag = false, int nbrEigenstates = 0);

  // copy constructor 
  //
  // operation = reference on operation to copy
  HamiltonianFullDiagonalizeOperation(const HamiltonianFullDiagonalizeOperation& operation);

  // constructor from a master node information
  //
  // hamiltonian = pointer to the hamiltonian to use
  // architecture = pointer to the distributed architecture to use for communications
  HamiltonianFullDiagonalizeOperation(AbstractHamiltonian* hamiltonian, SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~HamiltonianFullDiagonalizeOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();

  // get the diagonalized hamiltonian
  //
  // return value = diagonalized hamiltonian
  RealDiagonalMatrix GetDiagonalizedHamiltonian();

  // get the hamiltonian eigenstates
  //
  // eigenstates = reference on the matrix where the eigenstates will be stored
  void GetHamiltonianEigenstates(RealMatrix& eigenstates);

  // get the hamiltonian eigenstates
  //
  // eigenstates = reference on the matrix where the eigenstates will be stored
  void GetHamiltonianEigenstates(ComplexMatrix& eigenstates);

 protected:

  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
};

// get the diagonalized hamiltonian
//
// return value = diagonalized hamiltonian

inline RealDiagonalMatrix HamiltonianFullDiagonalizeOperation::GetDiagonalizedHamiltonian()
{
  return this->DiagonalizedMatrix;
}

// get the hamiltonian eigenstates
//
// eigenstates = reference on the matrix where the eigenstates will be stored

inline void HamiltonianFullDiagonalizeOperation::GetHamiltonianEigenstates(RealMatrix& eigenstates)
{
  eigenstates = this->RealEigenstates;
}

// get the hamiltonian eigenstates
//
// eigenstates = reference on the matrix where the eigenstates will be stored

inline void HamiltonianFullDiagonalizeOperation::GetHamiltonianEigenstates(ComplexMatrix& eigenstates)
{
  eigenstates = this->ComplexEigenstates;
}

#endif
