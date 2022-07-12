////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of matrix full diagonalization operation             //
//                                                                            //
//                        last modification : 23/06/2016                      //
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


#ifndef MATRIXFULLDIAGONALIZEOPERATION_H
#define MATRIXFULLDIAGONALIZEOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"



class MatrixFullDiagonalizeOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the matrix to diagonalize (real version)
  RealSymmetricMatrix* InitialRealMatrix;
  // pointer to the matrix to diagonalize (complex veersion)
  HermitianMatrix* InitialComplexMatrix;

  // execution time measured in RawApply
  double ExecutionTime;

  // true if the matrix is only stored on the master node
  bool DistributeFlag;

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
  // matrix = pointer to the matrix to diagonalize
  // distributeFlag = true if the matrix is only stored on the master node
  // eigenstateFlag = true if the eigenstates have to be computed
  // nbrEigenstates = number of eigenstates that have to be computed (<=0 if all eigenstates have to be computed)
  MatrixFullDiagonalizeOperation(RealSymmetricMatrix* matrix, bool distributeFlag = false, bool eigenstateFlag = false, int nbrEigenstates = 0);

  // constructor 
  //
  // matrix = pointer to the matrix to diagonalize
  // distributeFlag = true if the matrix is only stored on the master node
  // eigenstateFlag = true if the eigenstates have to be computed
  // nbrEigenstates = number of eigenstates that have to be computed (<=0 if all eigenstates have to be computed)
  MatrixFullDiagonalizeOperation(HermitianMatrix* matrix, bool distributeFlag = false, bool eigenstateFlag = false, int nbrEigenstates = 0);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MatrixFullDiagonalizeOperation(const MatrixFullDiagonalizeOperation& operation);

  // constructor from a master node information
  //
  // hamiltonian = pointer to the hamiltonian to use
  // architecture = pointer to the distributed architecture to use for communications
  MatrixFullDiagonalizeOperation(SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~MatrixFullDiagonalizeOperation();
  
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
  RealDiagonalMatrix GetDiagonalizedMatrix();

  // get the hamiltonian eigenstates
  //
  // eigenstates = reference on the matrix where the eigenstates will be stored
  void GetMatrixEigenstates(RealMatrix& eigenstates);

  // get the hamiltonian eigenstates
  //
  // eigenstates = reference on the matrix where the eigenstates will be stored
  void GetMatrixEigenstates(ComplexMatrix& eigenstates);

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

inline RealDiagonalMatrix MatrixFullDiagonalizeOperation::GetDiagonalizedMatrix()
{
  return this->DiagonalizedMatrix;
}

// get the hamiltonian eigenstates
//
// eigenstates = reference on the matrix where the eigenstates will be stored

inline void MatrixFullDiagonalizeOperation::GetMatrixEigenstates(RealMatrix& eigenstates)
{
  eigenstates = this->RealEigenstates;
}

// get the hamiltonian eigenstates
//
// eigenstates = reference on the matrix where the eigenstates will be stored

inline void MatrixFullDiagonalizeOperation::GetMatrixEigenstates(ComplexMatrix& eigenstates)
{
  eigenstates = this->ComplexEigenstates;
}

#endif
