////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                       class of state creation from a MPS	              //
//                                                                            //
//                        last modification : 08/10/2012                      //
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

#include "config.h"
#include "Architecture/ArchitectureOperation/FQHEMPSEvaluateCFTOperation.h"
#include "Vector/RealVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Matrix/SparseComplexMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSN1SuperconformalMatrix.h"


#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor to compute the CFT matrix elements
//
// mPSMatrix = pointer to the MPS matrix 
// u1BosonBasis = array that contains the Hilbert space for the partition at each level
// leftLevel = level for the left state
// rightLevel = level for the right state
// centralCharge12 =value of the central charge divided by 12
// weightLeft = conformal weight of the left state at level 0 
// weightRight = weight of the right state  at level 0
// WeightMatrixElement = weight of the field whose matrix elements have to be evaluated
// previousMatrixElements = array where the already computed matrix element are stored
// nbrLeftPreviousMatrixElements = number of entry of the PreviousMatrixElements first index
// nbrRightPreviousMatrixElements = number of entry of the PreviousMatrixElements second index
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode


FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(FQHEMPSClustered2RMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel, int rightLevel,
							 const LongRational& centralCharge12, const LongRational& weightLeft, 
							 const LongRational& weightRight, const LongRational& weightMatrixElement,
							 LongRationalMatrix** previousRationalMatrixElements, int nbrLeftPreviousMatrixElements, int nbrRightPreviousMatrixElements,
							 int nbrMPIStage, int nbrSMPStage)
{
  this->MPSMatrix = mPSMatrix;
  this->SupersymmetricFlag = false;
  this->U1BosonBasis = u1BosonBasis;
  this->LeftLevel = leftLevel;
  this->RightLevel = rightLevel;
  this->CentralCharge12 = centralCharge12;
  this->InvCentralCharge3 = 0l;
  this->UseRationalFlag = true;
  this->WeightLeft = weightLeft;
  this->WeightRight = weightRight;
  this->WeightMatrixElement = weightMatrixElement;
  this->CentralCharge12Numerical = 0.0;
  this->InvCentralCharge3Numerical = 0.0;
  this->WeightLeftNumerical = 0.0;
  this->WeightRightNumerical = 0.0;
  this->WeightMatrixElementNumerical = 0.0;
  this->PreviousRationalOverlapMatrices = 0;
  this->PreviousMatrixElements = 0;
  this->PreviousOverlapMatrices = 0;
  this->PreviousRationalMatrixElements = previousRationalMatrixElements;
  this->NbrLeftPreviousMatrixElements = nbrLeftPreviousMatrixElements;
  this->NbrRightPreviousMatrixElements = nbrRightPreviousMatrixElements;
  this->FirstComponent = 0;
  this->OverlapMatrixFlag = false;
  this->NbrComponent = this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() * this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
  this->RationalMatrixElements = LongRationalMatrix (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension(), U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension(), true);
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}

// constructor to compute the CFT overlap matrix
//
// mPSMatrix = pointer to the MPS matrix 
// u1BosonBasis = array that contains the Hilbert space for the partition at each level
// leftLevel = level for the left or right state
// centralCharge12 =value of the central charge divided by 12
// weightLeft = conformal weight of the left or right state at level 0 
// previousOverlapMatrices = array where the already computed overlap matrices are stored
// nbrPreviousOverlapMatrices = number of entry of the PreviousMatrixElements first index
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode

FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(FQHEMPSClustered2RMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel,
							 const LongRational& centralCharge12, const LongRational& weightLeft,
							 LongRationalMatrix* previousOverlapMatrices, int nbrPreviousOverlapMatrices,
							 int nbrMPIStage, int nbrSMPStage)
{
  this->MPSMatrix = mPSMatrix;
  this->SupersymmetricFlag = false;
  this->U1BosonBasis = u1BosonBasis;
  this->LeftLevel = leftLevel;
  this->RightLevel = leftLevel;
  this->CentralCharge12 = centralCharge12;
  this->InvCentralCharge3 = 0l;
  this->UseRationalFlag = true;
  this->WeightLeft = weightLeft;
  this->WeightRight = weightLeft;
  this->WeightMatrixElement = 0l;
  this->CentralCharge12Numerical = 0.0;
  this->InvCentralCharge3Numerical = 0.0;
  this->WeightLeftNumerical = 0.0;
  this->WeightRightNumerical = 0.0;
  this->WeightMatrixElementNumerical = 0.0;
  this->PreviousMatrixElements = 0;
  this->PreviousOverlapMatrices = 0;
  this->PreviousRationalMatrixElements = 0;
  this->PreviousRationalOverlapMatrices = previousOverlapMatrices;
  this->NbrLeftPreviousMatrixElements = nbrPreviousOverlapMatrices;
  this->NbrRightPreviousMatrixElements = 0;
  this->FirstComponent = 0;
  this->OverlapMatrixFlag = true;
  this->NbrComponent = (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() * (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() + 1)) / 2;
  this->RationalMatrixElements = LongRationalMatrix (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension(), U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension(), true);
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}

// constructor to compute the CFT matrix elements for the supersymmetric case
//
// mPSMatrix = pointer to the MPS matrix 
// u1BosonBasis = array that contains the Hilbert space for the partition at each level
// leftLevel = level for the left state
// rightLevel = level for the right state
// centralCharge12 =value of the central charge divided by 12
// invCentralCharge3 = 3 / (central charge)
// weightLeft = conformal weight of the left state at level 0 
// previousMatrixElements = array where the already computed matrix element are stored
// nbrLeftPreviousMatrixElements = number of entry of the PreviousMatrixElements first index
// nbrRightPreviousMatrixElements = number of entry of the PreviousMatrixElements second index
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode


FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(FQHEMPSN1SuperconformalMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel, int rightLevel,
							 const LongRational& centralCharge12, const LongRational& invCentralCharge3, const LongRational& weightLeft, 
							 LongRationalMatrix** previousRationalMatrixElements, int nbrLeftPreviousMatrixElements, 
							 int nbrRightPreviousMatrixElements,
							 int nbrMPIStage, int nbrSMPStage)
{
  this->MPSMatrix = mPSMatrix;
  this->SupersymmetricFlag = true;
  this->U1BosonBasis = u1BosonBasis;
  this->LeftLevel = leftLevel;
  this->RightLevel = rightLevel;
  this->CentralCharge12 = centralCharge12;
  this->InvCentralCharge3 = invCentralCharge3;
  this->UseRationalFlag = true;
  this->WeightLeft = weightLeft;
  this->WeightRight = weightLeft;
  this->WeightMatrixElement = 0l;
  this->CentralCharge12Numerical = 0.0;
  this->InvCentralCharge3Numerical = 0.0;
  this->WeightLeftNumerical = 0.0;
  this->WeightRightNumerical = 0.0;
  this->WeightMatrixElementNumerical = 0.0;
  this->PreviousRationalOverlapMatrices = 0;
  this->PreviousMatrixElements = 0;
  this->PreviousOverlapMatrices = 0;
  this->PreviousRationalMatrixElements = previousRationalMatrixElements;
  this->NbrLeftPreviousMatrixElements = nbrLeftPreviousMatrixElements;
  this->NbrRightPreviousMatrixElements = nbrRightPreviousMatrixElements;
  this->FirstComponent = 0;
  this->OverlapMatrixFlag = false;
  this->NbrComponent = this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() * this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
  this->RationalMatrixElements = LongRationalMatrix (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension(), U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension(), true);
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}

// constructor to compute the CFT overlap matrix for the supersymmetric case
//
// mPSMatrix = pointer to the MPS matrix 
// u1BosonBasis = array that contains the Hilbert space for the partition at each level
// leftLevel = level for the left or right state
// centralCharge12 = value of the central charge divided by 12
// invCentralCharge3 = 3 / (central charge)
// weightLeft = conformal weight of the left or right state at level 0 
// previousOverlapMatrices = array where the already computed overlap matrices are stored
// nbrPreviousOverlapMatrices = number of entry of the PreviousMatrixElements first index
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode

FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(FQHEMPSN1SuperconformalMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel,
							 const LongRational& centralCharge12, const LongRational& invCentralCharge3, const LongRational& weightLeft,
							 LongRationalMatrix* previousOverlapMatrices, int nbrPreviousOverlapMatrices,
							 int nbrMPIStage, int nbrSMPStage)
{
  this->MPSMatrix = mPSMatrix;
  this->SupersymmetricFlag = true;
  this->U1BosonBasis = u1BosonBasis;
  this->LeftLevel = leftLevel;
  this->RightLevel = leftLevel;
  this->CentralCharge12 = centralCharge12;
  this->InvCentralCharge3 = invCentralCharge3;
  this->UseRationalFlag = true;
  this->WeightLeft = weightLeft;
  this->WeightRight = weightLeft;
  this->WeightMatrixElement = 0l;
  this->CentralCharge12Numerical = 0.0;
  this->InvCentralCharge3Numerical = 0.0;
  this->WeightLeftNumerical = 0.0;
  this->WeightRightNumerical = 0.0;
  this->WeightMatrixElementNumerical = 0.0;
  this->PreviousMatrixElements = 0;
  this->PreviousOverlapMatrices = 0;
  this->PreviousRationalMatrixElements = 0;
  this->PreviousRationalOverlapMatrices = previousOverlapMatrices;
  this->NbrLeftPreviousMatrixElements = nbrPreviousOverlapMatrices;
  this->NbrRightPreviousMatrixElements = 0;
  this->FirstComponent = 0;
  this->OverlapMatrixFlag = true;
  this->NbrComponent = (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() * (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() + 1)) / 2;
  this->RationalMatrixElements = LongRationalMatrix (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension(), U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension(), true);
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}

// constructor to compute the CFT matrix elements, using double instead of rational numbers
//
// mPSMatrix = pointer to the MPS matrix 
// u1BosonBasis = array that contains the Hilbert space for the partition at each level
// leftLevel = level for the left state
// rightLevel = level for the right state
// centralCharge12 =value of the central charge divided by 12
// weightLeft = conformal weight of the left state at level 0 
// weightRight = weight of the right state  at level 0
// WeightMatrixElement = weight of the field whose matrix elements have to be evaluated
// previousMatrixElements = array where the already computed matrix element are stored
// nbrLeftPreviousMatrixElements = number of entry of the PreviousMatrixElements first index
// nbrRightPreviousMatrixElements = number of entry of the PreviousMatrixElements second index
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode


FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(FQHEMPSClustered2RMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel, int rightLevel,
							 double centralCharge12, double weightLeft, 
							 double weightRight, double weightMatrixElement,
							 RealMatrix** previousMatrixElements, int nbrLeftPreviousMatrixElements, int nbrRightPreviousMatrixElements,
							 int nbrMPIStage, int nbrSMPStage)
{
  this->MPSMatrix = mPSMatrix;
  this->SupersymmetricFlag = false;
  this->U1BosonBasis = u1BosonBasis;
  this->LeftLevel = leftLevel;
  this->RightLevel = rightLevel;
  this->CentralCharge12 = 0l;
  this->InvCentralCharge3 = 0l;
  this->UseRationalFlag = false;
  this->WeightLeft = 0l;
  this->WeightRight = 0l;
  this->WeightMatrixElement = 0l;
  this->WeightLeftNumerical = weightLeft;
  this->WeightRightNumerical = weightRight;
  this->WeightMatrixElementNumerical = weightMatrixElement;
  this->CentralCharge12Numerical = centralCharge12;  
  this->InvCentralCharge3Numerical = 0.0;
  this->PreviousRationalMatrixElements = 0;
  this->PreviousRationalOverlapMatrices = 0;
  this->PreviousMatrixElements = previousMatrixElements;
  this->PreviousOverlapMatrices = 0;
  this->NbrLeftPreviousMatrixElements = nbrLeftPreviousMatrixElements;
  this->NbrRightPreviousMatrixElements = nbrRightPreviousMatrixElements;
  this->FirstComponent = 0;
  this->OverlapMatrixFlag = false;
  this->NbrComponent = this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() * this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
  this->MatrixElements = RealMatrix (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension(), this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension(), true);
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}

// constructor to compute the CFT overlap matrix, using double instead of rational numbers
//
// mPSMatrix = pointer to the MPS matrix 
// u1BosonBasis = array that contains the Hilbert space for the partition at each level
// leftLevel = level for the left or right state
// centralCharge12 =value of the central charge divided by 12
// weightLeft = conformal weight of the left or right state at level 0 
// previousOverlapMatrices = array where the already computed overlap matrices are stored
// nbrPreviousOverlapMatrices = number of entry of the PreviousMatrixElements first index
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode


FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(FQHEMPSClustered2RMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel,
							 double centralCharge12, double weightLeft,
							 RealSymmetricMatrix* previousOverlapMatrices, int nbrPreviousOverlapMatrices,
							 int nbrMPIStage, int nbrSMPStage)
{
  this->MPSMatrix = mPSMatrix;
  this->SupersymmetricFlag = false;
  this->U1BosonBasis = u1BosonBasis;
  this->LeftLevel = leftLevel;
  this->RightLevel = leftLevel;
  this->CentralCharge12 = 0l;
  this->InvCentralCharge3 = 0l;
  this->UseRationalFlag = false;
  this->WeightLeft = 0l;
  this->WeightRight = 0l;
  this->WeightMatrixElement = 0l;
  this->WeightLeftNumerical = weightLeft;
  this->WeightRightNumerical = weightLeft;
  this->WeightMatrixElementNumerical = 0.0;
  this->CentralCharge12Numerical = centralCharge12;  
  this->InvCentralCharge3Numerical = 0.0;
  this->PreviousRationalMatrixElements = 0;
  this->PreviousRationalOverlapMatrices = 0;
  this->PreviousMatrixElements = 0;
  this->PreviousOverlapMatrices = previousOverlapMatrices;
  this->NbrLeftPreviousMatrixElements = nbrPreviousOverlapMatrices;
  this->NbrRightPreviousMatrixElements = 0;
  this->FirstComponent = 0;
  this->OverlapMatrixFlag = true;
  this->NbrComponent = (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() * (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() + 1)) / 2;
  this->OverlapMatrix = RealSymmetricMatrix (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension(), true);
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}


// constructor to compute the CFT matrix elements for the supersymmetric case, using double instead of rational numbers
//
// mPSMatrix = pointer to the MPS matrix 
// u1BosonBasis = array that contains the Hilbert space for the partition at each level
// leftLevel = level for the left state
// rightLevel = level for the right state
// centralCharge12 =value of the central charge divided by 12
// invCentralCharge3 = 3 / (central charge)
// weightLeft = conformal weight of the left state at level 0 
// previousMatrixElements = array where the already computed matrix element are stored
// nbrLeftPreviousMatrixElements = number of entry of the PreviousMatrixElements first index
// nbrRightPreviousMatrixElements = number of entry of the PreviousMatrixElements second index
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode


FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(FQHEMPSN1SuperconformalMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel, int rightLevel,
							 double centralCharge12, double invCentralCharge3, double weightLeft, 
							 RealMatrix** previousMatrixElements, int nbrLeftPreviousMatrixElements, 
							 int nbrRightPreviousMatrixElements,
							 int nbrMPIStage, int nbrSMPStage)
{
  this->MPSMatrix = mPSMatrix;
  this->SupersymmetricFlag = true;
  this->U1BosonBasis = u1BosonBasis;
  this->LeftLevel = leftLevel;
  this->RightLevel = rightLevel;
  this->CentralCharge12 = 0l;
  this->InvCentralCharge3 = 0l;
  this->UseRationalFlag = false;
  this->WeightLeft = 0l;
  this->WeightRight = 0l;
  this->WeightMatrixElement = 0l;
  this->CentralCharge12Numerical = centralCharge12;
  this->InvCentralCharge3Numerical = invCentralCharge3;
  this->WeightLeftNumerical = weightLeft;
  this->WeightRightNumerical = weightLeft;
  this->WeightMatrixElementNumerical = 0.0;
  this->PreviousRationalOverlapMatrices = 0;
  this->PreviousMatrixElements = previousMatrixElements;
  this->PreviousOverlapMatrices = 0;
  this->PreviousRationalMatrixElements = 0;
  this->NbrLeftPreviousMatrixElements = nbrLeftPreviousMatrixElements;
  this->NbrRightPreviousMatrixElements = nbrRightPreviousMatrixElements;
  this->FirstComponent = 0;
  this->OverlapMatrixFlag = false;
  this->NbrComponent = this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() * this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
  this->MatrixElements = RealMatrix (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension(), U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension(), true);
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}

// constructor to compute the CFT overlap matrix for the supersymmetric case, using double instead of rational numbers
//
// mPSMatrix = pointer to the MPS matrix 
// u1BosonBasis = array that contains the Hilbert space for the partition at each level
// leftLevel = level for the left or right state
// centralCharge12 = value of the central charge divided by 12
// invCentralCharge3 = 3 / (central charge)
// weightLeft = conformal weight of the left or right state at level 0 
// previousOverlapMatrices = array where the already computed overlap matrices are stored
// nbrPreviousOverlapMatrices = number of entry of the PreviousMatrixElements first index
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode

FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(FQHEMPSN1SuperconformalMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel,
							 double centralCharge12, double invCentralCharge3, double weightLeft,
							 RealSymmetricMatrix* previousOverlapMatrices, int nbrPreviousOverlapMatrices,
							 int nbrMPIStage, int nbrSMPStage)
{
  this->MPSMatrix = mPSMatrix;
  this->SupersymmetricFlag = true;
  this->U1BosonBasis = u1BosonBasis;
  this->LeftLevel = leftLevel;
  this->RightLevel = leftLevel;
  this->CentralCharge12 = 0l;
  this->InvCentralCharge3 = 0l;
  this->UseRationalFlag = false;
  this->WeightLeft = 0l;
  this->WeightRight = 0l;
  this->WeightMatrixElement = 0l;
  this->CentralCharge12Numerical = centralCharge12;
  this->InvCentralCharge3Numerical = invCentralCharge3;
  this->WeightLeftNumerical = weightLeft;
  this->WeightRightNumerical = weightLeft;
  this->WeightMatrixElementNumerical = 0.0;
  this->PreviousMatrixElements = 0;
  this->PreviousOverlapMatrices = previousOverlapMatrices;
  this->PreviousRationalMatrixElements = 0;
  this->PreviousRationalOverlapMatrices = 0;
  this->NbrLeftPreviousMatrixElements = nbrPreviousOverlapMatrices;
  this->NbrRightPreviousMatrixElements = 0;
  this->FirstComponent = 0;
  this->OverlapMatrixFlag = true;
  this->NbrComponent = (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() * (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() + 1)) / 2;
  this->OverlapMatrix = RealSymmetricMatrix (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension(), true);
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}

// copy constructor 
//
// operation = reference on operation to copy

FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(const FQHEMPSEvaluateCFTOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->SupersymmetricFlag = operation.SupersymmetricFlag;
  this->OverlapMatrixFlag = operation.OverlapMatrixFlag;
  this->MPSMatrix = operation.MPSMatrix;
  this->U1BosonBasis = operation.U1BosonBasis;
  this->LeftLevel = operation.LeftLevel;
  this->RightLevel = operation.RightLevel;
  this->UseRationalFlag = operation.UseRationalFlag;
  this->CentralCharge12 = operation.CentralCharge12;
  this->InvCentralCharge3 = operation.InvCentralCharge3;
  this->WeightLeft = operation.WeightLeft;
  this->WeightRight = operation.WeightRight;
  this->WeightMatrixElement = operation.WeightMatrixElement;
  this->CentralCharge12Numerical = operation.CentralCharge12Numerical;
  this->InvCentralCharge3Numerical = operation.InvCentralCharge3Numerical;
  this->WeightLeftNumerical = operation.WeightLeftNumerical;
  this->WeightRightNumerical = operation.WeightRightNumerical;
  this->WeightMatrixElementNumerical = operation.WeightMatrixElementNumerical;
  this->PreviousRationalOverlapMatrices = operation.PreviousRationalOverlapMatrices;
  this->PreviousRationalMatrixElements = operation.PreviousRationalMatrixElements;
  this->PreviousMatrixElements = operation.PreviousMatrixElements;
  this->PreviousOverlapMatrices = operation.PreviousOverlapMatrices;
  this->NbrLeftPreviousMatrixElements = operation.NbrLeftPreviousMatrixElements;
  this->NbrRightPreviousMatrixElements = operation.NbrRightPreviousMatrixElements;
  this->RationalMatrixElements = operation.RationalMatrixElements;
  this->MatrixElements = operation.MatrixElements;
  this->OverlapMatrix = operation.OverlapMatrix;
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = operation.NbrMPIStage;
  this->NbrSMPStage = operation.NbrSMPStage;
  this->SMPStages = operation.SMPStages;
}

// destructor
//

FQHEMPSEvaluateCFTOperation::~FQHEMPSEvaluateCFTOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHEMPSEvaluateCFTOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHEMPSEvaluateCFTOperation::Clone()
{
  return new FQHEMPSEvaluateCFTOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHEMPSEvaluateCFTOperation::RawApplyOperation()
{
  if (this->NbrComponent == 0)
    return true;
  int LastComponent = this->NbrComponent + this->FirstComponent;
  if (this->SupersymmetricFlag == false)
    {      
      if (this->OverlapMatrixFlag == false)
	{
	  long* Partition = new long[this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() + this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension() + 1];
	  unsigned long* TmpPartition = new unsigned long [this->MPSMatrix->GetTruncationLevel() + 2];
	  unsigned long* TemporaryOccupationNumber = new unsigned long [this->MPSMatrix->GetTruncationLevel() + 2];
	  if (this->UseRationalFlag == true)
	    {
	      for (int Index = this->FirstComponent; Index < LastComponent; ++Index)
		{
		  int n = Index / this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  int m = Index % this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  int PartitionLength = 0;
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(n, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			++PartitionLength;
		      }
		  int Position = PartitionLength;
		  PartitionLength = 0;
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[Position - PartitionLength - 1] = (long) k;
			++PartitionLength;
		      }
		  this->U1BosonBasis[this->RightLevel]->GetOccupationNumber(m, TmpPartition);	    
		  for (int k = 1; k <= this->RightLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[PartitionLength] = -(long) k;
			++PartitionLength;		  
		      }
		  LongRational Tmp = this->MPSMatrix->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, this->CentralCharge12, 
										      this->WeightLeft, this->WeightRight, this->WeightMatrixElement,
										      this->PreviousRationalMatrixElements, this->NbrLeftPreviousMatrixElements, this->NbrRightPreviousMatrixElements, 
										      this->U1BosonBasis, TemporaryOccupationNumber);
		  this->RationalMatrixElements.SetMatrixElement(n, m, Tmp);
		}
	    }
	  else
	    {
	      for (int Index = this->FirstComponent; Index < LastComponent; ++Index)
		{
		  int n = Index / this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  int m = Index % this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  int PartitionLength = 0;
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(n, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			++PartitionLength;
		      }
		  int Position = PartitionLength;
		  PartitionLength = 0;
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[Position - PartitionLength - 1] = (long) k;
			++PartitionLength;
		      }
		  this->U1BosonBasis[this->RightLevel]->GetOccupationNumber(m, TmpPartition);	    
		  for (int k = 1; k <= this->RightLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[PartitionLength] = -(long) k;
			++PartitionLength;		  
		      }
		  double Tmp = this->MPSMatrix->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, this->CentralCharge12Numerical, 
										this->WeightLeftNumerical, this->WeightRightNumerical, this->WeightMatrixElementNumerical,
										this->PreviousMatrixElements, this->NbrLeftPreviousMatrixElements, this->NbrRightPreviousMatrixElements, 
										this->U1BosonBasis, TemporaryOccupationNumber);
		  this->MatrixElements.SetMatrixElement(n, m, Tmp);
		}
	    }
	  delete[] TemporaryOccupationNumber;
	  delete[] TmpPartition;
	  delete[] Partition;
	  return true;
	}
      else
	{
	  long* Partition = new long[this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() + this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension() + 1];
	  unsigned long* TmpPartition = new unsigned long [this->MPSMatrix->GetTruncationLevel() + 2];
	  unsigned long* TemporaryOccupationNumber = new unsigned long [this->MPSMatrix->GetTruncationLevel() + 2];
	  if (this->UseRationalFlag == true)
	    {
	      for (int Index = this->FirstComponent; Index < LastComponent; ++Index)
		{
		  int n = 0;
		  int MaxN = this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  while ((n < MaxN) && (((n * (n +1)) / 2) <= Index))
		    ++n;
		  --n;
		  int m = Index - ((n * (n +1)) / 2);
		  int PartitionLength = 0;
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(n, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			++PartitionLength;
		      }
		  int Position = PartitionLength;
		  PartitionLength = 0;
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[Position - PartitionLength - 1] = (long) k;
			++PartitionLength;
		      }
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(m, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[PartitionLength] = -(long) k;
			++PartitionLength;		  
		      }
		  LongRational Tmp = this->MPSMatrix->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, this->WeightLeft,
											      this->PreviousRationalOverlapMatrices, this->NbrLeftPreviousMatrixElements, this->U1BosonBasis,
											      TemporaryOccupationNumber);
		  this->RationalMatrixElements.SetMatrixElement(m, n, Tmp);
		  if (n != m)
		    {
		      this->RationalMatrixElements.SetMatrixElement(n, m, Tmp);	      
		    }
		}
	    }
	  else
	    {
	      for (int Index = this->FirstComponent; Index < LastComponent; ++Index)
		{
		  int n = 0;
		  int MaxN = this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  while ((n < MaxN) && (((n * (n +1)) / 2) <= Index))
		    ++n;
		  --n;
		  int m = Index - ((n * (n +1)) / 2);
		  int PartitionLength = 0;
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(n, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
		    ++PartitionLength;
		      }
		  int Position = PartitionLength;
		  PartitionLength = 0;
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[Position - PartitionLength - 1] = (long) k;
			++PartitionLength;
		      }
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(m, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[PartitionLength] = -(long) k;
			++PartitionLength;		  
		      }
		  double Tmp = this->MPSMatrix->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12Numerical, this->WeightLeftNumerical,
											this->PreviousOverlapMatrices, this->NbrLeftPreviousMatrixElements, this->U1BosonBasis,
											TemporaryOccupationNumber);
		  this->OverlapMatrix.SetMatrixElement(m, n, Tmp);
		}
	    }
	  delete[] TemporaryOccupationNumber;
	  delete[] TmpPartition;
	  delete[] Partition;
	  return true;
	}
    }
  else
    {
      if (this->OverlapMatrixFlag == false)
	{
	  int EffectivePTruncation = 2 * this->MPSMatrix->GetTruncationLevel() + 3;
	  long* Partition = new long[this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() + this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension() + 1];
	  unsigned long* TmpPartition = new unsigned long [EffectivePTruncation + 2];
	  unsigned long* TemporaryOccupationNumber = new unsigned long [EffectivePTruncation + 2];
	  if (this->UseRationalFlag == true)
	    {
	      for (int Index = this->FirstComponent; Index < LastComponent; ++Index)
		{
		  int n = Index / this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  int m = Index % this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  int PartitionLength = 0;
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(n, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			++PartitionLength;
		      }
		  int Position = PartitionLength;
		  PartitionLength = 0;
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[Position - PartitionLength - 1] = (long) k;
			++PartitionLength;
		      }
		  long TmpWeight = this->RightLevel - this->LeftLevel;
		  Partition[Position] = TmpWeight;
		  if (TmpWeight > 0)
		    ++Position;
		  ++PartitionLength;
		  this->U1BosonBasis[this->RightLevel]->GetOccupationNumber(m, TmpPartition);	    
		  for (int k = 1; k <= this->RightLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[PartitionLength] = -(long) k;
			++PartitionLength;		  
		      }
		  LongRational Tmp = ((FQHEMPSN1SuperconformalMatrix*) this->MPSMatrix)->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, 
																 this->CentralCharge12, this->InvCentralCharge3, 
																 this->WeightLeft);
		  this->RationalMatrixElements.SetMatrixElement(n, m, Tmp);
		}
	    }
	  else
	    {
	      for (int Index = this->FirstComponent; Index < LastComponent; ++Index)
		{
		  int n = Index / this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  int m = Index % this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
		  int PartitionLength = 0;
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(n, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			++PartitionLength;
		      }
		  int Position = PartitionLength;
		  PartitionLength = 0;
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[Position - PartitionLength - 1] = (long) k;
			++PartitionLength;
		      }
		  long TmpWeight = this->RightLevel - this->LeftLevel;
		  Partition[Position] = TmpWeight;
		  if (TmpWeight > 0)
		    ++Position;
		  ++PartitionLength;
		  this->U1BosonBasis[this->RightLevel]->GetOccupationNumber(m, TmpPartition);	    
		  for (int k = 1; k <= this->RightLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[PartitionLength] = -(long) k;
			++PartitionLength;		  
		      }
		  double Tmp = ((FQHEMPSN1SuperconformalMatrix*) this->MPSMatrix)->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, 
															   this->CentralCharge12Numerical, this->InvCentralCharge3Numerical, 
															   this->WeightLeftNumerical);
		  this->MatrixElements.SetMatrixElement(n, m, Tmp);
		}
	    }

	}
      else
	{
	  int EffectivePTruncation = 2 * this->MPSMatrix->GetTruncationLevel() + 3;
	  long* Partition = new long[this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() + this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension() + 1];
	  unsigned long* TmpPartition = new unsigned long [EffectivePTruncation + 2];
	  unsigned long* TemporaryOccupationNumber = new unsigned long [EffectivePTruncation + 2];
	  if (this->UseRationalFlag == true)
	    {
	      for (int Index = this->FirstComponent; Index < LastComponent; ++Index)
		{
		  int n = 0;
		  int MaxN = this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension();
		  while ((n < MaxN) && (((n * (n +1)) / 2) <= Index))
		    ++n;
		  --n;
		  int m = Index - ((n * (n +1)) / 2);
		  int PartitionLength = 0;
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(n, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			++PartitionLength;
		      }
		  int Position = PartitionLength;
		  PartitionLength = 0;
		  for (int k = 2; k <= this->LeftLevel; k += 2)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[Position - PartitionLength - 1] = (long) k;
			++PartitionLength;
		      }
		    for (int k = 1; k <= this->LeftLevel; k += 2)
		      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			{
			  Partition[Position - PartitionLength - 1] = (long) k;
			  ++PartitionLength;
			}
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(m, TmpPartition);	    
		  for (int k = 2; k <= this->LeftLevel; k += 2)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[PartitionLength] = -(long) k;
			++PartitionLength;		  
		      }
		    for (int k = 1; k <= this->LeftLevel; k += 2)
		      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			{
			  Partition[PartitionLength] = -(long) k;
			  ++PartitionLength;		  
			}
		    LongRational Tmp = ((FQHEMPSN1SuperconformalMatrix*) this->MPSMatrix)->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, InvCentralCharge3, this->WeightLeft,
																   this->PreviousRationalOverlapMatrices, this->NbrLeftPreviousMatrixElements, this->U1BosonBasis,
																   TemporaryOccupationNumber);
		  this->RationalMatrixElements.SetMatrixElement(m, n, Tmp);
		  if (n != m)
		    {
		      this->RationalMatrixElements.SetMatrixElement(n, m, Tmp);	      
		    }
		}
	    }
	  else
	    {
	      for (int Index = this->FirstComponent; Index < LastComponent; ++Index)
		{
		  int n = 0;
		  int MaxN = this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension();
		  while ((n < MaxN) && (((n * (n +1)) / 2) <= Index))
		    ++n;
		  --n;
		  int m = Index - ((n * (n +1)) / 2);
		  int PartitionLength = 0;
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(n, TmpPartition);	    
		  for (int k = 1; k <= this->LeftLevel; ++k)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			++PartitionLength;
		      }
		  int Position = PartitionLength;
		  PartitionLength = 0;
		  for (int k = 2; k <= this->LeftLevel; k += 2)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[Position - PartitionLength - 1] = (long) k;
			++PartitionLength;
		      }
		    for (int k = 1; k <= this->LeftLevel; k += 2)
		      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			{
			  Partition[Position - PartitionLength - 1] = (long) k;
			  ++PartitionLength;
			}
		  this->U1BosonBasis[this->LeftLevel]->GetOccupationNumber(m, TmpPartition);	    
		  for (int k = 2; k <= this->LeftLevel; k += 2)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[PartitionLength] = -(long) k;
			++PartitionLength;		  
		      }
		  for (int k = 1; k <= this->LeftLevel; k += 2)
		    for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		      {
			Partition[PartitionLength] = -(long) k;
			++PartitionLength;		  
		      }
		  double Tmp = ((FQHEMPSN1SuperconformalMatrix*) this->MPSMatrix)->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, 
															   this->CentralCharge12Numerical, 
															   this->InvCentralCharge3Numerical, this->WeightLeftNumerical,
															   this->PreviousOverlapMatrices, this->NbrLeftPreviousMatrixElements, this->U1BosonBasis,
															   TemporaryOccupationNumber);
		  this->OverlapMatrix.SetMatrixElement(m, n, Tmp);
		  if (n != m)
		    {
		      this->OverlapMatrix.SetMatrixElement(n, m, Tmp);	      
		    }
		}
	    }
	}
    }
  return true;
}



// apply operation for SMP using round robin scheduling
//
//  architecture = instance of architecture class
// return value = true if no error occurs

bool FQHEMPSEvaluateCFTOperation::ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID)
{      
  int TmpNbrComponents = this->NbrComponent;
  int TmpFirstComponent = this->FirstComponent;
  
   int NbrStages = this->NbrSMPStage * architecture->GetNbrThreads();
   int StageIndex = 0;

   bool LockFlag = false;
   while (StageIndex < NbrStages) 
     {
       if (LockFlag == false)   
         {
           architecture->LockMutex();
           LockFlag = true;
         }
       StageIndex = this->SMPStages[0];
       if (StageIndex < NbrStages) 
         {         
           this->SMPStages[0]++;   
           architecture->UnLockMutex();
           LockFlag = false;
           this->SetIndicesRange(TmpFirstComponent + this->GetRankChunkStart(TmpNbrComponents, StageIndex,  NbrStages),  
				 this->GetRankChunkSize(TmpNbrComponents, StageIndex,  NbrStages));
           this->RawApplyOperation();
           ++StageIndex;
         }
     }
   if (LockFlag == true)
     {
       architecture->UnLockMutex();
       LockFlag = false;
     }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHEMPSEvaluateCFTOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  if (Step == 0)
    Step = this->NbrComponent;
  this->SMPStages[0] = 0;
  int TotalNbrComponent = this->FirstComponent + this->NbrComponent;
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHEMPSEvaluateCFTOperation** TmpOperations = new FQHEMPSEvaluateCFTOperation * [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHEMPSEvaluateCFTOperation *) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  architecture->SendJobsRoundRobin();

//   for( int i = 0; i < architecture->GetNbrThreads(); i++)
//     {
//       TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
//       TmpFirstComponent += Step;
//       if (TmpFirstComponent >= TotalNbrComponent)
// 	{
// 	  Step = 0;
// 	  TmpFirstComponent = TotalNbrComponent;
// 	}
//       else
// 	{
// 	  if ((TmpFirstComponent + Step) >= TotalNbrComponent)
// 	    {
// 	      Step = TotalNbrComponent - TmpFirstComponent;	      
// 	    }
// 	}
//     }
//  architecture->SendJobs();

  for (int i = 0; i < architecture->GetNbrThreads(); i++)
    {
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}

// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHEMPSEvaluateCFTOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__    
//   int Step = this->NbrComponent / architecture->GetNbrNodes();
//   int TmpFirstComponent = this->FirstComponent + (Step * architecture->GetNodeNbr());
//   int TmpNbrComponent = Step;
//   if ((architecture->GetNodeNbr() + 1) == architecture->GetNbrNodes())
//     {
//       TmpNbrComponent += this->NbrComponent % Step;
//     }
//   this->SetIndicesRange(TmpFirstComponent, TmpNbrComponent); 
//   switch (architecture->GetArchitectureID())
//     {	 
//     case AbstractArchitecture::MixedMPISMP:
//       this->ArchitectureDependentApplyOperation((SMPArchitecture*) (architecture->GetLocalArchitecture())); 
//       break;
//     default:
//       this->RawApplyOperation();
//       break;
//     }		
//   MPI::COMM_WORLD.Barrier();
//   if (this->OutputState != 0)
//     architecture->SumVector(*(this->OutputState));	
//   else
//     architecture->SumVector(*(this->ComplexOutputState));	
  return true;
#else
  return this->RawApplyOperation();
#endif
}
