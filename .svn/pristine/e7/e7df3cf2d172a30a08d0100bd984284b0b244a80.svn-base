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
#include "Architecture/ArchitectureOperation/FQHEMPSCreateStateOperation.h"
#include "Vector/RealVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Matrix/SparseComplexMatrix.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor 
//
// space = pointer to the Hilbert space
// bMatrices = array that gives the B matrices 
// state = pointer to the vector where the MPS state will be stored
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// blockSize = indicates the size of the block for precalculations

FQHEMPSCreateStateOperation::FQHEMPSCreateStateOperation(ParticleOnSphere* space, SparseRealMatrix* bMatrices, RealVector* state, 
							 int mPSRowIndex, int mPSColumnIndex, int blockSize)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->Space = (ParticleOnSphere*) space->Clone();
  this->OutputState = state;
  this->ComplexOutputState = 0;  
  this->BMatrices = bMatrices;
  this->ComplexBMatrices = 0;
  this->QuasiholeBMatrices = 0;
  this->NbrQuasiholes = 0;
  this->MPSRowIndex = mPSRowIndex;  
  this->MPSColumnIndex = mPSColumnIndex;
  this->TorusSpace = 0;
  this->TopologicalSectorIndices = 0;
  this->TopologicalSectorNbrIndices = 0;
  this->PrecalculationBlockSize = blockSize;
  this->OperationType = AbstractArchitectureOperation::FQHEMPSCreateStateOperation;
}

// constructor for MPS with quasiholes
//
// space = pointer to the Hilbert space
// bMatrices = array that gives the B matrices 
// quasiholeBMatrices = array that gives the B matrices for quasiholes 
// nbrQuasiholes = number of quasiholes 
// state = pointer to the vector where the MPS state will be stored
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// blockSize = indicates the size of the block for precalculations

FQHEMPSCreateStateOperation::FQHEMPSCreateStateOperation(ParticleOnSphere* space, SparseRealMatrix* bMatrices, SparseComplexMatrix* quasiholeBMatrices, int nbrQuasiholes,
							 ComplexVector* state, int mPSRowIndex, int mPSColumnIndex, int blockSize)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->Space = (ParticleOnSphere*) space->Clone();
  this->OutputState = 0;
  this->ComplexOutputState = state;  
  this->BMatrices = bMatrices;
  this->ComplexBMatrices = 0;
  this->QuasiholeBMatrices = quasiholeBMatrices;
  this->NbrQuasiholes = nbrQuasiholes;
  this->MPSRowIndex = mPSRowIndex;  
  this->MPSColumnIndex = mPSColumnIndex;
  this->PrecalculationBlockSize = blockSize;
  this->TorusSpace = 0;
  this->TopologicalSectorIndices = 0;
  this->TopologicalSectorNbrIndices = 0;
  this->OperationType = AbstractArchitectureOperation::FQHEMPSCreateStateOperation;
}

// constructor for the torus geometry
//
// space = pointer to the Hilbert space
// bMatrices = array that gives the B matrices 
// state = pointer to the vector where the MPS state will be stored
// stringMatrix = matrix that takes into account the Jordan Wigner string on the torus geometry
// topologicalSectorIndices = array that contains the auxiliary space indices related to the selected topological sector
// topologicalSectorNbrIndices = number of indices in TopologicalSectorIndices
// blockSize = indicates the size of the block for precalculations

FQHEMPSCreateStateOperation::FQHEMPSCreateStateOperation(ParticleOnTorus* space, SparseRealMatrix* bMatrices, SparseRealMatrix& stringMatrix, 
							 RealVector* state, int* topologicalSectorIndices, int topologicalSectorNbrIndices, int blockSize)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->TorusSpace = (ParticleOnTorus*) space->Clone();
  this->Space = 0;
  this->OutputState = state;
  this->ComplexOutputState = 0;  
  this->BMatrices = bMatrices;
  this->ComplexBMatrices = 0;
  this->NbrQuasiholes = 0;
  this->MPSRowIndex = 0;  
  this->MPSColumnIndex = 0;
  this->PrecalculationBlockSize = blockSize;
  this->TorusStringMatrix = stringMatrix;
  this->TopologicalSectorIndices = topologicalSectorIndices;
  this->TopologicalSectorNbrIndices = topologicalSectorNbrIndices;
  this->OperationType = AbstractArchitectureOperation::FQHEMPSCreateStateOperation;
}

// constructor for the torus geometry
//
// space = pointer to the Hilbert space
// bMatrices = array that gives the B matrices 
// state = pointer to the vector where the MPS state will be stored
// stringMatrix = matrix that takes into account the Jordan Wigner string on the torus geometry
// topologicalSectorIndices = array that contains the auxiliary space indices related to the selected topological sector
// topologicalSectorNbrIndices = number of indices in TopologicalSectorIndices
// blockSize = indicates the size of the block for precalculations

FQHEMPSCreateStateOperation::FQHEMPSCreateStateOperation(ParticleOnTorus* space, SparseComplexMatrix* bMatrices, SparseRealMatrix& stringMatrix, 
							ComplexVector* state, int* topologicalSectorIndices, int topologicalSectorNbrIndices, int blockSize)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->TorusSpace = (ParticleOnTorus*) space->Clone();
  this->Space = 0;
  this->ComplexOutputState = state;
  this->OutputState = 0;  
  this->ComplexBMatrices = bMatrices;
  this->BMatrices = 0;
  this->NbrQuasiholes = 0;
  this->MPSRowIndex = 0;  
  this->MPSColumnIndex = 0;
  this->PrecalculationBlockSize = blockSize;
  this->TorusStringMatrix = stringMatrix;
  this->TopologicalSectorIndices = topologicalSectorIndices;
  this->TopologicalSectorNbrIndices = topologicalSectorNbrIndices;
  this->OperationType = AbstractArchitectureOperation::FQHEMPSCreateStateOperation;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHEMPSCreateStateOperation::FQHEMPSCreateStateOperation(const FQHEMPSCreateStateOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;

  if (operation.Space != 0)
    {
      this->Space = (ParticleOnSphere*) operation.Space->Clone();
      this->TorusSpace = 0;
    }
  else
    {
      this->TorusSpace = (ParticleOnTorus*) operation.TorusSpace->Clone();
      this->Space = 0;
    }
  this->OutputState = operation.OutputState;
  this->ComplexOutputState = operation.ComplexOutputState;  
  this->BMatrices = operation.BMatrices;
  this->ComplexBMatrices = operation.ComplexBMatrices;
  this->QuasiholeBMatrices = operation.QuasiholeBMatrices;
  this->NbrQuasiholes = operation.NbrQuasiholes;
  this->MPSRowIndex = operation.MPSRowIndex;  
  this->MPSColumnIndex = operation.MPSColumnIndex;
  this->PrecalculationBlockSize = operation.PrecalculationBlockSize;
  this->TopologicalSectorIndices = operation.TopologicalSectorIndices;
  this->TopologicalSectorNbrIndices = operation.TopologicalSectorNbrIndices;
  this->TorusStringMatrix = operation.TorusStringMatrix;
  this->OperationType = AbstractArchitectureOperation::FQHEMPSCreateStateOperation;	
}

// destructor
//

FQHEMPSCreateStateOperation::~FQHEMPSCreateStateOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHEMPSCreateStateOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set the output state 
// 
// state = pointer to the output state

void FQHEMPSCreateStateOperation::SetOutputState (RealVector* state)
{
  this->OutputState = state;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHEMPSCreateStateOperation::Clone()
{
  return new FQHEMPSCreateStateOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHEMPSCreateStateOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  if (this->Space != 0)
    {
      if (this->NbrQuasiholes == 0)
	{
	  this->Space->CreateStateFromMPSDescription(this->BMatrices, *(this->OutputState), this->MPSRowIndex, this->MPSColumnIndex, (long) this->PrecalculationBlockSize, this->FirstComponent, this->NbrComponent);
	}
      else
	{
	  // FIXME: the third argument "1" is a hack for the case of a single matrix at edge
	  this->Space->CreateStateFromMPSDescription(this->BMatrices, this->QuasiholeBMatrices, 1, *(this->ComplexOutputState), 
						     this->MPSRowIndex, this->MPSColumnIndex, (long) this->PrecalculationBlockSize, this->FirstComponent, this->NbrComponent);
	}
    }
  else
    {
      if (this->BMatrices != 0)
	{
	  this->TorusSpace->CreateStateFromMPSDescription(this->BMatrices, this->TorusStringMatrix, *(this->OutputState), 
							  this->TopologicalSectorIndices, this->TopologicalSectorNbrIndices, 
							  (long) this->PrecalculationBlockSize, this->FirstComponent, this->NbrComponent);
	}
      else
	{
	  this->TorusSpace->CreateStateFromMPSDescription(this->ComplexBMatrices, this->TorusStringMatrix, *(this->ComplexOutputState), 
							  this->TopologicalSectorIndices, this->TopologicalSectorNbrIndices, 
							  (long) this->PrecalculationBlockSize, this->FirstComponent, this->NbrComponent);
	}
    }
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;
  return true;
}



// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHEMPSCreateStateOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHEMPSCreateStateOperation** TmpOperations = new FQHEMPSCreateStateOperation * [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHEMPSCreateStateOperation *) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  for( int i = 0; i <  ReducedNbrThreads ; i++)
    {
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpFirstComponent += Step;
    }

  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent + this->NbrComponent - TmpFirstComponent);
 
  architecture->SendJobs();

  for (int i = 1; i < architecture->GetNbrThreads(); i++)
    {
      delete TmpOperations[i];
    }

  delete TmpOperations[0];
  delete[] TmpOperations;
  return true;
}

// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHEMPSCreateStateOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__    
  int Step = this->NbrComponent / architecture->GetNbrNodes();
  int TmpFirstComponent = this->FirstComponent + (Step * architecture->GetNodeNbr());
  int TmpNbrComponent = Step;
  if ((architecture->GetNodeNbr() + 1) == architecture->GetNbrNodes())
    {
      TmpNbrComponent += this->NbrComponent % Step;
    }
  this->SetIndicesRange(TmpFirstComponent, TmpNbrComponent); 
  switch (architecture->GetArchitectureID())
    {	 
    case AbstractArchitecture::MixedMPISMP:
      this->ArchitectureDependentApplyOperation((SMPArchitecture*) (architecture->GetLocalArchitecture())); 
      break;
    default:
      this->RawApplyOperation();
      break;
    }		
  MPI::COMM_WORLD.Barrier();
  if (this->OutputState != 0)
    architecture->SumVector(*(this->OutputState));	
  else
    architecture->SumVector(*(this->ComplexOutputState));	
  return true;
#else
  return this->RawApplyOperation();
#endif
}
