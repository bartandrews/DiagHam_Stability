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


#ifndef FHQEMPSCREATESTATEOPERATION_H
#define FHQEMPSCREATESTATEOPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "Matrix/SparseRealMatrix.h"



class RealVector;
class SparseComplexMatrix;


class FQHEMPSCreateStateOperation: public AbstractArchitectureOperation
{

 protected:
  
  // pointer to the Hilbert space
  ParticleOnSphere* Space;
  // pointer to the Hilbert space for the torus geometry
  ParticleOnTorus* TorusSpace;

  // vector where the MPS state will be stored
  RealVector* OutputState;
  // vector where the MPS state will be stored (complex version)
  ComplexVector* ComplexOutputState;

  // array that gives the B matrices 
  SparseRealMatrix* BMatrices;
  // array that gives the B matrices 
  SparseComplexMatrix* ComplexBMatrices;

  // array that gives the B matrices for quasiholes 
  SparseComplexMatrix* QuasiholeBMatrices;
  // number of quasiholes 
  int NbrQuasiholes;

  // row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  int MPSRowIndex;
  // column index of the MPS element that has to be evaluated
  int MPSColumnIndex;

  // matrix that takes into account the Jordan Wigner string on the torus geometry
  SparseRealMatrix TorusStringMatrix;
  // array that contains the auxiliary space indices related to the selected topological sector
  int* TopologicalSectorIndices;
  // number of indices in TopologicalSectorIndices
  int TopologicalSectorNbrIndices;

  // indicates the size of the block for precalculations
  int PrecalculationBlockSize;

  // index of the first component
  long FirstComponent;
  // number of component to compute 
  long NbrComponent;
  
 public:
  
  // constructor 
  //
  // space = pointer to the Hilbert space
  // bMatrices = array that gives the B matrices 
  // state = pointer to the vector where the MPS state will be stored
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // blockSize = indicates the size of the block for precalculations
  FQHEMPSCreateStateOperation(ParticleOnSphere* space, SparseRealMatrix* bMatrices, RealVector* state, 
			      int mPSRowIndex, int mPSColumnIndex, int blockSize);
  
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
  FQHEMPSCreateStateOperation(ParticleOnSphere* space, SparseRealMatrix* bMatrices, SparseComplexMatrix* quasiholeBMatrices, int nbrQuasiholes,
			      ComplexVector* state, int mPSRowIndex, int mPSColumnIndex, int blockSize);
  
  // constructor for the torus geometry
  //
  // space = pointer to the Hilbert space
  // bMatrices = array that gives the B matrices 
  // state = pointer to the vector where the MPS state will be stored
  // stringMatrix = matrix that takes into account the Jordan Wigner string on the torus geometry
  // topologicalSectorIndices = array that contains the auxiliary space indices related to the selected topological sector
  // topologicalSectorNbrIndices = number of indices in TopologicalSectorIndices
  // blockSize = indicates the size of the block for precalculations
  FQHEMPSCreateStateOperation(ParticleOnTorus* space, SparseRealMatrix* bMatrices, SparseRealMatrix& stringMatrix, 
			      RealVector* state, int* topologicalSectorIndices, int topologicalSectorNbrIndices, int blockSize);

  // constructor for the torus geometry
  //
  // space = pointer to the Hilbert space
  // bMatrices = array that gives the B matrices 
  // state = pointer to the vector where the MPS state will be stored
  // stringMatrix = matrix that takes into account the Jordan Wigner string on the torus geometry
  // topologicalSectorIndices = array that contains the auxiliary space indices related to the selected topological sector
  // topologicalSectorNbrIndices = number of indices in TopologicalSectorIndices
  // blockSize = indicates the size of the block for precalculations
  FQHEMPSCreateStateOperation(ParticleOnTorus* space, SparseComplexMatrix* bMatrices, SparseRealMatrix& stringMatrix, 
			      ComplexVector* state, int* topologicalSectorIndices, int topologicalSectorNbrIndices, int blockSize);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHEMPSCreateStateOperation(const FQHEMPSCreateStateOperation & operation);
  
  // constructor from a master node information
  //
  // space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  FQHEMPSCreateStateOperation(ParticleOnSphere* space,  SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~FQHEMPSCreateStateOperation();  
  
  // set the output state 
  // 
  // state = pointer to the output state
  virtual void SetOutputState (RealVector* state);
  
  // get the output state 
  // 
  // return value = pointer to the output state 
  Vector* GetOutputState ();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  virtual void SetIndicesRange (const long& firstComponent, const long& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  virtual AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  virtual bool RawApplyOperation();
  
 protected:
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
};

// get the output state 
// 
// return value = pointer to the output state 

inline Vector* FQHEMPSCreateStateOperation::GetOutputState ()
{
  return this->OutputState;
}

#endif
