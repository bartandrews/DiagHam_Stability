////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of operations that compute the matrix elements            //
//                  of the many-body interaction on the torus                 // 
//                                                                            //
//                        last modification : 06/03/2015                      //
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


#ifndef FQHETORUSCOMPUTEMATRIXELEMENTOPERATION_H
#define FQHETORUSCOMPUTEMATRIXELEMENTOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Hamiltonian/ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian.h"


class FQHETorusComputeMatrixElementOperation : public AbstractArchitectureOperation
{

 protected:
  
  // pointer to the hamiltonian
  ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian* Hamiltonian;
  // pointer to the hamiltonian (twisted torus version)
  ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian* TwistedHamiltonian;

  // number of unique matrix elements
  long NbrUniqueMatrixElements;

  // array that contains the momentum sector of each matrix element
  int* MomentumSectorIndices;
  // array that contains the creation indices of each matrix element
  int* J1Indices;
  // array that contains the annihilation indices of each matrix element
  int* J2Indices;
  // array that contains the matrix elements
  Complex* MatrixElements;

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;
  // index of the first component (if long numbers are required)
  long LargeFirstComponent;
  // number of component  (if long numbers are required)
  long LargeNbrComponent;
  
 // number of part in the CFT calculation will be separated in MPI mode
  int NbrMPIStage;
  // number of part in the CFT calculation will be separated in SMP mode
  int NbrSMPStage;
  // array with size of SMP stages used to distribute work
  int* SMPStages; 

 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the generic n-body Hamiltonian
  // nbrUniqueMatrixElements = number of unique matrix elements
  // momentumSectorIndices = array that contains the momentum sector of each matrix element
  // j1Indices = array that contains the creation indices of each matrix element
  // j2Indices = array that contains the annihilation indices of each matrix element
  // matrixElements = array where the matrix elements will be stored
  // nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
  // nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode
  FQHETorusComputeMatrixElementOperation(ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian* hamiltonian, long nbrUniqueMatrixElements,
					 int* momentumSectorIndices, int* j1Indices, int* j2Indices, Complex* matrixElements,
					 int nbrMPIStage = 2, int nbrSMPStage = 2);
    
  // constructor 
  //
  // hamiltonian = pointer to the generic n-body Hamiltonian
  // nbrUniqueMatrixElements = number of unique matrix elements
  // momentumSectorIndices = array that contains the momentum sector of each matrix element
  // j1Indices = array that contains the creation indices of each matrix element
  // j2Indices = array that contains the annihilation indices of each matrix element
  // matrixElements = array where the matrix elements will be stored
  // nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
  // nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode
  FQHETorusComputeMatrixElementOperation(ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian* hamiltonian, long nbrUniqueMatrixElements,
					 int* momentumSectorIndices, int* j1Indices, int* j2Indices, Complex* matrixElements,
					 int nbrMPIStage = 2, int nbrSMPStage = 2);
    
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHETorusComputeMatrixElementOperation(const FQHETorusComputeMatrixElementOperation& operation);
  
  // destructor
  //
  ~FQHETorusComputeMatrixElementOperation();
  
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
  
  // apply operation for SMP using round robin scheduling
  //
  //  architecture = instance of architecture class
  // return value = true if no error occurs
  virtual bool ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID);

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

#endif
