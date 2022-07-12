////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                 Copyright (C) 2001-2002 Antoine Sterdyniak                 //
//                                                                            //
//                                                                            //
//                       class of computation operation                       //
//                                                                            //
//                        last modification : 09/05/2010                      //
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


#ifndef FQHESPHEREBOSONICSTATETIMESPOLARIZEDSLATERPROJECTIONOPERATION_H
#define FQHESPHEREBOSONICSTATETIMESPOLARIZEDSLATERPROJECTIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Architecture/MixedMPISMPArchitecture.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereShort.h"
//#include "GeneralTools/Permutations.h"



class FQHESphereBosonsWithSpinLandauLevelLiftOperation : public AbstractArchitectureOperation
{
 protected:
  
  // index of the first component
  int FirstComponent;
  
  // number of component 
  int NbrComponent;
  
  // pointer to the HilbertSpace after product and projection if there is one
  BosonOnSphereWithSpin * FinalSpace;
  
  // pointer to the initial Bosonic HilbertSpace
  BosonOnSphereWithSpin * InitialSpace;
  
  // pointer to the initial Fermionic HilbertSpace
  BosonOnSphereShort * PolarizedSpace;
	
  // RealVector where the result will be store
  RealVector * OutputVector;
  
  // RealVector where the spinful bosonic input state is stored
  RealVector * InitialVector;

  // RealVector where the polarized product state is stored
  RealVector * PolarizedVector;
  
  // number of part in which the initial bosonic vector will be separated
  int NbrMPIStage;
  
  // number of part in which the initial bosonic vector will be separated
  int NbrSMPStage;
  	
  
  // number of MPI process. If not using MPI this is 0
  int MPINodeNbr;
  
  // execution time measured in RawApply
  double ExecutionTime;
  
  // array with size of SMP stages used to distribute work
  int *SMPStages; 
    
  // Index to start from
  int ResumeIdx;
  
 public:
  
  // constructor 
  //
  // initialSpace = pointer to the inital spinful HilbertSpace
  // polarizedSpace = pointer to a spinless bosonic HilbertSpace
  // finalSpace = pointer to the final spinful HilbertSpace
  // initialVector = pointer to the inital spinful state vector
  // polarizedVector = pointer to a spinless bosonic state vector
  // outputVector = pointer to the final spinful state vector
  // nbrMPIStage = number of chunks for MPI parallelism
  // nbrSMPStage = number of chunks for SMP parallelism
  // resumeIdx = index where to resume operation
  FQHESphereBosonsWithSpinLandauLevelLiftOperation(BosonOnSphereWithSpin * initialSpace, BosonOnSphereShort * polarizedSpace, BosonOnSphereWithSpin * finalSpace, RealVector* initialVector, RealVector *polarizedVector, RealVector* outputVector, int nbrMPIStage = 20, int nbrSMPStage = 20, int resumeIdx = 0);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereBosonsWithSpinLandauLevelLiftOperation(const FQHESphereBosonsWithSpinLandauLevelLiftOperation & operation);
  
  // destructor
  //
  ~FQHESphereBosonsWithSpinLandauLevelLiftOperation();
  
   // This function calculates the size for a process when using multiple processes. 
  //
  // n = total size
  // rank = rank of process
  // size = number of processes
  // return value = size for process
  int GetRankChunkSize(int n, int rank, int size);
  
  // This function calculates the starting index for a process when using multiple processes. 
  //
  // n = total size
  // rank = rank of process
  // size = number of processes
  // return value = starting index for process
  int GetRankChunkStart(int n, int rank, int size);
    
 protected:
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetOutputVector (RealVector* outputVector);
  
  // get destination vector 
  // 
  // return value = pointer to destination vector
  RealVector* GetDestinationVector ();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
  // apply operation for SMP using round robin scheduling
  //
  //  architecture = instance of architecture class
  // return value = true if no error occurs
  bool ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID);
  
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
  
//   // apply operation for MixedMPISMP architecture
//   //
//   // architecture = pointer to the architecture
//   // return value = true if no error occurs
//   bool ArchitectureDependentApplyOperation(MixedMPISMPArchitecture* architecture);
  
};


// This function calculates the size for a process when using multiple processes. 
//
// n = total size
// rank = rank of process
// size = number of processes
// return value = size for process

inline int FQHESphereBosonsWithSpinLandauLevelLiftOperation::GetRankChunkSize(int n, int rank, int size) 
{
  int  MaxChunk, MinChunk, CutOff; 
  
  MinChunk = (int)floor((double)n / (double)size); 
  MaxChunk = (int)ceil((double)n / (double)size);
  CutOff = n - (size * MinChunk); 
  
  if ( rank < CutOff ) 
    {
      return MaxChunk;
    } 
  else 
    {
      return MinChunk;
    }
}


// This function calculates the starting index for a process when using multiple processes. 
//
// n = total size
// rank = rank of process
// size = number of processes
// return value = starting index for process

inline int FQHESphereBosonsWithSpinLandauLevelLiftOperation::GetRankChunkStart(int n, int rank, int size)
{
  int  MaxChunk, MinChunk, CutOff; 
  
  MinChunk = (int)floor((double)n / (double)size); 
  MaxChunk =(int)ceil((double)n / (double)size);
  CutOff = n - (size * MinChunk); 

  if ( rank < CutOff ) 
    {
      return rank * MaxChunk;
    } 
  else 
    {
      return (CutOff * MaxChunk) + ((rank - CutOff) * MinChunk);
    }
}

#endif	
