////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//               class of operation to generate the product of a              //
//                       bosonic state and fermionic state	              //
//                                                                            //
//                        last modification : 28/04/2017                      //
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


#ifndef FQHESPHEREBOSONICSTATETIMESFERMIONICSTATEOPERATION_H
#define FQHESPHEREBOSONICSTATETIMESFERMIONICSTATEOPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphere.h"



class FQHESphereBosonicStateTimesFermionicStateOperation : public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  long FirstComponent;
  // number of component to compute 
  long NbrComponent;
  
  // bosonic state
  RealVector BosonicState;
  // fermionic state
  RealVector FermionicState;
  // vector where the resulting product will be stored
  RealVector OutputState;
			     
  // pointer to the Hilbert Space associated to the bosonic state
  BosonOnSphereShort* BosonicSpace;
  // pointer to the Hilbert Space associated to the fermionic state
  FermionOnSphere* FermionicSpace;
  // pointer to the Hilbert Space associated to the resulting state
  FermionOnSphere* OutputSpace; 
  
  // pointer to the Hilbert Space associated to the second bosonic state
  BosonOnSphereShort* BosonicSpace2;
  // pointer to the Hilbert Space associated to the resulting state if bosonic
  BosonOnSphereShort* BosonicOutputSpace; 

  // true if the state should be written in the unnormalized basis
  bool UnnormalizedFlag;

  // number of part in the CFT calculation will be separated in MPI mode
  int NbrMPIStage;
  // number of part in the CFT calculation will be separated in SMP mode
  int NbrSMPStage;
  // array with size of SMP stages used to distribute work
  int* SMPStages; 

  // if non 0, use it as a prefix to store the partial results
  char* StorePartialProductPrefix;

 public:
  
  // constructor 
  //
  // bosonicState = reference on the bosonic state
  // fermionicState = reference on the fermionic state
  // bosonicSpace = pointer to the Hilbert Space associated to the bosonic state
  // fermionicSpace = pointer to the Hilbert Space associated to the fermionic state
  // outputSpace = pointer to the Hilbert Space associated to the resulting state
  // unnormalizedFlag = true if the state should be written in the unnormalized basis
  // storePartialProducts = if non 0, use it as a prefix to store the partial results
  // nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
  // nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode
  FQHESphereBosonicStateTimesFermionicStateOperation(RealVector& bosonicState, RealVector& fermionicState,
						     BosonOnSphereShort* bosonicSpace, FermionOnSphere* fermionicSpace, FermionOnSphere* outputSpace, 
						     bool unnormalizedFlag, char* storePartialProducts = 0, int nbrMPIStage = 10, int nbrSMPStage = 10);
  
  // constructor 
  //
  // bosonicState = reference on the bosonic state
  // bosonicState2 = reference on the second bosonic state
  // bosonicSpace = pointer to the Hilbert Space associated to the bosonic state
  // bosonicSpace = pointer to the Hilbert Space associated to the second bosonic state
  // outputSpace = pointer to the Hilbert Space associated to the resulting state
  // unnormalizedFlag = true if the state should be written in the unnormalized basis
  // storePartialProducts = if non 0, use it as a prefix to store the partial results
  // nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
  // nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode
  FQHESphereBosonicStateTimesFermionicStateOperation(RealVector& bosonicState, RealVector& bosonicState2,
						     BosonOnSphereShort* bosonicSpace, BosonOnSphereShort* bosonicSpace2, BosonOnSphereShort* outputSpace, 
						     bool unnormalizedFlag, char* storePartialProducts = 0, int nbrMPIStage = 10, int nbrSMPStage = 10);
  

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereBosonicStateTimesFermionicStateOperation(const FQHESphereBosonicStateTimesFermionicStateOperation & operation);
  
  // constructor from a master node information
  //
  // space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  //  FQHESphereBosonicStateTimesFermionicStateOperation(ParticleOnSphere* space,  SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~FQHESphereBosonicStateTimesFermionicStateOperation();  
  
  // get the state obtained from the product of the bosonic state and the fermionic state
  // 
  // return value = output vector
  virtual RealVector GetState ();
  
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

// get the state obtained from the product of the bosonic state and the fermionic state
// 
// return value = output vector

inline RealVector FQHESphereBosonicStateTimesFermionicStateOperation::GetState ()
{
  return this->OutputState;
}

#endif
