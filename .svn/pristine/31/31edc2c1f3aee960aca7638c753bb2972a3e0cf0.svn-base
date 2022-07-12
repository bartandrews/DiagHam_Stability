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


#include "config.h"
#include "Architecture/ArchitectureOperation/FQHESphereBosonicStateTimesFermionicStateOperation.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"


#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor to compute the CFT matrix elements
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


FQHESphereBosonicStateTimesFermionicStateOperation::FQHESphereBosonicStateTimesFermionicStateOperation(RealVector& bosonicState, RealVector& fermionicState,
												       BosonOnSphereShort* bosonicSpace, FermionOnSphere* fermionicSpace, FermionOnSphere* outputSpace, 
												       bool unnormalizedFlag, char* storePartialProducts, int nbrMPIStage, int nbrSMPStage)
{
  this->BosonicSpace = (BosonOnSphereShort*) bosonicSpace->Clone();
  this->FermionicSpace = (FermionOnSphere*) fermionicSpace->Clone();
  this->OutputSpace = (FermionOnSphere*) outputSpace->Clone(); 
  this->BosonicSpace2 = 0;
  this->BosonicOutputSpace = 0;

  this->BosonicState = bosonicState;
  this->FermionicState = fermionicState;
  this->OutputState = RealVector(this->OutputSpace->GetHilbertSpaceDimension(), true);

  this->UnnormalizedFlag = unnormalizedFlag;

  this->FirstComponent = 0;
  this->NbrComponent = this->BosonicSpace->GetHilbertSpaceDimension();

  this->OperationType = AbstractArchitectureOperation::FQHESphereBosonicStateTimesFermionicState;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 

  if (storePartialProducts != 0)
    {
      this->StorePartialProductPrefix = new char[strlen(storePartialProducts) + 1];
      strcpy (this->StorePartialProductPrefix, storePartialProducts);
    }
  else
    {
      this->StorePartialProductPrefix = 0;
    }
}

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

FQHESphereBosonicStateTimesFermionicStateOperation::FQHESphereBosonicStateTimesFermionicStateOperation(RealVector& bosonicState, RealVector& bosonicState2,
												       BosonOnSphereShort* bosonicSpace, BosonOnSphereShort* bosonicSpace2, BosonOnSphereShort* outputSpace, 
												       bool unnormalizedFlag, char* storePartialProducts, int nbrMPIStage, int nbrSMPStage)
{
  this->BosonicSpace = (BosonOnSphereShort*) bosonicSpace->Clone();
  this->FermionicSpace = 0;
  this->OutputSpace = 0; 
  this->BosonicSpace2 = (BosonOnSphereShort*) bosonicSpace2->Clone();
  this->BosonicOutputSpace = (BosonOnSphereShort*) outputSpace->Clone();

  this->BosonicState = bosonicState;
  this->FermionicState = bosonicState2;
  this->OutputState = RealVector(this->BosonicOutputSpace->GetHilbertSpaceDimension(), true);

  this->UnnormalizedFlag = unnormalizedFlag;

  this->FirstComponent = 0;
  this->NbrComponent = this->BosonicSpace->GetHilbertSpaceDimension();

  this->OperationType = AbstractArchitectureOperation::FQHESphereBosonicStateTimesFermionicState;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 

  if (storePartialProducts != 0)
    {
      this->StorePartialProductPrefix = new char[strlen(storePartialProducts) + 1];
      strcpy (this->StorePartialProductPrefix, storePartialProducts);
    }
  else
    {
      this->StorePartialProductPrefix = 0;
    }
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESphereBosonicStateTimesFermionicStateOperation::FQHESphereBosonicStateTimesFermionicStateOperation(const FQHESphereBosonicStateTimesFermionicStateOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->OperationType = AbstractArchitectureOperation::FQHESphereBosonicStateTimesFermionicState;

  this->BosonicSpace = (BosonOnSphereShort*) operation.BosonicSpace->Clone();
  if (operation.FermionicSpace != 0)
    {
      this->FermionicSpace = (FermionOnSphere*) operation.FermionicSpace->Clone();
      this->OutputSpace = (FermionOnSphere*) operation.OutputSpace->Clone(); 
      this->BosonicSpace2 = 0;
      this->BosonicOutputSpace = 0;
      this->OutputState = RealVector(this->OutputSpace->GetHilbertSpaceDimension(), true);
    }
  else
    {
      this->FermionicSpace = 0;
      this->OutputSpace = 0;
      this->BosonicSpace2 = (BosonOnSphereShort*) operation.BosonicSpace2->Clone();
      this->BosonicOutputSpace = (BosonOnSphereShort*) operation.BosonicOutputSpace->Clone(); 
      this->OutputState = RealVector(this->BosonicOutputSpace->GetHilbertSpaceDimension(), true);
    }
  
  this->BosonicState = operation.BosonicState;
  this->FermionicState = operation.FermionicState;

  this->UnnormalizedFlag = operation.UnnormalizedFlag;

  this->NbrMPIStage = operation.NbrMPIStage;
  this->NbrSMPStage = operation.NbrSMPStage;
  this->SMPStages = operation.SMPStages;

  if (operation.StorePartialProductPrefix != 0)
    {
      this->StorePartialProductPrefix = new char[strlen(operation.StorePartialProductPrefix) + 1];
      strcpy (this->StorePartialProductPrefix, operation.StorePartialProductPrefix);
    }
  else
    {
      this->StorePartialProductPrefix = 0;
    }
}

// destructor
//

FQHESphereBosonicStateTimesFermionicStateOperation::~FQHESphereBosonicStateTimesFermionicStateOperation()
{						
  delete this->BosonicSpace;
  if (this->FermionicSpace != 0)
    {
      delete this->FermionicSpace;
      delete this->OutputSpace; 
    }
  else
    {
      delete this->BosonicSpace2;
      delete this->BosonicOutputSpace;
    }
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHESphereBosonicStateTimesFermionicStateOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereBosonicStateTimesFermionicStateOperation::Clone()
{
  return new FQHESphereBosonicStateTimesFermionicStateOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesFermionicStateOperation::RawApplyOperation()
{
  if (this->NbrComponent == 0)
    return true;

  cout << "processing " << this->FirstComponent << " " << this->NbrComponent << endl;
  if (this->StorePartialProductPrefix == 0)
    {
      if (this->FermionicSpace != 0)
	{
	  this->OutputSpace->BosonicStateTimeFermionicState(this->BosonicState, this->FermionicState, this->OutputState, 
							    this->BosonicSpace, this->FermionicSpace, this->FirstComponent, 
							    this->NbrComponent, this->UnnormalizedFlag, 0);
	}
      else
	{
	  this->BosonicOutputSpace->BosonicStateTimeBosonicState(this->BosonicState, this->FermionicState, this->OutputState, 
								 this->BosonicSpace, this->BosonicSpace2, this->FirstComponent, 
								 this->NbrComponent, this->UnnormalizedFlag, 0);
	}
    }
  else
    {
      RealVector TmpVector (this->BosonicState.GetVectorDimension());      
      char* TmpFileName =  new char[strlen(this->StorePartialProductPrefix) + 100];
      for (int i = 0; i < this->NbrComponent; ++i)
	{	  
	  TmpVector[this->FirstComponent + i] = 1.0;
	  this->OutputState.ClearVector();
	  if (this->FermionicSpace != 0)
	    {
	      this->OutputSpace->BosonicStateTimeFermionicState(TmpVector, this->FermionicState, this->OutputState, 
								this->BosonicSpace, this->FermionicSpace, this->FirstComponent + i, 
								1, this->UnnormalizedFlag, 0);
	    }
	  else
	    {
	      this->BosonicOutputSpace->BosonicStateTimeBosonicState(TmpVector, this->FermionicState, this->OutputState, 
								     this->BosonicSpace, this->BosonicSpace2, this->FirstComponent + i, 
								     1, this->UnnormalizedFlag, 0);
	    }
	  sprintf (TmpFileName, "%s.%d.vec", this->StorePartialProductPrefix, (int) (this->FirstComponent + i));
	  this->OutputState.WriteVector(TmpFileName);
	}
      delete[] TmpFileName;
    }
  return true;
}



// apply operation for SMP using round robin scheduling
//
//  architecture = instance of architecture class
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesFermionicStateOperation::ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID)
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

bool FQHESphereBosonicStateTimesFermionicStateOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  if (Step == 0)
    Step = this->NbrComponent;
  this->SMPStages[0] = 0;
  int TotalNbrComponent = this->FirstComponent + this->NbrComponent;
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHESphereBosonicStateTimesFermionicStateOperation** TmpOperations = new FQHESphereBosonicStateTimesFermionicStateOperation * [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHESphereBosonicStateTimesFermionicStateOperation *) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  architecture->SendJobsRoundRobin();

  for (int i = 0; i < architecture->GetNbrThreads(); i++)
    {
      this->OutputState += TmpOperations[i]->OutputState;
    }
  delete[] TmpOperations;
  return true;
}

// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesFermionicStateOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
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
