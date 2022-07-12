////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//               class of mulitple Jack generator basis to create an	      //
//                     orthonomal basis on the cylinder geometry              //
//                                                                            //
//                        last modification : 21/07/2016                      //
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
#include "Architecture/ArchitectureOperation/FQHECylinderMultipleJackGeneratorOperation.h"
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
// targetSpace = pointer to the Hilbert where the quasihole states should be expressed
// rootConfigurations = array that contains all the root configurations
// nbrRootConfigurations = number of root configurations
// kValue = k value of the clustered (k,r) principle
// rValue = r value of the clustered (k,r) principle
// nbrParticles = number of particles
// lzMax = number of flux quantum
// totalLz = system total momentum
// ratio = cylinder aspect ratio
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode


FQHECylinderMultipleJackGeneratorOperation::FQHECylinderMultipleJackGeneratorOperation(ParticleOnSphere* targetSpace, unsigned long** rootConfigurations, int nbrRootConfigurations, 
										       int kValue, int rValue, int nbrParticles, int lzMax, int totalLz, double ratio, int nbrMPIStage, int nbrSMPStage)
{
  this->FirstComponent = 0;
  this->NbrComponent = nbrRootConfigurations;
  this->NbrRootConfigurations = nbrRootConfigurations;
  this->TargetSpace = (ParticleOnSphere*) targetSpace->Clone();
  this->NbrParticles = nbrParticles;
  this->NbrFluxQuanta = lzMax;
  this->TotalLz = totalLz;
  this->KValue = kValue;
  this->RValue = rValue;
  this->CylinderRatio = ratio;
  this->RootConfigurations = new unsigned long*[this->NbrRootConfigurations];
  for (int i = 0; i < this->NbrRootConfigurations; ++i)
    {
      this->RootConfigurations[i] = new unsigned long[this->NbrFluxQuanta + 1];
      for (int j = 0; j <= this->NbrFluxQuanta; ++j)
	{
	  this->RootConfigurations[i][j] = rootConfigurations[i][j];
	}
    }
  this->QuasiholeVectors = new RealVector[this->NbrRootConfigurations];
  this->OperationType = AbstractArchitectureOperation::FQHECylinderMultipleJackGeneratorOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}

// copy constructor 
//
// operation = reference on operation to copy

FQHECylinderMultipleJackGeneratorOperation::FQHECylinderMultipleJackGeneratorOperation(const FQHECylinderMultipleJackGeneratorOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->OperationType = AbstractArchitectureOperation::FQHECylinderMultipleJackGeneratorOperation;
  this->NbrRootConfigurations = operation.NbrRootConfigurations;
  this->TargetSpace = (ParticleOnSphere*) operation.TargetSpace->Clone();
  this->NbrParticles = operation.NbrParticles;
  this->NbrFluxQuanta = operation.NbrFluxQuanta;
  this->TotalLz = operation.TotalLz;
  this->KValue = operation.KValue;
  this->RValue = operation.RValue;
  this->CylinderRatio = operation.CylinderRatio;
  this->RootConfigurations = new unsigned long*[this->NbrRootConfigurations];
  for (int i = 0; i < this->NbrRootConfigurations; ++i)
    {
      this->RootConfigurations[i] = new unsigned long[this->NbrFluxQuanta + 1];
      for (int j = 0; j <= this->NbrFluxQuanta; ++j)
	{
	  this->RootConfigurations[i][j] = operation.RootConfigurations[i][j];
	}
    }
  this->QuasiholeVectors = new RealVector[this->NbrRootConfigurations];
  this->NbrMPIStage = operation.NbrMPIStage;
  this->NbrSMPStage = operation.NbrSMPStage;
  this->SMPStages = operation.SMPStages;
}

// destructor
//

FQHECylinderMultipleJackGeneratorOperation::~FQHECylinderMultipleJackGeneratorOperation()
{						
  for (int i = 0; i < this->NbrRootConfigurations; ++i)
    {
      delete[] this->RootConfigurations[i];
    }
  delete[] this->RootConfigurations;
  delete[] this->QuasiholeVectors;
  delete this->TargetSpace;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHECylinderMultipleJackGeneratorOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHECylinderMultipleJackGeneratorOperation::Clone()
{
  return new FQHECylinderMultipleJackGeneratorOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHECylinderMultipleJackGeneratorOperation::RawApplyOperation()
{
  if (this->NbrComponent == 0)
    return true;
  int LastComponent = this->NbrComponent + this->FirstComponent;

   int* ReferenceState = new int[this->NbrFluxQuanta + 1];
 if (this->TargetSpace->GetParticleStatistic() == AbstractQHEParticle::BosonicStatistic)
    {
      for (int i = this->FirstComponent; i < LastComponent; ++i)
	{	  
	  for (int j = 0; j <= this->NbrFluxQuanta; ++j)
	    {
	      ReferenceState[j] = (int) this->RootConfigurations[i][j];
	    }
	  double Alpha = ((double) -this->KValue) / ((double) (this->RValue + this->KValue));
	  if (this->NbrParticles > 1)
	    {
	      BosonOnSphereHaldaneBasisShort SqueezedSpace (this->NbrParticles, this->TotalLz, this->NbrFluxQuanta, ReferenceState);
	      RealVector TmpState(SqueezedSpace.GetLargeHilbertSpaceDimension(), true);
	      SqueezedSpace.GenerateJackPolynomial(TmpState, Alpha);
	      SqueezedSpace.NormalizeJackToCylinder(TmpState, this->CylinderRatio);
	      this->QuasiholeVectors[i] = SqueezedSpace.ConvertToNbodyBasis(TmpState, *((BosonOnSphereShort*) this->TargetSpace));
	    }
	  else
	    {
	      this->QuasiholeVectors[i] = RealVector(1);
	      this->QuasiholeVectors[i][0] = 1.0;
	    }
	}
    }
  else
    {
      double Alpha = ((double) -(this->KValue + 1)) / ((double) (this->RValue - 1));
      for (int i = this->FirstComponent; i < LastComponent; ++i)
	{	  
	  for (int j = 0; j <= this->NbrFluxQuanta; ++j)
	    {
	      ReferenceState[j] = (int) this->RootConfigurations[i][j];
	    }
	  if (this->NbrParticles > 1)
	    {
	      FermionOnSphereHaldaneBasis SqueezedSpace (this->NbrParticles, this->TotalLz, this->NbrFluxQuanta, ReferenceState);
	      RealVector TmpState(SqueezedSpace.GetLargeHilbertSpaceDimension(), true);
	      SqueezedSpace.GenerateJackPolynomial(TmpState, Alpha);
	      SqueezedSpace.NormalizeJackToCylinder(TmpState, this->CylinderRatio);
	      this->QuasiholeVectors[i] = SqueezedSpace.ConvertToNbodyBasis(TmpState, *((FermionOnSphere*) this->TargetSpace));
	    }
	  else
	    {
	      this->QuasiholeVectors[i] = RealVector(1);
	      this->QuasiholeVectors[i][0] = 1.0;
	    }
	}
    }
  delete[] ReferenceState;
  return true;
}



// apply operation for SMP using round robin scheduling
//
//  architecture = instance of architecture class
// return value = true if no error occurs

bool FQHECylinderMultipleJackGeneratorOperation::ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID)
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

bool FQHECylinderMultipleJackGeneratorOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  if (Step == 0)
    Step = this->NbrComponent;
  this->SMPStages[0] = 0;
  int TotalNbrComponent = this->FirstComponent + this->NbrComponent;
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHECylinderMultipleJackGeneratorOperation** TmpOperations = new FQHECylinderMultipleJackGeneratorOperation * [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHECylinderMultipleJackGeneratorOperation *) this->Clone();
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
      for (int j = 0; j < this->NbrRootConfigurations; ++j)
	{
	  if (TmpOperations[i]->QuasiholeVectors[j].GetVectorDimension() > 0)
	    {
	      this->QuasiholeVectors[j] = TmpOperations[i]->QuasiholeVectors[j];
	    }
	}
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}

// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHECylinderMultipleJackGeneratorOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
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
