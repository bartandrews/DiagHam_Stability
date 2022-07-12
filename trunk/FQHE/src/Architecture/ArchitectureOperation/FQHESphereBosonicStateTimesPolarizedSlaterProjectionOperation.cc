////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of operator Kostka numbers computation operation         //
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


#include "config.h"
#include "Architecture/ArchitectureOperation/FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation.h"

#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereThreeLandauLevels.h"
#include "HilbertSpace/FermionOnSphereFourLandauLevels.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;

FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation(ParticleOnSphere * initialSpace, FermionOnSphere * fermionSpace, FermionOnSphereWithSpin * finalSpace, RealVector* bosonicVector, RealVector* outputVector,bool twoLandauLevels, bool twoLandauLevelLz, bool twoLandauLevelSz, int nbrMPIStage, int nbrSMPStage, int resumeIdx)
{
  this->Cloned = false;
  this->FirstComponent = resumeIdx;	
  this->NbrComponent = initialSpace->GetHilbertSpaceDimension();
  this->FermionSpace = fermionSpace;
  this->FermionSpaceDown = 0;
  this->InitialSpace = initialSpace;
  this->FinalSpace = finalSpace;
  this->OperationType = AbstractArchitectureOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation;
  this->OutputVector =  outputVector;
  this->BosonicVector = bosonicVector;
  this->IndexUp = 0;
  this->IndexDown = 0;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->TwoLandauLevels = twoLandauLevels;
  this->TwoLandauLevelLz = twoLandauLevelLz;
  this->TwoLandauLevelSz = twoLandauLevelSz;
  this->MPINodeNbr = 0;       
  this->ResumeIdx = resumeIdx;
  unsigned long* Slater = new unsigned long[this->FermionSpace->NbrFermions];       
  fermionSpace->ConvertToMonomial(this->FermionSpace->StateDescription[0], Slater);
  EvaluateMonomialPermutations(fermionSpace->NbrFermions, Slater, this->NbrSlaterPermutations, this->SlaterPermutations, this->SlaterSigns );
  delete[] Slater;
  this->SMPStages = new int[1]; 
}


FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation(ParticleOnSphere * initialSpace, FermionOnSphere * fermionSpace,FermionOnSphere * fermionSpaceDown, FermionOnSphereWithSpin * finalSpace, RealVector* bosonicVector, RealVector* outputVector, int indexUp, int indexDown, bool twoLandauLevels, bool twoLandauLevelLz, bool twoLandauLevelSz, int nbrMPIStage, int nbrSMPStage, int resumeIdx)
{
  this->Cloned = false;
  this->FirstComponent = resumeIdx;	
  this->NbrComponent = initialSpace->GetHilbertSpaceDimension();
  this->FermionSpace = fermionSpace;
  this->FermionSpaceDown = fermionSpaceDown;
  this->IndexUp = indexUp;
  this->IndexDown = indexDown;
  this->InitialSpace = initialSpace;
  this->FinalSpace = finalSpace;
  this->OperationType = AbstractArchitectureOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation;
  this->OutputVector =  outputVector;
  this->BosonicVector = bosonicVector;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->TwoLandauLevels = twoLandauLevels;
  this->TwoLandauLevelLz = twoLandauLevelLz;
  this->TwoLandauLevelSz = twoLandauLevelSz;
  this->MPINodeNbr = 0;       
  this->ResumeIdx = resumeIdx;
  this->SlaterPermutations = 0;
  this->SlaterSigns = 0;
  this->SMPStages = new int[1]; 
}

// copy constructor 
//
// operation = reference on operation to copy
FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation(const FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation & operation)	
{
  this->Cloned = true;
  this->FirstComponent = operation.FirstComponent;	
  this->NbrComponent 	= operation.NbrComponent;
  this->FermionSpace = operation.FermionSpace;		
  this->FermionSpaceDown = operation.FermionSpaceDown;
  this->InitialSpace = (BosonOnSphereTwoLandauLevels *)operation.InitialSpace->Clone();	
  this->FinalSpace = (FermionOnSphereWithSpin *) operation.FinalSpace->Clone();
  this->OperationType = AbstractArchitectureOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation;
  this->OutputVector =  operation.OutputVector;
  this->BosonicVector = operation.BosonicVector;
  this->NbrMPIStage = operation.NbrMPIStage;
  this->NbrSMPStage = operation.NbrSMPStage;
  this->IndexDown = operation.IndexDown;
  this->IndexUp = operation.IndexUp;
  this->TwoLandauLevels = operation.TwoLandauLevels;
  this->TwoLandauLevelLz= operation.TwoLandauLevelLz;
  this->TwoLandauLevelSz= operation.TwoLandauLevelSz;
  this->MPINodeNbr = operation.MPINodeNbr;
  this->SMPStages = operation.SMPStages;
  this->NbrSlaterPermutations = operation.NbrSlaterPermutations;
  this->SlaterPermutations = operation.SlaterPermutations;
  this->SlaterSigns = operation.SlaterSigns;
  this->ResumeIdx = operation.ResumeIdx;
}

// destructor
//

FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::~FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation()
{
  if ( ! this->Cloned )
    {
      if (this->FermionSpaceDown == 0)
	{
	  for ( int i = 0 ; i < this->NbrSlaterPermutations ; i++ )
	    {
	      delete[] this->SlaterPermutations[i];
	    }
	  delete [] this->SlaterPermutations;
	  delete [] this->SlaterSigns;
	}
      delete [] this->SMPStages;
    }
}

// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::SetOutputVector (RealVector* outputVector)
{
  this->OutputVector = outputVector;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation*FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::Clone()
{
  return new FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation (*this);
}


// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
		if (this->TwoLandauLevels)
		{
			if (this->FermionSpaceDown != 0)
			{
				((BosonOnSphereTwoLandauLevels * )this->InitialSpace)->BosonicStateTimePolarizedSlaters( *this->BosonicVector, *this->OutputVector,this->FermionSpace,this->FermionSpaceDown , this->FinalSpace,this->IndexUp, this->IndexDown, this->FirstComponent ,this->NbrComponent);
			}
			else
			{
      if ( this->TwoLandauLevelLz ) 
	{
	  if ( this->TwoLandauLevelSz ) 
	    {
	      ((BosonOnSphereTwoLandauLevels * )this->InitialSpace)->BosonicStateTimePolarizedSlatersLzSzSymmetry(*this->BosonicVector, *this->OutputVector, this->FermionSpace , this->FinalSpace, this->FirstComponent ,this->NbrComponent, this->SlaterPermutations, this->SlaterSigns, this->NbrSlaterPermutations);
	    }
	  else
	    {
	      ((BosonOnSphereTwoLandauLevels * )this->InitialSpace)->BosonicStateTimePolarizedSlatersLzSymmetry(*this->BosonicVector, *this->OutputVector, this->FermionSpace , this->FinalSpace, this->FirstComponent ,this->NbrComponent, this->SlaterPermutations, this->SlaterSigns, this->NbrSlaterPermutations);
	    }
	}
      else
	{
	  ((BosonOnSphereTwoLandauLevels * )this->InitialSpace)->BosonicStateTimePolarizedSlaters( *this->BosonicVector, *this->OutputVector,this->FermionSpace , this->FinalSpace, this->FirstComponent ,this->NbrComponent, this->SlaterPermutations, this->SlaterSigns, this->NbrSlaterPermutations);
	}
    }
		}
  else
    ((BosonOnSphereShort *)this->InitialSpace)->BosonicStateTimePolarizedSlaters( *this->BosonicVector, *this->OutputVector,this->FermionSpace , this->FinalSpace, this->FirstComponent ,this->NbrComponent);
	
	
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  this->ExecutionTime   = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  //cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;
  return true;
}

// apply operation for SMP using round robin scheduling
//
//  architecture = instance of architecture class
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID)
{      
  int NbrComponents = this->NbrComponent;
  int FirstComponent = this->FirstComponent;
  
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  int NbrStages = this->NbrSMPStage*architecture->GetNbrThreads();
  int StageIdx=0;
  bool locked = false;
  char TmpString[512];
  sprintf(TmpString,"");
  while ( StageIdx < NbrStages ) 
    {
      if ( ! locked )	
        {
          architecture->LockMutex();
          locked = true;
        }
      StageIdx = this->SMPStages[0];
      if ( StageIdx < NbrStages ) 
	{	  
	  this->SMPStages[0]++;	  
	  if (architecture->VerboseMode() == true && strcmp(TmpString,"") != 0 )
	    {
	      architecture->AddToLog(TmpString);	      
	    }
	  architecture->UnLockMutex();
	  locked = false;
	  timeval TotalStartingTime2;
	  gettimeofday (&TotalStartingTime2, 0);
	  this->SetIndicesRange(FirstComponent + this->GetRankChunkStart(NbrComponents, StageIdx,  NbrStages),  this->GetRankChunkSize(NbrComponents, StageIdx,  NbrStages));
	  this->RawApplyOperation();
	  timeval TotalEndingTime2;
	  gettimeofday (&TotalEndingTime2, 0);
	  sprintf (TmpString, "FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation core operation on SMP id %d finished stage %d with size %d in %.4f",  threadID, StageIdx, this->GetRankChunkSize(NbrComponents, StageIdx,  NbrStages),
		   (((double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec)) +(((double) (TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec)) / 1000000.0)) );	      
	  StageIdx++;
	}
    }
  if ( locked )
    {
      architecture->UnLockMutex();
      locked = false;
    }
    
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  this->ExecutionTime   = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));  
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{  
  char SaveFileName[200];
  sprintf(SaveFileName, "fermions_su2_slater_sym_tmp%d.vec",this->MPINodeNbr);        
  int NbrComponent = this->NbrComponent;
	
  /*this->SMPStages = new bool[this->NbrSMPStage*architecture->GetNbrThreads()];
  for ( int i = 0 ; i < this->NbrSMPStage*architecture->GetNbrThreads() ; i++ )
    {
      this->SMPStages[i] = false;
    }*/
  this->SMPStages[0] = 0;
  
  FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation** TmpOperations = new FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation* [architecture->GetNbrThreads()];
  for(int i = 0 ;i < architecture->GetNbrThreads() ;i++)
    {
      TmpOperations[i] = (FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      TmpOperations[i]->SetOutputVector((RealVector*)this->OutputVector->EmptyClone(true));
    }      
    
  architecture->SendJobsRoundRobin();
  if (architecture->VerboseMode() == true)
    {
      char TmpString[512];
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	{
	  sprintf (TmpString, "FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation core operation on SMP id %d done in %.3f seconds",  i, TmpOperations[i]->ExecutionTime);
	  architecture->AddToLog(TmpString);
	}
    }
    
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      (*(this->OutputVector)) += (*(TmpOperations[i]->OutputVector));	
    }
  
//       int StageStart = this->FirstComponent + this->GetRankChunkStart(NbrComponent, p,  this->NbrSMPStage);
//       int StageSize = this->GetRankChunkSize(NbrComponent, p,  this->NbrSMPStage);
//       for (int i = 0; i < architecture->GetNbrThreads() ; ++i)
// 	{	
// 	  TmpOperations[i]->SetIndicesRange(StageStart + this->GetRankChunkStart(StageSize, i, architecture->GetNbrThreads()) ,
// 								    this->GetRankChunkSize(StageSize, i, architecture->GetNbrThreads()) );	    
// 	}	
//       for (int i = 1; i < architecture->GetNbrThreads(); ++i)
// 	{
// 	  TmpOperations[i]->OutputVector->ClearVector();
// 	}
//   
      

  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {	
      delete TmpOperations[i]->OutputVector;
      delete TmpOperations[i];
    }
  delete TmpOperations[0];
  delete[] TmpOperations;  
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

/*bool FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{  
  char SaveFileName[200];
  sprintf(SaveFileName, "fermions_su2_slater_sym_tmp%d.vec",this->MPINodeNbr);        
  int NbrComponent = this->NbrComponent;
	
  FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation** TmpOperations = new FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation* [architecture->GetNbrThreads()];
  for(int i = 0 ;i < architecture->GetNbrThreads() ;i++)
    {
      TmpOperations[i] = (FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      TmpOperations[i]->SetOutputVector((RealVector*)this->OutputVector->EmptyClone(true));
    }
     
  for(int p = 0 ; p < this->NbrSMPStage ; p++)
    {
      int StageStart = this->FirstComponent + this->GetRankChunkStart(NbrComponent, p,  this->NbrSMPStage);
      int StageSize = this->GetRankChunkSize(NbrComponent, p,  this->NbrSMPStage);
      for (int i = 0; i < architecture->GetNbrThreads() ; ++i)
	{	
	  TmpOperations[i]->SetIndicesRange(StageStart + this->GetRankChunkStart(StageSize, i, architecture->GetNbrThreads()) ,
								    this->GetRankChunkSize(StageSize, i, architecture->GetNbrThreads()) );	    
	}	
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  TmpOperations[i]->OutputVector->ClearVector();
	}
      architecture->SendJobs();
      if (architecture->VerboseMode() == true)
	{
	  char TmpString[512];
	  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	    {
	      sprintf (TmpString, "FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation core operation stage %d on SMP id %d done in %.3f seconds", p,  i, TmpOperations[i]->ExecutionTime);
	      architecture->AddToLog(TmpString);
	    }
	}
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  (*(this->OutputVector)) += (*(TmpOperations[i]->OutputVector));	
	}
      //cout << TmpFirstComponent << " /  " << this->NbrComponent << " (" << ((TmpFirstComponent * 100) / this->NbrComponent) << "%)                   \r";
      //cout << TmpFirstComponent << " /  " << this->NbrComponent << " (" << ((TmpFirstComponent * 100) / this->NbrComponent) << "%)                   \n";
      //cout.flush();      
  
      //sthis->OutputVector->WriteVector(SaveFileName);
    }

  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {	
      delete TmpOperations[i]->OutputVector;
      delete TmpOperations[i];
    }
  delete TmpOperations[0];
  delete[] TmpOperations;
  return true;
}*/


// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__    
  this->MPINodeNbr = architecture->GetNodeNbr();  
  int NbrComponents = this->NbrComponent;
  int FirstComponent = this->FirstComponent;
  int StartingStage = 0;
  for ( int Stage = StartingStage; Stage < this->NbrMPIStage ; Stage++ )  
    {
      if ( !architecture->IsMasterNode() )
	{
	  this->OutputVector->ClearVector();
	}
      int StageDimension = this->GetRankChunkSize(NbrComponents, Stage,  this->NbrMPIStage);
      int StageStart = this->GetRankChunkStart(NbrComponents, Stage,  this->NbrMPIStage);
      if ( this->ResumeIdx >= (StageStart  + StageDimension) )
	{
	  continue;
	}
      else if ( this->ResumeIdx >= StageStart) 
	{      
	  StageDimension =  StageDimension - (this->ResumeIdx - StageStart);
	  StageStart = this->ResumeIdx; 
	}	        	
      this->SetIndicesRange(StageStart +this->GetRankChunkStart(StageDimension,  this->MPINodeNbr,  architecture->GetNbrNodes()), 
				       this->GetRankChunkSize(StageDimension,  this->MPINodeNbr,  architecture->GetNbrNodes())); 
      timeval TotalStartingTime;
      gettimeofday (&TotalStartingTime, 0);				       
      switch (architecture->GetArchitectureID())
	{	 
	  case AbstractArchitecture::MixedMPISMP:
	    this->ArchitectureDependentApplyOperation((SMPArchitecture*)architecture->LocalArchitecture); 
	    break;
	  default:
	    this->RawApplyOperation();
	    break;
	}		
	timeval TotalEndingTime;
	gettimeofday (&TotalEndingTime, 0);
	this->ExecutionTime   = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
	if (architecture->VerboseMode() == true)	
	  {
	    char TmpString[512];      
	    sprintf (TmpString, "FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation process operation stage %d on MPI id %d done in %.3f seconds", Stage, this->MPINodeNbr, this->ExecutionTime);
	    architecture->AddToLog(TmpString);
	  }
	MPI::COMM_WORLD.Barrier();
	architecture->SumVector(*(this->OutputVector));	
	char SaveFileName[200];	
	sprintf(SaveFileName, "fermions_su2_slater_sym_tmp_stage_%d_element_%d.vec", Stage, StageStart + StageDimension);  
	architecture->WriteVector(*(this->OutputVector), SaveFileName);
    }    
      
  return true;
#else
  return this->RawApplyOperation();
#endif
}


// // apply operation for MixedMPISMP architecture
// //
// // architecture = pointer to the architecture
// // return value = true if no error occurs
// 
// bool FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::ArchitectureDependentApplyOperation(MixedMPISMPArchitecture* architecture)
// {
// #ifdef __MPI__  
//   this->MPINodeNbr = architecture->GetNodeNbr();
//   this->OutputVector->ClearVector();
//   for ( int Stage = 0; Stage < this->NbrStage ; Stage++ )  
//     {
//       int StageDimension = this->GetRankChunkSize(this->InitialSpace->GetHilbertSpaceDimension(), Stage,  NbrStage);
//       int StageStart = this->GetRankChunkStart(this->InitialSpace->GetHilbertSpaceDimension(), Stage,  NbrStage);
//       this->SetIndicesRange(StageStart +this->GetRankChunkStart(StageDimension,  this->MPINodeNbr,  architecture->GetNbrNodes()), 
// 				       this->GetRankChunkSize(StageDimension,  this->MPINodeNbr,  architecture->GetNbrNodes())); 				       				             
//       this->ArchitectureDependentApplyOperation((SMPArchitecture*)architecture->LocalArchitecture);      
//     }    
//   MPI_Barrier(MPI_COMM_WORLD);
//   architecture->SumVector(*(this->OutputVector));
//   return true;
// #else
//   return this->RawApplyOperation();
// #endif
// }
