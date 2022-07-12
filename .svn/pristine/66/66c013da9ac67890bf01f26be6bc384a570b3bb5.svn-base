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
#include "Architecture/ArchitectureOperation/FQHESphereMultipleMonomialsTimesSlaterProjectionOperation.h"

#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereThreeLandauLevels.h"
#include "HilbertSpace/FermionOnSphereFourLandauLevels.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>

using std::cout;
using std::endl;
using std::ios;
class BosonOnSphereShort;

// convert a state such that its first nonzero component is equal to 1
//
// state = reference to the state to convert
// return value = converted state

LongRationalVector PutFirstComponentToOne(LongRationalVector& state);
RealVector PutFirstComponentToOne(RealVector& state);

// constructor 
//
// space = pointer to the HilbertSpace
// fileName = file where the kostka Numbers will be stored
// nbrLL = number of Landau levels

FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::FQHESphereMultipleMonomialsTimesSlaterProjectionOperation(ParticleOnSphere* fermionSpace, ParticleOnSphere* lllSpace, ParticleOnSphere* finalSpace, int * matchingConditionsIndex, RealVector* lllVector, int nbrLL, int nbrStates, bool projection, bool symmetry,bool normalize, char * outputFileName, int resumingIndex, bool reverseFluxFlag)
{
  this->FermionSpace = fermionSpace;
  this->LLLSpace = lllSpace;
  this->FinalSpace = finalSpace;
  this->OperationType = AbstractArchitectureOperation::FQHESphereMultipleMonomialsTimesSlaterProjection;
  this->MatchingConditionsIndex = matchingConditionsIndex;
  this->FermionLongRationalVector = 0;
  this->FermionRealVector = new RealVector(this->FermionSpace->GetHilbertSpaceDimension(),true);
  this->LLLRealVector = lllVector;
  this->LLLLongRationalVector = 0; 
  this->Projection = projection;
  this->NbrLL = nbrLL;
  this->Symmetry = symmetry;
  this->Normalize = normalize;
  this->OutputRealVector = new RealVector(this->FinalSpace->GetHilbertSpaceDimension(),true);
  if(lllSpace->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic)
    this->BosonFlag = true;
  else
    this->BosonFlag = false;
  this->OutputFileName = outputFileName;
  this->NbrStates = nbrStates;
  this->Index = -1;
  this->ResumingIndex = resumingIndex;
  this->Rational = false;
  this->ReverseFluxFlag = reverseFluxFlag;
}

// constructor 
//


FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::FQHESphereMultipleMonomialsTimesSlaterProjectionOperation(ParticleOnSphere* fermionSpace, ParticleOnSphere* lllSpace, ParticleOnSphere* finalSpace, int * matchingConditionsIndex, LongRationalVector* lllVector, int nbrLL, int nbrStates, bool projection, bool symmetry,bool normalize,char * outputFileName,int resumingIndex, bool reverseFluxFlag)
{
  this->FermionSpace = fermionSpace;
  this->LLLSpace = lllSpace;
  this->FinalSpace = finalSpace;
  this->OperationType = AbstractArchitectureOperation::FQHESphereMultipleMonomialsTimesSlaterProjection;
  this->MatchingConditionsIndex = matchingConditionsIndex;
  this->FermionLongRationalVector = new LongRationalVector (this->FermionSpace->GetHilbertSpaceDimension(),true);
  this->FermionRealVector = 0;
  this->LLLRealVector = 0;
  this->LLLLongRationalVector = lllVector; 
  this->Projection = projection;
  this->NbrLL = nbrLL;
  this->Symmetry = symmetry;
  this->Normalize = normalize;
  this->OutputLongRationalVector = new LongRationalVector (this->FinalSpace->GetHilbertSpaceDimension(),true);
  this->OutputRealVector = 0;
  if(lllSpace->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic)
    this->BosonFlag = true;
  else
    this->BosonFlag = false;
  this->OutputFileName = outputFileName;
  this->NbrStates = nbrStates;
  this->Index = -1;
  this->ResumingIndex = resumingIndex;
  this->Rational = true;
  this->ReverseFluxFlag = false;
}


// copy constructor 
//
// operation = reference on operation to copy

FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::FQHESphereMultipleMonomialsTimesSlaterProjectionOperation(const FQHESphereMultipleMonomialsTimesSlaterProjectionOperation& operation)
{
  this->FermionSpace = operation.FermionSpace;
  this->LLLSpace = operation.LLLSpace;
  this->FinalSpace = (FermionOnSphere *) operation.FinalSpace->Clone();
  this->OperationType = AbstractArchitectureOperation::FQHESphereMonomialsTimesSlaterProjection;
  this->MatchingConditionsIndex = operation.MatchingConditionsIndex;
  if (operation.Rational == true)
    {
      this->FermionLongRationalVector = new LongRationalVector (this->FermionSpace->GetHilbertSpaceDimension(),true);
      this->OutputLongRationalVector = new LongRationalVector (operation.FinalSpace->GetHilbertSpaceDimension(),true);
    }
  else
    {
      this->FermionRealVector = new RealVector(this->FermionSpace->GetHilbertSpaceDimension(),true);
      this->OutputRealVector = new RealVector(operation.FinalSpace->GetHilbertSpaceDimension(),true);
    }
  this->LLLLongRationalVector = operation.LLLLongRationalVector;
  this->LLLRealVector = operation.LLLRealVector;
  this->Projection = operation.Projection;
  this->NbrLL = operation.NbrLL;
  this->Symmetry = operation.Symmetry;
  this->BosonFlag = operation.BosonFlag;
  this->Normalize = operation.Normalize;
  this->OutputFileName = new char [64];
  strncpy (this->OutputFileName,operation.OutputFileName,strlen(operation.OutputFileName));
  this->NbrStates = operation.NbrStates;
  this->Index = operation.Index;
  this->ResumingIndex = operation.ResumingIndex;
  this->Rational = operation.Rational;
  this->ReverseFluxFlag = operation.ReverseFluxFlag;
}


// destructor
//

FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::~FQHESphereMultipleMonomialsTimesSlaterProjectionOperation()
{
  if (this->Rational == true)
    {
      delete this->FermionLongRationalVector;
      delete this->OutputLongRationalVector;
    }
  else
    {
      delete this->FermionRealVector;
      delete this->OutputRealVector;
    }
}


void FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::SetIndex(int index)
{
  this->Index = index;
  if (this->Index != -1 )
    {
      char *Tmp;
      Tmp = strstr(this->OutputFileName,".");
      if(Tmp == 0)
	cout<<"c'est pas bon"<<endl;
      sprintf (Tmp,".%d.vec",index);
      if (this->Rational == true)
	{
	  (*(this->FermionLongRationalVector))[this->MatchingConditionsIndex[this->Index]] = 1l;
	}
      else
	{
	  (*(this->FermionRealVector))[this->MatchingConditionsIndex[this->Index]] = 1l;
	}
    }
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::Clone()
{
  return new FQHESphereMultipleMonomialsTimesSlaterProjectionOperation (*this);
}

// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::RawApplyOperation()
{
  double Error = MACHINE_PRECISION;
  if(this->Index != -1)
    {
      timeval TotalStartingTime;
      gettimeofday (&TotalStartingTime, 0);
      
      if(this->Rational == true)
	{
	  if(this->BosonFlag == true)
	    {
	      if (this->NbrLL > 2)
		{
		  switch (this->NbrLL)
		    {
		    case 3:
		      ((FermionOnSphereThreeLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicState(*(this->LLLLongRationalVector), *(this->FermionLongRationalVector) , *(this->OutputLongRationalVector), (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphere*)this->FinalSpace,0, this->LLLSpace->GetHilbertSpaceDimension());
		      break;
		    case 4:
		      ((FermionOnSphereFourLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicState(*(this->LLLLongRationalVector), *(this->FermionLongRationalVector) , *(this->OutputLongRationalVector), (BosonOnSphereShort*) this->LLLSpace,(FermionOnSphere*)this->FinalSpace,0, this->LLLSpace->GetHilbertSpaceDimension());
		      break;
		    }
		}
	      else
		{
		  if( this->Projection == true )
		    {
		      if( this->Symmetry == true)
			((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicStateSymmetric( *this->LLLLongRationalVector , *this->FermionLongRationalVector, *this->OutputLongRationalVector, (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphere*) this->FinalSpace, 0 , this->LLLSpace->GetHilbertSpaceDimension());
		      else
			((FermionOnSphereTwoLandauLevels *) this->FermionSpace)->BosonicStateTimeFermionicState( *this->LLLLongRationalVector, *this->FermionLongRationalVector , *this->OutputLongRationalVector, (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphere*) this->FinalSpace, 0, this->LLLSpace->GetHilbertSpaceDimension());
		    }
		  else
		    ((FermionOnSphereTwoLandauLevels *) this->FermionSpace)->BosonicStateTimeFermionicState( *this->LLLLongRationalVector, *this->FermionLongRationalVector , *this->OutputLongRationalVector, (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphereTwoLandauLevels*) this->FinalSpace, 0, this->LLLSpace->GetHilbertSpaceDimension());
		}
	    }
	  else
	    {
	      if (this->NbrLL > 2)
		{
		  switch (this->NbrLL)
		    {
		    case 3:
		      ((FermionOnSphereThreeLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicState(*(this->LLLLongRationalVector), *(this->FermionLongRationalVector) , *(this->OutputLongRationalVector), (FermionOnSphere*)this->LLLSpace, (BosonOnSphereShort*)this->FinalSpace, 0, this->LLLSpace->GetHilbertSpaceDimension());
		      break;
		    case 4:
		      ((FermionOnSphereFourLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicState(*(this->LLLLongRationalVector), *(this->FermionLongRationalVector) , *(this->OutputLongRationalVector), (FermionOnSphere*)this->LLLSpace, (BosonOnSphereShort*)this->FinalSpace,0, this->LLLSpace->GetHilbertSpaceDimension());
		      break;
		    }
		}
	      else
		{
		  if( this->Symmetry == true)
		    ((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicStateSymmetric( *this->LLLLongRationalVector, *this->FermionLongRationalVector, *this->OutputLongRationalVector, (FermionOnSphere*) this->LLLSpace, (BosonOnSphereShort*) this->FinalSpace, 0, this->LLLSpace->GetHilbertSpaceDimension());
		  else
		    ((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicState(*(this->LLLLongRationalVector), *(this->FermionLongRationalVector) , *(this->OutputLongRationalVector), (FermionOnSphere*)this->LLLSpace, (BosonOnSphereShort*)this->FinalSpace,0, this->LLLSpace->GetHilbertSpaceDimension());
		}
	    }
	  timeval TotalEndingTime;
	  gettimeofday (&TotalEndingTime, 0);
	  
	  PutFirstComponentToOne(*(this->OutputLongRationalVector));
	  
	  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
	  cout <<this->Index<<" : "<< Dt << "s" << endl;
	  
	  this->OutputLongRationalVector->WriteVector(this->OutputFileName);
	  this->OutputLongRationalVector->ClearVector();
	  this->FermionLongRationalVector->ClearVector();
	}
      else
	{
	  if(this->BosonFlag == true)
	    {
	      if (this->NbrLL > 2)
		{
		  switch (this->NbrLL)
		    {
		    case 3:
		      ((FermionOnSphereThreeLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicState(*(this->LLLRealVector), *(this->FermionRealVector) , *(this->OutputRealVector), (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphere*)this->FinalSpace,0, this->LLLSpace->GetHilbertSpaceDimension(),this->ReverseFluxFlag);
		      break;
		    case 4:
		      ((FermionOnSphereFourLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicState(*(this->LLLRealVector), *(this->FermionRealVector) , *(this->OutputRealVector), (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphere*)this->FinalSpace,0, this->LLLSpace->GetHilbertSpaceDimension(),this->ReverseFluxFlag);
		      break;
		    }
		}
	      else
		{
		  if( this->Projection == true )
		    {
		      if( this->Symmetry == true )
			((FermionOnSphereTwoLandauLevels *) this->FermionSpace)->BosonicStateTimeFermionicStateSymmetric( *this->LLLRealVector, *this->FermionRealVector, *this->OutputRealVector, (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphere*) this->FinalSpace, 0, this->LLLSpace->GetHilbertSpaceDimension(),this->ReverseFluxFlag);
		      else
			((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicState( *this->LLLRealVector, *this->FermionRealVector, *this->OutputRealVector, (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphere*) this->FinalSpace, 0, this->LLLSpace->GetHilbertSpaceDimension(),this->ReverseFluxFlag);
		    }
		  else
		    ((FermionOnSphereTwoLandauLevels *) this->FermionSpace)->BosonicStateTimeFermionicState( *this->LLLRealVector, *this->FermionRealVector, *this->OutputRealVector, (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphereTwoLandauLevels*) this->FinalSpace, 0, this->LLLSpace->GetHilbertSpaceDimension());
		}
	    }
	  else
	    {
	      if (this->NbrLL > 2)
		{
		  switch (this->NbrLL)
		    {
		    case 3:
		      ((FermionOnSphereThreeLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicState(*(this->LLLRealVector), *(this->FermionRealVector) , *(this->OutputRealVector), (FermionOnSphere*)this->LLLSpace, (BosonOnSphereShort*)this->FinalSpace,0, this->LLLSpace->GetHilbertSpaceDimension());
		      break;
		    case 4:
		      ((FermionOnSphereFourLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicState(*(this->LLLRealVector), *(this->FermionRealVector) , *(this->OutputRealVector), (FermionOnSphere*)this->LLLSpace,  (BosonOnSphereShort*)this->FinalSpace,0, this->LLLSpace->GetHilbertSpaceDimension());
		      break;
		    }
		}
	      else
		{
		  if (this->Projection == true)
		    {
		      if(this->Symmetry)
			((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicStateSymmetric(*(this->LLLRealVector), *(this->FermionRealVector) , *(this->OutputRealVector), (FermionOnSphere*)this->LLLSpace, (BosonOnSphereShort*) this->FinalSpace, 0, this->LLLSpace->GetHilbertSpaceDimension());
		      else
			((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicState(*(this->LLLRealVector), *(this->FermionRealVector) , *(this->OutputRealVector), (FermionOnSphere*)this->LLLSpace, (BosonOnSphereShort*)this->FinalSpace,0, this->LLLSpace->GetHilbertSpaceDimension());
		    }
		  else
		    ((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicState(*(this->LLLRealVector), *(this->FermionRealVector) , *(this->OutputRealVector), (FermionOnSphere*) this->LLLSpace, (BosonOnSphereTwoLandauLevels *) this->FinalSpace , 0 , this->LLLSpace->GetHilbertSpaceDimension());
		}
	    }
	  timeval TotalEndingTime;
	  gettimeofday (&TotalEndingTime, 0);
	  if (this->Normalize == true)
	    {
	      if(this->Projection == true)
		{
		  if(this->BosonFlag == true)
		    ((FermionOnSphere*)FinalSpace)->ConvertFromUnnormalizedMonomial(*(this->OutputRealVector),0l,true);
		  else
		    ((BosonOnSphereShort*)FinalSpace)->ConvertFromUnnormalizedMonomial(*(this->OutputRealVector),0l,true);
		  this->OutputRealVector->WriteVector(this->OutputFileName);
		}
	      else
		{
		  if(this->BosonFlag == false)
		    {
		      ((BosonOnSphereTwoLandauLevels*)FinalSpace)->ConvertFromUnnormalizedMonomial(*(this->OutputRealVector),0l,true);
		      
		      this->OutputRealVector->WriteVector(this->OutputFileName);
		    }
		}
	    }
	  else
	    {
	      if(this->Projection == true)
		{
		  PutFirstComponentToOne(*(this->OutputRealVector));
		}
	      this->OutputRealVector->WriteVector(this->OutputFileName);
	    }
	  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
	  cout <<this->Index<<" : "<< Dt << "s" << endl;
	  this->OutputRealVector->ClearVector();
	  this->FermionRealVector->ClearVector();
	}
    }
  return true;
}


bool FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture)
{
  for(int i = this->ResumingIndex ;i < this->NbrStates ;i++)
    {
      this->SetIndex(i);
      this->RawApplyOperation();
    }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int CurrentIndex = this->ResumingIndex;
  FQHESphereMultipleMonomialsTimesSlaterProjectionOperation** TmpOperations = new FQHESphereMultipleMonomialsTimesSlaterProjectionOperation* [architecture->GetNbrThreads()];
  for(int i = 0 ;i < architecture->GetNbrThreads() ;i++)
    {
      TmpOperations[i] = (FQHESphereMultipleMonomialsTimesSlaterProjectionOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  while(CurrentIndex < this->NbrStates - architecture->GetNbrThreads())
    {
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	{	
	  TmpOperations[i]->SetIndex(CurrentIndex);
	  CurrentIndex++;
	}
      architecture->SendJobs();
    }
  
  int RemainingStates = this->NbrStates - CurrentIndex;
  if(RemainingStates != 0)
    {
      for (int i = 0; i < RemainingStates; ++i)
	{	
	  TmpOperations[i]->SetIndex(CurrentIndex);
	  CurrentIndex++;
	}

      for (int i = RemainingStates; i <  architecture->GetNbrThreads(); i++)
	 TmpOperations[i]->SetIndex(-1);

      architecture->SendJobs();
    }
  
  for (int i = 1; i < architecture->GetNbrThreads() ; ++i)
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

bool FQHESphereMultipleMonomialsTimesSlaterProjectionOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
  cout <<"Not yet implemented"<<endl;
  return false;
  /*
#ifdef __MPI__
   if (architecture->IsMasterNode())
     {
       this->OutputVector->ClearVector();
       RealVector* TmpVector = (RealVector*) this->OutputVector->EmptyClone(true);
       int TmpRange[2];
       int TmpNbrSlaves = architecture->GetNbrSlaveNodes();
       int TmpSlaveID = 0;
       for (int i = 0; (i < TmpNbrSlaves) ; ++i)
	 {
	   TmpRange[0] = TmpFirstComponent;
	   TmpRange[1] = Step;
	   TmpFirstComponent += Step;
	   if (TmpFirstComponent >= (this->FirstComponent +this->NbrComponent))
	     Step = this->FirstComponent + this->NbrComponent - TmpFirstComponent;
	 }
       while ((Step > 0) && ((TmpSlaveID = architecture->WaitAnySlave())))
	 {
	   architecture->SumVector(*(this->OutputVector));
	 }
       while ((TmpNbrSlaves > 0) && ((TmpSlaveID = architecture->WaitAnySlave())))
	 {
	   --TmpNbrSlaves;
	   architecture->SumVector(*(this->OutputVector));
	 }
     }
   else
     {
       int TmpRange[2];
       int TmpNbrElement = 0;
       while ((architecture->ReceiveFromMaster(TmpRange, TmpNbrElement) == true) && (TmpRange[1] > 0))
	 {
	   this->OutputVector->ClearVector();
	   this->SetIndicesRange(TmpRange[0], TmpRange[1]); 
	   if (architecture->GetLocalArchitecture()->GetArchitectureID() == AbstractArchitecture::SMP)
	     this->ArchitectureDependentApplyOperation((SMPArchitecture*) architecture->GetLocalArchitecture());
	   else
	     this->RawApplyOperation();
	   architecture->SendDone();
	   architecture->SumVector(*(this->OutputVector));
	 }
     }
  return true;
#else
  return this->RawApplyOperation();
#endif
  */
}
