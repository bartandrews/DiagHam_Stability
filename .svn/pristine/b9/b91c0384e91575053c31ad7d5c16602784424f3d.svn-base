////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of monomials product operation 	              //
//                                                                            //
//                        last modification : 23/10/2002                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereMonomialsProductOperation.h"
#include "Vector/RealVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor 
//
// space = pointer to the HilbertSpace
// sourceVector = vector to be multiplied by the operator
// destinationVector = vector where the result has to be stored

MonomialsProductOperation::MonomialsProductOperation (ParticleOnSphere* space,ParticleOnSphere * finalSpace, RealVector* sourceVector1,RealVector* sourceVector2, RealVector* destinationVector, bool normalize)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->Space1 = space;
  this->FinalSpace = finalSpace;
  this->SourceVector1 = sourceVector1;
  this->SourceVector2 = sourceVector2;
  this->DestinationVector = destinationVector;
  this->OperationType = AbstractArchitectureOperation::FQHESphereMonomialsProductOperation;
  this->Squaring = true;
  if(space->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic)
    this->FirstStateFlag = false;
  else
    this->FirstStateFlag = true;
	this->FinalStateFlag = false;
  this->NormalizeFlag = normalize;
}

MonomialsProductOperation::MonomialsProductOperation(ParticleOnSphere* space1,ParticleOnSphere* space2,ParticleOnSphere * finalSpace,RealVector* sourceVector1, RealVector* sourceVector2, RealVector* destinationVector, bool normalize)
{
  this->FirstComponent = 0;
  this->NbrComponent = space1->GetHilbertSpaceDimension();
  this->Space1 = space1;
  this->Space2 = space2;
  this->FinalSpace = finalSpace;
  this->SourceVector1 = sourceVector1;
  this->SourceVector2 = sourceVector2;
  this->DestinationVector = destinationVector;
  this->OperationType = AbstractArchitectureOperation::FQHESphereMonomialsProductOperation;
  this->Squaring = false;
  if(space1->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic)
    this->FirstStateFlag = false;
  else
    this->FirstStateFlag = true;
  if(finalSpace->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic)
    this->FinalStateFlag = false;
  else
    this->FinalStateFlag = true;
  this->NormalizeFlag = normalize;
}


// copy constructor 
//
// operation = reference on operation to copy

MonomialsProductOperation::MonomialsProductOperation(const MonomialsProductOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Space1 = operation.Space1;
  this->Space2 = operation.Space2;
  this->FinalSpace = operation.FinalSpace;
  this->SourceVector1 = operation.SourceVector1;
  this->SourceVector2 = operation.SourceVector2;
  this->DestinationVector = operation.DestinationVector;
  this->OperationType = AbstractArchitectureOperation::FQHESphereMonomialsProductOperation;
  this->Squaring = operation.Squaring;
  this->FinalStateFlag = operation.FinalStateFlag;
  this->FirstStateFlag = operation.FirstStateFlag;
  this->NormalizeFlag = operation.NormalizeFlag;
}
  
// constructor from a master node information
//
// operator = pointer to the operator to use
// architecture = pointer to the distributed architecture to use for communications

MonomialsProductOperation::MonomialsProductOperation(ParticleOnSphere* space,ParticleOnSphere * finalSpace, bool normalize, SimpleMPIArchitecture* architecture)
{
  this->Space1 = space;
  this->FinalSpace = finalSpace;
  this->OperationType = AbstractArchitectureOperation::FQHESphereMonomialsProductOperation;
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  this->FirstComponent = (int) TmpMinimumIndex;
  this->NbrComponent = (int) (TmpMaximumIndex - TmpMinimumIndex + 1l);
  this->SourceVector1 = (RealVector*) architecture->ScatterVector();
  this->SourceVector2 = (RealVector*) architecture->ScatterVector();;
  this->DestinationVector =   (RealVector*) architecture->BroadcastVectorType();
  this->Squaring = true;
  if(space->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic)
    this->FirstStateFlag = false;
  else
    this->FirstStateFlag = true;
  this->NormalizeFlag = normalize;
}

// constructor from a master node information
//
// operator = pointer to the operator to use
// architecture = pointer to the distributed architecture to use for communications

MonomialsProductOperation::MonomialsProductOperation(ParticleOnSphere* space1,ParticleOnSphere* space2,ParticleOnSphere * finalSpace, bool normalize, SimpleMPIArchitecture* architecture)
{
  this->Space1 = space1;
  this->Space2 = space2;
  this->FinalSpace = finalSpace;
  this->OperationType = AbstractArchitectureOperation::FQHESphereMonomialsProductOperation;
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  this->FirstComponent = (int) TmpMinimumIndex;
  this->NbrComponent = (int) (TmpMaximumIndex - TmpMinimumIndex + 1l);
  this->SourceVector1 = (RealVector*) architecture->ScatterVector();
  this->SourceVector2 = (RealVector*) architecture->ScatterVector();
  this->DestinationVector =   (RealVector*) architecture->BroadcastVectorType();
  this->Squaring = false;
  if(space1->GetParticleStatistic()==ParticleOnSphere::BosonicStatistic)
    this->FirstStateFlag = false;
  else
    this->FirstStateFlag = true;
  if(finalSpace->GetParticleStatistic()==ParticleOnSphere::BosonicStatistic)
    this->FinalStateFlag = false;
  else
    this->FinalStateFlag = true;
  this->NormalizeFlag = normalize;
}

// destructor
//

MonomialsProductOperation::~MonomialsProductOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void MonomialsProductOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void MonomialsProductOperation::SetDestinationVector (RealVector* DestinationVector)
{
  this->DestinationVector = DestinationVector;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* MonomialsProductOperation::Clone()
{
  return new MonomialsProductOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool MonomialsProductOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  if(this->FinalStateFlag == true)
    {
      ((BosonOnSphereShort *) this->Space1)->BosonicStateTimeFermionicState(*(this->SourceVector1),*(this->SourceVector2),*(this->DestinationVector),((FermionOnSphere *) this->Space2),this->FirstComponent,this->NbrComponent,((FermionOnSphere *) this->FinalSpace));
      if(this->NormalizeFlag == true)
	{
	  ((FermionOnSphere *) this->FinalSpace)->ConvertFromUnnormalizedMonomial(*(this->DestinationVector),0);
	}
    }
  else
    {
      if(this->FirstStateFlag == true)
	{
	  if(this->Squaring == true)
	    {
	      ((BosonOnSphereShort *) this->FinalSpace)->FermionicStateTimeFermionicState(*(this->SourceVector1),*(this->SourceVector1),*(this->DestinationVector),((FermionOnSphere *) this->Space1),this->FirstComponent,this->NbrComponent);
	    }
	  else
	    {
	      ((BosonOnSphereShort *) this->FinalSpace)->FermionicStateTimeFermionicState(*(this->SourceVector1),*(this->SourceVector2),*(this->DestinationVector),((FermionOnSphere *) this->Space1),((FermionOnSphere *) this->Space2),this->FirstComponent,this->NbrComponent);
	    }
	  
	}
      else
	{
	  if(this->Squaring == true)
	    {
	      ((BosonOnSphereShort *) this->Space1)->BosonicStateTimeBosonicState(*(this->SourceVector1), *(this->SourceVector2), *(this->DestinationVector),this->FirstComponent,this->NbrComponent,((BosonOnSphereShort *) this->FinalSpace));
	    }
	  else
	    {	
	      ((BosonOnSphereShort *)this->Space1)->BosonicStateTimeBosonicState(*(this->SourceVector1), *(this->SourceVector2),*(this->DestinationVector),((BosonOnSphereShort *)this->Space2),this->FirstComponent,this->NbrComponent,((BosonOnSphereShort *)this->FinalSpace));
	    }
	}
      if(this->NormalizeFlag == true)
	{
	  ((BosonOnSphereShort *) this->FinalSpace)->ConvertFromUnnormalizedMonomial(*(this->DestinationVector),0);
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

bool MonomialsProductOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  char **TmpFileNames = new char* [100];
  int NBR = (int) (this->NbrComponent/100);
  this->DestinationVector->ClearVector();
  if (NBR>=8)
    {
      for(int j = 0; j < 100; j++)
	{
	  if(j != 99)
	    {
	      timeval TotalStartingTime;
	      gettimeofday (&TotalStartingTime, 0);
	      int Step =  NBR/ architecture->GetNbrThreads();
	      int TmpFirstComponent = this->FirstComponent+j*NBR;
	      int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
	      MonomialsProductOperation** TmpOperations = new MonomialsProductOperation* [architecture->GetNbrThreads()];
	      for (int i = 0; i < ReducedNbrThreads; ++i)
		{
		  //cout << "Step : " <<Step<<endl;
		  //cout <<"TmpFirstComponent = "<<TmpFirstComponent<<endl;
		  TmpOperations[i] = (MonomialsProductOperation*) this->Clone();
		  architecture->SetThreadOperation(TmpOperations[i], i);
		  TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
		  TmpFirstComponent += Step;
		}
	      TmpOperations[ReducedNbrThreads] = (MonomialsProductOperation*) this->Clone();
	      architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
	      TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent,this->FirstComponent+(j+1)*NBR - TmpFirstComponent);
	      //cout<<"TmpFirstComponent= "<<TmpFirstComponent<<endl;
	      //cout <<"Step = "<<FirstComponent+(j+1)*NBR - TmpFirstComponent<<endl;
	      
	      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
		{
		  TmpOperations[i]->SetDestinationVector((RealVector*)this->DestinationVector->EmptyClone(true));
		}
	      architecture->SendJobs();
	      TmpFileNames[j] = new char [20];
	      sprintf(TmpFileNames[j],"vector_%d.vec",j);
	      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
		{	
		  (*(this->DestinationVector)) += (*(TmpOperations[i]->DestinationVector));
		  delete TmpOperations[i]->DestinationVector;
		  delete TmpOperations[i];
		}
	      this->DestinationVector->WriteVector(TmpFileNames[j]);
	      delete TmpOperations[0];
	      delete[] TmpOperations;
	      timeval TotalEndingTime;
	      gettimeofday (&TotalEndingTime, 0);
	      double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
			    (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
	      cout << this->FirstComponent+j*NBR << " " <<  NBR << " : " << Dt << "s" << endl;
	    }
	  else
	    {
	      timeval TotalStartingTime;
	      gettimeofday (&TotalStartingTime, 0);
	      int Step = (this->NbrComponent-(this->FirstComponent+99*NBR))/ architecture->GetNbrThreads();
	      int TmpFirstComponent = this->FirstComponent+j*NBR;
	      int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
	      MonomialsProductOperation** TmpOperations = new MonomialsProductOperation* [architecture->GetNbrThreads()];
	      for (int i = 0; i < ReducedNbrThreads; ++i)
		{
		  //cout << "Step : " <<Step<<endl;
		  //cout <<"TmpFirstComponent = "<<TmpFirstComponent<<endl;
		  TmpOperations[i] = (MonomialsProductOperation*) this->Clone();
		  architecture->SetThreadOperation(TmpOperations[i], i);
		  TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
		  TmpFirstComponent += Step;
		}
	      TmpOperations[ReducedNbrThreads] = (MonomialsProductOperation*) this->Clone();
	      architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
	      TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent,this->FirstComponent+this->NbrComponent - TmpFirstComponent);
	      //cout <<"step = "<<this->FirstComponent+this->NbrComponent - TmpFirstComponent<<endl;
	      //cout<<"TmpFirstComponent= "<<TmpFirstComponent<<endl;
	      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
		{
		  TmpOperations[i]->SetDestinationVector((RealVector*)this->DestinationVector->EmptyClone(true));
		}
	      architecture->SendJobs();
	      TmpFileNames[j] = new char [20];
	      sprintf(TmpFileNames[j],"vector_%d.vec",j);
	      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
		{
		  (*(this->DestinationVector)) += (*(TmpOperations[i]->DestinationVector));
		  delete TmpOperations[i]->DestinationVector;
		  delete TmpOperations[i];
		}
	      this->DestinationVector->WriteVector(TmpFileNames[j]);
	      delete TmpOperations[0];
	      delete[] TmpOperations;
	      timeval TotalEndingTime;
	      gettimeofday (&TotalEndingTime, 0);
	      double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
			    (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
	      cout << this->FirstComponent+j*NBR << " " << (this->NbrComponent-(this->FirstComponent+19*NBR))<< " : " << Dt << "s" << endl;
	    }
	}
      return true;
    }
  else
    {
      int Step;
      int TmpFirstComponent = this->FirstComponent;
      int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
      MonomialsProductOperation** TmpOperations = new MonomialsProductOperation* [architecture->GetNbrThreads()];
      unsigned long Sum;
      unsigned long Sigma = 0x0ul;
      for (int i= this->FirstComponent;i<this->NbrComponent;i++)
	{
	  Sigma += this->TimeEvaluationFonction(i);
	}
      Sigma = Sigma/architecture->GetNbrThreads();
      
      cout <<"Procesors Numbers: "<< architecture->GetNbrThreads()<<endl;
      //cout << "Sigma = " <<Sigma<<endl;
      //cout <<"this->NbrComponent = "<<this->NbrComponent<<endl;
      //cout <<"dimension"<< FinalSpace->GetHilbertSpaceDimension()<<endl;
      
      for (int i = 0; i < ReducedNbrThreads; ++i)
	{
	  Step = 0;
	  Sum = 0;
	  while(Sum < Sigma)
	    {
	      Sum += this->TimeEvaluationFonction(TmpFirstComponent+Step);
	      Step++;
	    }
	  //cout << "Step : " <<Step<<endl;
	  //cout <<"TmpFirstComponent = "<<TmpFirstComponent<<endl;
	  TmpOperations[i] = (MonomialsProductOperation*) this->Clone();
	  architecture->SetThreadOperation(TmpOperations[i], i);
	  TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
	  TmpFirstComponent += Step;
	}
      TmpOperations[ReducedNbrThreads] = (MonomialsProductOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
      TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent + this->NbrComponent - TmpFirstComponent);
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  TmpOperations[i]->SetDestinationVector((RealVector*)this->DestinationVector->EmptyClone(true));
	}
      architecture->SendJobs();
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  (*(this->DestinationVector)) += (*(TmpOperations[i]->DestinationVector));
	  delete TmpOperations[i]->DestinationVector;
	  delete TmpOperations[i];
	}
      delete TmpOperations[0];
      delete[] TmpOperations;
      return true;
    }
  return false;
}


// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool MonomialsProductOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
/*#ifdef __MPI__
  if (architecture->IsMasterNode())
    {
      if (architecture->RequestOperation(this->OperationType) == false)
	{
	  return false;
	}
      architecture->ScatterVector(this->SourceVector);
      architecture->BroadcastVectorType(this->DestinationVector);  
    }
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  this->FirstComponent = (int) TmpMinimumIndex;  
  this->NbrComponent = (int) (TmpMaximumIndex - TmpMinimumIndex + 1l);
  timeval TotalStartingTime;
  if (architecture->VerboseMode())
    gettimeofday (&TotalStartingTime, 0);
  if (architecture->GetLocalArchitecture()->GetArchitectureID() == AbstractArchitecture::SMP)
    this->ArchitectureDependentApplyOperation((SMPArchitecture*) architecture->GetLocalArchitecture());
  else
    this->RawApplyOperation();
  if (architecture->VerboseMode())
    {
      timeval TotalEndingTime;
      gettimeofday (&TotalEndingTime, 0);
      double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		    (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
      char TmpString[256];
      sprintf (TmpString, "SchurDecompositionOperation core operation done in %.3f seconds", Dt);
      architecture->AddToLog(TmpString);
    }
  if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
    gettimeofday (&TotalStartingTime, 0);
  architecture->SumVector(*(this->DestinationVector));
  if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
    {
      timeval TotalEndingTime;
      gettimeofday (&TotalEndingTime, 0);
      double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		    (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
      char TmpString[256];
      sprintf (TmpString, "SchurDecompositionOperation sum operation done in %.3f seconds", Dt);
    architecture->AddToLog(TmpString, true);
    }
  if (architecture->IsMasterNode() == false)
    {
      delete this->DestinationVector;
      delete this->SourceVector;
    }
  
  return true;
#else*/
  return this->RawApplyOperation();
//#endif
}	
