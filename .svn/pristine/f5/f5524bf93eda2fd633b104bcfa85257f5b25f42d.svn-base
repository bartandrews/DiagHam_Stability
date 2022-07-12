////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of SMP Architecture                         //
//                                                                            //
//                        last modification : 30/04/2002                      //
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
#include "Architecture/SMPArchitecture.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MatrixMatrixMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#ifdef __SMP__
#include <pthread.h>
#endif
#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


struct ThreadMultiplyParameter
{
  int ThreadID;
  int* Flag;
#ifdef __SMP__
  pthread_mutex_t* mut;
#endif

  int FirstComponent;
  int NbrComponent;
  AbstractHamiltonian* Hamiltonian;
  Vector* SourceVector;
  Vector* DestinationVector;
};


// function used by a thread for multiplication
//
// param = pointer to additional parameters, has to be cast into ThreadMultiplyParameter pointer
// return value = unused pointer (null)
void* ThreadMultiply(void* param);


// thread for multiplicationmain function
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)
void* ThreadMain(void* param);

// thread for multiplicationmain function
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)
void* ThreadExecuteOperation(void* param);


// constructor
//
// nbrProcesses = number of processes to run simultaneously (in principle, the number of processors that can be allocated)

SMPArchitecture::SMPArchitecture(int nbrProcesses)
{
  this->NbrProcesses = nbrProcesses;
#ifdef __SMP__
  this->ThreadParameters = new ThreadMainParameter [this->NbrProcesses];
  this->Threads = new pthread_t [this->NbrProcesses];
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      this->ThreadParameters[i].ThreadState = SMPArchitecture::Wait;
      this->ThreadParameters[i].ThreadID = i;
 /*     if (pthread_create (&(this->Threads[i]), 0, ThreadMain, (void*) &(this->ThreadParameters[i])))
	{
	  cout << "error, cannot create thread" << endl;
	  exit(1);
	}*/
    }
#endif
}
  
// destructor
//

SMPArchitecture::~SMPArchitecture()
{
#ifdef __SMP__
/*  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      if (this->ThreadParameters[i].ThreadState != SMPArchitecture::Dead)
	{
	  pthread_mutex_lock(this->ThreadParameters[i].LocalMutex);
	  this->ThreadParameters[i].ThreadState = SMPArchitecture::Exit;
	  pthread_mutex_unlock(this->ThreadParameters[i].LocalMutex);
	  pthread_cond_signal(this->ThreadParameters[i].LocalCondition);       
	}
    }
  void* ReturnValue;
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      (void) pthread_join (this->Threads[i], &ReturnValue);
    }*/
#endif
}
  
// execute an architecture-dependent vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (VectorHamiltonianMultiplyOperation* operation)
{
  int Step = operation->GetDestinationVector()->GetVectorDimension() / this->NbrProcesses;
  int FirstComponent = 0;
  int ReducedNbrProcesses = this->NbrProcesses - 1;
  operation->GetDestinationVector()->ClearVector();
  for (int i = 0; i < ReducedNbrProcesses; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrProcesses].Operation = operation->Clone();
  ((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[ReducedNbrProcesses].Operation))->SetIndicesRange(FirstComponent, 
												  operation->GetDestinationVector()->GetVectorDimension() - FirstComponent);  
  for (int i = 1; i < this->NbrProcesses; ++i)
    {
      ((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->SetDestinationVector(operation->GetDestinationVector()->EmptyClone(true));
    }
  this->SendJobs();
  for (int i = 1; i < this->NbrProcesses; ++i)
    {
      (*(operation->GetDestinationVector())) += (*(((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->GetDestinationVector()));
      delete ((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->GetDestinationVector();
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}
  
// execute an architecture-dependent add real linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (AddRealLinearCombinationOperation* operation)
{
  int Step = operation->GetDestinationVector()->GetVectorDimension() / this->NbrProcesses;
  int FirstComponent = 0;
  int ReducedNbrProcesses = this->NbrProcesses - 1;
  for (int i = 0; i < ReducedNbrProcesses; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((AddRealLinearCombinationOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrProcesses].Operation = operation->Clone();
  ((AddRealLinearCombinationOperation*) (this->ThreadParameters[ReducedNbrProcesses].Operation))->SetIndicesRange(FirstComponent, 
														  operation->GetDestinationVector()
														  ->GetVectorDimension() - FirstComponent);  
  this->SendJobs();
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}  

// execute an architecture-dependent multiple real scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (MultipleRealScalarProductOperation* operation)
{
  int Step = operation->GetNbrScalarProduct() / this->NbrProcesses;
  int FirstComponent = 0;
  int ReducedNbrProcesses = this->NbrProcesses - 1;
  for (int i = 0; i < ReducedNbrProcesses; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((AddRealLinearCombinationOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrProcesses].Operation = operation->Clone();
  ((AddRealLinearCombinationOperation*) (this->ThreadParameters[ReducedNbrProcesses].Operation))->SetIndicesRange(FirstComponent, 
														  operation->GetNbrScalarProduct() - FirstComponent);  
  this->SendJobs();
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}  

// execute an architecture-dependent matrix matrix multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (MatrixMatrixMultiplyOperation* operation)
{
  int Step = operation->GetDestinationMatrix()->GetNbrRow() / this->NbrProcesses;
  int FirstComponent = 0;
  int ReducedNbrProcesses = this->NbrProcesses - 1;
  for (int i = 0; i < ReducedNbrProcesses; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((MatrixMatrixMultiplyOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrProcesses].Operation = operation->Clone();
  ((MatrixMatrixMultiplyOperation*) (this->ThreadParameters[ReducedNbrProcesses].Operation))->SetIndicesRange(FirstComponent, 
													      operation->GetDestinationMatrix()->GetNbrRow() - 
													      FirstComponent);  
  this->SendJobs();
  for (int i = 1; i < this->NbrProcesses; ++i)
    {
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}
    
// multiply a vector by an hamiltonian and store the result in another vector
//
// hamiltonian = pointer to the hamiltonian to use
// vSource = vector to multiply 
// vDestination = vector where result has to be stored 

void SMPArchitecture::Multiply (AbstractHamiltonian* hamiltonian, Vector& vSource, Vector& vDestination)
{  
  return;
}

// execute an architecture-dependent QHE particle hamiltonian precalculation operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (QHEParticlePrecalculationOperation* operation)
{
  int Step = operation->GetHilbertSpaceDimension() / this->NbrProcesses;
  int FirstComponent = 0;
  int ReducedNbrProcesses = this->NbrProcesses - 1;
  for (int i = 0; i < ReducedNbrProcesses; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((QHEParticlePrecalculationOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrProcesses].Operation = operation->Clone();
  ((QHEParticlePrecalculationOperation*) (this->ThreadParameters[ReducedNbrProcesses].Operation))->
    SetIndicesRange(FirstComponent, operation->GetHilbertSpaceDimension() - FirstComponent);  
  this->SendJobs();
  return true;
}
    
// send jobs to threads
//

void SMPArchitecture::SendJobs ()
{
  int Flag = 0;
  void* ret;
#ifdef __SMP__
  pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#endif
  int ReducedNbrProcesses =  this->NbrProcesses - 1;
  for (int i = 0; i < ReducedNbrProcesses; ++i)
    {
      this->ThreadParameters[i].ThreadID = i;
      this->ThreadParameters[i].Flag = &Flag;
#ifdef __SMP__
      this->ThreadParameters[i].mut = &mut;
#endif
    }
  this->ThreadParameters[ReducedNbrProcesses].ThreadID = ReducedNbrProcesses;
  this->ThreadParameters[ReducedNbrProcesses].Flag = &Flag;
#ifdef __SMP__
  this->ThreadParameters[ReducedNbrProcesses].mut = &mut;
#endif

#ifdef __SMP__
  pthread_t* Threads2 = new pthread_t [this->NbrProcesses];
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      if (pthread_create (&(Threads2[i]), 0, ThreadExecuteOperation, (void*) &(this->ThreadParameters[i])) )
	{
	  cout << "error, cannot create thread" << endl;
	  exit(1);
	}
    }
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      (void) pthread_join (Threads2[i], &ret);
    }
#endif

/*  for (int i = 0; i < this->NbrProcesses; ++i)
    {
//      pthread_mutex_lock(this->ThreadParameters[i].LocalMutex);
      this->ThreadParameters[i].ThreadState = SMPArchitecture::Execute;
//      pthread_mutex_unlock(this->ThreadParameters[i].LocalMutex);
      pthread_cond_signal(this->ThreadParameters[i].LocalCondition);       
    }
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      while (this->ThreadParameters[i].ThreadState != SMPArchitecture::Accomplished)
	pthread_cond_wait (this->ThreadParameters[i].LocalCondition, this->ThreadParameters[i].LocalMutex);      
//      pthread_mutex_lock(this->ThreadParameters[i].LocalMutex);
      this->ThreadParameters[i].ThreadState = SMPArchitecture::Wait;
//      cout << "stop waiting thread " << i << endl;
//      pthread_mutex_unlock(this->ThreadParameters[i].LocalMutex);
    }*/
}
  
// function used by a thread for multiplication
//
// param = pointer to additional parameters, has to be cast into ThreadMultiplyParameter pointer
// return value = unused pointer (null)

void* ThreadMultiply(void* param)
{
#ifdef __SMP__
  ThreadMultiplyParameter* TmpParam = (ThreadMultiplyParameter*) param;
  TmpParam->Hamiltonian->Multiply((*(TmpParam->SourceVector)), (*(TmpParam->DestinationVector)), TmpParam->FirstComponent, 
				  TmpParam->NbrComponent);
  pthread_mutex_lock(TmpParam->mut);
  (*(TmpParam->Flag)) = TmpParam->ThreadID;
  pthread_mutex_unlock(TmpParam->mut);
#endif
  return 0;
}

// thread for multiplicationmain function
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)

void* ThreadExecuteOperation(void* param)
{
#ifdef __SMP__
  ThreadMainParameter* LocalThreadParamater = (ThreadMainParameter*) param;
  LocalThreadParamater->Operation->ApplyOperation();
  pthread_mutex_lock(LocalThreadParamater->mut);
  (*(LocalThreadParamater->Flag)) = LocalThreadParamater->ThreadID;
  pthread_mutex_unlock(LocalThreadParamater->mut);
#endif
  return 0;
}

// main function for thread
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)

void* ThreadMain(void* param)
{
  ThreadMainParameter* LocalThreadParamater = (ThreadMainParameter*) param;
#ifdef __SMP__
  pthread_mutex_t LaunchJobMutex = PTHREAD_MUTEX_INITIALIZER;
  pthread_cond_t LaunchJobCondition = PTHREAD_COND_INITIALIZER;
  LocalThreadParamater->LocalMutex = &LaunchJobMutex;
  LocalThreadParamater->LocalCondition = &LaunchJobCondition;
  while (LocalThreadParamater->ThreadState != SMPArchitecture::Exit)
    {
      if (LocalThreadParamater->ThreadState == SMPArchitecture::Execute)
	{
//	  cout << "thread " << LocalThreadParamater->ThreadID << ": receiving new job " << LocalThreadParamater->Operation->GetOperationType() << endl;
	  LocalThreadParamater->Operation->ApplyOperation();
	  LocalThreadParamater->ThreadState = SMPArchitecture::Accomplished;
//	  cout << "thread " << LocalThreadParamater->ThreadID << ": job accomplished" << endl;	    
	}
      pthread_cond_signal(LocalThreadParamater->LocalCondition);
      pthread_mutex_lock(LocalThreadParamater->LocalMutex);
      while ((LocalThreadParamater->ThreadState == SMPArchitecture::Wait) || (LocalThreadParamater->ThreadState == SMPArchitecture::Accomplished))
	pthread_cond_wait (LocalThreadParamater->LocalCondition, LocalThreadParamater->LocalMutex);      
      pthread_mutex_unlock(LocalThreadParamater->LocalMutex);
    }  
#endif
//  cout << "exit thread " << LocalThreadParamater->ThreadID << endl;
  LocalThreadParamater->ThreadState = SMPArchitecture::Dead;
  return 0;
}
