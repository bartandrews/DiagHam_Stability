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

#ifdef __SMP__
#include <pthread.h>
#endif
#include <stdlib.h>


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
}
  
// destructor
//

SMPArchitecture::~SMPArchitecture()
{
}
  
// execute an architecture-dependent vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (VectorHamiltonianMultiplyOperation* operation)
{
  int Flag = 0;
  void* ret;
#ifdef __SMP__
  pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#endif
  ThreadMainParameter* param = new ThreadMainParameter [this->NbrProcesses];
  int Step = operation->GetDestinationVector()->GetVectorDimension() / this->NbrProcesses;
  int FirstComponent = 0;
  int ReducedNbrProcesses = this->NbrProcesses - 1;
  operation->GetDestinationVector()->ClearVector();
  for (int i = 0; i < ReducedNbrProcesses; ++i)
    {
      param[i].ThreadID = i;
      param[i].Flag = &Flag;
#ifdef __SMP__
      param[i].mut = &mut;
#endif
      param[i].Operation = operation->Clone();
      ((VectorHamiltonianMultiplyOperation*) (param[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  param[ReducedNbrProcesses].ThreadID = ReducedNbrProcesses;
  param[ReducedNbrProcesses].Flag = &Flag;
#ifdef __SMP__
  param[ReducedNbrProcesses].mut = &mut;
#endif
  param[ReducedNbrProcesses].Operation = operation->Clone();
  ((VectorHamiltonianMultiplyOperation*) (param[ReducedNbrProcesses].Operation))->SetIndicesRange(FirstComponent, 
												  operation->GetDestinationVector()->GetVectorDimension() - FirstComponent);  
  for (int i = 1; i < this->NbrProcesses; ++i)
    {
      ((VectorHamiltonianMultiplyOperation*) (param[i].Operation))->SetDestinationVector(operation->GetDestinationVector()->EmptyClone(true));
    }
#ifdef __SMP__
  pthread_t* Threads = new pthread_t [this->NbrProcesses];
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      if (pthread_create (&(Threads[i]), 0, ThreadExecuteOperation, (void*) &(param[i])) )
	{
	  cout << "error, cannot create thread" << endl;
	  exit(1);
	}
    }
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      (void) pthread_join (Threads[i], &ret);
    }
  for (int i = 1; i < this->NbrProcesses; ++i)
    {
      (*(operation->GetDestinationVector())) += (*(((VectorHamiltonianMultiplyOperation*) (param[i].Operation))->GetDestinationVector()));
      delete ((VectorHamiltonianMultiplyOperation*) (param[i].Operation))->GetDestinationVector();
    }
#endif
  delete[] param;
  return true;
}
  
// execute an architecture-dependent add real linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (AddRealLinearCombinationOperation* operation)
{
  int Flag = 0;
  void* ret;
#ifdef __SMP__
  pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#endif
  ThreadMainParameter* param = new ThreadMainParameter [this->NbrProcesses];
  int Step = operation->GetDestinationVector()->GetVectorDimension() / this->NbrProcesses;
  int FirstComponent = 0;
  int ReducedNbrProcesses = this->NbrProcesses - 1;
  operation->GetDestinationVector()->ClearVector();
  for (int i = 0; i < ReducedNbrProcesses; ++i)
    {
      param[i].ThreadID = i;
      param[i].Flag = &Flag;
#ifdef __SMP__
      param[i].mut = &mut;
#endif
      param[i].Operation = operation->Clone();
      ((AddRealLinearCombinationOperation*) (param[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  param[ReducedNbrProcesses].ThreadID = ReducedNbrProcesses;
  param[ReducedNbrProcesses].Flag = &Flag;
#ifdef __SMP__
  param[ReducedNbrProcesses].mut = &mut;
#endif
  param[ReducedNbrProcesses].Operation = operation->Clone();
  ((AddRealLinearCombinationOperation*) (param[ReducedNbrProcesses].Operation))->SetIndicesRange(FirstComponent, 
												 operation->GetDestinationVector()->GetVectorDimension() - FirstComponent);  
#ifdef __SMP__
  pthread_t* Threads = new pthread_t [this->NbrProcesses];
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      if (pthread_create (&(Threads[i]), 0, ThreadExecuteOperation, (void*) &(param[i])) )
	{
	  cout << "error, cannot create thread" << endl;
	  exit(1);
	}
    }
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      (void) pthread_join (Threads[i], &ret);
    }
#endif
  delete[] param;
  return true;
}  

// multiply a vector by an hamiltonian and store the result in another vector
//
// hamiltonian = pointer to the hamiltonian to use
// vSource = vector to multiply 
// vDestination = vector where result has to be stored 

void SMPArchitecture::Multiply (AbstractHamiltonian* hamiltonian, Vector& vSource, Vector& vDestination)
{  
  int Flag = 0;
  void* ret;
#ifdef __SMP__
  pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#endif
  ThreadMultiplyParameter* param = new ThreadMultiplyParameter [this->NbrProcesses];
  int Step = hamiltonian->GetHilbertSpaceDimension() / this->NbrProcesses;
  int FirstComponent = 0;
  int ReducedNbrProcesses = this->NbrProcesses - 1;
  vDestination.ClearVector();
  for (int i = 0; i < ReducedNbrProcesses; ++i)
    {
      param[i].ThreadID = i;
      param[i].Flag = &Flag;
#ifdef __SMP__
      param[i].mut = &mut;
#endif
      param[i].FirstComponent = FirstComponent;
      param[i].NbrComponent = Step;
      FirstComponent += Step;
      param[i].Hamiltonian = hamiltonian;
      param[i].SourceVector = &vSource;
    }
  param[ReducedNbrProcesses].ThreadID = ReducedNbrProcesses;
  param[ReducedNbrProcesses].Flag = &Flag;
#ifdef __SMP__
  param[ReducedNbrProcesses].mut = &mut;
#endif
  param[ReducedNbrProcesses].FirstComponent = FirstComponent;
  param[ReducedNbrProcesses].NbrComponent = hamiltonian->GetHilbertSpaceDimension() - FirstComponent;
  param[ReducedNbrProcesses].Hamiltonian = hamiltonian;
  param[ReducedNbrProcesses].SourceVector = &vSource;  
  param[0].DestinationVector = &vDestination;
  for (int i = 1; i < this->NbrProcesses; ++i)
    {
      param[i].DestinationVector = vDestination.EmptyClone(true);
    }
#ifdef __SMP__
  pthread_t* Threads = new pthread_t [this->NbrProcesses];
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      if (pthread_create (&(Threads[i]), 0, ThreadMultiply, (void*) &(param[i])) )
	{
	  cout << "error, cannot create thread" << endl;
	  exit(1);
	}
    }
  for (int i = 0; i < this->NbrProcesses; ++i)
    {
      (void) pthread_join (Threads[i], &ret);
    }
  for (int i = 1; i < this->NbrProcesses; ++i)
    {
      vDestination += (*(param[i].DestinationVector));
      delete param[i].DestinationVector;
    }
#endif
  delete[] param;
  return;
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
  
#endif
  LocalThreadParamater->ThreadState = SMPArchitecture::Dead;
  return 0;
}
