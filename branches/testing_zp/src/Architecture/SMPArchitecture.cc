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
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"

#ifdef __SMP__
#include <pthread.h>
#endif
#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;



// thread for multiplicationmain function
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)
void* ThreadExecuteOperation(void* param);


// constructor
//
// nbrThreads = number of threads to run simultaneously (in principle, the number of processors that can be allocated)

SMPArchitecture::SMPArchitecture(int nbrThreads)
{
  this->ArchitectureID = AbstractArchitecture::SMP;
  this->NbrThreads = nbrThreads;
#ifdef __SMP__
  this->ThreadParameters = new ThreadMainParameter [this->NbrThreads];
  this->Threads = new pthread_t [this->NbrThreads];
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      this->ThreadParameters[i].ThreadState = SMPArchitecture::Wait;
      this->ThreadParameters[i].ThreadID = i;
    }
#endif
}
  
// destructor
//

SMPArchitecture::~SMPArchitecture()
{
#ifdef __SMP__
  delete[] this->Threads;
#endif
}
  
// set the operation that has to be executed by a given thread
//
// operation = pointer to the operation
// index = thread index

void SMPArchitecture::SetThreadOperation(AbstractArchitectureOperation* operation, int index)
{
  this->ThreadParameters[index].Operation = operation;
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
  int ReducedNbrThreads =  this->NbrThreads - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].ThreadID = i;
      this->ThreadParameters[i].Flag = &Flag;
#ifdef __SMP__
      this->ThreadParameters[i].mut = &mut;
#endif
    }
  this->ThreadParameters[ReducedNbrThreads].ThreadID = ReducedNbrThreads;
  this->ThreadParameters[ReducedNbrThreads].Flag = &Flag;
#ifdef __SMP__
  this->ThreadParameters[ReducedNbrThreads].mut = &mut;
#endif

#ifdef __SMP__
  pthread_t* Threads2 = new pthread_t [this->NbrThreads];
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      int code;
      if ( (code = pthread_create (&(Threads2[i]), (const pthread_attr_t *)NULL, ThreadExecuteOperation, (void*) &(this->ThreadParameters[i])))!=0 )
	{
	  cout << "error, cannot create thread" << endl;
	  cout << "pthread_create exit code: "<<code<<endl;
	  exit(1);
	}
    }
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      (void) pthread_join (Threads2[i], &ret);
    }
  delete[] Threads2;
#endif

}
  
// thread for multiplicationmain function
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)

void* ThreadExecuteOperation(void* param)
{
#ifdef __SMP__
  ThreadMainParameter* LocalThreadParamater = (ThreadMainParameter*) param;
  LocalThreadParamater->Operation->RawApplyOperation();
  pthread_mutex_lock(LocalThreadParamater->mut);
  (*(LocalThreadParamater->Flag)) = LocalThreadParamater->ThreadID;
  pthread_mutex_unlock(LocalThreadParamater->mut);
#endif
  return 0;
}

