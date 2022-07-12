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
#include "GeneralTools/StringTools.h"

#ifdef __SMP__
#include <pthread.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>


using std::ofstream;
using std::ios;
using std::cout;
using std::endl;



// thread for multiplicationmain function
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)
void* ThreadExecuteOperation(void* param);

// thread for operation using round robin scheduling
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)
void* ThreadExecuteOperationRoundRobin(void* param);


// constructor
//
// nbrThreads = number of threads to run simultaneously (in principle, the number of processors that can be allocated)
// logFile = name of the optional log file to allow code profiling on SMP architecture

SMPArchitecture::SMPArchitecture(int nbrThreads, char* logFile)
{
  this->ArchitectureID = AbstractArchitecture::SMP;
  this->NbrThreads = nbrThreads;
  this->HilbertSpaceDimension = 0;
  if (logFile != 0)
    {
      this->VerboseModeFlag = true;
      this->LogFile = new char [strlen(logFile) + 8];
      strcpy (this->LogFile, logFile);
      ofstream File;
      File.open(this->LogFile, ios::out);
      if (!File.is_open())
	{
	  cout << "ERROR : cannot write log file " << this->LogFile << endl;
	  this->VerboseModeFlag = false;
	}
      File.close();
    }
  else
    {
      this->VerboseModeFlag = false;
      this->LogFile = 0;
    }

#ifdef __SMP__
  this->ThreadParameters = new ThreadMainParameter [this->NbrThreads];
  this->Threads = new pthread_t [this->NbrThreads];
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      this->ThreadParameters[i].ThreadState = SMPArchitecture::Wait;
      this->ThreadParameters[i].ThreadID = i;
      this->ThreadParameters[i].Architecture = this;
    }
#endif
}
  
// destructor
//

SMPArchitecture::~SMPArchitecture()
{
#ifdef __SMP__
  delete[] this->Threads;
  delete[] this->ThreadParameters;
#endif
  if (this->LogFile != 0)
    delete[] this->LogFile;
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
// nbrJobs = optional parameter that indicates how many jobs have to be sent (0 if all)

void SMPArchitecture::SendJobs (int nbrJobs)
{
  int Flag = 0;
  void* ret;
#ifdef __SMP__
  pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#endif
  if ((nbrJobs <= 0) || (nbrJobs > this->NbrThreads))
    {
      nbrJobs = this->NbrThreads;
    }
  int ReducedNbrThreads =  nbrJobs - 1;
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
  pthread_t* Threads2 = new pthread_t [nbrJobs];
  for (int i = 0; i < nbrJobs; ++i)
    {
      int code;
      if ( (code = pthread_create (&(Threads2[i]), (const pthread_attr_t *)NULL, ThreadExecuteOperation, (void*) &(this->ThreadParameters[i])))!=0 )
	{
	  cout << "error, cannot create thread" << endl;
	  cout << "pthread_create exit code: "<<code<<endl;
	  exit(1);
	}
    }
  for (int i = 0; i < nbrJobs; ++i)
    {
      (void) pthread_join (Threads2[i], &ret);
    }
  delete[] Threads2;
#endif

}
  
 
// send jobs to threads
//
// nbrJobs = optional parameter that indicates how many jobs have to be sent (0 if all)

void SMPArchitecture::SendJobsRoundRobin (int nbrJobs)
{
  int Flag = 0;
  void* ret;
#ifdef __SMP__
  pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#endif
  if ((nbrJobs <= 0) || (nbrJobs > this->NbrThreads))
    {
      nbrJobs = this->NbrThreads;
    }
  int ReducedNbrThreads =  nbrJobs - 1;
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
  pthread_t* Threads2 = new pthread_t [nbrJobs];
  for (int i = 0; i < nbrJobs; ++i)
    {
      int code;
      if ( (code = pthread_create (&(Threads2[i]), (const pthread_attr_t *)NULL, ThreadExecuteOperationRoundRobin, (void*) &(this->ThreadParameters[i])))!=0 )
	{
	  cout << "error, cannot create thread" << endl;
	  cout << "pthread_create exit code: "<<code<<endl;
	  exit(1);
	}
    }
  for (int i = 0; i < nbrJobs; ++i)
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
  pthread_exit(0);
#endif
  return 0;
}

// thread for operation using round robin scheduling
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)

void* ThreadExecuteOperationRoundRobin(void* param)
{
#ifdef __SMP__
  ThreadMainParameter* LocalThreadParamater = (ThreadMainParameter*) param;
  LocalThreadParamater->Operation->ApplyOperationSMPRoundRobin((SMPArchitecture*)LocalThreadParamater->Architecture, LocalThreadParamater->ThreadID);
  pthread_mutex_lock(LocalThreadParamater->mut);
  (*(LocalThreadParamater->Flag)) = LocalThreadParamater->ThreadID;
  pthread_mutex_unlock(LocalThreadParamater->mut);
#endif
  return 0;
}

// add an entry to the log file
//
// message = string corresponding to entry to add to the log file
// masterFlag = true if only the master node should add the entry
// return value = true if no error occured

bool SMPArchitecture::AddToLog(const char * message, bool masterFlag)
{
  if ((this->VerboseModeFlag == true) && (this->LogFile != 0))
    {
      ofstream File;
      File.open(this->LogFile, ios::out | ios::app);
      if (!File.is_open())
	{
	  cout << "ERROR : cannot write log file " << this->LogFile << endl;
	  return false;
	}
      File << message << endl;
      File.close();
      return true;
    }
  return false;
}

// dump the log file into a string
//
// header = optional header to add before the log file
// footer = optional footer to add at the end of the log file
// return value = string or 0 if an error occured or log is not available

char* SMPArchitecture::DumpLog(const char* header, const char* footer)
{
  if ((this->VerboseModeFlag == false) || (this->LogFile == 0))
    return 0;
  return DumpTextFile(this->LogFile, header, footer);
}
