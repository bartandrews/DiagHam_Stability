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


#ifndef SMPARCHITECTURE_H
#define SMPARCHITECTURE_H


#include "config.h"
#include "Architecture/AbstractArchitecture.h"

#ifdef __SMP__
#include <pthread.h>
#endif


class AbstractArchitectureOperation;


struct ThreadMainParameter
{
  int ThreadID;
  int* Flag;
#ifdef __SMP__
  // pointer to the mutex condition
  pthread_mutex_t* mut; 
#endif
  
  int ThreadState;

  AbstractArchitectureOperation* Operation;
  AbstractArchitecture* Architecture;
  
};



class SMPArchitecture : public AbstractArchitecture
{

 private:

  // number of threads to run simultaneously (in principle, the number of processors that can be allocated)
  int NbrThreads;

  //  parameters for each thread
  ThreadMainParameter* ThreadParameters;

  // thread handles
#ifdef __SMP__
  pthread_t* Threads;
#endif
  
  // code for each thread state
  enum ThreadStateCode
    {
      Wait = 0x0,
      Exit = 0x1,
      Dead = 0x2,
      Execute = 0x4,
      Accomplished = 0x8
    };

  // name of the optional log file to allow code profiling on SMP architecture
  char* LogFile;
  // flag to indicate if the log file option is activated
  bool VerboseModeFlag;

public:
  
  // constructor
  //
  // nbrThreads = number of threads to run simultaneously (in principle, the number of processors that can be allocated)
  // logFile = name of the optional log file to allow code profiling on SMP architecture
  SMPArchitecture(int nbrThreads, char* logFile = 0);
  
  // destructor
  //
  ~SMPArchitecture();
  
  // get the  number of threads that run simultaneously
  //
  // return value = number of threads
  int GetNbrThreads();

  // set the operation that has to be executed by a given thread
  //
  // operation = pointer to the operation
  // index = thread index
  void SetThreadOperation(AbstractArchitectureOperation* operation, int index);

  // send jobs to threads
  //
  // nbrJobs = optional parameter that indicates how many jobs have to be sent (0 if all)
  void SendJobs (int nbrJobs = 0);
  
  // send jobs to threads
  //
  // nbrJobs = optional parameter that indicates how many jobs have to be sent (0 if all)
  void SendJobsRoundRobin (int nbrJobs = 0);

  // indicate if the log file option is activated
  //
  // return value = true if the option is activated
  bool VerboseMode();

  // add an entry to the log file
  //
  // message = string corresponding to entry to add to the log file
  // masterFlag = true if only the master node should add the entry
  // return value = true if no error occured
  bool AddToLog(const char * message, bool masterFlag = false);

  // dump the log file into a string
  //
  // header = optional header to add before the log file
  // footer = optional footer to add at the end of the log file
  // return value = string or 0 if an error occured or log is not available
  char* DumpLog(const char* header = 0, const char* footer = 0);

  // attempt to lock the mutex and block until it is available
  //
  void LockMutex( );
  
  // Unlock mutex
  //
  void UnLockMutex( );  
  
  // get the id of the current threads
  //
  // return = thread id
  int GetThreadID();
  
};

// get the  number of threads that run simultaneously
//
// return value = number of threads

inline int SMPArchitecture::GetNbrThreads()
{
  return this->NbrThreads;
}

// indicate if the log file option is activated
//
// return value = true if the option is activated

inline bool SMPArchitecture::VerboseMode()
{
  return this->VerboseModeFlag;
}

// attempt to lock the mutex and block until it is available
//
  
inline void SMPArchitecture::LockMutex()
{
  pthread_mutex_lock(this->ThreadParameters->mut);
}
 
// Unlock mutex
// 
 
inline void SMPArchitecture::UnLockMutex()
{
  pthread_mutex_unlock(this->ThreadParameters->mut);
}

// get the id of the current threads
//
// return = thread id

inline int SMPArchitecture::GetThreadID()
{
  return this->ThreadParameters->ThreadID;
}


#endif
