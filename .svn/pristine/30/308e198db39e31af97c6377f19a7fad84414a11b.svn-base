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

 public:
  
  // constructor
  //
  // nbrThreads = number of threads to run simultaneously (in principle, the number of processors that can be allocated)
  SMPArchitecture(int nbrThreads);
  
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
  void SendJobs ();
  
};

// get the  number of threads that run simultaneously
//
// return value = number of threads

inline int SMPArchitecture::GetNbrThreads()
{
  return this->NbrThreads;
}

#endif
