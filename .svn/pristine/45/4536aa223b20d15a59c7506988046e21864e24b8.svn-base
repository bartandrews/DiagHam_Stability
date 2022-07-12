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
  // pointer to the global mutex used for the rendering queue (common to all threads)
  pthread_mutex_t* GlobalMutex;
  // pointer to the global condition used for the rendering queue (common to all threads)
  pthread_cond_t* GlobalCondition;

  // pointer to the local mutex used for the rendering queue
  pthread_mutex_t* LocalMutex;
  // pointer to the local condition used for the rendering queue
  pthread_cond_t* LocalCondition;

  pthread_mutex_t* mut;
#endif
  
  int ThreadState;

  AbstractArchitectureOperation* Operation;
  
};



class SMPArchitecture : public AbstractArchitecture
{

 private:

  // number of processes to run simultaneously (in principle, the number of processors that can be allocated)
  int NbrProcesses;

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
  // nbrProcesses = number of processes to run simultaneously (in principle, the number of processors that can be allocated)
  SMPArchitecture(int nbrProcesses);
  
  // destructor
  //
  ~SMPArchitecture();
  
  // multiply a vector by an hamiltonian and store the result in another vector
  //
  // hamiltonian = pointer to the hamiltonian to use
  // vSource = vector to multiply 
  // vDestination = vector where result has to be stored 
  void Multiply (AbstractHamiltonian* hamiltonian, Vector& vSource, Vector& vDestination);

  // execute an architecture-dependent vector hamiltonian multiplication operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  bool ExecuteOperation (VectorHamiltonianMultiplyOperation* operation);
  
  // main function for thread
  //
  // param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
  // return value = unused pointer (null)
  friend void* ThreadMain(void* param);

  // execute an architecture-dependent add real linear combination operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  bool ExecuteOperation (AddRealLinearCombinationOperation* operation);

  // execute an architecture-dependent multiple real scalar product operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  bool ExecuteOperation (MultipleRealScalarProductOperation* operation);
  
  // execute an architecture-dependent matrix matrix multiplication operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  bool ExecuteOperation (MatrixMatrixMultiplyOperation* operation);
    
  // execute an architecture-dependent QHE particle hamiltonian precalculation operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  bool ExecuteOperation (QHEParticlePrecalculationOperation* operation);

 protected:
  
  // send jobs to threads
  //
  void SendJobs ();
  
};

#endif
