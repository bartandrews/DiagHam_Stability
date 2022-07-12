////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of server for ditributed architecture                //
//                                                                            //
//                        last modification : 09/04/2003                      //
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


#ifndef SERVERDISTRIBUTEDARCHITECTURE_H
#define SERVERDISTRIBUTEDARCHITECTURE_H


#include "config.h"
#include "Architecture/AbstractArchitecture.h"

#ifdef __SMP__
#include <pthread.h>
#endif


class AbstractArchitectureOperation;



class ServerDistributedArchitecture : public AbstractArchitecture
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
  
  // code describing each communication chunk
  enum CommunicationCode
    {
      IsAlive = 0x00000001,
      ImAlive = 0x00000002,
      Exit = 0x1,
      Dead = 0x2,
      Execute = 0x4,
      Accomplished = 0x8
    };

 public:
  
  // constructor
  //
  // configurationFileName = pointer to the name of the file containing distributed computation configuration
  // port = port to be opened by the server 
  ServerDistributedArchitecture(char* configurationFileName, int port);
  
  // destructor
  //
  ~ServerDistributedArchitecture();
  
  // main function for thread
  //
  // param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
  // return value = unused pointer (null)
  friend void* ThreadMain(void* param);

  // execute an architecture-dependent vector hamiltonian multiplication operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  bool ExecuteOperation (VectorHamiltonianMultiplyOperation* operation);
  
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
