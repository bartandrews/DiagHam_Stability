////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of simple MPI Architecture                     //
//                                                                            //
//                        last modification : 17/05/2004                      //
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
#include "Architecture/SimpleMPIArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/SortedRealUniqueArray.h"
#include "GeneralTools/SortedComplexUniqueArray.h"


#include <sys/time.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstring>
#ifdef __MPI__
#include <mpi.h>
#endif


using std::ofstream;
using std::ios;
using std::cout;
using std::endl;


// default constructor
//

SimpleMPIArchitecture::SimpleMPIArchitecture()
{
  this->MinimumIndices = 0;
  this->MaximumIndices = 0;
  this->OldMinimumIndices = 0;
  this->OldMaximumIndices = 0;
  this->AutomaticLoadBalancing = false;
}

// constructor
//
// logFile = name of the optional log file to allow code profiling on MPI architecture
// automaticLoadBalancing = flag that indicates if automatic load balancing have to be done, overriding any manual load balancing

SimpleMPIArchitecture::SimpleMPIArchitecture(char* logFile, bool automaticLoadBalancing)
{
  this->PerformanceIndex = 1.0;
  this->ArchitectureID = AbstractArchitecture::SimpleMPI;
  this->AutomaticLoadBalancing = automaticLoadBalancing;
#ifdef __MPI__
  MPI::Init();
  this->NbrMPINodes = MPI::COMM_WORLD.Get_size();
  this->MPIRank = MPI::COMM_WORLD.Get_rank();
  this->ClusterPerformanceArray = new double [this->NbrMPINodes];
  this->MinimumIndices = 0;
  this->MaximumIndices = 0;
  this->OldMinimumIndices = 0;
  this->OldMaximumIndices = 0;
  if (this->MPIRank != 0)
    {
      this->MasterNodeFlag = false;
      if (logFile != 0)
	this->VerboseModeFlag = true;
      else
	this->VerboseModeFlag = false;
      this->LogFile = 0;
      MPI::COMM_WORLD.Send(&this->PerformanceIndex, 1, MPI::DOUBLE, 0, 1);
    }
  else
    {
      this->MasterNodeFlag = true;
      this->TotalPerformanceIndex = this->PerformanceIndex;
      this->ClusterPerformanceArray[0] = this->PerformanceIndex;
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  MPI::COMM_WORLD.Recv(&this->ClusterPerformanceArray[i], 1, MPI::DOUBLE, i, 1);	  
	  this->TotalPerformanceIndex += this->ClusterPerformanceArray[i];
	}
      for (int i = 0; i < this->NbrMPINodes; ++i)
	{
	  this->ClusterPerformanceArray[i] /= this->TotalPerformanceIndex;
	}
      if (logFile != 0)
	{
	  this->LogFile = new char [strlen(logFile) + 1];
	  strcpy (this->LogFile, logFile);
	  ofstream File;
	  File.open(this->LogFile, ios::out);
	  if (!File.is_open())
	    {
	      cout << "ERROR : cannot write log file " << this->LogFile << endl;
	      this->VerboseModeFlag = false;
	    }
	  else
	    this->VerboseModeFlag = true;
	  File.close();
	}
      else
	{
	  this->VerboseModeFlag = false;
	  this->LogFile = 0;
	}
    }
  MPI::COMM_WORLD.Bcast(this->ClusterPerformanceArray, this->NbrMPINodes, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&this->TotalPerformanceIndex, 1, MPI::DOUBLE, 0);
#else
  this->MasterNodeFlag = true;
  this->NbrMPINodes = 1;
  this->MPIRank = 0;
  this->ClusterPerformanceArray = 0;
  this->TotalPerformanceIndex = this->PerformanceIndex;
#endif
  if (this->MasterNodeFlag == true)
    {
      cout << this->NbrMPINodes << " " << this->TotalPerformanceIndex << endl;
    }
  this->LocalArchitecture = new MonoProcessorArchitecture;
}
  
// destructor
//

SimpleMPIArchitecture::~SimpleMPIArchitecture()
{
#ifdef __MPI__
  MPI::Finalize();
#endif
  if (this->MinimumIndices != 0)
    {
      delete [] this->MinimumIndices;
      delete [] this->MaximumIndices;
    }
  if (this->OldMinimumIndices != 0)
    {
      delete [] this->OldMinimumIndices;
      delete [] this->OldMaximumIndices;
    }
  delete this->LocalArchitecture;
  if (this->ClusterPerformanceArray != 0)
    delete[] this->ClusterPerformanceArray;
  if (this->LogFile != 0)
    delete[] this->LogFile;
}
  
// get typical range of indices on which the local architecture acts
//
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)

void SimpleMPIArchitecture::GetTypicalRange (long& minIndex, long& maxIndex)
{
  minIndex = this->MinimumIndex;
  maxIndex = this->MaximumIndex;
}
  
// get typical range of indices on which the local architecture acts, providing the number of calculations that have to be performed per index
//
// mbrOperationPerIndex = reference on the number of calculations per index. If the return value is true, a new array will be allocated
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)
// return value = true if the range has been optimized

bool SimpleMPIArchitecture::GetOptimizedTypicalRange (int*& nbrOperationPerIndex, long& minIndex, long& maxIndex)
{
  if (this->AutomaticLoadBalancing == false)
    {
      this->GetTypicalRange(minIndex, maxIndex);
      return false;
    }
  if (this->MasterNodeFlag == true)
    { 
      cout << "Optimizing typical range"<<endl;
      cout << "Before: ";
      this->PrintLoadBalancing() << endl;
      int* TmpNbrOperationPerIndex = new int [this->MaximumIndices[this->NbrMPINodes - 1] - this->MinimumIndices[0] + 1l];
      long TmpIndex = 0;
      long EffectiveDimension = this->MaximumIndices[0] - this->MinimumIndices[0] + 1l;
      long TotalEffectiveDimension = EffectiveDimension;
      for (long i = 0l; i < EffectiveDimension; ++i)
	{
	  TmpNbrOperationPerIndex[this->MinimumIndices[0] + i] = nbrOperationPerIndex[i];
	}  
      TmpIndex += EffectiveDimension;      
      for (int i = 1; i < NbrMPINodes; ++i)
	{
	  long TmpNbrElements = this->MaximumIndices[i] - this->MinimumIndices[i] + 1l;
	  this->ReceiveFromSlave(i - 1, TmpNbrOperationPerIndex + this->MinimumIndices[i], TmpNbrElements);
	  TotalEffectiveDimension += TmpNbrElements;
	}
      long TotalNbrOperations = 0l;
      for (long i = 0l; i < TotalEffectiveDimension; ++i)
	{
	  TotalNbrOperations += (long) TmpNbrOperationPerIndex[i];
	}
      // write load profile to disk      
      ofstream Store("load-profile.dat", std::ios::binary | std::ios::out);
      long TotalDimension = this->MaximumIndices[this->NbrMPINodes - 1] - this->MinimumIndices[0] + 1l;
      Store.write((char*) &(TotalDimension), sizeof(long));
      Store.write((char*) &(TotalNbrOperations), sizeof(long));
      Store.write((char*) TmpNbrOperationPerIndex, (TotalDimension)*sizeof(int));
      long CheckMark=(long)(1e12*M_PI);
      Store.write((char*) &(CheckMark), sizeof(long));
      // done writing load profile

      long *TmpMinimumIndices, *TmpMaximumIndices;
      this->DeduceLoadDistribution(TotalNbrOperations, TmpNbrOperationPerIndex, TmpMinimumIndices, TmpMaximumIndices);

      if (this->OldMinimumIndices != NULL)
	{
	  delete[] this->OldMinimumIndices;
	  delete[] this->OldMaximumIndices;
	}
      this->OldMinimumIndices = this->MinimumIndices;
      this->OldMaximumIndices = this->MaximumIndices;
      this->MinimumIndices = TmpMinimumIndices;
      this->MaximumIndices = TmpMaximumIndices;
#ifdef __MPI__      
      MPI::COMM_WORLD.Bcast(this->MinimumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
      MPI::COMM_WORLD.Bcast(this->MaximumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
#endif
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  this->SendToSlaves(i - 1, (int*) (TmpNbrOperationPerIndex + this->MinimumIndices[i]), 
			     (int) (this->MaximumIndices[i] - this->MinimumIndices[i] + 1l));      
	}
      delete[] nbrOperationPerIndex;
      nbrOperationPerIndex = new int [this->MaximumIndices[0] - this->MinimumIndices[0] + 1l];      
      for (long i = this->MinimumIndices[0]; i <= this->MaximumIndices[0]; ++i)
	nbrOperationPerIndex[i - this->MinimumIndices[0]] = TmpNbrOperationPerIndex[i];
      cout << "After: ";
      this->PrintLoadBalancing() << endl;

      delete[] TmpNbrOperationPerIndex;
    }
  else
    {
      this->SendToMaster(nbrOperationPerIndex, this->MaximumIndices[this->MPIRank] - this->MinimumIndices[this->MPIRank] + 1l);
      int TmpNbrValues = 0;
      delete[] nbrOperationPerIndex;
      if (this->OldMinimumIndices == NULL)
	{
	  this->OldMinimumIndices = new long[this->NbrMPINodes];
	  this->OldMaximumIndices = new long[this->NbrMPINodes];;
	}
      for (int i=0; i<this->NbrMPINodes; ++i)
	{
	  this->OldMinimumIndices[i] = this->MinimumIndices[i];
	  this->OldMaximumIndices[i] = this->MaximumIndices[i];
	}
#ifdef __MPI__      
      MPI::COMM_WORLD.Bcast(this->MinimumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
      MPI::COMM_WORLD.Bcast(this->MaximumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
#endif
      int TmpNbrElements = (int) (this->MaximumIndices[this->MPIRank] - this->MinimumIndices[this->MPIRank] + 1l);
      nbrOperationPerIndex = new int [TmpNbrElements];
      this->ReceiveFromMaster(nbrOperationPerIndex, TmpNbrElements);      
    }
  this->MinimumIndex = this->MinimumIndices[this->MPIRank];
  this->MaximumIndex = this->MaximumIndices[this->MPIRank];
  minIndex = this->MinimumIndex;
  maxIndex = this->MaximumIndex;
  return true;
}
  
// get typical range of indices on which the local architecture acts, providing the number of calculations that have to be performed per index
//
// mbrOperationPerIndex = reference on the number of calculations per index. If the return value is true, a new array will be allocated
// memoryPerOperation = memory required per operation (in bytes)
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)
// return value = true if the range has been optimized

bool SimpleMPIArchitecture::GetOptimizedTypicalRange (int*& nbrOperationPerIndex, int memoryPerOperation, long& minIndex, long& maxIndex)
{
  cout << "Warning, SimpleMPIArchitecture::GetOptimizedTypicalRange with memory calculation is not fully supported" << endl;
  return this->GetOptimizedTypicalRange (nbrOperationPerIndex, minIndex, maxIndex);
}

// get typical range of indices on which the local architecture acts, providing the number of calculations that have to be performed per index
// and merge lists of unique matrix elements on different nodes.
//
// mbrOperationPerIndex = reference on the number of calculations per index. If the return value is true, a new array will be allocated
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)
// realEntries = array of real interaction coefficients
// complexEntries = array of complex interaction coefficients
// return value = true if the range has been optimized
bool SimpleMPIArchitecture::GetOptimizedTypicalRange (int*& nbrOperationPerIndex, long& minIndex, long& maxIndex, 
				       SortedRealUniqueArray &realEntries, SortedComplexUniqueArray &complexEntries)
{
  bool Flag = true;
  if (this->AutomaticLoadBalancing == true)
    {
#ifdef __MPI__      
      Flag &= realEntries.MergeAcrossNodes(MPI::COMM_WORLD);
      Flag &= complexEntries.MergeAcrossNodes(MPI::COMM_WORLD);
#endif
      if (Flag) 
	{
	  if (this->MasterNodeFlag)
	    {
	      realEntries.WriteArray("real-entries.dat");
	      complexEntries.WriteArray("complex-entries.dat");
	    }
	  return this->GetOptimizedTypicalRange (nbrOperationPerIndex, minIndex, maxIndex);
	}
      else
	this->GetTypicalRange(minIndex, maxIndex);
      return false;
    }
  else 
    return this->GetOptimizedTypicalRange (nbrOperationPerIndex, minIndex, maxIndex);
}

// load a typical range of indices and the corresponding operations from a previous run of the calculation
//
// nbrOperationPerIndex = reference on the number of calculations per index. If the return value is true, a new array will be allocated
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)
// return value = true if the range has been optimized
bool SimpleMPIArchitecture::LoadOptimizedTypicalRange (int*& nbrOperationPerIndex, long& minIndex, long& maxIndex, const char* filename)
{
  if (this->MasterNodeFlag == true)
    {
      ifstream Store;
      Store.open(filename, std::ios::binary | std::ios::in);
      int Acknowledge=0;
      if (!Store.is_open())
	{
	  this->BroadcastToSlaves(Acknowledge);
	  return false; // no load pattern found.
	}
      
      // read load balancing information from disk     
      long TotalDimension;
      Store.read((char*) &(TotalDimension), sizeof(long));
      
      if (TotalDimension != this->MaximumIndices[this->NbrMPINodes - 1] - this->MinimumIndices[0] + 1l)
	{
	  cout << "Attention, load profile does not match current dimension - ignoring data"<<endl;
	  this->BroadcastToSlaves(Acknowledge);
	  return false;
	}
      long TotalNbrOperations;
      Store.read((char*) &(TotalNbrOperations), sizeof(long));
      nbrOperationPerIndex = new int[TotalDimension];
      Store.read((char*) (nbrOperationPerIndex), (TotalDimension)*sizeof(int));
      long CheckMark;
      Store.read((char*) &(CheckMark), sizeof(long));
      if (CheckMark!=(long)(1e12*M_PI))
	{
	  cout << "Error reading load profile"<<endl;
	  this->BroadcastToSlaves(Acknowledge);
	  return false;
	}
      long Sum=0;
      for (long i=0; i<TotalDimension; ++i)
	Sum += nbrOperationPerIndex[i];
      if (TotalNbrOperations!=Sum)
	{
	  cout << "Error recovering TotalNbrOperations from file "<<filename<<endl;
	  this->BroadcastToSlaves(Acknowledge);
	  return false;
	}
      Store.close();
      // done reading - signal success
      Acknowledge=1;
      this->BroadcastToSlaves(Acknowledge);

      long *TmpMinimumIndices, *TmpMaximumIndices;
      this->DeduceLoadDistribution(TotalNbrOperations, nbrOperationPerIndex, TmpMinimumIndices, TmpMaximumIndices);

      if (this->OldMinimumIndices != NULL)
	{
	  delete[] this->OldMinimumIndices;
	  delete[] this->OldMaximumIndices;
	}
      this->OldMinimumIndices = this->MinimumIndices;
      this->OldMaximumIndices = this->MaximumIndices;
      this->MinimumIndices = TmpMinimumIndices;
      this->MaximumIndices = TmpMaximumIndices;
#ifdef __MPI__      
      MPI::COMM_WORLD.Bcast(this->MinimumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
      MPI::COMM_WORLD.Bcast(this->MaximumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
#endif
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  this->SendToSlaves(i - 1, (int*) (nbrOperationPerIndex + this->MinimumIndices[i]), 
			     (int) (this->MaximumIndices[i] - this->MinimumIndices[i] + 1l));      
	}
      int *TmpNbrOperationPerIndex = new int [this->MaximumIndices[0] - this->MinimumIndices[0] + 1l];      
      for (long i = this->MinimumIndices[0]; i <= this->MaximumIndices[0]; ++i)
	TmpNbrOperationPerIndex[i - this->MinimumIndices[0]] = nbrOperationPerIndex[i];
      cout << "Deduced load-balancing from saved profile: ";
      delete[] nbrOperationPerIndex;
      nbrOperationPerIndex = TmpNbrOperationPerIndex;
      this->PrintLoadBalancing() << endl;
    }
  else
    {
      int Acknowledge=0;
      this->BroadcastToSlaves(Acknowledge);
      if (Acknowledge==0)
	return false; // data could not be loaded successfully.
      
      if (this->OldMinimumIndices == NULL)
	{
	  this->OldMinimumIndices = new long[this->NbrMPINodes];
	  this->OldMaximumIndices = new long[this->NbrMPINodes];;
	}
      for (int i=0; i<this->NbrMPINodes; ++i)
	{
	  this->OldMinimumIndices[i] = this->MinimumIndices[i];
	  this->OldMaximumIndices[i] = this->MaximumIndices[i];
	}
#ifdef __MPI__      
      MPI::COMM_WORLD.Bcast(this->MinimumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
      MPI::COMM_WORLD.Bcast(this->MaximumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
#endif
      int TmpNbrElements = (int) (this->MaximumIndices[this->MPIRank] - this->MinimumIndices[this->MPIRank] + 1l);
      nbrOperationPerIndex = new int [TmpNbrElements];
      this->ReceiveFromMaster(nbrOperationPerIndex, TmpNbrElements);      
    }
  this->MinimumIndex = this->MinimumIndices[this->MPIRank];
  this->MaximumIndex = this->MaximumIndices[this->MPIRank];
  minIndex = this->MinimumIndex;
  maxIndex = this->MaximumIndex;
  return true;
}



// balance an array differently across nodes (currently supporting T=int, long)
//
// array = reference on a local array holding one entry per local state prior to last call to GetOptimizedTypicalRange; on return - local array holding data for the new range of the MPIArchitecture
// filename = name of a file to which the data should be saved (on master node)
// return value = true if the array has been rebalanced
bool SimpleMPIArchitecture::RebalanceArray (int*& array, const char* filename) 
{
  return RebalanceArrayImplementation(array, filename);
}

bool SimpleMPIArchitecture::RebalanceArray (long*& array, const char* filename) {
  return RebalanceArrayImplementation(array, filename);
}


// Load an array from disk
bool SimpleMPIArchitecture::LoadArray (const char* filename, int*& array)
{
  return LoadArrayImplementation(array, filename);
}

bool SimpleMPIArchitecture::LoadArray (const char* filename, long*& array)
{
  return LoadArrayImplementation(array, filename);
}


// get the ID of the node that handles a given index
//
// index = index to check
// return value = corresponding node ID

int SimpleMPIArchitecture::GetNodeIDFromIndex(long index)
{
  int TmpID = 0;
  while (index > this->MaximumIndices[TmpID]) 
    ++TmpID;
  return TmpID;
}

// set dimension of the Hilbert space on which the architecture has to work
// 
// dimension = dimension of the Hilbert space

void SimpleMPIArchitecture::SetDimension (long dimension)
{
  this->HilbertSpaceDimension = dimension;
  double Tmp = 0.0;
  for (int i = 0; i < this->MPIRank; ++i)
    Tmp += this->ClusterPerformanceArray[i];
  if (this->MPIRank == 0)
    {
      this->MinimumIndex = (long) 0;
    }
  else
    {
      this->MinimumIndex = (long) (Tmp * ((double) dimension));
    }
  if (this->MPIRank == (this->NbrMPINodes - 1))
    {
      this->MaximumIndex = dimension - 1;
    }
  else
    {
      Tmp += this->ClusterPerformanceArray[this->MPIRank];      
      this->MaximumIndex = (long) (Tmp * ((double) dimension)) - (long) 1;
    }
  this->MinimumIndices = new long[this->NbrMPINodes];
  this->MaximumIndices = new long[this->NbrMPINodes];
  if (this->MasterNodeFlag)
    {
      this->MinimumIndices[0] = this->MinimumIndex;
      this->MaximumIndices[0] = this->MaximumIndex;
#ifdef __MPI__
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  MPI::COMM_WORLD.Recv(&(this->MinimumIndices[i]), 2, MPI::INT, i, 1);      
	  MPI::COMM_WORLD.Recv(&(this->MaximumIndices[i]), 2, MPI::INT, i, 1);      
	}
#endif 
    }
  else
    {
#ifdef __MPI__      
      MPI::COMM_WORLD.Send(&this->MinimumIndex, 2, MPI::INT, 0, 1); 
      MPI::COMM_WORLD.Send(&this->MaximumIndex, 2, MPI::INT, 0, 1);
#endif
    }
#ifdef __MPI__      
  MPI::COMM_WORLD.Bcast(this->MinimumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(this->MaximumIndices, 2 * this->NbrMPINodes, MPI::INT, 0);
#endif
}

// request an operation to the slave nodes and wait till they are ready to get operation parameters
//
// operationType = operation ID
// return value = true if no error occured

bool SimpleMPIArchitecture::RequestOperation (int operationType)
{
#ifdef __MPI__
  if (this->MasterNodeFlag)
    {
      MPI::COMM_WORLD.Bcast(&operationType, 1, MPI::INT, 0);
      int NbrMPINodes = MPI::COMM_WORLD.Get_size();
      bool Flag = true;
      int Acknowledge = 0;
      for (int i = 1; i < NbrMPINodes; ++i)
	{
	  MPI::COMM_WORLD.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	  if ((Flag == true) && (Acknowledge == 0))
	    Flag = false;
	}
      return Flag;
    }
#endif
  return false;
}

// wait an operation request from the master node  (without sending acknowledge)
//
// operationType = reference on the integer where the operation ID will be stored
// return value = true until the free slave signal is sent or an error occurs

bool SimpleMPIArchitecture::WaitOperation (int& operationType)
{
#ifdef __MPI__
  if (this->MasterNodeFlag == false)
    {
      MPI::COMM_WORLD.Bcast(&operationType, 1, MPI::INT, 0);
      if (operationType == SimpleMPIArchitecture::FreeSlaveSignal)
	{
	  return false;
	}
      else
	{
	  return true;
	}
    }
#endif
  return false;
}


// send acknowledge to the master node 
//
// acknowledge = true to send a positive answer
// return value = true if no error occured

bool SimpleMPIArchitecture::SendAcknowledge (bool acknowledge)
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 0;
      if (acknowledge == true)
	 Acknowledge = 1;
      MPI::COMM_WORLD.Send(&Acknowledge, 1, MPI::INT, 0, 1); 
      return true;
    }
#endif
  return false;
}

// wait for a slave to send the done signal
//
// return value = id of the slave that send the done signal

int SimpleMPIArchitecture::WaitAnySlave ()
{
#ifdef __MPI__
  if (this->MasterNodeFlag)
    {
      return 0;
    }
#endif
  return -1;
}

// send done signal to the master node 
//
// return value = true if no error occured
  
bool SimpleMPIArchitecture::SendDone ()
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 1;
      MPI::COMM_WORLD.Send(&Acknowledge, 1, MPI::INT, 0, 1); 
      return true;
    }
#endif
  return false;
}

// broadcast an integer from master node to slave nodes
// 
// value = integer to broadcast
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(int& value)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(&value, 1, MPI::INT, 0); 
  return true;
#else
  return false;
#endif  
}

// broadcast an integer from master node to slave nodes
// 
// value = integer to broadcast
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(long& value)
{
#ifdef __MPI__
#ifdef __64_BITS__
  MPI::COMM_WORLD.Bcast(&value, 2, MPI::INT, 0); 
#else
  MPI::COMM_WORLD.Bcast(&value, 1, MPI::INT, 0); 
#endif
  return true;
#else
  return false;
#endif  
}

// broadcast an integer array from master node to slave nodes
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(int* values, int nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(values, nbrValues, MPI::INT, 0); 
  return true;
#else
  return false;
#endif  
}

// broadcast an integer array from master node to slave nodes
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(long* values, int nbrValues)
{
#ifdef __MPI__
#ifdef __64_BITS__
  MPI::COMM_WORLD.Bcast(values, 2 * nbrValues, MPI::INT, 0); 
#else
  MPI::COMM_WORLD.Bcast(values, nbrValues, MPI::INT, 0); 
#endif
  return true;
#else
  return false;
#endif  
}

// send an integer array from master node to a given slave node
// 
// slaveID = slave ID
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::SendToSlaves(int slaveID, int* values, int nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Send(values, nbrValues, MPI::INT, slaveID + 1, 1);
  return true;
#else
  return false;
#endif  
}

// send an integer array from master node to a given slave node
// 
// slaveID = slave ID
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::SendToSlaves(int slaveID, long* values, int nbrValues)
{
#ifdef __MPI__
#ifdef __64_BITS__
  MPI::COMM_WORLD.Send(values, 2 * nbrValues, MPI::INT, slaveID + 1, 1);
#else
  MPI::COMM_WORLD.Send(values, nbrValues, MPI::INT, slaveID + 1, 1);
#endif  
  return true;
#else
  return false;
#endif  
}

// receive an integer array from master node to the given slave node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromMaster(int* values, int& nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::INT, 0, 1);
  return true;
#else
  return false;
#endif  
}

// receive an integer array from master node to the given slave node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromMaster(long* values, int& nbrValues)
{
#ifdef __MPI__
#ifdef __64_BITS__
  MPI::COMM_WORLD.Recv(values, 2 * nbrValues, MPI::INT, 0, 1);
#else
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::INT, 0, 1);
#endif  
  return true;
#else
  return false;
#endif  
}

// receive an integer array from master node to the given slave node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromMaster(int* values, long& nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::INT, 0, 1);
  return true;
#else
  return false;
#endif  
}

// send an integer array from the current slave node to master node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured
  
bool SimpleMPIArchitecture::SendToMaster(int* values, int nbrValues)
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 1;
      MPI::COMM_WORLD.Send(values, nbrValues, MPI::INT, 0, 1); 
      return true;
    }
#endif
  return false;
}

// send an integer array from the current slave node to master node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured
  
bool SimpleMPIArchitecture::SendToMaster(long* values, int nbrValues)
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 1;
#ifdef __64_BITS__
      MPI::COMM_WORLD.Send(values, 2 * nbrValues, MPI::INT, 0, 1); 
#else
      MPI::COMM_WORLD.Send(values, nbrValues, MPI::INT, 0, 1); 
#endif
      return true;
    }
#endif
  return false;
}

// send an integer array from the current slave node to master node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured
  
bool SimpleMPIArchitecture::SendToMaster(int* values, long nbrValues)
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 1;
      MPI::COMM_WORLD.Send(values, nbrValues, MPI::INT, 0, 1); 
      return true;
    }
#endif
  return false;
}

// receive an integer array from master node to the current slave node
// 
// slaveID = slave ID
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromSlave(int slaveID, int* values, int& nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::INT, slaveID + 1, 1);
  return true;
#else
  return false;
#endif  
}

// receive an integer array from master node to the current slave node
// 
// slaveID = slave ID
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromSlave(int slaveID, long* values, int& nbrValues)
{
#ifdef __MPI__
#ifdef __64_BITS__
  MPI::COMM_WORLD.Recv(values,  2 * nbrValues, MPI::INT, slaveID + 1, 1);
#else
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::INT, slaveID + 1, 1);
#endif
  return true;
#else
  return false;
#endif  
}

// receive an integer array from master node to the current slave node
// 
// slaveID = slave ID
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromSlave(int slaveID, int* values, long& nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::INT, slaveID + 1, 1);
  return true;
#else
  return false;
#endif  
}

// send a double array from the current slave node to master node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured
  
bool SimpleMPIArchitecture::SendToMaster(double* values, int nbrValues)
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 1;
      MPI::COMM_WORLD.Send(values, nbrValues, MPI::DOUBLE, 0, 1); 
      return true;
    }
#endif
  return false;
}

// receive a double array from master node to the current slave node
// 
// slaveID = slave ID
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromSlave(int slaveID, double* values, int& nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::DOUBLE, slaveID + 1, 1);
  return true;
#else
  return false;
#endif  
}

// send a double array from the current slave node to master node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured
  
bool SimpleMPIArchitecture::SendToMaster(double* values, long nbrValues)
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 1;
      MPI::COMM_WORLD.Send(values, nbrValues, MPI::DOUBLE, 0, 1); 
      return true;
    }
#endif
  return false;
}

// receive a double array from master node to the current slave node
// 
// slaveID = slave ID
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromSlave(int slaveID, double* values, long& nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::DOUBLE, slaveID + 1, 1);
  return true;
#else
  return false;
#endif  
}

#ifdef __LAPACK__

// send a double complex array from the current slave node to master node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured
  
bool SimpleMPIArchitecture::SendToMaster(doublecomplex* values, int nbrValues)
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 1;
      MPI::COMM_WORLD.Send(values, nbrValues, MPI::DOUBLE_COMPLEX, 0, 1); 
      return true;
    }
#endif
  return false;
}

// receive a double complex array from master node to the current slave node
// 
// slaveID = slave ID
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromSlave(int slaveID, doublecomplex* values, int& nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::DOUBLE_COMPLEX, slaveID + 1, 1);
  return true;
#else
  return false;
#endif  
}

// send a double complex array from the current slave node to master node
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured
  
bool SimpleMPIArchitecture::SendToMaster(doublecomplex* values, long nbrValues)
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 1;
      MPI::COMM_WORLD.Send(values, nbrValues, MPI::DOUBLE_COMPLEX, 0, 1); 
      return true;
    }
#endif
  return false;
}

// receive a double complex array from master node to the current slave node
// 
// slaveID = slave ID
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::ReceiveFromSlave(int slaveID, doublecomplex* values, long& nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Recv(values, nbrValues, MPI::DOUBLE_COMPLEX, slaveID + 1, 1);
  return true;
#else
  return false;
#endif  
}

#endif

// broadcast a double from master node to slave nodes
// 
// value = integer to broadcast
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(double& value)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(&value, 1, MPI::DOUBLE, 0); 
  return true;
#else
  return false;
#endif  
}

// broadcast a double array from master node to slave nodes
// 
// values = array of integers to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(double* values, int nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(values, nbrValues, MPI::DOUBLE, 0); 
  return true;
#else
  return false;
#endif  
}
  
// broadcast a vector on each slave node
//
// vector = pointer to the vector tobroadcast  (only useful for the master node)
// return value = pointer to the broadcasted vector or null pointer if an error occured

Vector* SimpleMPIArchitecture::BroadcastVector(Vector* vector)
{
#ifdef __MPI__
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      vector->BroadcastClone(MPI::COMM_WORLD, this->MPIRank);
      return vector;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	Vector TmpVector;
	return TmpVector.BroadcastClone(MPI::COMM_WORLD, 0);
      }
#endif  
  return 0;
}

// broadcast a vector from a node to the others 
//
// nodeID = id of the mode that broadcasts its vector
// vector = pointer to the vector to broadcast or to the vector where the content will be stored

void SimpleMPIArchitecture::BroadcastVector(int nodeID, Vector& vector)
{
#ifdef __MPI__
  vector.BroadcastVector(MPI::COMM_WORLD, nodeID);
#endif  
  return;
}

// scatter a vector upon each slave node
//
// vector = pointer to the vector to scatter  (only useful for the master node)
// return value = pointer to the broadcasted vector or null pointer if an error occured

Vector* SimpleMPIArchitecture::ScatterVector(Vector* vector)
{
#ifdef __MPI__
  int TmpIndices[2];
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  MPI::COMM_WORLD.Recv(TmpIndices, 2, MPI::INT, i, 1); 	      
	  vector->SendPartialClone(MPI::COMM_WORLD, i, TmpIndices[0], TmpIndices[1]);
	}
      return vector;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	TmpIndices[0] = (int) this->MinimumIndex;
	TmpIndices[1] = ((int) this->MaximumIndex) - TmpIndices[0] + 1;
	MPI::COMM_WORLD.Send(TmpIndices, 2, MPI::INT, 0, 1); 	      
	Vector TmpVector;
	return TmpVector.ReceivePartialClone(MPI::COMM_WORLD, 0);
      }
#endif  
  return 0;
}


// scatter a vector onto slave nodes, using the scatterV function
//
// vector = pointer to the vector to scatter  (only useful for the master node)
// return value = pointer to the broadcasted vector or null pointer if an error occured

Vector* SimpleMPIArchitecture::ScatterVectorNew(Vector* vector)
{
#ifdef __MPI__
  int TmpIndices[2];
  if ((this->MasterNodeFlag) && (vector != 0))    
    {
      vector->ScatterPartialClones(MPI::COMM_WORLD, this->MinimumIndices, this->MaximumIndices, 0);
      return vector;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	Vector TmpVector;
	return TmpVector.ReceiveScatteredClone(MPI::COMM_WORLD, 0);
      }
#endif  
  return 0;
}

// broadcast a vector type and allocate a vector based on it on each slave node
//
// vector = pointer to the vector to be used as reference (only useful for the master node)
// return value = pointer to the cloned vector or null pointer if an error occured

Vector* SimpleMPIArchitecture::BroadcastVectorType(Vector* vector)
{
#ifdef __MPI__
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      vector->BroadcastEmptyClone(MPI::COMM_WORLD, this->MPIRank);
      return vector;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	Vector TmpVector;
	return TmpVector.BroadcastEmptyClone(MPI::COMM_WORLD, 0);
      }
#endif  
  return 0;
}

// broadcast an array of vectors on each slave node
//
// nbrVectors = reference on the number of vectors to broadcast or get
// vector = pointer to the vector tobroadcast  (only useful for the master node)
// return value =  pointer to the array of broadcasted vectors or null pointer if an error occured null pointer if an error occured

Vector** SimpleMPIArchitecture::BroadcastVectorArray(int& nbrVectors, Vector* vector)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(&nbrVectors, 1, MPI::INT, 0); 
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      switch (vector->GetVectorType())
	{
	case Vector::RealDatas:
	  for (int i = 0; i < nbrVectors; ++i)
	    ((RealVector*) vector)[i].BroadcastClone(MPI::COMM_WORLD, this->MPIRank);
	  break;
	case Vector::ComplexDatas:	    
	  for (int i = 0; i < nbrVectors; ++i)
	    ((ComplexVector*) vector)[i].BroadcastClone(MPI::COMM_WORLD, this->MPIRank);
	  break;
	default :
	  {
	    cout << "warning : SimpleMPIArchitecture::BroadcastVectorArray can't guess vector type" << endl;	    
	  }
	  break;
	}
      return 0;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	
	Vector** TmpVectorArray = new Vector*[nbrVectors];
	Vector TmpVector;
	for (int i = 0; i < nbrVectors; ++i)
	  TmpVectorArray[i] = TmpVector.BroadcastClone(MPI::COMM_WORLD, 0);
	return TmpVectorArray;
      }
#endif  
  return 0;
}

// broadcast vector type and allocate an array of vectors based on it on each slave node
//
// nbrVectors = reference on the number of vectors to broadcast or get
// vector = pointer to the vector to be used as reference (only useful for the master node)
// return value =  pointer to the array of cloned vector or null pointer if an error occurednull pointer if an error occured

Vector** SimpleMPIArchitecture::BroadcastVectorTypeArray(int& nbrVectors, Vector* vector)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(&nbrVectors, 1, MPI::INT, 0); 
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      vector[0].BroadcastEmptyClone(MPI::COMM_WORLD, this->MPIRank);
      return 0;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	Vector** TmpVectorArray = new Vector*[nbrVectors];
	Vector TmpVector;
	TmpVectorArray[0] = TmpVector.BroadcastEmptyClone(MPI::COMM_WORLD, 0);
	for (int i = 1; i < nbrVectors; ++i)
	  TmpVectorArray[i] = TmpVectorArray[0]->EmptyClone();
	return TmpVectorArray;
      }
#endif  
  return 0;
}

// scatter an array of vectors upon each slave node
//
// nbrVectors = reference on the number of vectors to broadcast or get
// vector = pointer to the vector to be used as reference (only useful for the master node)
// return value = pointer to the broadcasted vector or null pointer if an error occured

Vector** SimpleMPIArchitecture::ScatterVectorArray(int& nbrVectors, Vector* vector)
{
#ifdef __MPI__
  int TmpIndices[2];
  MPI::COMM_WORLD.Bcast(&nbrVectors, 1, MPI::INT, 0); 
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      switch (vector->GetVectorType() & Vector::DataTypeMask)
	{
	case Vector::RealDatas:
	  {
	    for (int i = 1; i < this->NbrMPINodes; ++i)
	      {
		MPI::COMM_WORLD.Recv(TmpIndices, 2, MPI::INT, i, 1); 	      
		for (int j = 0; j < nbrVectors; ++j)		 
		  ((RealVector*) vector)[j].SendPartialClone(MPI::COMM_WORLD, i, TmpIndices[0], TmpIndices[1]);
	      }
	  }
	  break;
	case Vector::ComplexDatas:
	  {
	    for (int i = 1; i < this->NbrMPINodes; ++i)
	      {
		MPI::COMM_WORLD.Recv(TmpIndices, 2, MPI::INT, i, 1); 	      
		for (int j = 0; j < nbrVectors; ++j)		 
		  ((ComplexVector*) vector)[j].SendPartialClone(MPI::COMM_WORLD, i, TmpIndices[0], TmpIndices[1]);
	      }
	  }
	  break;
	default :
	  {
	    cout << "warning : ScatterVectorArray can't guess vector type" << endl;	    
	  }
	  break;
	}
      return 0;
    }
  else
    {
      if (this->MasterNodeFlag == false)
	{
	  TmpIndices[0] = (int) this->MinimumIndex;
	  TmpIndices[1] = ((int) this->MaximumIndex) - TmpIndices[0] + 1;
	  MPI::COMM_WORLD.Send(TmpIndices, 2, MPI::INT, 0, 1); 	      
	  Vector** TmpVectorArray = new Vector*[nbrVectors];
	  for (int j = 0; j < nbrVectors; ++j)		 
	    {	    
	      Vector TmpVector;
	      TmpVectorArray[j] = TmpVector.ReceivePartialClone(MPI::COMM_WORLD, 0);
	    }
	  return TmpVectorArray;
	}
      else
	{
	  cout << "warning : master node uses SimpleMPIArchitecture::ScatterVectorArray with an empty array" << endl;
	}
    }
#endif  
  return 0;
}

// broadcast a matrix on each slave node
//
// matrix = atrix to broadcast or to the matrix where the content will be stored

void SimpleMPIArchitecture::BroadcastMatrix(Matrix& matrix)
{
#ifdef __MPI__
  matrix.BroadcastMatrix(MPI::COMM_WORLD, 0);
#endif  
}

// broadcast a matrix on each slave node
//
// matrix = pointer to the matrix tobroadcast  (only useful for the master node)
// return value = pointer to the broadcasted matrix or null pointer if an error occured

Matrix* SimpleMPIArchitecture::BroadcastMatrix(Matrix* matrix)
{
#ifdef __MPI__
  if ((this->MasterNodeFlag) && (matrix != 0))
    {
      matrix->BroadcastClone(MPI::COMM_WORLD, this->MPIRank);
      return matrix;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	Matrix TmpMatrix;
	return TmpMatrix.BroadcastClone(MPI::COMM_WORLD, 0);
      }
#endif  
  return 0;
}

// broadcast a matrix from a node to the others 
//
// nodeID = id of the mode that broadcasts its matrix
// matrix = matrix to broadcast or to the matrix where the content will be stored

void SimpleMPIArchitecture::BroadcastMatrix(int nodeID, Matrix& matrix)
{
#ifdef __MPI__
  matrix.BroadcastMatrix(MPI::COMM_WORLD, nodeID);
#endif  
}

// broadcast an array of matrix on each slave node
//
// nbrMatrices = reference on the number of matrices to broadcast or get
// matrix = pointer to the matrix to broadcast  (only useful for the master node)
// return value =  pointer to the array of broadcasted matrices or null pointer if an error occured null pointer if an error occured

Matrix** SimpleMPIArchitecture::BroadcastMatrixArray(int& nbrMatrices, Matrix* matrix)
{
#ifdef __MPI__
  int TmpIndices[2];
  MPI::COMM_WORLD.Bcast(&nbrMatrices, 1, MPI::INT, 0); 
  if ((this->MasterNodeFlag) && (matrix != 0))
    {
      switch (matrix->GetMatrixType())
	{
	case Matrix::RealElements:
	  {
	    for (int i = 1; i < this->NbrMPINodes; ++i)
	      {
		MPI::COMM_WORLD.Recv(TmpIndices, 2, MPI::INT, i, 1); 	      
		for (int j = 0; j < nbrMatrices; ++j)		 
		  ((RealMatrix*) matrix)[j].SendPartialClone(MPI::COMM_WORLD, i, TmpIndices[0], TmpIndices[1]);
	      }
	  }
	  break;
	case Matrix::ComplexElements:
	  {
	    for (int i = 1; i < this->NbrMPINodes; ++i)
	      {
		MPI::COMM_WORLD.Recv(TmpIndices, 2, MPI::INT, i, 1); 	      
		for (int j = 0; j < nbrMatrices; ++j)		 
		  ((ComplexMatrix*) matrix)[j].SendPartialClone(MPI::COMM_WORLD, i, TmpIndices[0], TmpIndices[1]);
	      }
	  }
	  break;
	default :
	  {
	    cout << "warning : ScatterMatrixArray can't guess matrix type" << endl;	    
	  }
	  break;
	}
      return 0;
    }
  else
    {
      if (this->MasterNodeFlag == false)
	{
	  TmpIndices[0] = (int) this->MinimumIndex;
	  TmpIndices[1] = ((int) this->MaximumIndex) - TmpIndices[0] + 1;
	  MPI::COMM_WORLD.Send(TmpIndices, 2, MPI::INT, 0, 1); 	      
	  Matrix** TmpMatrixArray = new Matrix*[nbrMatrices];
	  for (int j = 0; j < nbrMatrices; ++j)		 
	    {	    
	      Matrix TmpMatrix;
	      TmpMatrixArray[j] = TmpMatrix.ReceivePartialClone(MPI::COMM_WORLD, 0);
	    }
	  return TmpMatrixArray;
	}
      else
	{
	  cout << "warning : master node uses SimpleMPIArchitecture::ScatterMatrixArray with an empty array" << endl;
	}
    }
#endif  
  return 0;
}

// get a temporary file name
//
// return value = string corresponding to a temporary file name

char* SimpleMPIArchitecture::GetTemporaryFileName()
{
  timeval Time;
  gettimeofday (&Time, 0);
  char* TmpString = new char [64];
  sprintf (TmpString, "diagam%d%d_node%d.tmp",(int)  Time.tv_sec, (int)  Time.tv_usec, this->MPIRank);
  return TmpString;
}
  
// add an entry to the log file 
//
// message = string corresponding to entry to add to the log file
// masterFlag = true if only the master node should add the entry
// return value = true if no error occured

bool SimpleMPIArchitecture::AddToLog(const char * message, bool masterFlag)
{
#ifdef __MPI__
  if (this->MasterNodeFlag == false)
    {
      if (masterFlag == false)
	{
	  int TmpMessageLength = strlen(message);
	  MPI::COMM_WORLD.Send(&TmpMessageLength, 1, MPI::INT, 0, 1);
	  MPI::COMM_WORLD.Send(message, TmpMessageLength, MPI::CHAR, 0, 1);
	  return true;
	}
      else
	return false;
    }
  else
    {
      ofstream File;
      File.open(this->LogFile, ios::out | ios::app);
      if (!File.is_open())
	cout << "ERROR : cannot write log file " << this->LogFile << endl;
      int TmpMessageLength = strlen(message);
      int TmpInc;
      char* TmpMessage = new char[TmpMessageLength + 256];
      sprintf (TmpMessage, "node 0: %s", message);
      File << TmpMessage << endl;
      delete [] TmpMessage;
      if (masterFlag == false)
	for (int i = 1; i < NbrMPINodes; ++i)
	  {
	    MPI::COMM_WORLD.Recv(&TmpMessageLength, 1, MPI::INT, i, 1);      
	    TmpMessage = new char[TmpMessageLength + 256];
	    sprintf (TmpMessage, "node %d: ", i);
	    TmpInc = strlen(TmpMessage);
	    TmpMessage[TmpInc + TmpMessageLength] = '\0';
	    MPI::COMM_WORLD.Recv(TmpMessage + strlen(TmpMessage), TmpMessageLength, MPI::CHAR, i, 1);  
	    File << TmpMessage << endl;
	    delete [] TmpMessage;
	  }
      File.close();
    }
  return true;
#else
  return false;
#endif  
}

// dump the log file into a string
//
// header = optional header to add before the log file
// footer = optional footer to add at the end of the log file
// return value = string or 0 if an error occured or log is not available

char* SimpleMPIArchitecture::DumpLog(const char* header, const char* footer)
{
  if ((this->VerboseModeFlag == false) || (this->LogFile == 0))
    return 0;
  return DumpTextFile(this->LogFile, header, footer);
}


// write vector in a file 
//
// vector = vector to write
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs
  
bool SimpleMPIArchitecture::WriteVector(RealVector& vector, const char* fileName)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Barrier();
  if ( this->IsMasterNode() )
    {
      return vector.WriteVector(fileName);
    }
  else
    {
      return true;
    }
#else
  return vector.WriteVector(fileName);
#endif
}


// read a vector from a file but only on master
//
// vector = vector to read
// fileName = name of the file where the vector is read from
// return value = true if no error occurs
  
bool SimpleMPIArchitecture::ReadVector(RealVector& vector, const char* fileName)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Barrier();
  int Value = 0;
  if ( this->IsMasterNode() )
    {
      if ( vector.ReadVector(fileName) )
	{
	  Value = 1;
	}
      else
	{
	  Value = 0;
	}
      this->BroadcastToSlaves(Value);
    }
  else
    {
      this->BroadcastToSlaves(Value);
    }
  if ( Value == 1 ) 
    {
      return true;
    }
  else 
    return false;
#else
  return vector.WriteVector(fileName);
#endif
}

// indicate if the current architecture allows to write on disk
//
// return value = true if the current architecture allows to write on disk

bool SimpleMPIArchitecture::CanWriteOnDisk()
{
  return this->IsMasterNode();
}


// deduce the load distribution for the given load profile (running on Master node)
// totalNbrOperations = total number of operations
// nbrOperationPerIndex = number of operations per index
// minimumIndices = minimum indices to be determined
// maximumIndices = maximum indices to be determined
 void SimpleMPIArchitecture::DeduceLoadDistribution(long totalNbrOperations, int *nbrOperationPerIndex, long *&minimumIndices, long *&maximumIndices)
{
  char* TmpString = new char [512];
  sprintf (TmpString, "total number of operations = %ld", totalNbrOperations);
  if (this->LogFile != 0)
    this->AddToLog(TmpString, true);

  minimumIndices = new long[this->NbrMPINodes];
  maximumIndices = new long[this->NbrMPINodes];
  minimumIndices[0] = this->MinimumIndices[0];
  maximumIndices[this->NbrMPINodes - 1] = this->MaximumIndices[this->NbrMPINodes - 1];
  long TmpMininumIndex = minimumIndices[0];
  long TmpMaximumIndex = maximumIndices[this->NbrMPINodes - 1];
  for (int i = 0; i < (this->NbrMPINodes - 1); ++i)
    {      
      long LocalNbrOperations = (long) ((((double) totalNbrOperations) * this->ClusterPerformanceArray[i]));
      if (LocalNbrOperations == 0l)
	LocalNbrOperations = 1l;
      long TrueLocalNbrOperations = 0l;
      while ((TmpMininumIndex <= TmpMaximumIndex) && (TrueLocalNbrOperations < LocalNbrOperations))
	{
	  TrueLocalNbrOperations += (long) nbrOperationPerIndex[TmpMininumIndex];
	  ++TmpMininumIndex;
	}
      maximumIndices[i] = TmpMininumIndex - 1l;
      minimumIndices[i + 1] = TmpMininumIndex;
      sprintf (TmpString, "number of operations on node %d = %ld", i, TrueLocalNbrOperations);	  
      if (this->LogFile != 0)
	this->AddToLog(TmpString, true);     	  
    }
  long TrueLocalNbrOperations = 0;
  while (TmpMininumIndex <= TmpMaximumIndex)
    {
      TrueLocalNbrOperations += nbrOperationPerIndex[TmpMininumIndex];
      ++TmpMininumIndex;
    }
  sprintf (TmpString, "number of operations on node %d = %ld", (this->NbrMPINodes - 1), TrueLocalNbrOperations);	  
  if (this->LogFile != 0)
    this->AddToLog(TmpString, true);  
  delete[] TmpString;
}

// balance an array differently across nodes (currently supporting T=int, long)
//
// nbrOperationPerIndex = reference on a local array holding one entry per local state prior to last call to GetOptimizedTypicalRange; on return - local array holding data for the new range of the MPIArchitecture
// return value = true if the array has been rebalanced
template<typename T>
bool SimpleMPIArchitecture::RebalanceArrayImplementation(T*& array, const char* filename)
{
  // don't do anything if the old min/max indices are not set, or are identical to the current ones
  if (this->OldMaximumIndices == NULL)
    return false;
  bool Flag=false;
  for (int i=0; i<this->NbrMPINodes; ++i)
    if ( (this->OldMaximumIndices[i] != this->MaximumIndices[i]) || (this->OldMinimumIndices[i] != this->MinimumIndices[i]))
      {
	Flag=true;
	break;
      }
  if (Flag == false)
    return false;
  // else, rebalance the arrays:
  if (this->MasterNodeFlag == true)
    { 
      // gather full array on Master node
      T* TmpArray = new T [this->MaximumIndices[this->NbrMPINodes - 1] - this->MinimumIndices[0] + 1l];
      long TmpIndex = 0;
      long EffectiveDimension = this->OldMaximumIndices[0] - this->OldMinimumIndices[0] + 1l;
      long TotalEffectiveDimension = EffectiveDimension;
      for (long i = 0l; i < EffectiveDimension; ++i)
	{
	  TmpArray[this->OldMinimumIndices[0] + i] = array[i];
	}
      TmpIndex += EffectiveDimension;      
      for (int i = 1; i < NbrMPINodes; ++i)
	{
	  int TmpNbrElements = (int) (this->OldMaximumIndices[i] - this->OldMinimumIndices[i] + 1);
	  if (TmpNbrElements < (T) (this->OldMaximumIndices[i] - this->OldMinimumIndices[i] + 1))
	    {
	      char errMsg[100] = "Count exceeds int maximum in RebalanceArrayImplementation.";
	      this->AddToLog(errMsg);	      
	    }
	  this->ReceiveFromSlave(i - 1, TmpArray + this->OldMinimumIndices[i], TmpNbrElements);
	  TotalEffectiveDimension += TmpNbrElements;
	}

      if (filename != NULL)
	{
	  // write array to disk

	  ofstream Store(filename, std::ios::binary | std::ios::out);
	  long TotalDimension = this->MaximumIndices[this->NbrMPINodes - 1] - this->MinimumIndices[0] + 1l;
	  long Sum=0;
	  for (long i=0; i<TotalDimension; ++i)
	    Sum += TmpArray[i];
	  Store.write((char*) &(TotalDimension), sizeof(long));
	  Store.write((char*) &(Sum), sizeof(long));
	  Store.write((char*) TmpArray, (TotalDimension)*sizeof(T));
	  long CheckMark=(long)(1e12*M_PI);
	  Store.write((char*) &(CheckMark), sizeof(long));
	  // done writing array
	}

      // redistribute array according to new load balancing
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  cout << "On node "<<i<<": "<<this->MaximumIndices[i] - this->MinimumIndices[i] + 1l<<endl;
      
	  this->SendToSlaves(i - 1, (T*) (TmpArray + this->MinimumIndices[i]), 
			     (int) (this->MaximumIndices[i] - this->MinimumIndices[i] + 1l));      
	}
      delete[] array;
      array = new T [this->MaximumIndices[0] - this->MinimumIndices[0] + 1l];
      cout << "On Master node: "<<this->MaximumIndices[0] - this->MinimumIndices[0] + 1l<<endl;
      for (long i = this->MinimumIndices[0]; i <= this->MaximumIndices[0]; ++i)
	array[i - this->MinimumIndices[0]] = TmpArray[i];
      delete[] TmpArray;
    }
  else
    {
      this->SendToMaster(array, (int)(this->OldMaximumIndices[this->MPIRank] - this->OldMinimumIndices[this->MPIRank] + 1));
      int TmpNbrElements = (int) (this->MaximumIndices[this->MPIRank] - this->MinimumIndices[this->MPIRank] + 1l);
      array = new T [TmpNbrElements];
      this->ReceiveFromMaster(array, TmpNbrElements);
    }
  return true;
}



// load an array on Master node and distribute among nodes according to current balancing
//
// array = reference on a local array, to be allocated and filled with the loaded data
// return value = true if the array has been rebalanced
template<typename T>
bool SimpleMPIArchitecture::LoadArrayImplementation(T*& array, const char* filename)
{
  if (this->MasterNodeFlag == true)
    {
      ifstream Store;
      Store.open(filename, std::ios::binary | std::ios::in);
      int Acknowledge=0;
      if (!Store.is_open())
	{
	  this->BroadcastToSlaves(Acknowledge);
	  return false; // no load pattern found.
	}
      
      // read load balancing information from disk     
      long TotalDimension;
      Store.read((char*) &(TotalDimension), sizeof(long));
      
      if (TotalDimension != this->MaximumIndices[this->NbrMPINodes - 1] - this->MinimumIndices[0] + 1l)
	{
	  cout << "Attention, load profile does not match current dimension - ignoring data"<<endl;
	  this->BroadcastToSlaves(Acknowledge);
	  return false;
	}
      long CheckSum = 0;
      Store.read((char*) &(CheckSum), sizeof(long));
      T *TmpArray = new T[TotalDimension];
      Store.read((char*) (TmpArray), (TotalDimension)*sizeof(T));
      long CheckMark;
      Store.read((char*) &(CheckMark), sizeof(long));
      if (CheckMark!=((long)(1e12*M_PI)))
	{
	  cout << "Error reading load profile"<<endl;
	  this->BroadcastToSlaves(Acknowledge);
	  return false;
	}
      long Sum=0;
      for (long i=0; i<TotalDimension; ++i)
	Sum += TmpArray[i];
      if (CheckSum!=Sum)
	{
	  cout << "Error recovering checksum from file "<<filename<<endl;
	  this->BroadcastToSlaves(Acknowledge);
	  return false;
	}
      // done reading - signal success
      Acknowledge=1;
      this->BroadcastToSlaves(Acknowledge);

      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  this->SendToSlaves(i - 1, (int*) (TmpArray + this->MinimumIndices[i]), 
			     (int) (this->MaximumIndices[i] - this->MinimumIndices[i] + 1l));      
	}
      array = new T [this->MaximumIndices[0] - this->MinimumIndices[0] + 1l];
      for (long i = this->MinimumIndices[0]; i <= this->MaximumIndices[0]; ++i)
	array[i - this->MinimumIndices[0]] = TmpArray[i];
      delete[] TmpArray;
    }
  else
    {
      int Acknowledge=0;
      this->BroadcastToSlaves(Acknowledge);
      if (Acknowledge==0)
	return false; // data could not be loaded successfully.
      
      int TmpNbrElements = (int) (this->MaximumIndices[this->MPIRank] - this->MinimumIndices[this->MPIRank] + 1l);
      array = new T [TmpNbrElements];
      this->ReceiveFromMaster(array, TmpNbrElements);
    }
  return true;
}

// print load balancing
ostream& SimpleMPIArchitecture::PrintLoadBalancing(ostream &Str)
{
  Str << "Nodes operating on ranges: ";
  Str << "["<<this->MinimumIndices[0]<<"->"<<this->MaximumIndices[0]<<"]";
  for (int i=1; i<this->NbrMPINodes; ++i)
    Str << ", ["<<this->MinimumIndices[i]<<"->"<<this->MaximumIndices[i]<<"]";
  return Str;
}
