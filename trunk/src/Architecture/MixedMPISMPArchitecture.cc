////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of mixed MPI - SMP Architecture                  //
//                                                                            //
//                        last modification : 23/04/2007                      //
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
#include "Architecture/MixedMPISMPArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/SMPArchitecture.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <fstream>
#include <cstring>
#include <sys/time.h>
#include <string.h>
#include <iostream>
#include <unistd.h>
#ifdef __MPI__
#include <mpi.h>
#endif


using std::ofstream;
using std::ios;
using std::cout;
using std::endl;


// constructor
//
// clusterFileName = name of the file that describes the cluster, if none assume one cpu per MPI node. The file should be at least accessible by the master mode
// logFile = name of the optional log file to allow code profiling on MPI architecture
// automaticLoadBalancing = flag that indicates if automatic load balancing have to be done, overriding any manual load balancing

MixedMPISMPArchitecture::MixedMPISMPArchitecture(char* clusterFileName, char* logFile, bool automaticLoadBalancing)
{
  this->PerformanceIndex = 1.0;
  this->ArchitectureID = AbstractArchitecture::MixedMPISMP;
  this->AutomaticLoadBalancing = automaticLoadBalancing;
#ifdef __MPI__
  MPI::Init();
  this->NbrMPINodes = MPI::COMM_WORLD.Get_size();
  this->MPIRank = MPI::COMM_WORLD.Get_rank();
  this->NbrCPUPerNode = new int [this->NbrMPINodes];
  this->ClusterMemoryArray = new long [this->NbrMPINodes];
  this->ClusterPerformanceArray = new double [this->NbrMPINodes];

  char* TmpLocalHostname = new char [512];
  gethostname(TmpLocalHostname, 511);

  if (this->MPIRank != 0)
    {
      this->NodeHostnames = 0;
      int HostnameStringSize = strlen(TmpLocalHostname) + 1;
      MPI::COMM_WORLD.Send(&HostnameStringSize, 1, MPI::INT, 0, 1);      
      MPI::COMM_WORLD.Send(TmpLocalHostname, HostnameStringSize, MPI::CHAR, 0, 1);      
    }
  else
    {
      this->NodeHostnames = new char*[this->NbrMPINodes];
      int HostnameStringSize = strlen(TmpLocalHostname) + 1;
      this->NodeHostnames[0] = new char[HostnameStringSize];
      strncpy (this->NodeHostnames[0], TmpLocalHostname, HostnameStringSize);
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  MPI::COMM_WORLD.Recv(&HostnameStringSize, 1, MPI::INT, i, 1);
	  this->NodeHostnames[i] = new char[HostnameStringSize];
	  MPI::COMM_WORLD.Recv(this->NodeHostnames[i], HostnameStringSize, MPI::CHAR, i, 1);	  
	}
    }
  delete[] TmpLocalHostname;

  if (this->MPIRank == 0)
    {
      this->MasterNodeFlag = true;
      if (clusterFileName != 0)
	{
	  MultiColumnASCIIFile ClusterFile;
	  bool ErrorFlag = false;
	  if (ClusterFile.Parse(clusterFileName) == true)
	    {
	      if (ClusterFile.GetNbrColumns() < 2)
		ErrorFlag = true;
	      else
		{
		  double* TmpClusterPerformanceArray = 0;
		  long* TmpClusterMemoryArray = 0;
		  if (ClusterFile.GetNbrColumns() > 2)
		    TmpClusterPerformanceArray = ClusterFile.GetAsDoubleArray(2);
		  if (ClusterFile.GetNbrColumns() > 3)
		    TmpClusterMemoryArray = ClusterFile.GetAsLongArray(3);
		  int* TmpNbrCPUNode = ClusterFile.GetAsIntegerArray(1);
		  if (TmpNbrCPUNode != 0)
		    {
		      int DefaultCPUPerNode=1;
		      double DefaultPerNodePerformance=1.0;
		      long DefaultPerNodeMemory=-1l;
		      for (int j = 0; j < ClusterFile.GetNbrLines(); ++j)
			{
			  if (strcmp("default", ClusterFile(0,j)) == 0) 
			    {
			      DefaultCPUPerNode = TmpNbrCPUNode[j];
			      DefaultPerNodePerformance = (double)DefaultCPUPerNode;
			      if (TmpClusterPerformanceArray !=0)
				DefaultPerNodePerformance *= TmpClusterPerformanceArray[j];
			      if (TmpClusterMemoryArray != 0)
				DefaultPerNodeMemory = (TmpClusterMemoryArray[j]) << 20;
			      cout << "Default node configuration: " << DefaultCPUPerNode <<" CPUs, Performance: "
				   <<DefaultPerNodePerformance<<", Memory: "<<DefaultPerNodeMemory<<endl;
			      j = ClusterFile.GetNbrLines();
			    }
			}
		      int DefaultMasterCPUs=DefaultCPUPerNode;
		      double DefaultMasterPerformance=1.0;
		      long DefaultMasterMemory=-1l;
		      bool OverrideMaster=false;
		      for (int j = 0; j < ClusterFile.GetNbrLines(); ++j)
			{
			  if (strcmp("master", ClusterFile(0,j)) == 0) 
			    {
			      DefaultMasterCPUs = TmpNbrCPUNode[j];
			      DefaultMasterPerformance  = (double)DefaultMasterCPUs;
			      if (TmpClusterPerformanceArray !=0)
				DefaultMasterPerformance *= TmpClusterPerformanceArray[j];
			      DefaultMasterMemory = (TmpClusterMemoryArray[j]) << 20;
			      cout << "Master node configuration: "<<DefaultMasterCPUs<<" CPUs, Performance: "
				   <<DefaultMasterPerformance<<", Memory: "<<DefaultMasterMemory<<endl;

			      OverrideMaster = true;
			      j = ClusterFile.GetNbrLines();
			    }
			}
		      for (int i = 0; i < this->NbrMPINodes; ++i)
			{
			  this->NbrCPUPerNode[i] = DefaultCPUPerNode;
			  this->ClusterPerformanceArray[i] = DefaultCPUPerNode;
			  this->ClusterMemoryArray[i] = DefaultPerNodeMemory;
			  for (int j = 0; j < ClusterFile.GetNbrLines(); ++j)
			    {
			      if ((strcmp(this->NodeHostnames[i], ClusterFile(0,j)) == 0) || 
				  ((strncmp(this->NodeHostnames[i], ClusterFile(0,j), strlen(ClusterFile(0,j))) == 0) && 
				   (this->NodeHostnames[i][strlen(ClusterFile(0,j))] == '.')) ||
				  ((strncmp(this->NodeHostnames[i], ClusterFile(0,j), strlen(ClusterFile(0,j)) - 1) == 0) && 
				   (ClusterFile(0,j)[strlen(ClusterFile(0,j)) - 1] == '*')))
				{
				  this->NbrCPUPerNode[i] = TmpNbrCPUNode[j];
				  this->ClusterPerformanceArray[i] = (double) this->NbrCPUPerNode[i];
				  if (TmpClusterPerformanceArray !=0)
				    this->ClusterPerformanceArray[i] *= TmpClusterPerformanceArray[j];
				  if (TmpClusterMemoryArray != 0)
				    this->ClusterMemoryArray[i] = (TmpClusterMemoryArray[j]) << 20;
				  else
				    this->ClusterMemoryArray[i] = -1l;
				  j = ClusterFile.GetNbrLines();
				}
			    }
			}
		      if (OverrideMaster)
			{
			  this->NbrCPUPerNode[0] = DefaultMasterCPUs;
			  this->ClusterPerformanceArray[0] = DefaultMasterPerformance;
			  ClusterMemoryArray[0] = DefaultMasterMemory;
			}
		      delete[] TmpNbrCPUNode;
		    }
		  else
		    ErrorFlag = true;
		  // clean up temporary arrays
		  if (ClusterFile.GetNbrColumns() > 2)
		    delete [] TmpClusterPerformanceArray;
		  if (ClusterFile.GetNbrColumns() > 3)
		    delete [] TmpClusterMemoryArray;
		}
	    }
	  else
	    {
	      cout << "Errors parsing cluster configuration file " << clusterFileName << endl;
	      ClusterFile.DumpErrors(cout);
	      cout << "Aborting calculation."<<endl;
	      exit(1);
	    }
	  if (ErrorFlag == true)
	    {
	      cout << "an error occured while opening " << clusterFileName << endl;
	      ClusterFile.DumpErrors(cout);
	      cout << "switching to one cpu per node mode" << endl;
	      clusterFileName = 0;
	    }
	}
      if (clusterFileName == 0)
	for (int i = 0; i < this->NbrMPINodes; ++i)
	  {
	    this->NbrCPUPerNode[i] = 1;
	    this->ClusterPerformanceArray[i] = 1.0;
	  }
      this->PerformanceIndex = this->ClusterPerformanceArray[this->MPIRank];
      this->TotalPerformanceIndex = this->PerformanceIndex;
      for (int i = 1; i < this->NbrMPINodes; ++i)
	this->TotalPerformanceIndex += this->ClusterPerformanceArray[i];
      for (int i = 0; i < this->NbrMPINodes; ++i)
	this->ClusterPerformanceArray[i] /= this->TotalPerformanceIndex;      
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
	    {
	      File << "number of nodes = " << this->TotalPerformanceIndex << endl;
	      File << "cluster description : " << endl;
	      for (int i = 0; i < this->NbrMPINodes; ++i)
		File << "    node " << i << " :  hostname=" << this->NodeHostnames[i] << "  cpu=" << this->NbrCPUPerNode[i] << "  perf. index=" << this->ClusterPerformanceArray[i] 
		     <<  " mem=" << (this->ClusterMemoryArray[i] >> 20) << "Mb" << endl;
	      File << " ---------------------------------------------" << endl
		   << "                    profiling                 " << endl
		   << " ---------------------------------------------" << endl;
	      
	      this->VerboseModeFlag = true;
	    }
	  File.close();
	}
      else
	{
	  this->VerboseModeFlag = false;
	  this->LogFile = 0;
	}
    }
  else
    {
      this->MasterNodeFlag = false;
      if (logFile != 0)
	this->VerboseModeFlag = true;
      else
	this->VerboseModeFlag = false;
      this->LogFile = 0;
    }
  MPI::COMM_WORLD.Bcast(this->NbrCPUPerNode, this->NbrMPINodes, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(this->ClusterPerformanceArray, this->NbrMPINodes, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&this->TotalPerformanceIndex, 1, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(this->ClusterMemoryArray, this->NbrMPINodes, MPI::LONG, 0);
#else
  this->MasterNodeFlag = true;
  this->NbrMPINodes = 1;
  this->MPIRank = 0;
  this->ClusterPerformanceArray = 0;
  this->TotalPerformanceIndex = this->PerformanceIndex;
  this->NbrCPUPerNode = new int[1];
  this->NbrCPUPerNode[0] = 1;
  this->ClusterMemoryArray = 0;
#endif
  if (this->MasterNodeFlag == true)
    {
      cout << "number of nodes = " << this->NbrMPINodes << endl << "total performance index = " << this->TotalPerformanceIndex << endl;
      for (int i = 0; i < this->NbrMPINodes; ++i)
	cout << "node " << i << " using " << this->NbrCPUPerNode[i] << " CPU/cores, performance index = " << this->ClusterPerformanceArray[i] << endl;
    }
  if (this->NbrCPUPerNode[this->MPIRank] > 1)
    {
      //      delete this->LocalArchitecture;
      if (logFile == 0)
	{
	  this->LocalArchitecture = new SMPArchitecture(this->NbrCPUPerNode[this->MPIRank]);
	}
      else
	{
	  char* TmpLogFileName = new char [strlen(logFile) + 128];
	  sprintf (TmpLogFileName, "%s.local.%i.dat", logFile, this->MPIRank);
	  this->LocalArchitecture = new SMPArchitecture(this->NbrCPUPerNode[this->MPIRank], TmpLogFileName);
	  delete[] TmpLogFileName;
	}
    }
  else
    {
      this->LocalArchitecture = new MonoProcessorArchitecture;
    }
}
  
// destructor
//

MixedMPISMPArchitecture::~MixedMPISMPArchitecture()
{
  if (this->VerboseModeFlag == true)
    {
      if (this->NbrCPUPerNode[this->MPIRank] > 1)
	{
	  char* TmpString = this->LocalArchitecture->DumpLog("detailed log :\n", "-------------------------");
	  this->AddToLog(TmpString);
	  delete[] TmpString;
	}
      else
	{
	  this->AddToLog("detailed log : none\n-------------------------");
	}
    }
  delete[] this->NbrCPUPerNode;
  delete[] this->ClusterMemoryArray;
  if (this->NodeHostnames != 0)
    {
      for (int i = 0; i < this->NbrMPINodes; ++i)
	delete[] this->NodeHostnames[i];
      delete[] this->NodeHostnames;
    }
}
  
