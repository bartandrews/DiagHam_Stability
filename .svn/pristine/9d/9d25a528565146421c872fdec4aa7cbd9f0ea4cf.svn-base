////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of Architecture Manager                      //
//                                                                            //
//                        last modification : 28/05/2004                      //
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
#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Architecture/MixedMPISMPArchitecture.h"
#include <iostream>

// default constructor
//

ArchitectureManager::ArchitectureManager()
{
  this->Architecture = 0;
  this->Options = 0;
}

// destructor
//

ArchitectureManager::~ArchitectureManager()
{
  if (this->Architecture != 0)
    {
      delete this->Architecture;
    }
}

// add an option group containing all options related to the architecture 
//
// manager = pointer to the option manager

void ArchitectureManager::AddOptionGroup(OptionManager* manager)
{
  this->Options = manager;
  OptionGroup* ParallelizationGroup  = new OptionGroup ("parallelization options");
  (*(this->Options)) += ParallelizationGroup;
#ifdef __SMP__
  (*ParallelizationGroup) += new BooleanOption  ('S', "SMP", "enable SMP mode");
  (*ParallelizationGroup) += new SingleIntegerOption  ('\n', "processors", "number of processors to use in SMP mode", 2);
  (*ParallelizationGroup) += new SingleStringOption  ('\n', "smp-profil", "enable SMP profiling, the name of the log file  has to be passed as argument");
#endif
#ifdef __MPI__
  (*ParallelizationGroup) += new BooleanOption  ('\n', "mpi", "enable MPI mode");  
#ifdef __SMP__
  (*ParallelizationGroup) += new SingleStringOption ('\n', "mpi-smp", "enable both MPI and SMP mode, the name file describing the cluster has to be passed as argument");
#endif
  (*ParallelizationGroup) += new SingleStringOption  ('\n', "cluster-profil", "enable cluster profiling, the name of the log file  has to be passed as argument");  
  (*ParallelizationGroup) += new BooleanOption ('\n', "auto-loadbalancing", "use automatic load balancing, overriding any manual load balancing");
#endif
  
}

// get the best avalaible architecture in agrement with the option constraints
//
// return value = pointer to the architecture

AbstractArchitecture* ArchitectureManager::GetArchitecture()
{
  if ((this->Options != 0) && (this->Architecture == 0))
    {
#ifdef __SMP__
      bool SMPFlag = this->Options->GetBoolean("SMP");
#else 
      bool SMPFlag = false;
#endif
#ifdef __MPI__
      bool MPIFlag = this->Options->GetBoolean("mpi");
      char* MPILogFile = this->Options->GetString("cluster-profil");
      bool AutomaticLoadBalancing = this->Options->GetBoolean("auto-loadbalancing");
#ifdef __SMP__
      if (this->Options->GetString("mpi-smp") != 0)
	{
	  SMPFlag = true;
	  MPIFlag = true;
	}
      else
	if (MPIFlag == true)
	  SMPFlag = false;
#endif      
#else      
      bool MPIFlag = false;
      char* MPILogFile = 0;
      bool AutomaticLoadBalancing = false;
#endif
      int NbrProcessor = (this->Options->GetInteger("processors"));
      if (SMPFlag == false)
	if (MPIFlag == false)
	  this->Architecture = new MonoProcessorArchitecture;
	else
	  this->Architecture = new SimpleMPIArchitecture(MPILogFile, AutomaticLoadBalancing);
      else
	if (MPIFlag == false)
	  this->Architecture = new SMPArchitecture(NbrProcessor, this->Options->GetString("smp-profil"));
	else
	  this->Architecture = new MixedMPISMPArchitecture(this->Options->GetString("mpi-smp"), MPILogFile, AutomaticLoadBalancing);
    }
  return this->Architecture;
}
