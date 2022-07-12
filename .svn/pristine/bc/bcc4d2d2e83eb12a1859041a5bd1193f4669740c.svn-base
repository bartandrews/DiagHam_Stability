////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of Abstract architecture operation                 //
//                                                                            //
//                        last modification : 23/10/2002                      //
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
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/MixedMPISMPArchitecture.h"

#include <iostream>


using std::cout;
using std::endl;



// destructor
//

AbstractArchitectureOperation::~AbstractArchitectureOperation()
{
}

// apply operation for a given architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool AbstractArchitectureOperation::ApplyOperation(AbstractArchitecture* architecture)
{
  switch (architecture->GetArchitectureID())
    {
    case AbstractArchitecture::MonoProcessor:
      return this->ArchitectureDependentApplyOperation((MonoProcessorArchitecture*) architecture);
    case AbstractArchitecture::SMP:
      return this->ArchitectureDependentApplyOperation((SMPArchitecture*) architecture);
    case AbstractArchitecture::SimpleMPI:
      return this->ArchitectureDependentApplyOperation((SimpleMPIArchitecture*) architecture);
    case AbstractArchitecture::MixedMPISMP:      
      return this->ArchitectureDependentApplyOperation((SimpleMPIArchitecture*) architecture);
    default:
      return this->RawApplyOperation();
    }
}
  
// apply operation for mono processor architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool  AbstractArchitectureOperation::ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture)
{
  return this->RawApplyOperation();
}
  
// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool AbstractArchitectureOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  return this->RawApplyOperation();
}
 
// apply operation for simple MPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool AbstractArchitectureOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
  return this->RawApplyOperation();
}

// apply an SMP round robin operation 
//
// return value = true if no error occurs
  
bool AbstractArchitectureOperation::ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID)
{
  return this->RawApplyOperation();
}
  
