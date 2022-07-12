////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Antoine Sterdyniak                //
//                                                                            //
//                                                                            //
//                   class of multiplication by Jastrow Operation	      //
//                                                                            //
//                        last modification : 03/03/2010                      //
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


#ifndef JASTROWMULTIPLICATIONOPERATION_H
#define JASTROWMULTIPLICATIONOPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/BosonOnSphereShort.h"


class BosonOnSphereShort;
class RealVector;


class FQHESphereJastrowMultiplicationOperation: public AbstractArchitectureOperation
{

 protected:
  
  // pointer to the HilbertSpace
  BosonOnSphereShort * Space;
  
  // vector where the component are stored
  RealVector* SourceVector;
  
  
  // vector where the result has to be stored
  RealVector* DestinationVector;
  
  // vector where the monomials decomposition of the Slater determinant that correspond to the FirstComponent is to be stored
  RealVector* JackVector;
  
  
  // index of the first component
  long FirstComponent;
  
  int NbrState;
		
 public:
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // sourceVector = vector where the component are stored
  // destinationVector = vector where the result has to be stored
  FQHESphereJastrowMultiplicationOperation(BosonOnSphereShort* space, RealVector* sourceVector, RealVector* destinationVector,int nbrState);
  
  
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereJastrowMultiplicationOperation(const FQHESphereJastrowMultiplicationOperation & operation);
  
  // constructor from a master node information
  //
  // Space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  FQHESphereJastrowMultiplicationOperation(BosonOnSphereShort* space,  SimpleMPIArchitecture* architecture,int nbrState);
  
  // destructor
  //
  ~FQHESphereJastrowMultiplicationOperation();
  
  
  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetDestinationVector (RealVector* DestinationVector);
  
  // set Jack Vector
  //
  // vectors where the Jack Polynomial will be stored
  void SetJackVector (RealVector* jackVector);
  
  // get destination vector 
  // 
  // return value = pointer to destination vector
  Vector* GetDestinationVector ();
  
  //set FirstCoimponent
  //
  // firstComponent  
  void SetFirstComponent (const long& firstComponent);
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
 protected:
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation for MonoProcessor architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture);
  
};

// get destination vector 
// 
// return value = pointer to destination vector

inline Vector* FQHESphereJastrowMultiplicationOperation::GetDestinationVector()
{
  return this->DestinationVector;
}

#endif
