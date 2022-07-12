////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Antoine Sterdyniak                //
//                                                                            //
//                                                                            //
//     class of U(1) states symmetrization operation for the torus geometry   //
//                                                                            //
//                        last modification : 11/11/2015                      //
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


#ifndef FQHETORUSSYMMETRIZEU1U1STATEOPERATION_H
#define FQHETORUSSYMMETRIZEU1U1STATEOPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/ParticleOnTorus.h"




class FQHETorusSymmetrizeU1U1StateOperation: public AbstractArchitectureOperation
{

 protected:
  
  // pointer to the target Hilbert space
  ParticleOnTorus * FinalSpace;
	
  // pointer to the Hilbert space for the left state
  ParticleOnTorus * LeftSpace;
  // pointer to the Hilbert space for the right state	
  ParticleOnTorus * RightSpace;
  
  // left vector 
  RealVector* LeftVector;
  
  // right vector
  RealVector* RightVector;
  
  // vector where the result has to be stored
  RealVector* DestinationVector;
  
  // left vector  (complex version)
  ComplexVector* ComplexLeftVector;
  
  // right vector (complex version)
  ComplexVector* ComplexRightVector;
  
  // vector where the result has to be stored (complex version)
  ComplexVector* ComplexDestinationVector;
  
  // index of the first component
  long FirstComponent;
  // number of components to process  
  long NbrComponent;
  
  
 public:
  
  // constructor 
  //
  // finalSpace = pointer to the Hilbert space of the target space
  // leftSpace = pointer to the Hilbert space of the first state
  // rightSpace = pointer to the Hilbert space of the second state
  // destinationVector = vector where the result has to be stored
  // leftVector = vector that contains the first state
  // rightVector = vector that contains the second state
  FQHETorusSymmetrizeU1U1StateOperation(ParticleOnTorus* finalSpace, ParticleOnTorus* leftSpace , ParticleOnTorus* rightSpace , 
					RealVector* destinationVector, RealVector* leftVector, RealVector* rightVector);
  
  
  // constructor for complex vector input
  //
  // finalSpace = pointer to the Hilbert space of the target space
  // leftSpace = pointer to the Hilbert space of the first state
  // rightSpace = pointer to the Hilbert space of the second state
  // destinationVector = vector where the result has to be stored
  // leftVector = vector that contains the first state
  // rightVector = vector that contains the second state
  FQHETorusSymmetrizeU1U1StateOperation(ParticleOnTorus* finalSpace, ParticleOnTorus* leftSpace , ParticleOnTorus* rightSpace, 
					ComplexVector* destinationVector, ComplexVector* leftVector, ComplexVector* rightVector);
  
  
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHETorusSymmetrizeU1U1StateOperation(const FQHETorusSymmetrizeU1U1StateOperation & operation);
  
  // constructor from a master node information
  //
  // Space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  FQHETorusSymmetrizeU1U1StateOperation(ParticleOnTorus* space,  SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~FQHETorusSymmetrizeU1U1StateOperation();
  
  
  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetDestinationVector (RealVector* DestinationVector);
  
  // get destination vector 
  // 
  // return value = pointer to destination vector
  Vector* GetDestinationVector ();
  
  //set FirstCoimponent
  //
  // firstComponent  
  void SetIndicesRange (const long& firstComponent, const long& nbrComponent);
  
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
  
  
};

// get destination vector 
// 
// return value = pointer to destination vector

inline Vector* FQHETorusSymmetrizeU1U1StateOperation::GetDestinationVector()
{
  return this->DestinationVector;
}

#endif
