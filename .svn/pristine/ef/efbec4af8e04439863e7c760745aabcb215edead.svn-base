////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of operations that apply a one body               //
//                        transformation to a many body state                 //
//                                                                            //
//                        last modification : 24/11/2015                      //
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


#ifndef FQHESPHEREWITHSPINAPPLYONEBODYTRANSFORMATIONOPERATION_H
#define FQHESPHEREWITHSPINAPPLYONEBODYTRANSFORMATIONOPERATION_H


#include "config.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"


class FQHESphereWithSpinApplyOneBodyTransformationOperation : public AbstractArchitectureOperation
{

 protected:
  
  // pointer to the Hilbert space of the initial state
  ParticleOnSphereWithSpin* InputSpace;
  
  // vector where the initial state is stored
  ComplexVector* InputState;
  // vector where the rotated state is stored
  ComplexVector* OutputState;
  
  // matrices describing the one-body tranformation per orbital
  ComplexMatrix* RotationMatrices;

  // vector where the initial state is stored (real case)
  RealVector* RealInputState;
  // vector where the rotated state is stored (real case)
  RealVector* RealOutputState;
  
  // matrices describing the one-body tranformation per orbital (real case)
  RealMatrix* RealRotationMatrices;

  // index of the first component
  long FirstComponent;
  // number of component 
  long NbrComponent;
  
 public:
  
  // constructor 
  //
  // inputState = vector where the initial state is stored
  // outputState = vector where the rotated state is stored
  // rotationMatrices =  matrices describing the one-body tranformation per orbital
  // inputSpace = pointer to the Hilbert space
  FQHESphereWithSpinApplyOneBodyTransformationOperation(RealVector* inputState, RealVector* outputState, RealMatrix* rotationMatrices,
							ParticleOnSphereWithSpin* inputSpace);
    
  // constructor 
  //
  // inputState = vector where the initial state is stored
  // outputState = vector where the rotated state is stored
  // rotationMatrices =  matrices describing the one-body tranformation per orbital
  // inputSpace = pointer to the Hilbert space
  FQHESphereWithSpinApplyOneBodyTransformationOperation(ComplexVector* inputState, ComplexVector* outputState, ComplexMatrix* rotationMatrices,
							ParticleOnSphereWithSpin* inputSpace);
    
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereWithSpinApplyOneBodyTransformationOperation(const FQHESphereWithSpinApplyOneBodyTransformationOperation & operation);
  
  // destructor
  //
  ~FQHESphereWithSpinApplyOneBodyTransformationOperation();
  
  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetDestinationVector (ComplexVector* DestinationVector);
  
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

inline Vector* FQHESphereWithSpinApplyOneBodyTransformationOperation::GetDestinationVector()
{
  return this->OutputState;
}

#endif
