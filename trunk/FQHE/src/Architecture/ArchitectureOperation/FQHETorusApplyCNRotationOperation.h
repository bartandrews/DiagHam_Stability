////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of operations that apply a Cn rotation             //
//                                                                            //
//                        last modification : 18/05/2014                      //
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


#ifndef FQHETORUSAPPLYCNROTATIONOPERATION_H
#define FQHETORUSAPPLYCNROTATIONOPERATION_H


#include "config.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"


class FQHETorusApplyCNRotationOperation : public AbstractArchitectureOperation
{

 protected:
  
  // N Value of the rotation (negative if clockwise)
  int NValue;

  // pointer to the Hilbert space of the initial state
  ParticleOnTorus* InputSpace;
  // pointer to the Hilbert space of the rotated state	
  ParticleOnTorus* OutputSpace;

  // pointer to the Hilbert space of the initial state
  ParticleOnTorusWithMagneticTranslations* InputSpaceWithMagneticTranslations;
  // pointer to the Hilbert space of the rotated state	
  ParticleOnTorusWithMagneticTranslations* OutputSpaceWithMagneticTranslations;

  // vector where the initial state is stored
  ComplexVector* InputState;
  // vector where the rotated state is stored
  ComplexVector* OutputState;
  
  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;
  
 public:
  
  // constructor 
  //
  // nValue = N Value of the rotation (negative if clockwise)
  // inputState = pointer to the Hilbert space of the initial state
  // outputState = pointer to the Hilbert space of the rotated state
  // inputSpace = vector where the initial state is stored
  // outputSpace = vector where the rotated state is stored
  FQHETorusApplyCNRotationOperation(int nValue, ComplexVector* inputState, ComplexVector* outputState, ParticleOnTorus* inputSpace, ParticleOnTorus* outputSpace);
    
  // constructor 
  //
  // nValue = N Value of the rotation (negative if clockwise)
  // inputState = pointer to the Hilbert space of the initial state
  // outputState = pointer to the Hilbert space of the rotated state
  // inputSpace = vector where the initial state is stored
  // outputSpace = vector where the rotated state is stored
  FQHETorusApplyCNRotationOperation(int nValue, ComplexVector* inputState, ComplexVector* outputState, ParticleOnTorusWithMagneticTranslations* inputSpace, ParticleOnTorusWithMagneticTranslations* outputSpace);
    
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHETorusApplyCNRotationOperation(const FQHETorusApplyCNRotationOperation & operation);
  
  // destructor
  //
  ~FQHETorusApplyCNRotationOperation();
  
  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetDestinationVector (ComplexVector* DestinationVector);
  
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

inline Vector* FQHETorusApplyCNRotationOperation::GetDestinationVector()
{
  return this->OutputState;
}

#endif
