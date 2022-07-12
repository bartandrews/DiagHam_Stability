////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Antoine Sterdyniak                //
//                                                                            //
//                                                                            //
//                  class of division by a Jastrow operation                  //
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


#ifndef FQHELATTICEFOURIERTRANSFORMOPERATION_H
#define FQHELATTICEFOURIERTRANSFORMOPERATION_H


#include "config.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/BosonOnLatticeKy.h"



class FQHELatticeFourierTransformOperation : public AbstractArchitectureOperation
{
	
 protected:
  
  // pointer to the HilbertSpace
  BosonOnLattice * FullSpace;
  BosonOnLatticeKy * KySpace;
  // vector where the component are stored
  ComplexVector * SourceVector;
  
  // vector where the result has to be stored
  ComplexVector * DestinationVector;
  
  int NbrVectors;
  // number of component
  long NbrComponent;
  
  // index of the first component
  long FirstComponent;
  
  // number of states to handle
  //int NbrStates;
  
 public:
  
  // constructor 
  //
  // space = pointer to the HilbertSpace to use
  // sourceVector = array of vectors describing the fermionic states
  // destinationVector = array of vectors where the resulting bosonic states have to be stored
  // nbrStates = number of states to handle
  FQHELatticeFourierTransformOperation(BosonOnLattice * fullSpace, BosonOnLatticeKy * kySpace, ComplexVector* sourceVector,int nbrComponent, ComplexVector* destinationVector, int nbrVectors);
    
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHELatticeFourierTransformOperation(const FQHELatticeFourierTransformOperation & operation);
  
  // destructor
  //
  ~FQHELatticeFourierTransformOperation();
  
  
  // set destination vector 
  //	 
  // vector where the result has to be stored
  void SetDestinationVector (ComplexVector* DestinationVector);
  
  // set range of indices	
  // 	
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const long& firstComponent, const long& nbrComponent);
  
  // get destination vector 
  // 
  // return value = pointer to destination vector
  Vector* GetDestinationVector();
  
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

inline Vector* FQHELatticeFourierTransformOperation::GetDestinationVector()
{
  return this->DestinationVector;
}


#endif
