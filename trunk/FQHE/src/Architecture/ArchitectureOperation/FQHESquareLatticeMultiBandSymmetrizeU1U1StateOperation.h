////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Antoine Sterdyniak                //
//                                                                            //
//                                                                            //
//                   class of U1U1 states symmetrization Operation	      //
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


#ifndef FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation_H
#define FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation_H


#include "config.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMultiBandMomentumSpace.h"


class ComplexVector;


class FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation: public AbstractArchitectureOperation
{

 protected:
  
  // pointer to the Hilbert spaces
  BosonOnSquareLatticeMultiBandMomentumSpace* FinalSpace;
  BosonOnSquareLatticeMultiBandMomentumSpace* LeftSpace;
  BosonOnSquareLatticeMultiBandMomentumSpace* RightSpace;
  
  // pointer to the Hilbert spaces (long version)
  BosonOnSquareLatticeMomentumSpaceLong* FinalSpaceLong;
  BosonOnSquareLatticeMomentumSpaceLong* LeftSpaceLong;
  BosonOnSquareLatticeMomentumSpaceLong* RightSpaceLong;
  
  // vector where the component are stored
  ComplexVector* LeftVector;
  
  // vector where the component are stored
  ComplexVector* RightVector;
  
  // vector where the result has to be stored
  ComplexVector* DestinationVector;
  
  // index of the first component
  long FirstComponent;
  
  long NbrComponent;
  
  bool UnnormalizedBasisFlag;
  
 public:
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // sourceVector = vector where the component are stored
  // destinationVector = vector where the result has to be stored
  FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation( BosonOnSquareLatticeMultiBandMomentumSpace * finalSpace, BosonOnSquareLatticeMultiBandMomentumSpace * leftSpace , BosonOnSquareLatticeMultiBandMomentumSpace * rightSpace , ComplexVector* destinationVector, ComplexVector* leftVector, ComplexVector* rightVector,  bool unnormalizedBasisFlag);
  
  
  // constructor using long Hilbert spaces
  //
  // Space = pointer to the HilbertSpace to use
  // sourceVector = vector where the component are stored
  // destinationVector = vector where the result has to be stored
  FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation( BosonOnSquareLatticeMomentumSpaceLong* finalSpace, BosonOnSquareLatticeMomentumSpaceLong* leftSpace , BosonOnSquareLatticeMomentumSpaceLong* rightSpace , ComplexVector* destinationVector, ComplexVector* leftVector, ComplexVector* rightVector,  bool unnormalizedBasisFlag);
  
  
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation(const FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation & operation);
  
  // constructor from a master node information
  //
  // Space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation( BosonOnSquareLatticeMultiBandMomentumSpace* space,  SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation();
  
  
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

inline Vector* FQHESquareLatticeMultiBandSymmetrizeU1U1StateOperation::GetDestinationVector()
{
  return this->DestinationVector;
}

#endif
