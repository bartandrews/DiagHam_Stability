////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of monomials product operation 	              //
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


#ifndef MONOMIALSPRODUCTOPERATION_H
#define MONOMIALSPRODUCTOPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/BosonOnSphereShort.h"


class BosonOnSphereShort;
class RealVector;


class MonomialsProductOperation: public AbstractArchitectureOperation
{
  
 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the first Hilbert Space
  ParticleOnSphere * Space1;
  
  //pointer to the eventual second Hilbert Space
  ParticleOnSphere * Space2;
  
  //pointer to the final Hilbert Space
  ParticleOnSphere * FinalSpace;

  // pointer to the first vector 
  RealVector* SourceVector1;

  // pointer to the second vector  
  RealVector* SourceVector2;

  // vector where the result has to be stored
  RealVector* DestinationVector;

  // true if only one Hilbert space  
  bool Squaring;
	
  // true if final state is fermionic
  bool FinalStateFlag;
	
  //true is first state is fermionic
  bool FirstStateFlag;
  
  //true if the result state will be normalized at the end of the calculation
  bool NormalizeFlag;
	
	bool ReverseFluxFlag;

 public:
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // sourceVector = vector where the component are stored
  // destinationVector = vector where the result has to be stored 
  MonomialsProductOperation(ParticleOnSphere * space,ParticleOnSphere * FinalSpace, RealVector* sourceVector1,RealVector* sourceVector2, RealVector* destinationVector, bool normalize);
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // sourceVector = vector where the component are stored
  // destinationVector = vector where the result has to be stored 
  MonomialsProductOperation(ParticleOnSphere * space1,ParticleOnSphere * space2,ParticleOnSphere * FinalSpace,RealVector* sourceVector1, RealVector* sourceVector2, RealVector* destinationVector, bool normalize );

  // copy constructor 
  //
  // operation = reference on operation to copy
  MonomialsProductOperation(const MonomialsProductOperation & operation);

  // constructor from a master node information
  //
  // Space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  MonomialsProductOperation(ParticleOnSphere * space,ParticleOnSphere * FinalSpace,  bool normalize, SimpleMPIArchitecture* architecture);
  
  
  MonomialsProductOperation(ParticleOnSphere * space1,ParticleOnSphere * space2,ParticleOnSphere * finalSpace,  bool normalize, SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~MonomialsProductOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetDestinationVector (RealVector* DestinationVector);

  // get destination vector 
  // 
  // return value = pointer to destination vector
  Vector* GetDestinationVector ();

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
  
  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
  // evaluate the time needed for the computation of the product for the monomial considered (assuming linear dependency)
  // 
  // index = index of the considered monomial
  // return value = time estimation
  int TimeEvaluationFonction(int index);
  
  //evaluate the time needed for the computation of the product for the monomial considered (assuming square dependency)
  //
  // index = index of the considered monomial
  // return value = time estimation
  unsigned long TimeEvaluationFonctionSquare(int index);
  
};

// get destination vector 
// 
// return value = pointer to destination vector

inline Vector* MonomialsProductOperation::GetDestinationVector()
{
  return this->DestinationVector;
}

inline int MonomialsProductOperation::TimeEvaluationFonction(int index)
{
  return  (this->NbrComponent-index)*(this->NbrComponent-1-index)/2;
}

inline unsigned long  MonomialsProductOperation::TimeEvaluationFonctionSquare(int index)
{
  return  ((unsigned long) (this->NbrComponent-index)*(this->NbrComponent-index-1)*(2*(this->NbrComponent-index)-1)/6);
}

#endif
