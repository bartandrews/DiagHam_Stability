////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of operations that perform the Van der Monde            //
//                      times spinful Slater multiplication                   //
//                                                                            //
//                        last modification : 01/12/2016                      //
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


#ifndef FQHESPHEREWITHSU2SPINVANDERMODETIMESSLATEROPERATION_H
#define FQHESPHEREWITHSU2SPINVANDERMODETIMESSLATEROPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"


class FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation : public AbstractArchitectureOperation
{

 protected:
  
  // use reverse flux attachment
  bool ReverseFluxAttachment;

  // vector where the produced state will be stored
  RealVector OutputState;

  // pointer to the Hilbert space
  BosonOnSphereWithSU2Spin* Space;

  // slaterUp = monomial representation of the Slater spin up part
  unsigned long* SlaterUp;
  // slaterDown = monomial representation of the Slater spin up part
  unsigned long* SlaterDown;

  // slaterUp = monomial representation of the second Landau part of the Slater spin up part
  unsigned long* Slater2LLUp;
  // slaterDown = monomial representation of the second Landau part of the Slater spin up part
  unsigned long* Slater2LLDown;

  // number of spin up bosons in the lowest Landau level 
  int NbrBosonsLLLUp;
  // number of spin down bosons in the lowest Landau level
  int NbrBosonsLLLDown;

  // array where the integrals of the three orbital product are stored
  double*** ThreeOrbitalOverlaps;

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;
  
 public:
  
  // constructor 
  //
  // space = pointer to the Hilbert space
  // reverseFluxAttachment = use reverse flux attachment
  // slaterUp = monomial representation of the Slater spin up part
  // slaterDown = monomial representation of the Slater spin up part
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
  FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation(BosonOnSphereWithSU2Spin* space, bool reverseFluxAttachment, unsigned long* slaterUp, unsigned long* slaterDown, 
						       double** threeOrbitalOverlaps);
    
  // constructor 
  //
  // space = pointer to the Hilbert space
  // reverseFluxAttachment = use reverse flux attachment
  // slaterLLLUp = monomial representation of the lowest Landau part of the Slater spin up part
  // slater2LLUp = monomial representation of the second Landau part of the Slater spin up part
  // slaterLLLDown = monomial representation of the lowest Landau part  of the Slater spin down part
  // slater2LLDown = monomial representation of the second Landau part of the Slater spin down part
  // nbrBosonsLLLUp = number of spin up bosons in the lowest Landau level
  // nbrBosonsLLLDown = number of spin down bosons in the lowest Landau level
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored  
  FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation(BosonOnSphereWithSU2Spin* space, bool reverseFluxAttachment, 
						       unsigned long* slaterLLLUp, unsigned long* slater2LLUp, 
						       unsigned long* slaterLLLDown, unsigned long* slater2LLDown, 
						       int nbrBosonsLLLUp, int nbrBosonsLLLDown, double*** threeOrbitalOverlaps);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation(const FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation & operation);
  
  // destructor
  //
  ~FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation();
  
  // set destination vector 
  // 
  // destinationVector = new destination vector
  // set destination vector 
  void SetDestinationVector (RealVector destinationVector);
  
  // get destination vector 
  // 
  // return value = pointer to destination vector
  RealVector GetDestinationVector ();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of components
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

inline RealVector FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::GetDestinationVector()
{
  return this->OutputState;
}

#endif
