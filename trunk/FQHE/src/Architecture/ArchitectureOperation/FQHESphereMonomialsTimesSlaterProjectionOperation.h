////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                 Copyright (C) 2001-2002 Antoine Sterdyniak                 //
//                                                                            //
//                                                                            //
//             class of Kostka number computation operation                   //
//                                                                            //
//                        last modification : 09/05/2010                      //
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


#ifndef FQHESPHEREMONOMIALSTIMESSLATERPROJECTIONOPERATION_H
#define FQHESPHEREMONOMIALSTIMESSLATERPROJECTIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphereThreeLandauLevels.h"




class FQHESphereMonomialsTimesSlaterProjectionOperation : public AbstractArchitectureOperation
{
  
 protected:
  
  // index of the first component
  int FirstComponent;
  
  // number of component 
  int NbrComponent;
  
  // pointer to the HilbertSpace after product and projection if there is one
  ParticleOnSphere * FinalSpace;
  
  // pointer to the initial Fermionic HilbertSpace
  ParticleOnSphere * FermionSpace;
  
  // pointer to the initial lowest Landau level HilbertSpace
  ParticleOnSphere * LLLSpace;
  
  // RealVector where the result will be store
  RealVector * OutputRealVector;
  
  // RealVector where the fermionic state is stored
  RealVector * FermionRealVector;
  
  // RealVector where the bosonic state is stored
  RealVector * LLLRealVector;
  
  // RealVector where the result will be store
  LongRationalVector * OutputLongRationalVector;
  
  // RealVector where the fermionic state is stored
  LongRationalVector * FermionLongRationalVector;
  
  // RealVector where the bosonic state is stored
  LongRationalVector * LLLLongRationalVector;
  
  //true if the result vector is to be projected
  bool Projection;
  
  //true if the result vector is to be normalized
  bool Normalize;
  
  // number of Landau levels
  int NbrLL;
  
  //true if the lz->-lz symmetry is used
  bool Symmetry;
  
  //true if the initial lll HilbertSpace is bosonic
  bool BosonFlag;
  
  // number of part in which the initial bosonic vector will be separated
  int NbrStage;
  
  // true if the mode rational is used
  bool Rational;

  // true if the flux are in the opposite direction to the magnetic field
  bool ReverseFluxFlag;

 public:
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // fileName = name of the file where the kostka number will be store
  // nbrLL = number of Landau levels
  FQHESphereMonomialsTimesSlaterProjectionOperation(ParticleOnSphere* fermionSpace, ParticleOnSphere* lllSpace, ParticleOnSphere* finalSpace, RealVector* fermionVector, RealVector* lllVector, RealVector* outputVector, int resume, int nbrComponent, bool projection, int step, int nbrLL, bool symmetry, bool reverseFluxFlag);
  
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // fileName = name of the file where the kostka number will be store
  // nbrLL = number of Landau levels
  FQHESphereMonomialsTimesSlaterProjectionOperation(ParticleOnSphere* fermionSpace, ParticleOnSphere* lllSpace, ParticleOnSphere* finalSpace, LongRationalVector* fermionVector, LongRationalVector* lllVector, LongRationalVector* outputVector, int resume, int nbrComponent, bool projection, int step, int nbrLL, bool symmetry, bool reverseFluxFlag);
  
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereMonomialsTimesSlaterProjectionOperation(const FQHESphereMonomialsTimesSlaterProjectionOperation & operation);
  
  // destructor
  //
  ~FQHESphereMonomialsTimesSlaterProjectionOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetOutputVector (RealVector* outputVector);
  void SetOutputVector (LongRationalVector* outputVector);
  
  // get destination vector 
  // 
  // return value = pointer to destination vector
  Vector* GetDestinationVector ();
  
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
  
};


#endif
