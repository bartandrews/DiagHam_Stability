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


#ifndef FQHESPHEREMULTIPLEMONOMIALSTIMESSLATERPROJECTIONOPERATION_H
#define FQHESPHEREMULTIPLEMONOMIALSTIMESSLATERPROJECTIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/ParticleOnSphere.h"


class FQHESphereMultipleMonomialsTimesSlaterProjectionOperation : public AbstractArchitectureOperation
{
  
 protected:
  
  // pointer to the HilbertSpace after product and projection if there is one
  ParticleOnSphere * FinalSpace;
  
  // pointer to the initial Fermionic HilbertSpace
  ParticleOnSphere * FermionSpace;
  
  //Array where the index of the non-zero components of the fermionic states are stored 
  int * MatchingConditionsIndex;
  
  // pointer to the initial lowest Landau level HilbertSpace
  ParticleOnSphere * LLLSpace;
  
  // array of vectors describing the fermionic states
  RealVector * FermionRealVector;
  
  // array of vectors describing the fermionic states
  LongRationalVector * FermionLongRationalVector;
  
  // array of vectors describing the fermionic states
  RealVector * LLLRealVector;
  
  // array of vectors describing the fermionic states
  LongRationalVector * LLLLongRationalVector;
  
  //true if the result vector is to be projected
  bool Projection;
  
  // number of Landau levels
  int NbrLL;
  
  //true if the lz->-lz symmetry is used
  bool Symmetry;
  
  //true if the initial lll HilbertSpace is bosonic
  bool BosonFlag;
  
  // Number of States 
  int NbrStates;
  
  //true if the final state should be normalize (only possible if projection)
  bool Normalize;
  
  //name of the file where the vectors will be stored
  char * OutputFileName;
  
  // index of the state in the n-body basis to be computed
  int Index;
  
  // RealVector where the result will be store
  RealVector * OutputRealVector;
	
  LongRationalVector * OutputLongRationalVector;
  
  // index from wich the computation will be resumed in case of interruption
  int ResumingIndex;
  
  // true if the calculation are made with rational
  bool Rational;
	
  // true if the flux attached are in the opposite direction to the magnetic field
  bool ReverseFluxFlag;
  
 public:
  
  
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // fileName = name of the file where the kostka number will be store
  // nbrLL = number of Landau levels
  FQHESphereMultipleMonomialsTimesSlaterProjectionOperation(ParticleOnSphere* fermionSpace, ParticleOnSphere* lllSpace, ParticleOnSphere* finalSpace, int * matchingConditionsIndex, RealVector* lllVector, int nbrLL, int nbrStates,bool projection, bool symmetry,bool normalize,char * outputFileName, int resumingIndex, bool reverseFluxFlag);
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // fileName = name of the file where the kostka number will be store
  // nbrLL = number of Landau levels
  FQHESphereMultipleMonomialsTimesSlaterProjectionOperation(ParticleOnSphere* fermionSpace, ParticleOnSphere* lllSpace, ParticleOnSphere* finalSpace, int * matchingConditionsIndex, LongRationalVector* lllVector, int nbrLL, int nbrStates,bool projection, bool symmetry,bool normalize,char * outputFileName,int resumingIndex,bool reverseFluxFlag);
  
  // destructor
  //
  ~FQHESphereMultipleMonomialsTimesSlaterProjectionOperation();													
 protected :
  
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereMultipleMonomialsTimesSlaterProjectionOperation(const FQHESphereMultipleMonomialsTimesSlaterProjectionOperation & operation);
  
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // get destination vector 
  // 
  // return value = pointer to destination vector
  Vector* GetDestinationVector ();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
  //change the value of Index
  //
  // index = new Index value
  void SetIndex(int index);
  
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
  
  // apply operation for MonoProcessorArchitecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool  ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture);
  
};


#endif
