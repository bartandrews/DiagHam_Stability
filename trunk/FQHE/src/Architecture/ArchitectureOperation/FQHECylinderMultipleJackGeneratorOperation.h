////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//               class of mulitple Jack generator basis to create an	      //
//                     orthonomal basis on the cylinder geometry              //
//                                                                            //
//                        last modification : 21/07/2016                      //
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


#ifndef FHQECYLINDERMULTIPLEJACKGENERATOROPERATION_H
#define FHQECYLINDERMULTIPLEJACKGENERATOROPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"

#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"



class RealVector;
class ParticleOnSphere;


class FQHECylinderMultipleJackGeneratorOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  long FirstComponent;
  // number of component to compute 
  long NbrComponent;

  //  the non-orthogonalized Jack polynomial basis
  RealVector* QuasiholeVectors;

 // number of part in the CFT calculation will be separated in MPI mode
  int NbrMPIStage;
  // number of part in the CFT calculation will be separated in SMP mode
  int NbrSMPStage;
  // array with size of SMP stages used to distribute work
  int* SMPStages; 


  // cylinder aspect ratio
  double CylinderRatio;

  // Hilbert space where all the Jack states have to be expressed
  ParticleOnSphere* TargetSpace;
  // number of particles 
  int NbrParticles;
  // number of flux quanta
  int NbrFluxQuanta;  
  // system total momentum
  int TotalLz;
  // array that contains all the root configurations
  unsigned long** RootConfigurations;
  // number of root configurations
  int NbrRootConfigurations;
  // k value of the clustered (k,r) principle
  int KValue;
  // r value of the clustered (k,r) principle
  int RValue;

 public:
  
  // constructor 
  //
  // targetSpace = pointer to the Hilbert where the quasihole states should be expressed
  // rootConfigurations = array that contains all the root configurations
  // nbrRootConfigurations = number of root configurations
  // kValue = k value of the clustered (k,r) principle
  // rValue = r value of the clustered (k,r) principle
  // nbrParticles = number of particles
  // lzMax = number of flux quantum
  // totalLz = system total momentum
  // ratio = cylinder aspect ratio
  // nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
  // nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode
  FQHECylinderMultipleJackGeneratorOperation(ParticleOnSphere* targetSpace, unsigned long** rootConfigurations, int nbrRootConfigurations, 
					     int kValue, int rValue, int nbrParticles, int lzMax, int totalLz, double ratio, int nbrMPIStage = 10, int nbrSMPStage = 10);
  

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHECylinderMultipleJackGeneratorOperation(const FQHECylinderMultipleJackGeneratorOperation & operation);
  
  // constructor from a master node information
  //
  // space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  //  FQHECylinderMultipleJackGeneratorOperation(ParticleOnSphere* space,  SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~FQHECylinderMultipleJackGeneratorOperation();  
  
  // get the rthogonalized Jack polynomial basis
  // 
  // return value = reference on the matrix that contains the basis
  virtual RealMatrix GetBasis ();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  virtual void SetIndicesRange (const long& firstComponent, const long& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  virtual AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  virtual bool RawApplyOperation();
  
  // apply operation for SMP using round robin scheduling
  //
  //  architecture = instance of architecture class
  // return value = true if no error occurs
  virtual bool ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID);

 protected:
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  

};

// get the orthogonalized Jack polynomial basis
// 
// return value = reference on the matrix that contains the basis

inline RealMatrix FQHECylinderMultipleJackGeneratorOperation::GetBasis ()
{
  RealVector* TmpVectors = new RealVector[this->NbrRootConfigurations];	
  for (int i =0; i < this->NbrRootConfigurations; ++i)
    {
      TmpVectors[i] = this->QuasiholeVectors[i];
    }
  RealMatrix TmpMatrix(TmpVectors, this->NbrRootConfigurations); 
  TmpMatrix.OrthoNormalizeColumns();
  return TmpMatrix;
}

#endif
