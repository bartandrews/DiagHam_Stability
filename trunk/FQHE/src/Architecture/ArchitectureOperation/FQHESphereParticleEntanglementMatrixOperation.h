////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of FQHE particle entanglement matrix parallelization operation    //
//                                                                            //
//                        last modification : 23/12/2016                      //
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


#ifndef FQHESPHEREPARTICLEENTANGLEMENTMATRIXOPERATION_H
#define FQHESPHEREPARTICLEENTANGLEMENTMATRIXOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Vector/ComplexVector.h"


class RealVector;
class RealMatrix;


class FQHESphereParticleEntanglementMatrixOperation: public AbstractPrecalculationOperation
{

 protected:

  // pointer to the full Hilbert space to use
  ParticleOnSphere* FullSpace;
  // pointer to the destination Hilbert space (i.e. part A)
  ParticleOnSphere* DestinationHilbertSpace;
  // pointer to the complementary Hilbert space (i.e. part B)
  ParticleOnSphere* ComplementaryHilbertSpace;

  // total system ground state
  RealVector GroundState;

  // entanglement matrix where result is stored
  RealMatrix EntanglementMatrix;

  // remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  bool RemoveBinomialCoefficientFlag;

  // upper bound on the number of non zero matrix element in the reduced density matrix
  long NbrNonZeroElements;

  // a temporary array to store copies of operations in SMP mode
  FQHESphereParticleEntanglementMatrixOperation** LocalOperations;
  // number of operation copies
  int NbrLocalOperations;

 public:
  
  // constructor 
  //
  // fullSpace = pointer to the full Hilbert space to use
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // groundState = reference on the total system ground state
  // entanglementMatrix = reference on the entanglement matrix where result has to stored
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  FQHESphereParticleEntanglementMatrixOperation(ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, 
						ParticleOnSphere* complementarySpace, RealVector& groundState, RealMatrix& entanglementMatrix,
						bool removeBinomialCoefficient = false);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereParticleEntanglementMatrixOperation(const FQHESphereParticleEntanglementMatrixOperation& operation);
  
  // destructor
  //
  ~FQHESphereParticleEntanglementMatrixOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // get hilbert space dimension
  // 
  // return value = hilbert space dimension  
  int GetHilbertSpaceDimension ();

  // upper bound on the number of non zero matrix element in the reduced density matrix
  //
  // return value = upper bound
  long GetNbrNonZeroMatrixElements();

  // get the reference on the reduced density matrix
  // 
  // return value = eference on the reduced density matrix
  HermitianMatrix& GetMatrix();


 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

// get hilbert space dimension
// 
// return value = hilbert space dimension

inline int FQHESphereParticleEntanglementMatrixOperation::GetHilbertSpaceDimension ()
{
  return this->ComplementaryHilbertSpace->GetHilbertSpaceDimension();
}

// upper bound on the number of non zero matrix element in the reduced density matrix
//
// return value = upper bound

inline long FQHESphereParticleEntanglementMatrixOperation::GetNbrNonZeroMatrixElements()
{
  return this->NbrNonZeroElements;
}


#endif
