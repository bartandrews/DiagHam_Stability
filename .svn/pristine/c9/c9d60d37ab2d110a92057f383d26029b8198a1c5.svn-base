////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of FQHE particle entanglement spectrum parallelization operation    //
//                                                                            //
//                        last modification : 15/12/2010                      //
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


#ifndef FQHESPHEREPARTICLEENTANGLEMENTSPECTRUMOPERATION_H
#define FQHESPHEREPARTICLEENTANGLEMENTSPECTRUMOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Vector/ComplexVector.h"

class AbstractFunctionBasis;
class RealVector;
class RealSymmetricMatrix;
class HermitianMatrix;


class FQHESphereParticleEntanglementSpectrumOperation: public AbstractPrecalculationOperation
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
  // reduced density matrix where result is stored
  RealSymmetricMatrix DensityMatrix;
  // total system ground state (complex version)
  ComplexVector ComplexGroundState;
  // reduced density matrix where result is stored (complex version)
  HermitianMatrix ComplexDensityMatrix;
  // upper bound on the number of non zero matrix element in the reduced density matrix
  long NbrNonZeroElements;
  // number of projectors
  int NbrGroundStates;
  // array of degenerate groundstates associated to each projector
  ComplexVector* ComplexGroundStates;
  // array of weights in front of each projector
  double* GroundStateWeights;


  // pointer to the array where the top part coefficients are stored
  double* IncompleteBetaThetaTop;
  // pointer on the pointer to the array where the bottom part coefficients are stored
  double* IncompleteBetaThetaBottom;
  // The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
  double PhiRange;

  // a temporary array to store copies of operations in SMP mode
  FQHESphereParticleEntanglementSpectrumOperation** LocalOperations;
  // number of operation copies
  int NbrLocalOperations;

 public:
  
  // constructor 
  //
  // fullSpace = pointer to the full Hilbert space to use
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  FQHESphereParticleEntanglementSpectrumOperation(ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, ParticleOnSphere* complementarySpace, RealVector& groundState, RealSymmetricMatrix& densityMatrix);

  // constructor 
  //
  // fullSpace = pointer to the full Hilbert space to use
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  FQHESphereParticleEntanglementSpectrumOperation(ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, ParticleOnSphere* complementarySpace, ComplexVector& groundState, HermitianMatrix& densityMatrix);

  // constructor when using a sum of projectors
  //
  // fullSpace = pointer to the full Hilbert space to use
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // densityMatrix = reference on the density matrix where result has to stored
  FQHESphereParticleEntanglementSpectrumOperation(ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, ParticleOnSphere* complementarySpace, int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix& densityMatrix);

  // constructor 
  //
  // fullSpace = pointer to the full Hilbert space to use
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // incompleteBetaThetaTop = pointer to the array where the top part coefficients are stored
  // incompleteBetaThetaBotton = pointer on the pointer to the array where the bottom part coefficients are stored
  // phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
  FQHESphereParticleEntanglementSpectrumOperation(ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, ParticleOnSphere* complementarySpace, RealVector& groundState, RealSymmetricMatrix& densityMatrix, double* incompleteBetaThetaBottom, double* incompleteBetaThetaTop, double phiRange);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereParticleEntanglementSpectrumOperation(const FQHESphereParticleEntanglementSpectrumOperation& operation);
  
  // destructor
  //
  ~FQHESphereParticleEntanglementSpectrumOperation();
  
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

inline int FQHESphereParticleEntanglementSpectrumOperation::GetHilbertSpaceDimension ()
{
  return this->ComplementaryHilbertSpace->GetHilbertSpaceDimension();
}

// upper bound on the number of non zero matrix element in the reduced density matrix
//
// return value = upper bound

inline long FQHESphereParticleEntanglementSpectrumOperation::GetNbrNonZeroMatrixElements()
{
  return this->NbrNonZeroElements;
}

#endif
