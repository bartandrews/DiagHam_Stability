////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of FQHE particle entanglement spectrum parallelization operation    //
//                    for torus with magnetic translations                    //
//                                                                            //
//                        last modification : 13/06/2011                      //
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


#ifndef FQHETORUSPARTICLEENTANGLEMENTSPECTRUMOPERATION_H
#define FQHETORUSPARTICLEENTANGLEMENTSPECTRUMOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Vector/ComplexVector.h"


class AbstractFunctionBasis;
class RealVector;
class RealSymmetricMatrix;
class HermitianMatrix;


class FQHETorusParticleEntanglementSpectrumOperation: public AbstractPrecalculationOperation
{

 protected:

  // pointer to the full Hilbert space to use
  ParticleOnTorusWithMagneticTranslations* FullSpace;
  // pointer to the destination Hilbert space (i.e. part A)
  ParticleOnTorusWithMagneticTranslations* DestinationHilbertSpace;
  // pointer to the complementary Hilbert space (i.e. part B)
  ParticleOnTorusWithMagneticTranslations* ComplementaryHilbertSpace;
  // pointer to the full spinful Hilbert space to use
  ParticleOnTorusWithSpinAndMagneticTranslations* SpinfulFullSpace;
  // pointer to the destination spinful Hilbert space (i.e. part A)
  ParticleOnTorusWithSpinAndMagneticTranslations* SpinfulDestinationHilbertSpace;
  // pointer to the complementary spinful Hilbert space (i.e. part B)
  ParticleOnTorusWithSpinAndMagneticTranslations* SpinfulComplementaryHilbertSpace;
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

  // a temporary array to store copies of operations in SMP mode
  FQHETorusParticleEntanglementSpectrumOperation** LocalOperations;
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
  FQHETorusParticleEntanglementSpectrumOperation(ParticleOnTorusWithMagneticTranslations* fullSpace, ParticleOnTorusWithMagneticTranslations* destinationSpace, 
						 ParticleOnTorusWithMagneticTranslations* complementarySpace, ComplexVector& groundState, HermitianMatrix& densityMatrix);

  // constructor for the spinful case
  //
  // fullSpace = pointer to the full Hilbert space to use
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  FQHETorusParticleEntanglementSpectrumOperation(ParticleOnTorusWithSpinAndMagneticTranslations* fullSpace, ParticleOnTorusWithSpinAndMagneticTranslations* destinationSpace, 
						 ParticleOnTorusWithSpinAndMagneticTranslations* complementarySpace, ComplexVector& groundState, HermitianMatrix& densityMatrix);

  // constructor  when using a sum of projectors
  //
  // fullSpace = pointer to the full Hilbert space to use
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // densityMatrix = reference on the density matrix where result has to stored
  FQHETorusParticleEntanglementSpectrumOperation(ParticleOnTorusWithMagneticTranslations* fullSpace, ParticleOnTorusWithMagneticTranslations* destinationSpace, 
						 ParticleOnTorusWithMagneticTranslations* complementarySpace, 
						 int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix& densityMatrix);

  // constructor for the spinful case when using a sum of projectors
  //
  // fullSpace = pointer to the full Hilbert space to use
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // densityMatrix = reference on the density matrix where result has to stored
  FQHETorusParticleEntanglementSpectrumOperation(ParticleOnTorusWithSpinAndMagneticTranslations* fullSpace, ParticleOnTorusWithSpinAndMagneticTranslations* destinationSpace, 
						 ParticleOnTorusWithSpinAndMagneticTranslations* complementarySpace, 
						 int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix& densityMatrix);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHETorusParticleEntanglementSpectrumOperation(const FQHETorusParticleEntanglementSpectrumOperation& operation);
  
  // destructor
  //
  ~FQHETorusParticleEntanglementSpectrumOperation();
  
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
  
  
  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

// get hilbert space dimension
// 
// return value = hilbert space dimension

inline int FQHETorusParticleEntanglementSpectrumOperation::GetHilbertSpaceDimension ()
{
  if (this->SpinfulComplementaryHilbertSpace == 0)
    return this->ComplementaryHilbertSpace->GetHilbertSpaceDimension();
  else
    return this->SpinfulComplementaryHilbertSpace->GetHilbertSpaceDimension();
}

// upper bound on the number of non zero matrix element in the reduced density matrix
//
// return value = upper bound

inline long FQHETorusParticleEntanglementSpectrumOperation::GetNbrNonZeroMatrixElements()
{
  return this->NbrNonZeroElements;
}

// get the reference on the reduced density matrix
// 
// return value = eference on the reduced density matrix

inline HermitianMatrix& FQHETorusParticleEntanglementSpectrumOperation::GetMatrix()
{
  return this->ComplexDensityMatrix;
}

#endif
