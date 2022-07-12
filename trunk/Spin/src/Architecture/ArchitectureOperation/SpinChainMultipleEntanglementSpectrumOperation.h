////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of spin chain multiple entanglement spectrum calculation        //
//                                                                            //
//                        last modification : 06/08/2016                      //
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


#ifndef SPINCHAINMULTIPLEENTANGLEMENTSPECTRUMOPERATION_H
#define SPINCHAINMULTIPLEENTANGLEMENTSPECTRUMOPERATION_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"


class AbstractOperator;


class SpinChainMultipleEntanglementSpectrumOperation: public AbstractArchitectureOperation
{

 protected:

  // pointer to the Hilbert space
  AbstractSpinChain* Space;

  // eigenstates (real version)
  RealMatrix RealEigenstates;
  // eigenstates (complex version)
  ComplexMatrix ComplexEigenstates;

  // number of sites in the subsystem
  int SubsystemSize;
  // array that describes the subsystem (if null, consider the subsystemSize first sites)
  int* SubsystemSites;
  // twice the total Sz value for the subsystem
  int SubsystemSz; 

  // index of the first state to evaluate 
  int FirstState;
  // number of states that have to be evaluated by the local architecture
  int NbrStates;

  // array where the entanglement spectra will be stored
  double** EntanglementSpectra;

  //  maximum number of levels in a entanglement spectrum
  int EntanglementSpectrumDimension;

 public:

  // constructor 
  //
  // space = pointer to the Hilbert space
  // eigenstates = matrix that contains the eigenstates
  // firstEigenstate = first eigenstate to consider
  // lastEigenstate = last eigenstate to consider
  // subsystemSize = number of sites in the subsystem
  // subsystemSz = twice the total Sz value for the subsystem
  // subsystemSites = array that describes the subsystem (if null, consider the subsystemSize first sites)
  SpinChainMultipleEntanglementSpectrumOperation(AbstractSpinChain* space, RealMatrix& eigenstates, int firstEigenstate, int lastEigenstate, 
						 int subsystemSize, int subsystemSz, int* subsystemSites = 0);

  // constructor 
  //
  // space = pointer to the Hilbert space
  // eigenstates = matrix that contains the eigenstates
  // firstEigenstate = first eigenstate to consider
  // lastEigenstate = last eigenstate to consider
  // subsystemSize = number of sites in the subsystem
  // subsystemSz = twice the total Sz value for the subsystem
  // subsystemSites = array that describes the subsystem (if null, consider the subsystemSize first sites)
  SpinChainMultipleEntanglementSpectrumOperation(AbstractSpinChain* space, ComplexMatrix& eigenstates, int firstEigenstate, int lastEigenstate, 
						 int subsystemSize, int subsystemSz, int* subsystemSites = 0);

  // copy constructor 
  //
  // operation = reference on operation to copy
  SpinChainMultipleEntanglementSpectrumOperation(const SpinChainMultipleEntanglementSpectrumOperation& operation);
  
  // destructor
  //
  ~SpinChainMultipleEntanglementSpectrumOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // get the maximum number of levels in a entanglement spectrum
  // 
  // return value = maximum number of levels in a entanglement spectrum
  int GetEntanglementSpectrumDimension();
  
  // get the entanglement spectra
  //
  // return value = array that contains the entanglement spectra (the first index is the entanglement spectrum index, the second index is the entanglement energy index)
  double** GetEntanglementSpectra();

  // set the number of states that have to be locally evaluated
  // 
  // firstState = index of the first eigenstate to evaluate 
  // nbrStates = number of eigenstates that have to be evaluated by the local architecture
  void SetStateRange(int firstState, int nbrStates);

 protected:

  // apply operation(architecture independent)
  //
  // return value = true if no error occurs
  virtual bool RawApplyOperation();
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
};

// get the maximum number of levels in a entanglement spectrum
// 
// return value = maximum number of levels in a entanglement spectrum

inline int SpinChainMultipleEntanglementSpectrumOperation::GetEntanglementSpectrumDimension()
{
  return this->EntanglementSpectrumDimension;
}  

// get the entanglement spectra
//
// return value = array that contains the entanglement spectra (the first index is the entanglement spectrum index, the second index is the entanglement energy index)

inline double** SpinChainMultipleEntanglementSpectrumOperation::GetEntanglementSpectra()
{
  return this->EntanglementSpectra;
}

#endif
