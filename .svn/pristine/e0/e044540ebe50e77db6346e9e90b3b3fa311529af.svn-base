////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of QHE particle hamiltonian precalculation operation        //
//                                                                            //
//                        last modification : 11/03/2003                      //
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


#ifndef QHEPARTICLEPRECALCULATIONOPERATION_H
#define QHEPARTICLEPRECALCULATIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"


class ParticleOnSphereDeltaHamiltonian;


class QHEParticlePrecalculationOperation: public AbstractPrecalculationOperation
{

 protected:

  // pointer to the hamiltonian
  AbstractQHEHamiltonian* Hamiltonian;

  // flag to indicate if the operation has to be applied to the first pass of the precalculations
  bool FirstPass;

  // variable to store memory requirements - non-zero matrix elements
  long RequiredMemory;

 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the hamiltonian to use
  // firstPass = flag to indicate if the operation has to be applied to the first pass of the precalculations
  QHEParticlePrecalculationOperation(AbstractQHEHamiltonian* hamiltonian, bool firstPass = true);

  // copy constructor 
  //
  // operation = reference on operation to copy
  QHEParticlePrecalculationOperation(const QHEParticlePrecalculationOperation& operation);
  
  // destructor
  //
  ~QHEParticlePrecalculationOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // get hilbert space dimension
  // 
  // return value = hilbert space dimension  
  int GetHilbertSpaceDimension ();

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // mpiNodeNbr = provide the additional MPI node ID
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture, int mpiNodeNbr);


  
};

// get hilbert space dimension
// 
// return value = hilbert space dimension

inline int QHEParticlePrecalculationOperation::GetHilbertSpaceDimension ()
{
  return this->Hamiltonian->GetHilbertSpaceDimension();
}

#endif
