////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of NDMAP hamiltonian precalculation operation           //
//                                                                            //
//                        last modification : 03/11/2003                      //
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


#ifndef NDMAPPRECALCULATIONOPERATION_H
#define NDMAPPRECALCULATIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "Hamiltonian/NDMAPSpinChainHamiltonian.h"


class NDMAPPrecalculationOperation: public AbstractPrecalculationOperation
{

 protected:

  // pointer to the hamiltonian
  NDMAPSpinChainHamiltonian* Hamiltonian;

  // flag to indicate if the operation has to be applied to the first pass of the precalculations
  bool FirstPass;

 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the hamiltonian to use
  // firstPass = flag to indicate if the operation has to be applied to the first pass of the precalculations
  NDMAPPrecalculationOperation(NDMAPSpinChainHamiltonian* hamiltonian, bool firstPass = true);

  // copy constructor 
  //
  // operation = reference on operation to copy
  NDMAPPrecalculationOperation(const NDMAPPrecalculationOperation& operation);
  
  // destructor
  //
  ~NDMAPPrecalculationOperation();
  
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
  
 protected:

  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

// get hilbert space dimension
// 
// return value = hilbert space dimension

inline int NDMAPPrecalculationOperation::GetHilbertSpaceDimension ()
{
  return this->Hamiltonian->GetHilbertSpaceDimension();
}

#endif
