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


#include "config.h"
#include "Architecture/ArchitectureOperation/SUNSpinPrecalculationOperation.h"


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// firstPass = flag to indicate if the operation has to be applied to the first pass of the precalculations

SUNSpinPrecalculationOperation::SUNSpinPrecalculationOperation (AbstractSUNSpinHamiltonian* hamiltonian, bool firstPass)
{
  this->FirstComponent = 0;
  this->NbrComponent = hamiltonian->GetHilbertSpaceDimension();
  this->Hamiltonian = hamiltonian;
  this->OperationType = AbstractArchitectureOperation::QHEParticlePrecalculation;
  this->FirstPass = firstPass;
}

// copy constructor 
//
// operation = reference on operation to copy

SUNSpinPrecalculationOperation::SUNSpinPrecalculationOperation(const SUNSpinPrecalculationOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Hamiltonian = operation.Hamiltonian;
  this->OperationType = AbstractArchitectureOperation::QHEParticlePrecalculation;
  this->FirstPass = operation.FirstPass;
}
  
// destructor
//

SUNSpinPrecalculationOperation::~SUNSpinPrecalculationOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void SUNSpinPrecalculationOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* SUNSpinPrecalculationOperation::Clone()
{
  return new SUNSpinPrecalculationOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool SUNSpinPrecalculationOperation::RawApplyOperation()
{
  if (this->FirstPass ==  true)
    {
      this->Hamiltonian->PartialFastMultiplicationMemory(this->FirstComponent, this->NbrComponent);
    }
  else
    {
      this->Hamiltonian->PartialEnableFastMultiplication(this->FirstComponent, this->NbrComponent);
    }
  return true;
}


